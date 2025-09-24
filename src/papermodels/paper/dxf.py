from __future__ import annotations
from copy import deepcopy
from typing import Optional, Any, Union
import ezdxf
import re
from decimal import Decimal
import math
import ezdxf.entities
from papermodels.datatypes.annotation import Annotation
from papermodels.paper.annotations import scale_annotations
from papermodels.geometry import geom_ops
import pathlib
import parse
import numpy as np


def load_dxf_directory(
    directory_path: pathlib.Path | str,
    directory_page_idx: Optional[int] = None
) -> list[Annotation]:
    """
    Returns a list of Annotations representing the annotations in all of the 
    .dxf files within the 'directory_path'. 

    If 'directory_page_idx' is not None, then its value will be applied to the
    .page_idx attribute for all annotations in all files within the directory.
    Otherwise, a page idx will be generated based on the "glob order" of the
    files in the directory and applied to the annotations originating from
    that file.
    """
    dir_path = pathlib.Path(directory_path)
    annotations = []
    for page_idx, dxf_path in enumerate(dir_path.glob("*.dxf")):
        if directory_page_idx is not None:
            page_idx = directory_page_idx
        file_annotations = load_dxf_annotations(dxf_path, page_idx)
        annotations += file_annotations
    return annotations


def load_dxf_annotations(
        dxf_path: pathlib.Path | str,
        page_idx: int = 0
) -> list[Annotation]:
    """
    Returns a lists of pdf annotations keyed by page index.

    'dxf_path': Path-like object representing the path to the PDF file to open.
    'dxf_dir': If provided, a list of paths which to find DXF files comprising a 
        single model.
        The order of the paths is important and will be used to create the order
        of spatial planes in the model, in ascending order (the first path will
        be the lowest plane in the model).
    'annotations_layer': if None, all entities will be treated as annotation entities
    'pages_layer': if None, will treat all entities in the modelspace as
        belonging to the same spatial plane.
        If provided, then entities existing within each page rectangle will be treated
        as belonging to the same plane. Each set of annotation entities within the
        page must also include an "origin" annotation. The origin must be in the same
        place for each "page".
    'pages_layer_order': Optional[str], one of {"ltr", "rtl", "ttb", "btt"}
        left-to-right, right-to-left, top-to-bottom, bottom-to-top
        Which will be used to order the pages in ascending order.
    """
    dxf_path = pathlib.Path(dxf_path)
    doc = ezdxf.readfile(dxf_path)
    layers = doc.layers.entries
    msp = doc.modelspace()
    lines = msp.query("LINE")
    lwpolylines = msp.query("LWPOLYLINE")
    blocks = msp.query("INSERT")
    all_entities = list(lines) + list(lwpolylines) + list(blocks)
    annotations = []
    for local_idx, entity in enumerate(all_entities):
        annotation = dxf_entity_to_annotation(entity, page_idx, local_idx)
        if annotation.vertices:
            annotations.append(annotation)
    return annotations


def dxf_entity_to_annotation(entity: ezdxf.entities.DXFGraphic, page_idx: int, local_idx: int) -> Annotation:
    """
    Converts the entity into an Annotation
    """
    dxf_type = entity.dxftype()
    layer = entity.dxf.layer
    if dxf_type == "LINE":
        object_type = "Line"
        coords = parse_line_coords(entity)
        text = layer
    elif dxf_type == "LWPOLYLINE":
        object_type = "Polygon"
        coords = parse_polyline_coords(entity)
        text = layer
    elif dxf_type == "INSERT":
        object_type = "Polygon"
        block_name = entity.get_dxf_attrib('name', default='')
        coords = parse_block_coordinates(entity)
        text = layer
    else:
        print(f"{dxf_type=}")
    line_color = (0, 0, 0)#entity.dxf.color # convert to RBG tuple
    line_type = None #entity.dxf.linetype
    line_weight = 1.0#entity.dxf.thickness
    # transparency = 1.0 #entity.dxf.transparency or 1.0
    opacity = 1.0 #- transparency
    vertices = coords_to_vertices_list(coords)
    return Annotation(
        page=page_idx,
        object_type=object_type,
        text=text,
        vertices=vertices,
        line_color=line_color,
        fill_color = None,
        line_type=line_type,
        line_weight=line_weight,
        line_opacity=opacity,
        fill_opacity=opacity,
        matrix=(1, 0, 0, 1, 0, 0),
        local_id=local_idx

    )


def parse_block_coordinates(entity):
    """
    Extract sequenced coordinates from the block.
    """
    geoms = []
    for e in entity.virtual_entities():
        coords = None
        if e.dxftype() == 'LINE':
            coords = parse_line_coords(e)
        elif e.dxftype() == 'LWPOLYLINE':
            coords = parse_polyline_coords(e)
        elif e.dxftype() == 'ARC':
            coords = parse_arc_coords(e)
        
        if coords is not None:
            geoms += coords
    return geoms


def parse_line_coords(entity: ezdxf.entities.DXFGraphic):
    coords = [
        (entity.dxf.start[0], entity.dxf.start[1]), 
        (entity.dxf.end[0], entity.dxf.end[1])
    ]
    return coords


def parse_polyline_coords(entity: ezdxf.entities.DXFGraphic):
    # point_pairs = zip(entity.get_points(), entity.get_points()[1:])
    # coords = [[(p[0][0], p[0][1]), (p[1][0], p[1][1])] for p in point_pairs]
    coords = [(p[0], p[1]) for p in entity.get_points()]

    return coords


def parse_arc_coords(arc: ezdxf.entities.DXFGraphic, num_segments=12):
    """
    Approximate an ARC entity as a series of linear segments
    """
    center = arc.dxf.center
    radius = arc.dxf.radius
    start_angle = math.radians(arc.dxf.start_angle)
    end_angle = math.radians(arc.dxf.end_angle)
    if end_angle < start_angle:
        end_angle += 2 * math.pi
    
    angles = [start_angle + (end_angle - start_angle) * i / (num_segments - 1) for i in range(num_segments)]
    points = [(center[0] + radius * math.cos(a), center[1] + radius * math.sin(a)) for a in angles]
    return points


def coords_to_vertices_list(coords: list[tuple[float, float]]) -> list:
    """
    Formats the ordered pairs of coordinates in 'coords' as a flattened
    list of coordinates in this order [x0, y0, x1, y1, x2, y2, ..., xn, yn]
    """
    vertices = []
    for coord in coords:
        x, y = coord[0], coord[1]
        vertices.append(Decimal(x))
        vertices.append(Decimal(y))
    return tuple(vertices)