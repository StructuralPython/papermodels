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
import shapely as shp
import pathlib
import parse
import numpy as np


def load_dxf_directory(
    directory_path: pathlib.Path | str, directory_page_idx: Optional[int] = None
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
    page_idx: int = 0,
    polygonize_layers: Optional[list[str]] = None,
) -> list[Annotation]:
    dxf_entities = dxf_file_to_dxf_entities(
        dxf_path, page_idx, polygonize_layers=polygonize_layers
    )
    annotations = dxf_entities_to_annotations(dxf_entities)
    return annotations


def dxf_file_to_dxf_entities(
    dxf_path: pathlib.Path | str,
    page_idx: int = 0,
    polygonize_layers: Optional[list[str]] = None,
) -> dict:
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
    'polygonize_layers': A list of layers if passed will polygonize the individual
        linestrings on these layers.
    """
    dxf_path = pathlib.Path(dxf_path)
    doc = ezdxf.readfile(dxf_path)
    layers = doc.layers.entries
    msp = doc.modelspace()
    lines = msp.query("LINE")
    lwpolylines = msp.query("LWPOLYLINE")
    blocks = msp.query("INSERT")
    all_entities = list(lines) + list(lwpolylines) + list(blocks)
    dxf_entities = {}
    filtered_entities = {}
    for entity in all_entities:
        dxf_type = entity.dxftype()
        layer = entity.dxf.layer

        if dxf_type == "LINE":
            object_type = "Line"
            geom = parse_line_coords(entity)
            text = layer
        elif dxf_type == "LWPOLYLINE":
            object_type = "Polygon"
            geom = parse_polyline_coords(entity)
            text = layer
        elif dxf_type == "INSERT":
            object_type = "Polygon"
            block_name = entity.get_dxf_attrib("name", default="")
            geom = parse_block_coordinates(entity)
            text = layer
        else:
            print(f"{dxf_type=}")
        dxf_entities["page_idx"] = page_idx
        dxf_entities.setdefault(layer, [])
        dxf_entities[layer].append(geom)
    poly_layers = []
    if polygonize_layers is not None:
        poly_layers = [f"{poly_layer}".lower() for poly_layer in polygonize_layers]
    lower_layers = [l.lower() for l in layers]
    if polygonize_layers is not None and set(poly_layers) & set(lower_layers):
        filtered_entities = {}
        filtered_entities["page_idx"] = page_idx
        for layer, geoms in dxf_entities.items():
            if layer == "page_idx" or layer.lower() not in poly_layers:
                continue
            ls = [geom for geom in geoms if geom.geom_type == "LineString"]
            others = [geom for geom in geoms if geom.geom_type != "LineString"]
            line_unions = shp.polygonize(ls)
            polys = shp.MultiPolygon(
                [lu for lu in line_unions.geoms if lu.geom_type == "Polygon"]
            )
            dxf_polys = [poly for poly in polys.geoms if geom.geom_type == "Polygon"]
            filtered_entities.update({layer: dxf_polys + others})

    dxf_entities.update(filtered_entities)
    return dxf_entities


def dxf_entity_to_shapely(
    entity: ezdxf.entities.DXFGraphic,
    page_idx: int,
    local_idx: int,
    polygonize_lines: bool = False,
) -> Annotation:
    """
    Converts the entity into an Annotation
    """
    dxf_type = entity.dxftype()
    layer = entity.dxf.layer
    dxf_entities = {}
    if dxf_type == "LINE":
        object_type = "Line"
        geom = parse_line_coords(entity)
        text = layer
    elif dxf_type == "LWPOLYLINE":
        object_type = "Polygon"
        geom = parse_polyline_coords(entity)
        text = layer
    elif dxf_type == "INSERT":
        object_type = "Polygon"
        block_name = entity.get_dxf_attrib("name", default="")
        geom = parse_block_coordinates(entity)
        text = layer
    else:
        print(f"{dxf_type=}")
    dxf_entities["page_idx"] = page_idx
    dxf_entities.setdefault(layer, [])
    dxf_entities[layer].append(geom)
    if polygonize_lines:
        filtered_entities = {}
        filtered_entities["page_idx"] = page_idx
        for layer, geoms in dxf_entities.items():
            ls = [geom for geom in geoms if geom.geom_type == "LineString"]
            others = [geom for geom in geoms if geom.geom_type != "LineString"]
            polys = shp.polygonize(ls)
            dxf_polys = [poly for poly in polys if geom.geom_type == "Polygon"]
            filtered_entities.update({layer: dxf_polys + others})
        return filtered_entities
    return dxf_entities


def dxf_entities_to_annotations(
    dxf_entities: dict[str, list | str],
) -> list[Annotation]:
    """
    Returns a list of annotation for the dxf-entities
    """
    # dxf_path = pathlib.Path(dxf_path)
    # doc = ezdxf.readfile(dxf_path)
    # layers = doc.layers.entries
    # msp = doc.modelspace()
    page_idx = dxf_entities.pop("page_idx")
    annotations = []
    counter = 0
    for layer, geoms in dxf_entities.items():
        for geom in geoms:
            if geom.geom_type == "Polygon":
                vertices = geom_ops.flatten_vertex_array(np.array(geom.exterior.coords))
            else:
                vertices = geom_ops.flatten_vertex_array(np.array(geom.coords))
            line_color = (0, 0, 0)  # entity.dxf.color # convert to RBG tuple
            line_type = None  # entity.dxf.linetype
            line_weight = 1.0  # entity.dxf.thickness
            transparency = 1.0  # entity.dxf.transparency or 1.0
            opacity = 0.5  # - transparency
            annot = Annotation(
                page=page_idx,
                object_type=geom.geom_type,
                text=layer,
                vertices=vertices,
                line_color=line_color,
                fill_color=None,
                line_type=line_type,
                line_weight=line_weight,
                line_opacity=opacity,
                fill_opacity=opacity,
                matrix=(1, 0, 0, 1, 0, 0),
                local_id=counter,
            )
            annotations.append(annot)
            counter += 1
    return annotations


def parse_block_coordinates(entity):
    """
    Extract sequenced coordinates from the block.
    """
    geoms = []
    for e in entity.virtual_entities():
        geom = None
        if e.dxftype() == "LINE":
            geom = parse_line_coords(e)
        elif e.dxftype() == "LWPOLYLINE":
            geom = parse_polyline_coords(e)
        elif e.dxftype() == "ARC":
            geom = parse_arc_coords(e)
        geoms.append(geom)
    return geoms


def parse_line_coords(entity: ezdxf.entities.DXFGraphic):
    coords = [
        (entity.dxf.start[0], entity.dxf.start[1]),
        (entity.dxf.end[0], entity.dxf.end[1]),
    ]
    return shp.set_precision(shp.LineString(coords), grid_size=1e-3)


def parse_polyline_coords(entity: ezdxf.entities.DXFGraphic):
    # point_pairs = zip(entity.get_points(), entity.get_points()[1:])
    # coords = [[(p[0][0], p[0][1]), (p[1][0], p[1][1])] for p in point_pairs]
    coords = [(p[0], p[1]) for p in entity.get_points()]

    return shp.set_precision(shp.Polygon(coords), grid_size=1e-3)


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

    angles = [
        start_angle + (end_angle - start_angle) * i / (num_segments - 1)
        for i in range(num_segments)
    ]
    points = [
        (center[0] + radius * math.cos(a), center[1] + radius * math.sin(a))
        for a in angles
    ]
    return shp.set_precision(shp.LineString(points), grid_size=1e-3)


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
