from __future__ import annotations
from shapely.wkt import loads as wkt_loads
from shapely.geometry import GeometryCollection
from papermodels.datatypes.annotation import Annotation
from typing import Any, Optional
import numpy as np


def annotations_to_shapely(
    annots: list[Annotation], as_geometry_collection=False
) -> list[Any]:
    """
    Returns a WKT string representing the geometry in 'annot'
    'annots' - a list of Annotation objects
    'as_geometry_collection' - If
    """
    geoms = [annotation_to_shapely(annot) for annot in annots]
    if as_geometry_collection:
        return GeometryCollection(geoms)
    return geoms


def annotation_to_shapely(annot: Annotation) -> Any:
    """
    Returns an shapely geometry created from the annotation type and
    vertices in 'annot'.
    """
    return wkt_loads(_annotation_to_wkt(annot))


def _annotation_to_wkt(annot: Annotation) -> str:
    """
    Returns a WKT string representing the geometry in 'annot'
    """
    if annot.object_type == "PolyLine" or annot.object_type == "Line":
        grouped_vertices = _group_vertices_str(annot.vertices)
        return f"LINESTRING({grouped_vertices})"
    elif annot.object_type == "Polygon" or annot.object_type == "Rectangle":
        grouped_vertices = _group_vertices_str(annot.vertices, close=True)
        return f"POLYGON(({grouped_vertices}))"


def filter_annotations(annots: list[Annotation], properties: dict) -> list[Annotation]:
    """
    Returns a list of annotations from 'annots' that have properties that match
    the keywords in 'properties'.
    Note: The filtering process currently requires that both the keys AND values in 'properties'
    be hashable.
    """
    filtered = []
    for annot in annots:
        if (annot.__dict__.items() & properties.items()) == properties.items():
            filtered.append(annot)
    return filtered


def scale_annotations(
    annots: list[Annotation],
    scale: float,
    paper_origin: Optional[tuple[float, float]] = None,
) -> list[Annotation]:
    """
    Scale the annotations in 'annots'. Each vertex in each annotation in 'annots' will be multiplied
    by 'scale'.
    If 'paper_origin' is provided, then the annotation coordinates will have their origin reset
    to 'paper_origin'. Note that 'paper_origin' is the unscaled coordinate space (i.e. in points)
    """
    scaled_annotations = []
    for annot in annots:
        scaled_annotations.append(scale_annotation(annot, scale, paper_origin))
    return scaled_annotations


def scale_annotation(
    annot: Annotation, scale: float, paper_origin: Optional[tuple[float, float]] = None
) -> Annotation:
    """
    Scale the annotation. Each vertex in 'annot' will be multiplied
    by 'scale'.
    If 'paper_origin' is provided, then the annotation coordinates will have their origin reset
    to 'paper_origin'. Note that 'paper_origin' is the unscaled coordinate space (i.e. in points)
    """
    vertices = annot.vertices
    if paper_origin is not None:
        offset_x = paper_origin[0]
        offset_y = paper_origin[1]
        vertices = _translate_vertices(vertices, offset_x, offset_y)

    scaled_vertices = [vertex * scale for vertex in vertices]
    annot.vertices = scaled_vertices
    return annot


def _translate_vertices(
    vertices: list[float], offset_x: float, offset_y: float
) -> Annotation:
    """
    Returns a list of float representing 'verticies' translated by 'offset_x' and 'offset_y'.
    """
    coord_array = np.array(_group_vertices(vertices))
    offset_array = np.array([offset_x, offset_y])
    translated_array = coord_array + offset_array
    flattened_array = translated_array.flatten()
    return list(flattened_array)


def _group_vertices(vertices: str, close=False) -> list[tuple[float, float]]:
    """
    Returns a list of (x, y) tuples from a list of vertices in the format of:
    'x1 y1 x2 y2 x3 y3 ... xn yn'
    """
    grouped_vertices = []
    coordinates = []
    for idx, ordinate in enumerate(vertices):
        if idx % 2:
            coordinates.append(ordinate)
            grouped_vertices.append([coordinates])
            coordinates = []
        else:
            coordinates.append(ordinate)
    if close:
        grouped_vertices.append(grouped_vertices[0])
    return grouped_vertices


def _group_vertices_str(vertices: str, close=False) -> str:
    """
    Returns a list of (x, y) tuples from a list of vertices in the format of:
    'x1 y1 x2 y2 x3 y3 ... xn yn'
    """
    acc = []
    coordinates = []
    for idx, ordinate in enumerate(vertices):
        if idx % 2:
            coordinates.append(f"{ordinate}")
            acc.append(" ".join(coordinates))
            coordinates = []
        else:
            coordinates.append(f"{ordinate}")
    if close:
        acc.append(acc[0])
    return ", ".join(acc)
