from __future__ import annotations
from shapely.wkt import loads as wkt_loads
from shapely.geometry import GeometryCollection
from papermodels.db.data_model import Annotation
from typing import Any

def annotations_to_shapely(annots: list[Annotation], as_geometry_collection=False) -> list[Any]:
    """
    Returns a WKT string representing the geometry in 'annot'
    'annots' - a list of Annotation objects
    'as_geometry_collection' - If 
    """
    geoms = [wkt_loads(_annotation_to_wkt(annotation)) for annotation in annots]
    if as_geometry_collection:
        return GeometryCollection(geoms)
    return geoms

def _annotation_to_wkt(annot: Annotation) -> str:
    """
    Returns a WKT string representing the geometry in 'annot'
    """
    if annot.object_type == "PolyLine" or annot.object_type == "Line":
        grouped_vertices = _group_vertices(annot.vertices)
        return f"LINESTRING({grouped_vertices})"
    elif annot.object_type == "Polygon" or annot.object_type == "Rectangle":
        grouped_vertices = _group_vertices(annot.vertices, close = True)
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
        if (
            annot.__dict__.items() & properties.items()
        ) == properties.items():
            filtered.append(annot)
    return filtered


def scale_annotation(annot: Annotation, scale: float) -> str:
    """
    Scale the annotation. Each vertex in 'annot' will be multiplied
    by 'scale'
    """
    return

    
def _group_vertices(vertices: str, close = False) -> list[tuple[float, float]]:
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