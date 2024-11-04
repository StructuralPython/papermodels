from __future__ import annotations
from copy import deepcopy
from dataclasses import asdict
from shapely.wkt import loads as wkt_loads
from shapely import Geometry, GeometryCollection, Point
from papermodels.datatypes.annotation import Annotation
from papermodels.datatypes.element import Element
from papermodels.geometry.geom_ops import get_intersection
from papermodels.fileio.utils import str_to_int
from typing import Any, Optional
import numpy as np


def annotations_to_shapely(
    annots: list[Annotation], as_geometry_collection=False
) -> list[Any]:
    """
    Returns a shapely geometry representing the geometry in 'annot'
    'annots' - a list of Annotation objects
    'as_geometry_collection' - If True, the geometries in 'annots'
        will all be grouped together in a shapely.geometry.GeometryCollection
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


def get_annotation_geometry_pairs(
    annots: list[Annotation],
) -> dict[Annotation, Geometry]:
    """
    Returns a dict of annotation, shapely geometry pairs
    """
    return {annot: annotation_to_shapely(annot) for annot in annots}


def parse_annotations(
    annots: list[Annotation], legend: list[Annotation]
) -> dict[Annotation, dict]:
    """
    Returns a dictionary of annotations organized by their legend entry. If the annotation type is not
    in the 'legend', then it is excluded from the results.
    """
    # TODO: Make this more configurable by the user
    properties_to_match = [
        "object_type",
        "line_color",
        "fill_color",
        "line_type",
        "line_weight",
    ]
    parsed_annotations = {}
    for legend_item in legend:
        legend_properties = {
            prop: getattr(legend_item, prop) for prop in properties_to_match
        }
        matching_annots = filter_annotations(annots, legend_properties)
        legend_data = legend_item.text.replace("Legend\n", "").split("\n")
        annot_attributes = {
            legend_attr.split(": ")[0].lower(): legend_attr.split(": ")[1]
            for legend_attr in legend_data
        }
        for annot in matching_annots:
            annot_geom = annotation_to_shapely(annot)
            annot_attrs = {}
            annot_attrs["geometry"] = annot_geom
            for annot_key, annot_attr in annot_attributes.items():
                annot_attrs[annot_key] = str_to_int(
                    annot_attr.split("<")[0]
                )  # .split() to remove trailing HTML tags
            # annot_attrs["rank"] = int(annot_attributes["rank"])
            parsed_annotations.update({annot: annot_attrs})
    return parsed_annotations


def tag_parsed_annotations(
    parsed_annots: dict[Annotation, dict]
) -> dict[Annotation, dict]:
    """
    Adds an identifying tag to the annotation based on the page number of the annotation and
    its identified type.
    """
    counts = {}
    annots_to_tag = parsed_annots.copy()
    for annot, annot_attrs in annots_to_tag.items():
        type_initials = "".join(
            [label[0].upper() for label in annot_attrs["type"].split(" ")]
        )
        tag_prefix = f"{type_initials}{annot.page}"
        if tag_prefix not in counts:
            counts[tag_prefix] = 1
        count = counts[tag_prefix]
        tag = f"{tag_prefix}.{count}"
        annot_attrs["tag"] = tag
        counts[tag_prefix] += 1
    return annots_to_tag


def get_structural_elements(
    annots: list[Annotation], legend: list[Annotation]
) -> list[Element]:
    """
    Returns a list of Element generated from the annotations in 'annots' according to the element
    types described in the 'legend'. If an annotation is not described in the legend then it will
    not be included in the result list of Elements.
    """
    sorted_by_page_annotations = sorted(annots, key=lambda x: x.page, reverse=True)
    parsed_annotations = parse_annotations(sorted_by_page_annotations, legend)
    tagged_annotations = tag_parsed_annotations(parsed_annotations)
    intersecting_annotations = get_geometry_intersections(tagged_annotations)
    corresponding_annotations = get_geometry_correspondents(intersecting_annotations)
    elements = []
    for annot, annot_attrs in corresponding_annotations.items():
        element = Element(
            tag=annot_attrs["tag"],
            type=annot_attrs["type"],
            page=annot.page,
            geometry=annot_attrs["geometry"],
            intersections=tuple(annot_attrs["intersections"]),
            correspondents=tuple(annot_attrs["correspondents"]),
            page_label=annot_attrs.get("page_label", None),
        )
        elements.append(element)
    return elements


def get_geometry_intersections(
    tagged_annotations: dict[Annotation, dict]
) -> dict[Annotation, dict]:
    """
    Returns a dictionary of
    """
    annots = list(tagged_annotations.keys())
    intersected_annotations = tagged_annotations.copy()
    for i_annot in annots:
        i_attrs = intersected_annotations[i_annot]
        i_rank = i_attrs["rank"]
        i_page = i_annot.page
        intersections = []
        for j_annot in annots:
            j_attrs = intersected_annotations[j_annot]
            j_rank = j_attrs["rank"]
            j_page = j_annot.page
            if i_rank < j_rank and i_page == j_page:
                i_geom = i_attrs["geometry"]
                j_geom = j_attrs["geometry"]
                intersection = get_intersection(i_geom, j_geom, j_attrs["tag"])
                if intersection is not None:
                    intersections.append(intersection)
        i_attrs["intersections"] = intersections
    return intersected_annotations


def get_geometry_correspondents(
    tagged_annotations: dict[Annotation, dict]
) -> dict[Annotation, dict]:
    """
    Returns a copy of 'tagged_annotations' with a 'correspondents' field added to that
    attributes dictionary of each Annotation key.
    """
    annots_by_page = annotations_by_page(tagged_annotations)
    descending_pages = sorted(annots_by_page.keys(), reverse=True)
    last_page = descending_pages[-1]
    corresponding_annotations = {}
    for page in descending_pages:
        if page != last_page:
            next_page = page - 1
            annots_here = annots_by_page[page]
            annots_below = annots_by_page[next_page]
            for i_annot, i_attrs in annots_here.items():
                i_page = i_annot.page
                correspondents = []
                for j_annot, j_attrs in annots_below.items():
                    j_attrs = annots_below[j_annot]
                    j_page = j_annot.page
                    i_geom = i_attrs["geometry"]
                    j_geom = j_attrs["geometry"]
                    if j_page in (i_page + 1, i_page - 1) and i_geom.contains(j_geom):
                        correspondents.append(j_attrs["tag"])
                i_attrs["correspondents"] = correspondents
                corresponding_annotations.update({i_annot: i_attrs})
        else:
            annots_here = annots_by_page[page]
            for i_annot, i_attrs in annots_here.items():
                i_attrs["correspondents"] = []
                corresponding_annotations.update({i_annot: i_attrs})
    return corresponding_annotations


def _annotation_to_wkt(annot: Annotation) -> str:
    """
    Returns a WKT string representing the geometry in 'annot'. The WKT
    string can be loaded with shapely.wkt.loads (see shapely documentation)
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
    'properties' is a dictionary of annotation properties and their values, e.g.
        {'line_weight': 3.0, 'line_color': (1, 0, 0)}
        or
        {'text': "Slab Outline"}
    The returned annotations will only be annotations that match ALL of the properties
    described.
    """
    filtered = []
    for annot in annots:
        if (asdict(annot).items() & properties.items()) == properties.items():
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
    The purpose of setting 'paper_origin' if the annotation has a "datum" that is set somewhere in the
    file. Leave as None if there is no datum set.
    """
    scaled_annotations = []
    for annot in annots:
        annot_dict = asdict(annot).copy()
        scaled_vertices = scale_vertices(annot.vertices, scale)
        annot_dict["vertices"] = scaled_vertices
        scaled_annotations.append(Annotation(**annot_dict))
    return scaled_annotations


def scale_vertices(
    vertices: list[float],
    scale: float,
    paper_origin: Optional[tuple[float, float]] = None,
) -> Annotation:
    """
    Scale the annotation. Each vertex in 'annot' will be multiplied
    by 'scale'.
    If 'paper_origin' is provided, then the annotation coordinates will have their origin reset
    to 'paper_origin'. Note that 'paper_origin' is the unscaled coordinate space (i.e. in points)
    """
    if paper_origin is not None:
        offset_x = paper_origin[0]
        offset_y = paper_origin[1]
        vertices = _translate_vertices(vertices, offset_x, offset_y)

    scaled_vertices = [vertex * scale for vertex in vertices]
    return scaled_vertices


def annotations_by_page(
    annots: dict[Annotation, dict], ascending=False
) -> dict[int, dict[Annotation, dict]]:
    """
    Returns 'annots' in a dictionary keyed by page number
    """
    annots_by_page = {}
    for annot, annot_attrs in annots.items():
        annots_on_page = annots_by_page.get(annot.page, {})
        annots_on_page.update({annot: annot_attrs})
        annots_by_page[annot.page] = annots_on_page
    return annots_by_page


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
