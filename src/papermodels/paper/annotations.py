from __future__ import annotations
from copy import deepcopy
from decimal import Decimal
from dataclasses import asdict
from shapely.wkt import loads as wkt_loads
from shapely import Geometry, GeometryCollection, Point, Polygon
from papermodels.datatypes.annotation import Annotation
from papermodels.loads.load_distribution import LoadingGeometry
from papermodels.geometry import geom_ops
from papermodels.datatypes.exceptions import LegendError, GeometryError, AnnotationError
from typing import Any, Optional
import re
import numpy as np
from numpy.typing import ArrayLike



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


def parsed_annotations_to_loading_geometry(
    parsed_annots: dict[Annotation, dict],
) -> list[LoadingGeometry]:
    """
    Convert annotations representing loading areas into a list of LoadingArea
    """
    acc = []
    for annot, annot_attrs in parsed_annots.items():
        lg = LoadingGeometry(
            geometry=annot_attrs['geometry'],
            occupancy=annot_attrs.get("occupancy", None),
            load_components=annot_attrs.get("components", None),
            plane_id=annot_attrs['page_label'],
        )
        acc.append(lg)
    return acc


def parse_annotations(
    annots: list[Annotation], legend: list[Annotation], legend_identifier: str
) -> dict[Annotation, dict]:
    """
    Returns a dictionary of annotations organized by their legend entry. If the annotation type is not
    in the 'legend', then it is excluded from the results.
    """
    # TODO: Make this more configurable by the user
    properties_to_match = [
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
        legend_text = strip_html_tags(legend_item.text)
        legend_data = legend_text.lower().replace("\r\n", "\n").replace("\r", "\n").replace(f"{legend_identifier.lower()}\n", "").split("\n")
        legend_data = [elem for elem in legend_data if elem]
        annot_attributes = {}
        for legend_attr in legend_data:
            try:
                key, value = legend_attr.split(":")
            except ValueError as e:
                raise LegendError(f"Incorrect legend format on the following annotation: {legend_data}")
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()
            annot_attributes.update({key: value})
        for annot in matching_annots:
            if annot in legend: 
                continue
            existing_annot_tag = parse_existing_annot_tag(annot.text)
            annot_geom = annotation_to_shapely(annot)
            annot_attrs = {}
            annot_attrs["geometry"] = annot_geom
            annot_attrs['page_label'] = annot.page
            annot_attrs['tag'] = existing_annot_tag
            for annot_key, annot_attr in annot_attributes.items():
                annot_attrs[annot_key] = str_to_int(
                    annot_attr.split("<")[0]
                )  # .split() to remove trailing HTML tags
            # annot_attrs["rank"] = int(annot_attributes["rank"])
            annot_attrs.setdefault("reaction_type", "point")
            annot_attrs['reaction_type'] = annot_attrs['reaction_type'].lower()
            if annot_geom.geom_type == "Polygon" and annot_attrs['reaction_type'] == "linear":
                annot_attrs['length'] = geom_ops.get_rectangle_centerline(annot_geom).length
            parsed_annotations.update({annot: annot_attrs})
    return parsed_annotations


def tag_parsed_annotations(
    parsed_annots: dict[Annotation, dict],
) -> dict[Annotation, dict]:
    """
    Adds an identifying tag to the annotation based on the page number of the annotation and
    its identified type. Prioritizes manually named tags according to thsi format: 
    "{ABBR}{PAGE_ID}.{INDEX}". If an annotation is not already tagged accordingly then
    it will be auto-assigned an integer index after all tagged annotations have been accounted
    for.

    All tags are guaranteed to be unique.
    """
    counts = {}
    annots_to_tag = parsed_annots.copy()
    annots_to_enumerate = {}
    for annot, annot_attrs in annots_to_tag.items():
        tag = annot_attrs.get('tag')
        # There is an existing annotation
        if tag is not None:
            parsed = parse_tag_components(tag)
            if parsed is not None and len(parsed) == 3: # Tag is not in correct format so ignore
                type_initials, page_id, tag_idx = parsed
                tag_prefix = f"{type_initials}{page_id}"
                counts.setdefault(tag_prefix, set())
                counts[tag_prefix].add(int(tag_idx))
            else:
                annots_to_enumerate.update({annot: annot_attrs})
        else:
            annots_to_enumerate.update({annot: annot_attrs})


    for annot, annot_attrs in annots_to_enumerate.items():
        type_initials = "".join(
            [label[0].upper() for label in annot_attrs["type"].split(" ")]
        )
        tag_prefix = f"{type_initials}{annot.page}"
        counts.setdefault(tag_prefix, set())
        count = 0
        while count in counts[tag_prefix]:
            count += 1
        tag = f"{tag_prefix}.{count}"
        annot_attrs['tag'] = tag
        counts[tag_prefix].add(count)
        prev_count = count

    return annots_to_tag


def parse_tag_components(tag: str) -> tuple[str, int, int]:
    """
    Returns a tuple of the tag components: type abbrev., page_id, tag_idx
    """
    pattern = re.compile(r"([A-Za-z]+)([0-9]+).([0-9]+)")
    results = pattern.search(tag)
    if results is not None:
        return results.groups()


def _annotation_to_wkt(annot: Annotation) -> str:
    """
    Returns a WKT string representing the geometry in 'annot'. The WKT
    string can be loaded with shapely.wkt.loads (see shapely documentation)
    """
    if annot.object_type == "PolyLine" or annot.object_type == "Line":
        grouped_vertices = geom_ops._group_vertices_str(annot.vertices)
        return f"LINESTRING({grouped_vertices})"
    elif annot.object_type in ("Polygon", "Rectangle", "Square"):
        grouped_vertices = geom_ops._group_vertices_str(annot.vertices, close=True)
        return f"POLYGON(({grouped_vertices}))"
    
    
def assign_page_id_to_annotations(annots: list[Annotation], page_containers: list[Annotation], left_to_right: bool = True) -> list[list[Annotation]]:
    """
    Returns a list of Annotation representing 'annots' after they have been assigned to their page container
    in 'page_containers'.
    """
    page_annot_map = enumerate_page_annotations(page_containers, left_to_right)
    sorted_annots = sort_annotations_by_page_polygon(annots, page_annot_map)
    aligned_annotations = align_annotations_to_pages(sorted_annots)
    page_indexed_annots = []
    for page_id, annots_on_page in enumerate(aligned_annotations.values()):
        page_indexed_annots.append([])
        for annot in annots_on_page:
            new_annotation = Annotation(
                page=page_id,
                object_type=annot.object_type,
                text=annot.text,
                vertices=annot.vertices,
                line_color=annot.line_color,
                fill_color=annot.fill_color,
                line_type=annot.line_type,
                line_weight=annot.line_weight,
                line_opacity=annot.line_opacity,
                fill_opacity=annot.fill_opacity,
                matrix=annot.matrix,
                local_id=annot.local_id
            )
            page_indexed_annots[page_id].append(new_annotation)
    return page_indexed_annots


def enumerate_page_annotations(page_annots: list[Annotation], left_to_right: bool = True) -> dict[Polygon, int]:
    """
    Returns the annotations in 'annots' that correspond to page-demarcation
    polygons organized in a seqential fashion.
    """
    page_geoms_by_page = get_page_geom_by_page_index(page_annots)
    enumerated_page_annotations = {}
    counter = 0
    for page_geoms in page_geoms_by_page.values():
        reverse=False
        if not left_to_right:
            reverse=True
        sorted_page_geoms = sorted(page_geoms, key=lambda x: x.centroid.coords[0], reverse=reverse)
        for page_geom in sorted_page_geoms:
            enumerated_page_annotations.update({counter: page_geom})
            counter += 1
    return enumerated_page_annotations


def sort_annotations_by_page_polygon(annots: list[Annotation], page_geom_map: dict[int, Polygon]) -> dict[Polygon, list[Annotation]]:
    """
    Sorts each annotation into its own list[Annotation] coresponding to which page they are contained in.
    """
    acc = {}
    for annot in annots:
        annot_geom = annotation_to_shapely(annot)
        for page_poly in page_geom_map.values():
            acc.setdefault(page_poly, [])
            if page_poly.contains(annot_geom):
                acc[page_poly].append(annot)
    return acc


def align_annotations_to_pages(
    annotations_by_page: dict[Polygon, list[Annotation]],
) -> list[list[Annotation]]:
    """
    Returns a copy of 'annotations_by_page' but with all annotations aligned
    within their respective pages by their origin points. The distance between
    the page origin points and the page corner is set by the distance of the
    first origin on the first page.
    """
    acc = {}
    counter = 0
    for page_poly, page_annots in annotations_by_page.items():
        try:
            origin_annot = next((page_annot for page_annot in page_annots if page_annot.text == "origin"))
        except StopIteration:
            raise AnnotationError(f"The 'page' annotation with index={counter} appears to be missing an 'origin' annotation.")
        if counter == 0:
            global_offset = get_origin_offset(origin_annot, page_poly)
        origin_centroid = get_origin_centroid(origin_annot)
        page_origin = get_page_bottom_left_corner(page_poly)
        shifted_annots = reset_annotations_to_origin(page_annots, page_origin, origin_centroid, global_offset)
        acc.update({page_poly: shifted_annots})
        counter += 1
    return acc
        

def reset_annotations_to_origin(
    annots_on_page: list[Annotation],
    page_origin_xy: ArrayLike,
    origin_centroid_xy: ArrayLike,
    global_offset: ArrayLike
) -> list[Annotation]:
    """
    Translates all of the annotations in 'annots_on_page' so that the origin
    annotation (contained within 'annots_on_page') is located at 'xy_offset' 
    from the bottom-left corner.
    """
    updated_annots = []
    local_offset = origin_centroid_xy - page_origin_xy
    local_delta = local_offset - global_offset
    translation_vector = page_origin_xy - global_offset + local_delta
    for annot in annots_on_page:
        vertices = annot.vertices
        translated_vertices = geom_ops._translate_vertices(vertices, offset_x=-translation_vector[0], offset_y=-translation_vector[1])
        updated_annot = Annotation(
            annot.page,
            annot.object_type,
            text=annot.text,
            vertices=translated_vertices,
            line_color=annot.line_color,
            fill_color=annot.fill_color,
            line_type=annot.line_type,
            line_weight=annot.line_weight,
            line_opacity=annot.line_opacity,
            fill_opacity=annot.fill_opacity,
            matrix=annot.matrix,
            local_id=annot.local_id
        )
        updated_annots.append(updated_annot)
    return updated_annots
        

def get_origin_offset(origin_annot: Annotation, page_poly: Polygon) -> ArrayLike:
    """
    Returns the xy offset of the 'origin_annot' to the bottom-left corner of the 'page_annot'.
    """
    origin_centroid = get_origin_centroid(origin_annot)
    page_origin = get_page_bottom_left_corner(page_poly)
    return origin_centroid - page_origin


def get_origin_centroid(origin_annot: Annotation) -> ArrayLike:
    """
    Returns the centroid of the 'origin_annot' as an xy tuple.
    """
    poly = annotation_to_shapely(origin_annot)
    return np.array(poly.centroid.coords[0])


def get_page_bottom_left_corner(page_poly: Polygon) -> ArrayLike:
    """
    Returns the coordinate fo the bottom left hand corner of the 'page_annot'
    (a rectangle)
    """
    minx, miny, _, __ = page_poly.bounds
    page_point = Point(minx, miny)
    return np.array(page_point.coords[0])


def get_page_geom_by_page_index(page_annots: list[Annotation]) -> dict[int, list[Polygon]]:
    """
    Returns a dictionary that organizes the 'page_annots' by their page id.
    """
    page_geoms = {}
    for page_annot in page_annots:
        page_geom = annotation_to_shapely(page_annot)
        page_id = page_annot.page
        page_geoms.setdefault(page_id, [])
        page_geoms[page_id].append(page_geom)
    return page_geoms


def parse_existing_annot_tag(text_contents: str) -> Optional[str]:
    """
    Returns a tag if there is a tag field existing within the annotation's
    text field, 'text_contents'. Returns None, otherwise.

    A tag field looks like: "tag: <tag name>".
    """
    index = text_contents.lower().find("tag:")
    if index == -1:
        return None
    tag_text = text_contents[index:]
    tag_value = re.search(r"^tag:[\s]*([A-Za-z.0-9\-]+)", tag_text).groups()[0]
    return tag_value


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
    scale: Decimal,
    paper_origin: Optional[tuple[float, float]] = None,
    round_precision: int = 4
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
        scaled_vertices = geom_ops.scale_vertices(annot.vertices, scale, round_precision=round_precision, paper_origin=paper_origin)
        annot_dict["vertices"] = scaled_vertices
        scaled_annotations.append(Annotation(**annot_dict))
    return scaled_annotations


def strip_html_tags(s: str) -> str:
    """
    Removes but does not sanitize HTML tags from strings
    """
    return re.sub('<[^<]+?>', '\n', s)

def str_to_int(s: str) -> int | str:
    try:
        return int(s)
    except ValueError:
        return s
