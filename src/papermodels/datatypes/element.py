from dataclasses import dataclass
from typing import Optional
from shapely import Point, Geometry, LineString, Polygon
import math
import parse
from papermodels.datatypes.joist import Joist
from PyNite import FEModel3D


@dataclass(frozen=True)
class Element:
    """
    A dataclass to generically represent a structural element on a PDF page
    whether it is a beam, column, joist, or otherwise.

    Note: each Element has no "knowledge" of its intersections and correspondents.
    They are generated in a separate process and are not part of the capability
    of the class. These attributes are simply to keep track of pre-discovered
    intersections and correspondents.

    'tag': str, represent a unique name for this element, as per the designer's preference
    'type': str, describing what "type" of element it is. This is not an enumeration
        and can take any designer-defined value. It is for user-level categorization.
    'page': int, describing the page index of the PDF that this element is found on
    'geometry': the shapely geometry that represents this element in 2D space
    'intersections': list[tuple[str, Point]], where each 2-tuple represents the tag
        of the other Element that this element intersects with on plan and the point
        of intersection. Intersections exist only between elements that are on the same PDF page
        and thus are used to describe connections that occur on the horizontal plane.
    'correspondents': list[str], where each item in the list represents the tag of another
        Element that is approximately in the exact same position as this element but on the
        adjacent page. Correspondents exist only between elements that are on adjacent PDF
        pages and thus are used to describe connections that occur on the vertical plane
        such as the the top and bottom of a wall or column (the bottom may be on page 0 and
        the top may be on page 1)
    'page_label': A designer-defined str label for the page (e.g. "Ground Floor" or "L03", etc.)
    """

    tag: str
    type: str
    page: int
    geometry: Geometry
    intersections: list[tuple[str, Point]]
    correspondents: list[tuple[str, Geometry]]
    page_label: Optional[str] = None

    def get_intersection(self, tag: str) -> Point | None:
        """
        Returns the corresponding Point object for 'tag' if 'tag' is in the tuple of
        self.intersections. Returns None if not.
        """
        lookup = dict(self.intersections)
        return lookup.get(tag, None)


# Examples
## This example shows a beam that is connected to a joist and a column on the same page
## and with that column having a correspondent on the page below
E00 = Element(
    tag="FB1.1",
    type="Flush Beam",
    page=1,
    geometry=LineString([[101.5, 52.0], [101.5, 85.3]]),
    intersections=[("J1.1", Point([101.5, 65.2]))],
    correspondents=[],
)

E01 = Element(
    tag="C1.1",
    type="Column",
    page=1,
    geometry=Polygon([[100.0, 100.0], [100.0, 103.0], [103.0, 103.0], [103.0, 100.0]]),
    intersections=[("FB1.1", Point([101.5, 53.5]))],
    correspondents=["C0.1"],
    page_label="L02",
)

E02 = Element(
    tag="C0.1",
    type="Column",
    page=0,
    geometry=Polygon([[100.0, 100.0], [100.0, 103.0], [103.0, 103.0], [103.0, 100.0]]),
    intersections=[],
    correspondents=["C1.1"],
    page_label="L01",
)


def element_to_beam_model(element: Element) -> FEModel3D:
    """
    Returns an FEModel3D for the data in 'element'
    """
    model = FEModel3D()
    model.add_material("default", 1, 1, 1, 1)
    element_length = element.geometry.length
    start_point = Point(element.geometry.coords[0])
    add_start_node = True
    add_end_node = True
    start_node = "N0"
    model.add_node(start_node, 0, 0, 0)
    for node_name, point in element.intersections:
        x_coord = start_point.distance(point)
        if math.isclose(x_coord, element_length, abs_tol=1e-3):
            last_node = node_name
            add_end_node = False
        if math.isclose(x_coord, 0, abs_tol=1e-3):
            model.Nodes.pop(start_node)
            start_node = node_name
        model.add_node(node_name, x_coord, 0, 0)
        model.def_support(node_name, support_DY=True)
    model.def_support(node_name, 1, 1, 1, 1, 1, 0)

    if add_end_node:
        last_node = "Nn"
        model.add_node(last_node, element_length, 0, 0)
    model.add_member(element.tag, "N0", last_node, "default", 1, 1, 1, 1)
    model.add_load_combo("Pass", {"pass": 1.0})
    return model


def element_to_joist_model(element: Element, w: float = 0.0) -> Joist:
    """
    Returns a Joist object based on the data in 'element'
    """
    try:
        r1, r2 = element.intersections
    except ValueError:
        raise ValueError(
            f"Joists currently need to have two supports. {element.tag=} | {element.intersections=}"
        )
    elem_geometry = element.geometry
    i_end = Point(elem_geometry.coords[0])
    j_end = Point(elem_geometry.coords[1])
    length = elem_geometry.length

    # R1 should be closest to I-end
    r1_geom = r1[1]
    r2_geom = r2[1]
    if i_end.distance(r1_geom) > i_end.distance(r2_geom):
        r1, r2 = r2, r1

    span = r1_geom.distance(r2_geom)
    a_cantilever = round(abs(i_end.distance(r1_geom)), 6)
    b_cantilever = round(abs(r2_geom.distance(j_end)), 6)
    return Joist(span, a_cantilever, b_cantilever)


def get_tag_type(this_element_tag: str) -> str:
    """
    Returns the prefix portion of 'this_element_tag'. The prefix portion is the
    alphabetical portion of the tag at the beginning.
    """
    format = "{type_tag}{page_tag:d}.{enum_tag:d}"
    result = parse.parse(format, this_element_tag)
    return result.named["type_tag"]


def get_elements_by_page(elements: list[Element]) -> dict[int, list[Element]]:
    """
    Returns 'elements' sorted by page
    """
    elements_by_page = {}
    for element in elements:
        elements_on_page = elements_by_page.get(element.page, [])
        elements_on_page.append(element)
        elements_by_page[element.page] = elements_on_page
    return elements_by_page


def get_normalized_coordinate(element: Element, intersection_point: Point) -> float:
    """
    Returns a normalized x-coordinate for the given 'coord' as it is located on the geometry
    of 'element'. Returns None if the 'coord' provided is not colinear with
    the element geometry.
    """
    geom = element.geometry
    i_coord = Point(geom.coords[0])
    distance = i_coord.distance(intersection_point)
    return distance
