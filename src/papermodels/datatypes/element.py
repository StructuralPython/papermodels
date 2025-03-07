from dataclasses import dataclass
from typing import Optional
from shapely import Point, Geometry, LineString, Polygon
from .annotation import Annotation
from ..geometry import geom_ops
import parse


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

    geometry: Geometry
    tag: Optional[str] = None
    isjoist: bool = False
    annotation: Optional[Annotation] = None
    intersections: Optional[dict[str | int, Point]] = None
    correspondents: Optional[dict[str | int, Geometry]] = None

    @classmethod
    def from_geometries(
        cls,
        subject_geom: Geometry,
        supporting_geoms: list[Geometry] | dict[str, Geometry],
        subject_tag: Optional[str] = None,
        corresponding_geoms: Optional[list[Geometry] | dict[str, Geometry]] = None,
        isjoist: bool = False
    ):
        """
        Returns an Element instance created from the main 'subject' geometry
        and 'supporting' geometries which represent structural supporting
        members.

        subject_geom: The structural element that is the subject of study.
            Any singular shapely Geometry type.
        supporting: A list of any singular shapely Geometry types. If a dict
            is passed, then the keys will be used as identifying tags for the
            provided supporting geometries. The supporting geometries must
            intersect with the subject geometry _on the same plane_.
        subject_tag: a unique str label to identify the subject element.
        corresponding_geoms: A list of any singular shapely Geometry types.
            If a dict is passed, then the keys will be used as identifying
            tags for the provided corresponding geometries.
            The
        """
        if isinstance(supporting_geoms, list):
            supporting_geoms = {
                idx: supporting_geom
                for idx, supporting_geom in enumerate(supporting_geoms)
            }

        intersections = {}
        for support_tag, support_geom in supporting_geoms.items():
            intersection = geom_ops.get_intersection(
                subject_geom, support_geom
            )
            if intersection is not None:
                intersections.update({support_tag: intersection})

        return cls(
            tag=subject_tag,
            geometry=subject_geom,
            intersections=intersections,
            correspondents=corresponding_geoms,
        )
    
    @property
    def supporting_geometries(self):
        acc = []
        for intersection_tuple in self.intersections.values():
            acc.append(intersection_tuple[1])
        return acc


# Examples
## This example shows a beam that is connected to a joist and a column on the same page
## and with that column having a correspondent on the page below
E00 = Element(
    tag="FB1.1",
    # type="Flush Beam",
    # page=1,
    geometry=LineString([[101.5, 52.0], [101.5, 85.3]]),
    intersections=[("J1.1", Point([101.5, 65.2]))],
    # correspondents=[],
)

E01 = Element(
    tag="C1.1",
    # type="Column",
    # page=1,
    geometry=Polygon([[100.0, 100.0], [100.0, 103.0], [103.0, 103.0], [103.0, 100.0]]),
    intersections=[("FB1.1", Point([101.5, 53.5]))],
    correspondents=["C0.1"],
    # page_label="L02",
)

E02 = Element(
    tag="C0.1",
    # type="Column",
    # page=0,
    geometry=Polygon([[100.0, 100.0], [100.0, 103.0], [103.0, 103.0], [103.0, 100.0]]),
    intersections=[],
    correspondents=["C1.1"],
    # page_label="L01",
)


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
    Returns a normalized x-coordinate for the given 'intersection_point' as it is located on the geometry
    of 'element'.
    """
    geom = element.geometry
    i_coord = Point(geom.coords[0])
    distance = i_coord.distance(intersection_point)
    return distance


def get_structured_model_data(element: Element) -> dict:
    """ """
    pass
