from dataclasses import dataclass
from typing import Optional, Union, NamedTuple
from shapely import Point, LineString, Polygon
from .annotation import Annotation
from ..paper.annotations import (
    parse_annotations,
    tag_parsed_annotations,
    get_geometry_intersections,
    get_geometry_correspondents
)
from ..geometry import geom_ops
import parse

Geometry = Union[LineString, Polygon]


class Intersection(NamedTuple):
    """
    A class to represent an intersection of geometries
    """
    intersecting_region: Point
    other_geometry: Union[LineString, Polygon]
    other_tag: str
    other_index: Optional[int] = None


class Correspondent(NamedTuple):
    """
    A class to represent the correspondence of two Polygons
    """
    overlap_ratio: float
    other_geometry: Polygon
    other_tag: str


@dataclass
class Element:
    """
    A class to generically represent a connected 2D geometry within a 3D "geometry graph".

    The 2D geometry can exist independent of an interconnected 3D graph by not having
    any 'intersections' or 'correspondents'. The existence of 'intersections' and/or
    'correspondents' indicates that the 2D geometry is part of a graph.

    'intersections' describe interconnectivity on the same 2D plane.
    'correspondents' describe interconnectivity between adjacent 2D planes,
        whether above or below.

    geometry: Union[LineString, Polygon], The geometry for the Element
    tag: str | int,  An optional unique name or integer ID for the Element.
    intersections_above/below: a dict whose keys are the 'tag' of an intersecting
        Element and the values are the Point within self.geometry
        where the intersection occurs. _above represents geometries that occur
        "above" the element in the directed geoemtry graph and while _below
        represents those "below". 
    correspondents_above/below: a dict whose keys are the 'tag' of a corresponding
        geometry on an adjacent 2D plane and the values are the corresponding
        Geometry on the adjacent plane. _above represents geometries that occur
        "above" the element in the directed geoemtry graph and while _below
        represents those "below".
    plane_id: Optional[str | int] = None, An optional unique identifier for the 2D
        plane that this Element resides on
    """

    geometry: Geometry
    tag: Optional[str | int] = None
    intersections_above: Optional[list[tuple]] = None
    intersections_below: Optional[list[tuple]] = None
    correspondents_above: Optional[list[dict]] = None
    correspondents_below: Optional[list[dict]] = None
    plane_id: Optional[str | int] = None

    def __post_init__(self):
        if self.geometry.geom_type == "LineString" and len(self.geometry.coords) != 2:
            raise ValueError(
                "Element objects of LineStrings must have LineStrings containing only one segment.\n"
                f"Element with {self.tag} has {len(self.geometry.coords -1)} segments."
            )
 

    @classmethod
    def from_geometries(
        cls,
        elem_geom: Geometry,
        elem_tag: str | int,
        intersections_above: Optional[dict[str | int, Geometry]] = None,
        intersections_below: Optional[dict[str | int, Geometry]] = None,
        correspondents_above: Optional[dict[str | int, Geometry]] = None,
        correspondents_below: Optional[dict[str | int, Geometry]] = None,
        plane_id: Optional[str | int] = None
    ):
        """
        Generates an Element from provided geometries
        """
        inters_above = {
            above_tag: geom_ops.get_intersection(elem_geom, above_geom, above_tag)
            for above_tag, above_geom in intersections_above.items()

        } if intersections_above is not None else {}
        inters_below = {
            below_tag: geom_ops.get_intersection(elem_geom, below_geom, below_tag)
            for below_tag, below_geom in intersections_below.items()
        } if intersections_below is not None else {}

        return cls(
            tag=elem_tag,
            geometry=elem_geom,
            intersections_above=inters_above,
            intersections_below=inters_below,
            correspondents_above=correspondents_above or {},
            correspondents_below=correspondents_below or {},
        )

    @classmethod
    def from_annotations(
        cls,
        annots: list[Annotation], 
        legend: list[Annotation],
        correspond_with_like_only: bool = True,
    ) -> list["Element"]:
        """
        Returns a list of Element generated from the annotations in 'annots' according to the element
        types described in the 'legend'. If an annotation is not described in the legend then it will
        not be included in the result list of Elements.
        """
        sorted_by_page_annotations = sorted(annots, key=lambda x: x.page, reverse=True)
        parsed_annotations = parse_annotations(sorted_by_page_annotations, legend)
        tagged_annotations = tag_parsed_annotations(parsed_annotations)
        annotations_w_intersect = get_geometry_intersections(tagged_annotations)
        annotations_w_intersect_corrs = get_geometry_correspondents(annotations_w_intersect)
        # print(corresponding_annotations)

        elements = []
        for annot_attrs in annotations_w_intersect_corrs.values():
            if correspond_with_like_only:
                corrs_a = [cor for cor in annot_attrs['correspondents_above'] if annot_attrs['tag'][0] == cor.other_tag[0]]
                corrs_b = [cor for cor in annot_attrs['correspondents_below'] if annot_attrs['tag'][0] == cor.other_tag[0]]
            element = cls(
                tag=annot_attrs["tag"],
                geometry=annot_attrs["geometry"],
                intersections_above=annot_attrs["intersections_above"],
                intersections_below=annot_attrs["intersections_below"],
                correspondents_above=corrs_a,
                correspondents_below=corrs_b,
                plane_id=annot_attrs.get("page_label", None),
            )
            elements.append(element)
        return elements
    
    @property
    def supporting_geometries(self):
        acc = []
        for intersection_tuple in self.intersections.values():
            acc.append(intersection_tuple[1])
        return acc


# Examples
E00 = Element(
    tag="FB1.1",
    # type="Flush Beam",
    # page=1,
    geometry=LineString([[101.5, 52.0], [101.5, 85.3]]),
    intersections_above=[
        (Point([101.5, 65.2]), LineString([(84.2, 65.2), (120.0, 65.2)]), "J1.1")
    ],
    # correspondents=[],
)


# ## This example shows a beam that is connected to a joist and a column on the same page
# ## and with that column having a correspondent on the page below
# E01 = Element(
#     tag="C1.1",
#     # type="Column",
#     # page=1,
#     geometry=Polygon([[100.0, 100.0], [100.0, 103.0], [103.0, 103.0], [103.0, 100.0]]),
#     intersections_below=[("FB1.1", Point([101.5, 53.5]))],
#     correspondents_below=["C0.1"],
#     # page_label="L02",
# )

# E02 = Element(
#     tag="C0.1",
#     # type="Column",
#     # page=0,
#     geometry=Polygon([[100.0, 100.0], [100.0, 103.0], [103.0, 103.0], [103.0, 100.0]]),
#     intersections=[],
#     correspondents=["C1.1"],
#     # page_label="L01",
# )


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
