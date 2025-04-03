from dataclasses import dataclass
from typing import Optional, Union, NamedTuple
import numpy as np
import numpy.typing as npt
from shapely import Point, LineString, Polygon
import shapely
from .annotation import Annotation
from ..paper.annotations import (
    parse_annotations,
    tag_parsed_annotations,
)
from ..geometry import geom_ops
from ..loads import load_distribution as ld
import parse
import math
import tomli_w
import json

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
    element_type: One of {"collector", "transfer"} or None. Assigned within the
        GeometryGraph.
    subelements: list[Element] or None. Assigned within the GeometryGraph.
    trib_area: Optional[Polygon] or None. # Not sure if adding this here is the right
        thing to do. Currently in use for the creation of collector subelements and for storing
        their trib areas.
    """

    geometry: Geometry
    tag: Optional[str | int] = None
    intersections_above: Optional[list[tuple]] = None
    intersections_below: Optional[list[tuple]] = None
    correspondents_above: Optional[list[dict]] = None
    correspondents_below: Optional[list[dict]] = None
    plane_id: Optional[str | int] = None
    element_type: Optional[str] = None
    subelements: list["Element"] = None
    trib_area: Optional[Polygon] = None

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
            above_tag: Intersection(*geom_ops.get_intersection(elem_geom, above_geom, above_tag))
            for above_tag, above_geom in intersections_above.items()

        } if intersections_above is not None else {}
        inters_below = {
            below_tag: Intersection(*geom_ops.get_intersection(elem_geom, below_geom, below_tag))
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
    def from_parsed_annotations(
        cls,
        parsed_annotations: dict[Annotation, dict], 
        correspond_with_like_only: bool = True,
    ) -> list["Element"]:
        """
        Returns a list of Element generated from the annotations in 'annots' according to the element
        types described in the 'legend'. If an annotation is not described in the legend then it will
        not be included in the result list of Elements.
        """
        tagged_annotations = tag_parsed_annotations(parsed_annotations)
        annotations_w_intersect = get_geometry_intersections(tagged_annotations)
        annotations_w_intersect_corrs = get_geometry_correspondents(annotations_w_intersect)

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

@dataclass
class LoadedElement(Element):
    loading_geoms: Optional[ld.LoadingGeometry] = None
    applied_loading_areas: Optional[list[tuple[Polygon, npt.ArrayLike]]] = None
    model: Optional[dict] = None

    """
    'loading_areas' - A list of tuples. Each tuple consists of a Polygon and a dict of
        attributes associated with that Polygon. If no attributes are desired,
        pass a tuple with an empty dict.
    'applied_loading_areas' - A dict of Polygon/attributes that intersect with the
        trib area of this Element. The Polygons in the dict represent the intersecting
        area of the loading area and the trib area. These are computed from the provided
        'loading_areas' during initialization. The designer is not expected to populate
        this parameter.
    'model' - A dictionary that describes the LoadedElement in terms of a structural
        element. Populated during initialization. The designer is not expected to populate
        this parameter.


    TODO: # HERE: Need to apply loading to sub-elements 
    """
    def __post_init__(self):
        """
        Populates self.applied_loading_areas
        """
        self.applied_loading_areas = self._get_intersecting_loads()
        self.model = self._build_model()
                
        
    def _get_intersecting_loads(self) -> list[tuple[Polygon, dict]]:
        loading_array = np.array([loading_area.geometry for loading_area in self.loading_geoms])
        applied_loading_areas = []
        if self.trib_area is not None:
            intersecting_loads = self.trib_area.intersection(loading_array)
            for idx, intersecting_load in enumerate(intersecting_loads):
                if intersecting_load.is_empty or math.isclose(intersecting_load.area, 0): 
                    continue
                applied_loading_areas.append(
                    (intersecting_load, self.loading_geoms[idx])
                )
        return applied_loading_areas 

        
    def dump_analysis_model(self) -> dict:
        """
        Returns the structured beam data dict to go to analysis model
        """
        return {}
    
    def _build_model(self) -> dict:
        """
        Returns the structured beam dict for serialization
        """
        orientation = "unknown"
        if self.geometry.geom_type == "LineString": 
            orientation = "horizontal"
        elif self.geometry.geom_type == "Polygon":
            orientation = "vertical"

        support_locations = self._get_support_locations()
        transfer_loads = []
        if self.element_type == "transfer":
            transfer_loads = self._get_transfer_loads()
        distributed_loads = self._get_distributed_loads()

        model = {
            "element_attributes":
                {
                    "tag": self.tag,
                    "length": self.geometry.length if self.geometry.geom_type == "LineString" else {},
                    "orientation": orientation,
                },
            "element_geometry":
                {
                    "geometry": self.geometry.wkt,
                    "supports": support_locations,
                },
            "loads": {
                "point_loads": transfer_loads,
                "distributed_loads": distributed_loads or []
            }
        }
        return model


    def _get_support_locations(self):
        """
        Calculates the support locations from the intersections below
        """
        if self.geometry.geom_type == "LineString":
            coords_a, coords_b = self.geometry.coords
            coords_a, coords_b = Point(coords_a), Point(coords_b)
            ordered_coords = geom_ops.order_nodes_positive(coords_a, coords_b)
            start_coord = ordered_coords[0]
            support_locations = geom_ops.get_local_intersection_ordinates(
                start_coord,
                [intersection[0] for intersection in self.intersections_below]
            )
            supports_acc = []
            for idx, support_location in enumerate(support_locations):
                fixity = "roller"
                if idx == 0:
                    fixity = "pin"
                supports_acc.append({"location": support_location, "fixity": fixity})
            return supports_acc
        else:
            return []
    
    def _get_transfer_loads(self):
        """
        Calculates the transfer load locations from the intersections above
        """
        transfer_loads = []
        if self.geometry.geom_type == "LineString":
            coords_a, coords_b = self.geometry.coords
            coords_a, coords_b = Point(coords_a), Point(coords_b)
            ordered_coords = geom_ops.order_nodes_positive(coords_a, coords_b)
            start_coord = ordered_coords[0]
            transfer_locations = geom_ops.get_local_intersection_ordinates(
                start_coord,
                [intersection[0] for intersection in self.intersections_above]
            )
            for idx, transfer_location in enumerate(transfer_locations):
                intersection_data = self.intersections_above[idx]
                source_member = intersection_data.other_tag
                reaction_idx = intersection_data.other_index
                if reaction_idx is None:
                    raise ValueError(
                        "The .other_index attribute within the .intersections_above list"
                        " is not calculated. Generate LoadedElement objects through the GeometryGraph"
                        " interface in order to populate this necessary index."
                    )
                transfer_loads.append(
                    {
                        "location": transfer_location,
                        "magnitude": 0,
                        "transfer_source": f"{source_member}",
                        "transfer_reaction_idx": reaction_idx,
                        "direction": "gravity"
                    }
                )
        elif self.geometry.geom_type == "Polygon":
            for intersection in self.intersections_above:
                transfer_loads.append(
                    {
                        "location": [],
                        "magnitude": 0,
                        "transfer_source": intersection.other_tag,
                        "transfer_reaction_idx": intersection.other_index,
                        "direction": "gravity"
                    }
                )
        return transfer_loads
    
    def _get_distributed_loads(self):
        """
        Computes the resulting distributed loads from the applied
        loading areas
        """
        if self.geometry.geom_type == "LineString":
            raw_dist_loads = ld.get_distributed_loads_from_projected_polygons(
                self.geometry,
                self.applied_loading_areas
            )
            polygon_areas = geom_ops.calculate_trapezoid_area_sums(raw_dist_loads)
            distributed_loads = []
            for idx, dist_load_collection in enumerate(raw_dist_loads):
                total_polygon_area = polygon_areas[idx]
                for dist_load_element in dist_load_collection:
                    start_xy, end_xy = dist_load_element
                    start_x, start_y = start_xy
                    if math.isclose(start_x, 0, abs_tol=1e-6):
                        start_x = 0
                    end_x, end_y = end_xy
                    area_dist_load = geom_ops.trapezoid_area(
                        h=(end_x - start_x),
                        b1=start_y,
                        b2=end_y
                    )
                    trapezoid_ratio = area_dist_load / total_polygon_area
                    intersected_poly, applied_loading = self.applied_loading_areas[idx]
                    dist_load = {
                        "occupancy": applied_loading.occupancy,
                        "load_components": applied_loading.load_components or [],
                        "applied_area": intersected_poly.area * trapezoid_ratio,
                        "start_loc": start_x,
                        "start_magnitude": start_y,
                        "end_loc": end_x,
                        "end_magnitude":  end_y,
                    }
                    distributed_loads.append(dist_load)
            return distributed_loads


    def dump_toml(self, fp):
        """
        Dumps the .model attribute to a TOML file
        """
        tomli_w.dump(self.model, fp)
        return fp
        
    def dump_json(self, fp):
        """
        Dumps the .model attribute to a TOML file
        """
        json.dump(self.model, fp, indent=2)
        return fp
        
    
    @classmethod
    def from_element_with_loads(cls, elem: Element, loading_geoms: dict[Polygon, Union[str | npt.ArrayLike]], trib_area: Optional[Polygon] = None):
        """
        Returns a LoadedElement
        """
        return cls(
            elem.geometry,
            elem.tag,
            elem.intersections_above,
            elem.intersections_below,
            elem.correspondents_above,
            elem.correspondents_below,
            elem.plane_id,
            element_type=elem.element_type,
            subelements=elem.subelements,
            trib_area=elem.trib_area or trib_area,
            loading_geoms=loading_geoms,
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


def get_geometry_intersections(
    tagged_annotations: dict[Annotation, dict],
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
        intersections_above = []
        intersections_below = []
        for j_annot in annots:
            j_attrs = intersected_annotations[j_annot]
            j_rank = j_attrs["rank"]
            j_page = j_annot.page
            i_geom = i_attrs["geometry"]
            j_geom = j_attrs["geometry"]
            if i_page != j_page: 
                continue
            if i_rank < j_rank:
                intersection = geom_ops.get_intersection(i_geom, j_geom, j_attrs['tag'])
                if intersection is None: continue
                intersections_below.append(Intersection(*intersection))
            elif i_rank > j_rank:
                intersection = geom_ops.get_intersection(j_geom, i_geom, j_attrs['tag'])
                if intersection is None: continue
                intersections_above.append(Intersection(*intersection))
        i_attrs["intersections_above"] = intersections_above
        i_attrs["intersections_below"] = intersections_below
    return intersected_annotations


def get_geometry_correspondents(
    tagged_annotations: dict[Annotation, dict],
) -> dict[Annotation, dict]:
    """
    Returns a copy of 'tagged_annotations' with a 'correspondents' field added to that
    attributes dictionary of each Annotation key.
    """
    annots_by_page = annotations_by_page(tagged_annotations)
    descending_pages = sorted(annots_by_page.keys(), reverse=True)
    last_page = descending_pages[-1]
    corresponding_annotations = tagged_annotations.copy()
    prev_page = None
    for page in descending_pages:
        if page != last_page:
            next_page = page - 1
            annots_here = annots_by_page[page]
            annots_below = annots_by_page[next_page]
            correspondents_above = {j_attrs['tag']: [] for j_attrs in annots_below.values()}
            correspondents_below = []

            for i_annot, i_attrs in annots_here.items():
                i_page = i_annot.page
                correspondents_below = []
                for j_annot, j_attrs in annots_below.items():
                    j_attrs = annots_below[j_annot]
                    j_page = j_annot.page
                    i_geom = i_attrs["geometry"]
                    j_geom = j_attrs["geometry"]
                    i_tag = i_attrs['tag']
                    j_tag = j_attrs['tag']
                    correspondence_ratio = geom_ops.check_corresponds(i_geom, j_geom)
                    if correspondence_ratio:
                        correspondents_below.append(Correspondent(correspondence_ratio, j_geom, j_tag))
                        correspondents_above[j_attrs['tag']].append(Correspondent(correspondence_ratio, i_geom, i_attrs['tag']))
                corresponding_annotations[i_annot]["correspondents_above"] = correspondents_above.get(i_attrs['tag'], [])
                corresponding_annotations[i_annot]["correspondents_below"] = correspondents_below
                
        else:
            annots_here = annots_by_page[page]
            if len(descending_pages) == 1: 
                correspondents_above = {} # There are no correspondents above or below on a single page
            for i_annot, i_attrs in annots_here.items():
                corresponding_annotations[i_annot]["correspondents_below"] = []
                corresponding_annotations[i_annot]["correspondents_above"] = correspondents_above.get(i_attrs['tag'], [])
        if prev_page is None:
            prev_page = page
    return corresponding_annotations


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
