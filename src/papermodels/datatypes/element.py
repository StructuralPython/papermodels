from dataclasses import dataclass
from typing import Optional, Union, NamedTuple
import numpy as np
import numpy.typing as npt
from shapely import Point, LineString, Polygon
from shapely import wkt
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
    other_reaction_type: str = "point"
    other_extents: Optional[tuple] = None


class Correspondent(NamedTuple):
    """
    A class to represent the correspondence of two Polygons
    """
    overlap_ratio: float
    other_geometry: Polygon
    other_tag: str
    other_rank: int
    other_reaction_type: str = "point"
    other_extents: Optional[tuple] = None


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
    rank: Optional[int] = None
    intersections_above: Optional[list[tuple]] = None
    intersections_below: Optional[list[tuple]] = None
    correspondents_above: Optional[list[dict]] = None
    correspondents_below: Optional[list[dict]] = None
    plane_id: Optional[str | int] = None
    element_type: Optional[str] = None
    subelements: list["Element"] = None
    trib_area: Optional[Polygon] = None
    reaction_type: str = "point"

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
        inters_above = []
        inters_below = []
        if intersections_above is not None:
            inters_above = [
                Intersection(*geom_ops.get_intersection(elem_geom, above_geom, above_tag))
                for above_tag, above_geom in intersections_above.items()]
        
        if intersections_below is not None:
            inters_below = [
                Intersection(*geom_ops.get_intersection(elem_geom, below_geom, below_tag))
                for below_tag, below_geom in intersections_below.items() 
                if intersections_below is not None 
            ]

        return cls(
            tag=elem_tag,
            geometry=elem_geom,
            intersections_above=inters_above,
            intersections_below=inters_below,
            correspondents_above=correspondents_above or [],
            correspondents_below=correspondents_below or [],
        )

    @classmethod
    def from_parsed_annotations(
        cls,
        parsed_annotations: dict[Annotation, dict],
        trib_annotations: Optional[dict[Annotation, dict]] = None,
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
        trib_area_geoms = np.array([annot_attrs['geometry'] for annot_attrs in trib_annotations.values()])
        matching_trib_poly = None
        elements = []
        for annot_attrs in annotations_w_intersect_corrs.values():
            element_family = annot_attrs['tag'][0]
            geometry = annot_attrs['geometry']
            if geometry.geom_type == "LineString" and trib_annotations:
                intersection_mask = geometry.intersects(trib_area_geoms)
                intersection_lines = trib_area_geoms[intersection_mask]
                get_intersection_lengths = np.vectorize(lambda x: x.length)
                matching_trib_index = np.argmax(get_intersection_lengths) # The matching trib is the one that is mostly in the trib poly
                matching_trib_poly = trib_area_geoms[matching_trib_index]
            if correspond_with_like_only:
                corrs_a = prioritize_correspondents(annot_attrs['correspondents_above'], element_family)
                corrs_b = prioritize_correspondents(annot_attrs['correspondents_below'], element_family)

            element = cls(
                tag=annot_attrs["tag"],
                geometry=annot_attrs["geometry"],
                rank=annot_attrs['rank'],
                intersections_above=annot_attrs["intersections_above"],
                intersections_below=annot_attrs["intersections_below"],
                correspondents_above=corrs_a,
                correspondents_below=corrs_b,
                plane_id=annot_attrs.get("page_label", None),
                reaction_type=annot_attrs.get('reaction_type', 'point'),
                trib_area=matching_trib_poly
            )
            elements.append(element)
        return elements
    
    @property
    def supporting_geometries(self):
        acc = []
        for intersection_tuple in self.intersections.values():
            acc.append(intersection_tuple[1])
        return acc

def prioritize_correspondents(correspondents: list[Correspondent], family: str) -> list[Correspondent]:
    """
    Filters 'correspondents' such that:
        - Only one correspondent exists in the list
        - That correspondent is of the same "family"
        - The rank is either 0 or the same as 'rank'

    The purpose of this filtering is to ensure that the top of a vertical element only corresponds with
    one element below it, either as a continuation of the same type of element (e.g. a column continuing
    down the structure) or as the point load that results from transferring out the bottom of the
    element to some other element (e.g. the bottom of a column transferring out to a beam).
    """
    filtered_transfers = []
    filtered_same_family = []
    for correspondent in correspondents:
        corr_family = correspondent.other_tag[0]
        corr_rank = correspondent.other_rank
        if corr_rank == 0: # If the element transfers out
            filtered_transfers.append(correspondent)
        elif corr_family == family:
            filtered_same_family.append(correspondent)
    if filtered_transfers:
        largest_overlap = sorted(filtered_transfers, key=lambda x: x.overlap_ratio, reverse=True)[0]
        return [largest_overlap]
    if filtered_same_family:
        largest_overlap = sorted(filtered_same_family, key=lambda x: x.overlap_ratio, reverse=True)[0]
        return [largest_overlap]
    else:
        return []
    
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


    """
    def __post_init__(self):
        """
        Populates self.applied_loading_areas
        """
        self.applied_loading_areas = self._get_intersecting_loads()
                
        
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
    
    
    def model(self, precision=3) -> dict:
        """
        Returns the structured beam dict for serialization
        """
        orientation = "unknown"
        if self.geometry.geom_type == "LineString": 
            orientation = "horizontal"
        elif self.geometry.geom_type == "Polygon":
            orientation = "vertical"
        length = round(self.get_length(), precision) if self.get_length() else {}
        support_locations = self._get_support_locations(precision)
        transfer_loads = {}
        if self.element_type == "transfer":
            transfer_loads = self._get_transfer_loads(precision)
        distributed_loads = self._get_distributed_loads(precision)


        elem_model = {
            "element_attributes":
                {
                    "tag": self.tag,
                    "length": length,
                    "orientation": orientation,
                    "vert_correspondent_below": [corr.other_tag for corr in self.correspondents_below],
                    "vert_correspondent_above": [corr.other_tag for corr in self.correspondents_above],
                    "horz_intersects_above": [inter.other_tag for inter in self.intersections_above],
                    "horz_intersects_below": [inter.other_tag for inter in self.intersections_below]
                },
            "element_geometry":
                {
                    "geometry": wkt.dumps(self.geometry, trim=True, rounding_precision=precision),
                    "supports": support_locations,
                },
            "loads": {
                "point_loads": transfer_loads.get('point', []),
                "distributed_loads": transfer_loads.get('dist', []) + distributed_loads
            }
        }
        return elem_model


    def get_length(self):
        """
        Calculates the length fo the element, if applicable
        """
        if self.geometry.geom_type == "LineString":
            return self.geometry.length
        elif self.geometry.geom_type == "Polygon" and self.reaction_type == "linear":
            return geom_ops.get_rectangle_centerline(self.geometry).length
        else:
            return {}


    def _get_support_locations(self, precision: int):
        """
        Calculates the support locations from the intersections below
        """
        if self.geometry.geom_type == "LineString":
            coords_a, coords_b = self.geometry.coords
            coords_a, coords_b = Point(coords_a), Point(coords_b)
            ordered_coords = geom_ops.order_nodes_positive([coords_a, coords_b])
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
                supports_acc.append({"location": round(support_location, precision), "fixity": fixity})
            return supports_acc
        else:
            return []
    
    def _get_transfer_loads(self, precision: int):
        """
        Calculates the transfer load locations from the intersections above
        """
        transfer_loads = {"point": [], "dist": []}
        if self.geometry.geom_type == "LineString":
            coords_a, coords_b = self.geometry.coords
            coords_a, coords_b = Point(coords_a), Point(coords_b)
            ordered_coords = geom_ops.order_nodes_positive([coords_a, coords_b])
            start_coord = ordered_coords[0]
            transfer_locations = geom_ops.get_local_intersection_ordinates(
                start_coord,
                [intersection[0] for intersection in self.intersections_above]
            )
            for idx, transfer_location in enumerate(transfer_locations):
                intersection_above: Intersection = self.intersections_above[idx]
                transfer_type = intersection_above.other_reaction_type
                source_member = intersection_above.other_tag
                reaction_idx = intersection_above.other_index
                if reaction_idx is None:
                    raise ValueError(
                        "The .other_index attribute within the .intersections_above list"
                        " is not calculated. Generate LoadedElement objects through the GeometryGraph"
                        " interface in order to populate this necessary index."
                    )
                if transfer_type == "point":
                    point_load = self.create_point_load(
                        transfer_location=round(transfer_location, precision),
                        magnitude=0.0,
                        transfer_source=f"{source_member}",
                        transfer_reaction_index=reaction_idx,
                        direction="gravity",
                    )
                    transfer_loads['point'].append(point_load)
                elif transfer_type == "linear":
                    dist_load = self.create_distributed_load(
                        start_location=round(intersection_above.other_extents[0], precision),
                        start_magnitude=1.0,
                        end_location=round(intersection_above.other_extents[1], precision),
                        end_magnitude=1.0,
                        transfer_source=f"{source_member}",
                        transfer_reaction_index=intersection_above.other_index,
                        occupancy="",
                        load_components={},
                        applied_area=0.0,
                        direction="gravity"
                    )
                    transfer_loads['dist'].append(dist_load)

        elif self.geometry.geom_type == "Polygon":
            for correspondent in self.correspondents_above:
                if correspondent.other_reaction_type == "point":
                    point_load = self.create_point_load(
                        transfer_location=[],
                        magnitude=0.0,
                        transfer_source=correspondent.other_tag,
                        transfer_reaction_index=0,
                        direction="gravity"
                    )
                    transfer_loads['point'].append(point_load)
                elif correspondent.other_reaction_type == "linear":
                    dist_load = self.create_distributed_load(
                        start_location=0.0,
                        start_magnitude=1.0,
                        end_location=round(geom_ops.get_rectangle_centerline(self.geometry).length, precision),
                        end_magnitude=1.0,
                        transfer_source=correspondent.other_tag,
                        transfer_reaction_index=0,
                        occupancy="",
                        load_components={},
                        applied_area=0.0,
                        direction="gravity",
                    )
                    transfer_loads['dist'].append(dist_load)
            for intersection in self.intersections_above:
                if intersection.other_reaction_type == "point":
                    transfer_loads['point'].append(
                        {
                            "location": [],
                            "magnitude": 0,
                            "transfer_source": intersection.other_tag,
                            "transfer_reaction_idx": intersection.other_index,
                            "direction": "gravity"
                        }
                    )
                elif intersection.other_reaction_type == "linear":
                    source_member = intersection.other_tag
                    start_x, end_x = intersection.other_extents
                    transfer_loads['dist'].append(
                        self.create_distributed_load(
                            start_location=round(start_x, precision),
                            start_magnitude=1.0,
                            end_location=round(end_x, precision),
                            end_magnitude=1.0,
                            transfer_source=f"{source_member}",
                            transfer_reaction_index=intersection.other_index,
                            occupancy="",
                            load_components={},
                            applied_area=0.0,
                            direction="gravity"
                        )
                    )
                    
        return transfer_loads
    

    @staticmethod
    def create_point_load(
            transfer_location: str,
            magnitude: float,
            transfer_source: str,
            transfer_reaction_index: int,
            direction: str = "gravity"
    ):

        return {
            "location": transfer_location,
            "magnitude": magnitude,
            "transfer_source": transfer_source,
            "transfer_reaction_index": transfer_reaction_index,
            "direction": direction
        }
    
    @staticmethod
    def create_distributed_load(
        start_location: float,
        start_magnitude: float,
        end_location: float,
        end_magnitude: float,
        transfer_source: str,
        transfer_reaction_index: int,
        occupancy: str,
        load_components: dict,
        applied_area: float,
        direction: str,
    ):

        return {
            "transfer_source": transfer_source,
            "transfer_reaction_index": transfer_reaction_index,
            "occupancy": occupancy,
            "load_components": load_components,
            "applied_area": applied_area,
            "start_loc": start_location,
            "start_magnitude": start_magnitude,
            "end_loc": end_location,
            "end_magnitude": end_magnitude,
            "direction": direction,
        }


    def _get_distributed_loads(self, precision: int) -> list[dict]:
        """
        Computes the resulting distributed loads from the applied
        loading areas
        """            
        distributed_loads = []
        if self.geometry.geom_type == "LineString":
            raw_dist_loads = ld.get_distributed_loads_from_projected_polygons(
                self.geometry,
                self.applied_loading_areas
            )
            polygon_areas = geom_ops.calculate_trapezoid_area_sums(raw_dist_loads)
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
                        "transfer_source": "",
                        "transfer_reaction_index": "", 
                        "occupancy": applied_loading.occupancy,
                        "load_components": applied_loading.load_components or [],
                        "applied_area": round(intersected_poly.area * trapezoid_ratio, precision),
                        "start_loc": round(start_x, precision),
                        "start_magnitude": round(start_y, precision),
                        "end_loc": round(end_x, precision),
                        "end_magnitude":  round(end_y, precision),
                    }
                    distributed_loads.append(dist_load)
            return distributed_loads
        else:
            return []


    def dump_toml(self, fp, precision=3):
        """
        Dumps the .model attribute to a TOML file
        """
        tomli_w.dump(self.model(precision), fp)
        return fp
    
        
    def dump_json(self, fp, precision=3):
        """
        Dumps the .model attribute to a TOML file
        """
        json.dump(self.model(precision), fp, indent=2)
        return fp
        
    
    @classmethod
    def from_element_with_loads(
        cls, 
        elem: Element, 
        loading_geoms: dict[Polygon, Union[str | npt.ArrayLike]], 
        trib_area: Optional[Polygon] = None,
        predecessors: Optional[list[str]] = None,
        successors: Optional[list[str]] = None,
        ):
        """
        Returns a LoadedElement. Validates the intersections and correspondents against the
        supplied 'predecessors' and 'successors' from the graph. Intersections and correspondents
        that do not exist in the 'predecessors' or 'successors' are excluded.
        """
        cleaned_intersections_above = []
        for intersection in elem.intersections_above:
            if intersection.other_tag in predecessors:
                cleaned_intersections_above.append(intersection)
        cleaned_intersections_below = []
        for intersection in elem.intersections_below:
            if intersection.other_tag in successors:
                cleaned_intersections_below.append(intersection)
        cleaned_correspondents_above = []
        for correspondent in elem.correspondents_above:
            if correspondent.other_tag in predecessors:
                cleaned_correspondents_above.append(correspondent)
        cleaned_correspondents_below = []
        for correspondent in elem.correspondents_below:
            if correspondent.other_tag in successors:
                cleaned_correspondents_below.append(correspondent)


        return cls(
            elem.geometry,
            elem.tag,
            elem.rank,
            cleaned_intersections_above,
            cleaned_intersections_below,
            cleaned_correspondents_above,
            cleaned_correspondents_below,
            elem.plane_id,
            element_type=elem.element_type,
            subelements=elem.subelements,
            trib_area=elem.trib_area or trib_area,
            loading_geoms=loading_geoms,
            reaction_type=elem.reaction_type
        )


def get_collector_extents(
    collector_prototype: Element,
) -> dict[str, tuple]:
    """
    Returns a dict keyed by .tag attributes in collector_prototype.intersections_below
    and with values representing the (start_x, end_x) locations where the collector
    prototype would spread over the other_geometry. The (start_x, end_x) locations 
    refer to ordinates on other_geometry, not on the collector_prototype.
    """
    support_tags_by_geom = {
        geom_ops.clean_polygon_supports([ib.other_geometry], collector_prototype.geometry)[0]: ib.other_tag 
        for ib in collector_prototype.intersections_below
    }
    poly_support_geoms = list(support_tags_by_geom.keys())
    support_geoms = geom_ops.clean_polygon_supports(poly_support_geoms, collector_prototype.geometry)
    
    cleaned_supports_map = {}
    for idx, poly_support_geom in enumerate(poly_support_geoms):
        clean_support_geom = support_geoms[idx]
        cleaned_supports_map.update({clean_support_geom: poly_support_geom})
    ordered_support_geoms = geom_ops.sort_supports(collector_prototype.geometry, support_geoms)
    try:
        extents = geom_ops.get_joist_extents(collector_prototype.geometry, ordered_support_geoms)
    except AssertionError as e:
        raise AssertionError(f"No intersection within joist extents: {collector_prototype.tag=}")

    
    tagged_extents = {}
    for idx, extent in enumerate(extents):
        support_geom = ordered_support_geoms[idx]
        support_start, _ = geom_ops.get_start_end_nodes(support_geom)
        support_tag = support_tags_by_geom[support_geom]
        extent_start = extent[0].distance(support_start)
        extent_end = extent[1].distance(support_start)
        tagged_extents.update({support_tag: (extent_start, extent_end)})
    return tagged_extents


def get_transfer_extents(element: Element) -> tuple[str, dict]:
    """
    Returns a tuple of str, extents_dict

    e.g.
    ("FB0.1", {"A": 12, "B": 23.4})

    For the polygon element that has a linear reaction type.
    This element could have more than one intersection below
    if it overlaps multiple beams, for example. It can also 
    have correspondents that need to have a portion of teh
    linear load applied to it.
    """
    intersection_extents = {}
    for intersection_below in element.intersections_below:
        intersection_below: Intersection
        tag = intersection_below.other_tag
        other_geom = intersection_below.other_geometry
        if isinstance(other_geom, LineString):
            start_coord, _ = geom_ops.get_start_end_nodes(other_geom)
            overlapping_linestring = element.geometry.intersection(other_geom)
            overlap_start, overlap_end = geom_ops.get_start_end_nodes(overlapping_linestring)
            intersection_extents.update(
                {
                    tag: (
                     start_coord.distance(overlap_start), 
                     start_coord.distance(overlap_end)
                    )
                }
            )
        elif isinstance(other_geom, Polygon) and isinstance(element.geometry, Polygon):
            # The element will be a rank 0 element which means it is a load source
            # and the other_geom of the intersection below will be the physical
            # element of which the extents should be measured by.
            intersecting_region = intersection_below.intersecting_region
            other_geom = intersection_below.other_geometry
            other_geom_centerline = geom_ops.get_rectangle_centerline(other_geom)
            start_coord, _ = geom_ops.get_start_end_nodes(other_geom_centerline)
            intersecting_centerline = geom_ops.get_rectangle_centerline(intersecting_region)
            inter_start_coord, inter_end_coord = geom_ops.get_start_end_nodes(intersecting_centerline)
            intersection_extents.update(
                {
                    tag: (
                        start_coord.distance(inter_start_coord),
                        start_coord.distance(inter_end_coord),
                    )
                }
            )

    return intersection_extents



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
            if j_rank > i_rank:
                if i_geom.geom_type == j_geom.geom_type == "Polygon":
                    if not check_eligible_polygon_intersection(i_attrs['tag'], j_attrs['tag']):
                        continue
                intersection = geom_ops.get_intersection(i_geom, j_geom, j_attrs['tag'])
                if intersection is None: continue
                intersections_below.append(Intersection(*intersection))
            elif i_rank > j_rank:
                if i_geom.geom_type == j_geom.geom_type == "Polygon":
                    if not check_eligible_polygon_intersection(i_attrs['tag'], j_attrs['tag']):
                        continue
                intersection = geom_ops.get_intersection(j_geom, i_geom, j_attrs['tag'])
                # reaction_type
                if intersection is None: continue
                intersections_above.append(Intersection(*intersection))

        i_attrs["intersections_above"] = intersections_above
        i_attrs["intersections_below"] = intersections_below
    return intersected_annotations


def check_eligible_polygon_intersection(i_tag, j_tag) -> bool:
    """
    Returns True if 'i_tag' and 'j_tag' indicate the polygon elements
    are of the same 'family'
    """
    return i_tag[0] == j_tag[0]


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
    annots_prev = {}
    for page in descending_pages:
        if page != last_page:
            next_page = page - 1
            annots_here = annots_by_page[page]
            annots_below = annots_by_page[next_page]
            correspondents_above = {j_attrs['tag']: [] for j_attrs in annots_below.values()}
            correspondents_below = {}

            for i_annot, i_attrs in annots_here.items():
                i_page = i_annot.page
                i_tag = i_attrs['tag']
                i_geom = i_attrs["geometry"]
                i_rank = i_attrs["rank"]
                i_rxn_type = i_attrs.get('reaction_type', "point")
                for j_annot, j_attrs in annots_below.items():
                    j_attrs = annots_below[j_annot]
                    j_page = j_annot.page
                    j_geom = j_attrs["geometry"]
                    j_rank = j_attrs['rank']
                    j_tag = j_attrs['tag']
                    j_rxn_type = j_attrs.get('reaction_type', "point")
                    if i_geom.geom_type == "LineString" and j_geom.geom_type == "LineString":
                        continue # No correspondence between lines
                    correspondence_ratio = geom_ops.check_corresponds(i_geom, j_geom)
                    correspondents_below.setdefault(i_tag, [])
                    if correspondence_ratio and i_rank >= j_rank: # Same rank allowed to transfer in correspondents (e.g. column to column)
                        correspondents_below[i_tag].append(Correspondent(correspondence_ratio, j_geom, j_tag, other_rank=j_rank, other_reaction_type=j_rxn_type))
                        correspondents_above[j_tag].append(Correspondent(correspondence_ratio, i_geom, i_attrs['tag'], other_rank=i_rank, other_reaction_type=i_rxn_type))
                        corresponding_annotations[j_annot].setdefault("correspondents_above", [])
                        corresponding_annotations[i_annot].setdefault("correspondents_below", [])
                        corresponding_annotations[j_annot]["correspondents_above"] = correspondents_above[j_tag]
                        corresponding_annotations[i_annot]["correspondents_below"] = correspondents_below[i_tag]
                    else:
                        # Populate empty fields for annotations with no correspondents
                        corresponding_annotations[i_annot].setdefault("correspondents_below", [])
                        corresponding_annotations[i_annot].setdefault("correspondents_above", [])
        else:
            annots_last = annots_by_page[page]
            if len(descending_pages) == 1: 
                correspondents_above = {} # There are no correspondents above or below on a single page
            for i_annot, i_attrs in annots_last.items():
                i_tag = i_attrs['tag']
                corresponding_annotations[i_annot]['correspondents_above'] = correspondents_above.get(i_tag, [])
                corresponding_annotations[i_annot]['correspondents_below'] = []
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
