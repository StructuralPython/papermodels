from copy import deepcopy
from dataclasses import dataclass
from typing import Optional, Union, NamedTuple
import numpy as np
import numpy.typing as npt
from shapely import Point, LineString, Polygon, GeometryCollection
from shapely import wkt
from .annotation import Annotation
from ..paper.annotations import (
    parse_annotations,
    tag_parsed_annotations,
)
from ..geometry import geom_ops
import load_distribution as ld
import parse
import math
import tomli_w
import json

Geometry = Union[LineString, Polygon]

ELEMENT_ATTRS = {
    "tag",
    "geometry",
    "rank",
    "type",
    "length",
    "intersections_above",
    "intersections_below",
    "correspondents_above",
    "correspondents_below",
    "page_label",
    "reaction_type",
    "extent_line",
}


class Intersection(NamedTuple):
    """
    A class to represent an intersection of geometries
    """

    intersecting_region: Point | LineString
    other_geometry: Union[LineString, Polygon]
    other_tag: str
    other_overlap: Optional[LineString | Polygon] = None
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
    kwargs: Optional[dict]. Additional kwargs that are defined on the annotation are
        passed-through to the Element instance and stored here.
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
    kwargs: Optional[dict] = None
    extent_line: Optional[LineString] = None

    @property
    def extent_polygon(self):
        if self.extent_line is not None and self.geometry.geom_type == "LineString":
            return geom_ops.create_extent_polygon(self.geometry, self.extent_line)
        return None

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
        plane_id: Optional[str | int] = None,
    ):
        """
        Generates an Element from provided geometries
        """
        inters_above = []
        inters_below = []
        if intersections_above is not None:
            inters_above = [
                Intersection(
                    *geom_ops.get_intersection(elem_geom, above_geom, above_tag)
                )
                for above_tag, above_geom in intersections_above.items()
            ]

        if intersections_below is not None:
            inters_below = [
                Intersection(
                    *geom_ops.get_intersection(elem_geom, below_geom, below_tag)
                )
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
        annotations_w_intersect_corrs = get_geometry_correspondents(
            annotations_w_intersect
        )
        if correspond_with_like_only:
            filtered_annots = filter_correspondents(annotations_w_intersect_corrs)
        else:
            filtered_annots = annotations_w_intersect_corrs
        trib_area_geoms = np.array(
            [annot_attrs["geometry"] for annot_attrs in trib_annotations.values()]
        )
        matching_trib_poly = None
        elements = []
        for annot_attrs in filtered_annots.values():
            geometry = annot_attrs["geometry"]
            if geometry.geom_type == "LineString" and trib_annotations:
                intersection_mask = geometry.intersects(trib_area_geoms)
                intersection_lines = trib_area_geoms[intersection_mask]
                get_intersection_lengths = np.vectorize(lambda x: x.length)
                matching_trib_index = np.argmax(
                    get_intersection_lengths
                )  # The matching trib is the one that is mostly in the trib poly
                matching_trib_poly = trib_area_geoms[matching_trib_index]

            available_kwargs = {
                k: v for k, v in annot_attrs.items() if k not in ELEMENT_ATTRS
            }
            available_kwargs = available_kwargs or None
            element = cls(
                tag=annot_attrs["tag"],
                geometry=annot_attrs["geometry"],
                rank=annot_attrs["rank"],
                intersections_above=annot_attrs["intersections_above"],
                intersections_below=annot_attrs["intersections_below"],
                correspondents_above=annot_attrs["correspondents_above"],
                correspondents_below=annot_attrs["correspondents_below"],
                plane_id=annot_attrs.get("page_label", None),
                reaction_type=annot_attrs.get("reaction_type", "point"),
                trib_area=matching_trib_poly,
                kwargs=available_kwargs,
                extent_line=annot_attrs["extent_line"],
            )
            elements.append(element)
        return elements

    @property
    def supporting_geometries(self):
        acc = []
        for intersection_tuple in self.intersections.values():
            acc.append(intersection_tuple[1])
        return acc

    def get_collector_extents(self, relative: bool = True) -> dict[str, tuple]:
        """
        Returns a dict keyed by .tag attributes in self.intersections_below
        and with values representing the (start_x, end_x) locations where the collector
        prototype would spread over the other_geometry. The (start_x, end_x) locations
        refer to ordinates on other_geometry, not on the self.

        If 'relative' == False, then the values returned are absolute Point objects.
        """
        # When we have collector extents and an extent_polygon
        if self.trib_area is None and self.extent_polygon is not None:
            tagged_extents = {}
            for ib in self.intersections_below:
                region_start, region_end = geom_ops.get_start_end_nodes(
                    ib.intersecting_region
                )
                if ib.other_geometry.geom_type == "Polygon":
                    support_start, support_end = geom_ops.get_start_end_nodes(
                        geom_ops.get_rectangle_centerline(ib.other_geometry)
                    )
                else:  # LineString
                    support_start, support_end = geom_ops.get_start_end_nodes(
                        ib.other_geometry
                    )
                if relative:
                    extent_start = support_start.distance(region_start)
                    extent_end = support_start.distance(region_end)
                    tagged_extents.update({ib.other_tag: (extent_start, extent_end)})
                else:
                    tagged_extents.update({ib.other_tag: (support_start, support_end)})
        # When we have joist prototypes that have been drawn for all locations
        else:
            support_tags_by_geom = {
                geom_ops.clean_polygon_supports([ib.other_geometry], self.geometry)[
                    0
                ]: ib.other_tag
                for ib in self.intersections_below
            }
            support_geoms = list(support_tags_by_geom.keys())
            cleaned_supports_map = {}
            for idx, poly_support_geom in enumerate(support_geoms):
                clean_support_geom = support_geoms[idx]
                cleaned_supports_map.update({clean_support_geom: poly_support_geom})
                ordered_support_geoms = geom_ops.sort_supports(
                    self.geometry, support_geoms
                )
            try:
                extents = geom_ops.get_joist_extents(
                    self.geometry,
                    ordered_support_geoms,
                    self.trib_area,
                    extent_polygon=self.extent_polygon,
                )
            except (AssertionError, ValueError) as e:
                raise AssertionError(
                    f"No intersection within joist extents: {self.tag=}"
                )
            tagged_extents = {}
            for idx, extent in enumerate(extents):
                support_geom = ordered_support_geoms[idx]
                support_start, support_end = geom_ops.get_start_end_nodes(support_geom)
                support_tag = support_tags_by_geom[support_geom]
                extent_start = extent[0].distance(support_start)
                extent_end = extent[1].distance(support_start)
                if relative:
                    tagged_extents.update(
                        {support_tag: tuple(sorted((extent_start, extent_end)))}
                    )
                else:
                    tagged_extents.update({support_tag: (support_start, support_end)})

        return tagged_extents

    def get_transfer_extents(self) -> tuple[str, dict]:
        """
        Returns a tuple of str, extents_dict

        e.g.
        {"FB0.1": (2.5, 6.5, 0.3, 4.3)}

        Where the first two values describe the extents of the _transferred_ elemnt
        and the last two values describe the extents of the _transferring_ element.

        For the polygon element that has a linear reaction type.
        This element could have more than one intersection below
        if it overlaps multiple beams, for example. It can also
        have correspondents that need to have a portion of teh
        linear load applied to it.
        """
        intersection_extents = {}
        for intersection_below in self.intersections_below:
            intersection_below: Intersection
            tag = intersection_below.other_tag
            other_geom = intersection_below.other_geometry
            if isinstance(other_geom, LineString):
                below_start_coord, _ = geom_ops.get_start_end_nodes(
                    other_geom
                )  # extents are in reference to "below" geometry
                above_start_coord, _ = geom_ops.get_rectangle_centerline(
                    self.geometry
                ).coords
                above_start_coord = Point(above_start_coord)
                overlapping_linestring = self.geometry.intersection(other_geom)
                overlap_start, overlap_end = geom_ops.get_start_end_nodes(
                    overlapping_linestring
                )
                intersection_extents.update(
                    {
                        tag: (
                            below_start_coord.distance(overlap_start),
                            below_start_coord.distance(overlap_end),
                            above_start_coord.distance(
                                overlap_start
                            ),  # Extents in relation to transferring element
                            above_start_coord.distance(
                                overlap_end
                            ),  # Extents in relation to transferring element
                        )
                    }
                )
            elif isinstance(other_geom, Polygon) and isinstance(self.geometry, Polygon):
                # The element will be a rank 0 element which means it is a load source
                # and the other_geom of the intersection below will be the physical
                # element of which the extents should be measured by.
                intersecting_region = intersection_below.intersecting_region
                other_geom = intersection_below.other_geometry
                other_geom_centerline = geom_ops.get_rectangle_centerline(other_geom)
                below_start_coord, _ = geom_ops.get_start_end_nodes(
                    other_geom_centerline
                )
                above_start_coord, _ = geom_ops.get_rectangle_centerline(
                    self.geometry
                ).coords
                above_start_coord = Point(above_start_coord)
                intersecting_centerline = geom_ops.get_rectangle_centerline(
                    intersecting_region
                )
                inter_start_coord, inter_end_coord = geom_ops.get_start_end_nodes(
                    intersecting_centerline
                )
                intersection_extents.update(
                    {
                        tag: (
                            below_start_coord.distance(inter_start_coord),
                            below_start_coord.distance(inter_end_coord),
                            above_start_coord.distance(inter_start_coord),
                            above_start_coord.distance(inter_end_coord),
                        )
                    }
                )

        correspondent_extents = {}
        for correspondent_below in self.correspondents_below:
            tag = correspondent_below.other_tag
            other_geom = correspondent_below.other_geometry
            intersecting_region = self.geometry.intersection(other_geom)
            other_geom_centerline = geom_ops.get_rectangle_centerline(other_geom)
            below_start_coord, _ = geom_ops.get_start_end_nodes(other_geom_centerline)
            above_start_coord, _ = geom_ops.get_rectangle_centerline(
                self.geometry
            ).coords
            above_start_coord = Point(above_start_coord)
            intersecting_centerline = geom_ops.get_rectangle_centerline(
                intersecting_region
            )
            inter_start_coord, inter_end_coord = geom_ops.get_start_end_nodes(
                intersecting_centerline
            )
            correspondent_extents.update(
                {
                    tag: (
                        below_start_coord.distance(inter_start_coord),
                        below_start_coord.distance(inter_end_coord),
                        above_start_coord.distance(inter_start_coord),
                        above_start_coord.distance(inter_end_coord),
                    )
                }
            )

        return intersection_extents | correspondent_extents


def prioritize_correspondents(
    correspondents: list[Correspondent], family: str
) -> list[Correspondent]:
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
        if corr_rank == 0 and corr_family == family:  # If the element transfers out
            filtered_transfers.append(correspondent)
        elif corr_family == family:
            filtered_same_family.append(correspondent)
    if filtered_transfers:
        largest_overlap = sorted(
            filtered_transfers, key=lambda x: x.overlap_ratio, reverse=True
        )[0]
        return [largest_overlap]
    if filtered_same_family:
        largest_overlap = sorted(
            filtered_same_family, key=lambda x: x.overlap_ratio, reverse=True
        )[0]
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
        loading_array = np.array(
            [loading_area.geometry for loading_area in self.loading_geoms]
        )
        applied_loading_areas = []
        if self.trib_area is not None:
            intersecting_loads = self.trib_area.intersection(loading_array)
            for idx, intersecting_load in enumerate(intersecting_loads):
                if intersecting_load.is_empty or math.isclose(
                    intersecting_load.area, 0
                ):
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
        available_kwargs = self.kwargs or {}

        elem_model = {
            "element_attributes": {
                "tag": self.tag,
                "length": length,
                "orientation": orientation,
                "vert_correspondent_below": [
                    corr.other_tag for corr in self.correspondents_below
                ],
                "vert_correspondent_above": [
                    corr.other_tag for corr in self.correspondents_above
                ],
                "horz_intersects_above": [
                    inter.other_tag for inter in self.intersections_above
                ],
                "horz_intersects_below": [
                    inter.other_tag for inter in self.intersections_below
                ],
                "reaction_type": self.reaction_type,
                "rank": self.rank,
                "user_defined": available_kwargs,
            },
            "element_geometry": {
                "geometry": wkt.dumps(
                    self.geometry, trim=True, rounding_precision=precision
                ),
                "supports": support_locations,
            },
            "loads": {
                "point_loads": transfer_loads.get("point", []),
                "distributed_loads": transfer_loads.get("dist", []) + distributed_loads,
            },
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
            for ib in self.intersections_below:
                region = ib.intersecting_region
                if region.geom_type == "LineString":
                    coords_a, coords_b = region.coords
                    coords_a, coords_b = Point(coords_a), Point(coords_b)
                    # if ib.other_geometry
            support_locations = geom_ops.get_local_intersection_ordinates(
                start_coord,
                [intersection[0] for intersection in self.intersections_below],
            )
            terminus_supports = [
                any(
                    [
                        ib.other_geometry.contains(coords_a),
                        ib.other_geometry.contains(coords_b),
                    ]
                )
                for ib in self.intersections_below
            ]
            overlap_regions = [
                intersection.other_overlap for intersection in self.intersections_below
            ]
            supports_acc = []
            for idx, loc in enumerate(support_locations):
                overlap_region = overlap_regions[idx]
                terminus_support = (
                    terminus_supports[idx] if overlap_region is not None else None
                )
                overlap_length = 0.0
                if overlap_region is not None:
                    overlap_length = overlap_region.length

                supports_acc.append(
                    {
                        "location": round(loc, 3),
                        "overlap_length": round(overlap_length, 3),
                        "terminus": terminus_support,
                    }
                )

            # for idx, support_location in enumerate(support_locations):
            #     fixity = "roller"
            #     if idx == 0:
            #         fixity = "pin"
            #     supports_acc.append(
            #         {"location": round(support_location, precision), "fixity": fixity}
            #     )
            return sorted(supports_acc, key=lambda x: x["location"])
        else:
            return []

    def _get_transfer_loads(self, precision: int):
        """
        Calculates the transfer load locations from the intersections above
        """
        # It is possible to calculate the load eccentricity for columns based
        # on the getting the transfer_locations relative to the column centroid.
        # Perhaps a future feature.
        transfer_loads = {"point": [], "dist": []}
        if self.geometry.geom_type == "LineString":
            coords_a, coords_b = self.geometry.coords
        elif self.geometry.geom_type == "Polygon" and self.reaction_type == "linear":
            centerline = geom_ops.get_rectangle_centerline(self.geometry)
            coords_a, coords_b = centerline.coords
            coords_a, coords_b = Point(coords_a), Point(coords_b)
        elif self.geometry.geom_type == "Polygon" and self.reaction_type == "point":
            coords_a, coords_b = self.geometry.centroid, self.geometry.centroid

        # This applies to most scenarios
        if not self.geometry.geom_type == "Polygon" and self.reaction_type == "point":
            coords_a, coords_b = Point(coords_a), Point(coords_b)
            ordered_coords = geom_ops.order_nodes_positive([coords_a, coords_b])
            start_coord = ordered_coords[0]
            transfer_locations = geom_ops.get_local_intersection_ordinates(
                start_coord,
                [
                    intersection.intersecting_region
                    for intersection in self.intersections_above
                ],
            )
        elif self.geometry.geom_type == "Polygon" and self.reaction_type == "linear":
            coords_a, coords_b = Point(coords_a), Point(coords_b)
            ordered_coords = geom_ops.order_nodes_positive([coords_a, coords_b])
            start_coord = ordered_coords[0]
            transfer_locations = geom_ops.get_local_intersection_ordinates(
                start_coord,
                [
                    intersection.intersecting_region
                    for intersection in self.intersections_above
                ],
            )
        else:  # But not when it is a column
            start_coord = coords_a
            # This is where teh eccentricity can be calculated based on using the
            # intersection.intersecting_region instead of start_coord
            transfer_locations = geom_ops.get_local_intersection_ordinates(
                start_coord, [start_coord for intersection in self.intersections_above]
            )

        # Intersections
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
                transfer_loads["point"].append(point_load)
            elif transfer_type == "linear":
                if len(intersection_above.other_extents) == 2:
                    source_start_extent = []
                    source_end_extent = []
                elif len(intersection_above.other_extents) == 4:
                    source_start_extent = round(
                        intersection_above.other_extents[2], precision
                    )
                    source_end_extent = round(
                        intersection_above.other_extents[3], precision
                    )
                dist_load = self.create_distributed_load(
                    start_location=round(
                        intersection_above.other_extents[0], precision
                    ),
                    start_magnitude=1.0,
                    end_location=round(intersection_above.other_extents[1], precision),
                    end_magnitude=1.0,
                    transfer_source=f"{source_member}",
                    transfer_reaction_index=intersection_above.other_index,
                    transfer_source_start_extent=source_start_extent,
                    transfer_source_end_extent=source_end_extent,
                    occupancy="",
                    load_components={},
                    applied_area=0.0,
                    direction="gravity",
                )
                transfer_loads["dist"].append(dist_load)

        if self.geometry.geom_type == "Polygon":
            for correspondent_above in self.correspondents_above:
                if correspondent_above.other_reaction_type == "point":
                    point_load = self.create_point_load(
                        transfer_location=transfer_locations,
                        magnitude=0.0,
                        transfer_source=correspondent_above.other_tag,
                        transfer_reaction_index=0,
                        direction="gravity",
                    )
                    transfer_loads["point"].append(point_load)
                elif correspondent_above.other_reaction_type == "linear":
                    source_member = correspondent_above.other_tag
                    if correspondent_above.other_extents is None:
                        target_start_extent = []
                        target_end_extent = []
                        source_start_extent = []
                        source_end_extent = []
                    elif len(correspondent_above.other_extents) == 2:
                        target_start_extent = round(
                            correspondent_above.other_extents[0], precision
                        )
                        target_end_extent = round(
                            correspondent_above.other_extents[1], precision
                        )
                        source_start_extent = []
                        source_end_extent = []
                    elif len(correspondent_above.other_extents) == 4:
                        target_start_extent = round(
                            correspondent_above.other_extents[0], precision
                        )
                        target_end_extent = round(
                            correspondent_above.other_extents[1], precision
                        )
                        source_start_extent = round(
                            correspondent_above.other_extents[2], precision
                        )
                        source_end_extent = round(
                            correspondent_above.other_extents[3], precision
                        )

                    dist_load = self.create_distributed_load(
                        start_location=target_start_extent,
                        start_magnitude=1.0,
                        end_location=target_end_extent,
                        end_magnitude=1.0,
                        transfer_source=f"{source_member}",
                        transfer_reaction_index=0,
                        transfer_source_start_extent=source_start_extent,
                        transfer_source_end_extent=source_end_extent,
                        occupancy="",
                        load_components={},
                        applied_area=0.0,
                        direction="gravity",
                    )
                    transfer_loads["dist"].append(dist_load)

        return transfer_loads

    @staticmethod
    def create_point_load(
        transfer_location: str,
        magnitude: float,
        transfer_source: str,
        transfer_reaction_index: int,
        direction: str = "gravity",
    ):

        return {
            "location": transfer_location,
            "magnitude": magnitude,
            "transfer_source": transfer_source,
            "transfer_reaction_index": transfer_reaction_index,
            "direction": direction,
        }

    @staticmethod
    def create_distributed_load(
        start_location: float,
        start_magnitude: float,
        end_location: float,
        end_magnitude: float,
        transfer_source: str,
        transfer_reaction_index: int,
        transfer_source_start_extent: float,
        transfer_source_end_extent: float,
        occupancy: str,
        load_components: dict,
        applied_area: float,
        direction: str,
    ):

        return {
            "transfer_source": transfer_source,
            "transfer_reaction_index": transfer_reaction_index,
            "transfer_source_start_extent": transfer_source_start_extent,
            "transfer_source_end_extent": transfer_source_end_extent,
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
                self.geometry, self.applied_loading_areas
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
                        h=(end_x - start_x), b1=start_y, b2=end_y
                    )
                    try:
                        if self.reaction_type == "point":
                            # collectors from the JoistArrayModel have trib areas that reflect
                            # the actual trib area for the subelement. Thus, the area reatio
                            # should reflect the percentage of area from each distributed load
                            # sub-section and what percentage that reflects of the complete
                            # distributed load system which may consist of several trapezoids
                            area_ratio = area_dist_load / total_polygon_area
                        # THIS GIVES THE CORRECT TRAPEZOID RATIO FOR COLLECTORTRIBMODEL
                        elif self.reaction_type == "linear":
                            # In the CollectorTribModel the trib area reflects the size of a whole
                            # spread of joists which will be reduced down to a reaction over a unit length.
                            # In this case, we need an area ratio to reflect the percentage of area that
                            # an area load covers in relation to the total trib area.
                            area_ratio = (
                                self.applied_loading_areas[idx][0].area
                                / self.trib_area.area
                            )
                    except ZeroDivisionError:
                        continue  # Skip this dist load if there is no polygon area
                    intersected_poly, applied_loading = self.applied_loading_areas[idx]
                    if area_ratio == 0.0 and intersected_poly.area == 0.0:
                        continue  # Skip this dist load if there is no intersection area
                    dist_load = {
                        "transfer_source": "",
                        "transfer_reaction_index": "",
                        "occupancy": applied_loading.occupancy,
                        "load_components": applied_loading.load_components or [],
                        "applied_area": round(intersected_poly.area, precision),
                        "applied_area_ratio": area_ratio,
                        "start_loc": round(start_x, precision),
                        "start_magnitude": round(start_y, precision),
                        "end_loc": round(end_x, precision),
                        "end_magnitude": round(end_y, precision),
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
            other_tag = intersection.other_tag.split("-")[
                0
            ]  # Sub-elements have a hyphen in their name
            if other_tag in predecessors:
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
            reaction_type=elem.reaction_type,
            kwargs=elem.kwargs,
        )


def create_element_filter(
    page_idxs: Optional[list[int]] = None,
    tags: Optional[list[str]] = None,
    element_types: Optional[list[str]] = None,
    user_defined: Optional[dict[str, str]] = None,
    exclude_tags: Optional[list[str]] = None,
    exclude_element_types: Optional[list[str]] = None,
    exclude_user_defined: Optional[dict[str, str]] = None,
) -> callable:
    """
    Returns a function with this signature:
        func(element: Element) -> bool

    The returned function returns True if the supplied
    element matches the criteria defined in this function.

    'page_idxs' - A list of ints. Any collector elements on those page
        indexes will be included.
    'tags' - A list of joist tags to be included by the filter
    'element_types' - A list of element types that exist in the graph.
        In this context, an "element_type" refers to the tag prefix and
        not the designation of "collector" or "transfer". The element
        will return True if the element matches ANY of the element_types
        listed (e.g. ['FB', 'SJ', 'CT'])
    'user_defined' - A dict of values assigned manually to an element.
        If ALL the dict of values provided matches the kwarg values in an
        element, it will be included.
    'exclude_tags': A list of tags that will be excluded from the filter
    'exclude_element_types': A list of element types which will be excluded
    'exclude_user_defined': If ANY of the key/values in this dict match the
        kwarg values in an element, it will be excluded.
    """

    def filter_function(element: Element) -> bool:
        """
        Returns True if the element matches the filter criteria.
        Returns False otherwise.
        """
        include_by_page = element.plane_id in page_idxs if page_idxs else True
        include_by_tag = element.tag in tags if tags else True
        include_by_element_type = (
            any([element.tag.startswith(e_type) for e_type in element_types])
            if element_types
            else True
        )

        nonlocal user_defined
        if user_defined is None:
            include_by_user_defined = True
        else:
            user_defined = user_defined or {}
            user_defined_acc = []
            for k, v in user_defined.items():
                if element.kwargs is None:
                    user_defined_acc.append(False)
                    continue
                ev = element.kwargs.get(k)
                user_defined_acc.append(ev == v)
            include_by_user_defined = all(user_defined_acc)

        include_element = all(
            [
                include_by_page,
                include_by_tag,
                include_by_element_type,
                include_by_user_defined,
            ]
        )

        exclude_by_tag = element.tag not in exclude_tags if exclude_tags else True
        exclude_by_element_type = (
            not any(
                [element.tag.startswith(e_type) for e_type in exclude_element_types]
            )
            if exclude_element_types
            else True
        )

        nonlocal exclude_user_defined
        if exclude_user_defined is None:
            exclude_by_user_defined = True
        else:
            exclude_user_defined = exclude_user_defined or {}
            exclude_user_defined_acc = []
            for k, v in exclude_user_defined.items():
                if element.kwargs is None:
                    exclude_user_defined_acc.append(False)
                    continue
                ev = element.kwargs.get(k)
                exclude_user_defined_acc.append(ev == v)
            exclude_by_user_defined = not any(exclude_user_defined_acc)

        return all(
            [
                include_element,
                exclude_by_tag,
                exclude_by_element_type,
                exclude_by_user_defined,
            ]
        )

    return filter_function


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
        i_geom = i_attrs["geometry"]
        i_extent_line = i_attrs["extent_line"]
        i_extent_poly = geom_ops.create_extent_polygon(i_geom, i_extent_line)
        i_attrs.setdefault("intersections_below", [])
        i_attrs.setdefault("intersections_above", [])
        for j_annot in annots:
            j_attrs = intersected_annotations[j_annot]
            try:
                j_rank = j_attrs["rank"]
            except KeyError:
                print(j_annot, j_attrs)
                raise ValueError
            j_page = j_annot.page
            j_geom = j_attrs["geometry"]
            if i_geom.is_empty or j_geom.is_empty:
                continue
            i_tag = i_attrs["tag"]
            j_tag = j_attrs["tag"]
            if i_page != j_page:
                continue
            if j_rank > i_rank:  # When i transfers to j
                if i_geom.geom_type == j_geom.geom_type == "Polygon":
                    if not check_eligible_polygon_intersection(
                        i_attrs["tag"], j_attrs["tag"]
                    ):
                        continue
                # Use the extent polygon to find intersections (if it exists)
                extent_intersection = False
                if (
                    # i_rank == 0
                    # and i_extent_poly is not None
                    i_extent_poly is not None
                    and check_eligible_collector_extent_polygon_intersection(
                        j_geom.geom_type, j_attrs["reaction_type"]
                    )
                ):
                    intersection_below = geom_ops.get_intersection(
                        i_geom, j_geom, j_tag, i_extent_poly
                    )
                    j_intersection_above = geom_ops.get_intersection(
                        i_geom, j_geom, i_tag, i_extent_poly
                    )
                else:
                    intersection_below = geom_ops.get_intersection(
                        i_geom, j_geom, j_tag
                    )
                    j_intersection_above = geom_ops.get_intersection(
                        j_geom, i_geom, i_tag
                    )

                if intersection_below is None:
                    continue

                i_attrs["intersections_below"].append(
                    Intersection(
                        *intersection_below,
                        other_reaction_type=j_attrs["reaction_type"],
                    )
                )
                if j_intersection_above is not None:
                    j_attrs.setdefault("intersections_above", [])
                    j_attrs["intersections_above"].append(
                        Intersection(
                            *j_intersection_above,
                            other_reaction_type=i_attrs["reaction_type"],
                        )
                    )

    return intersected_annotations


def check_eligible_polygon_intersection(i_tag, j_tag) -> bool:
    """
    Returns True if 'i_tag' and 'j_tag' indicate the polygon elements
    are of the same 'family'
    """
    return i_tag[0] == j_tag[0]


def check_eligible_collector_extent_polygon_intersection(
    j_geomtype: str, j_reaction_type: str
) -> bool:
    """
    Returns True if the "j" annotation is eligible to receive collector reactions
    i.e. is not a column/post
    """
    if j_geomtype == "Polygon":
        if j_reaction_type == "linear":
            return True
        else:
            return False
    else:  # LineString
        return True


def get_geometry_correspondents(
    tagged_annotations: dict[Annotation, dict],
) -> dict[Annotation, dict]:
    """
    Returns a copy of 'tagged_annotations' with a 'correspondents' field added to that
    attribute's dictionary of each Annotation key.
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
            # correspondents_above = {
            #     j_attrs["tag"]: [] for j_attrs in annots_below.values()
            # }
            # correspondents_above = {}
            # correspondents_below = {}

            for i_annot, i_attrs in annots_here.items():
                corresponding_annotations[i_annot].setdefault(
                    "correspondents_below", []
                )
                i_page = i_annot.page
                i_tag = i_attrs["tag"]
                i_geom = i_attrs["geometry"]
                i_rank = i_attrs["rank"]
                i_rxn_type = i_attrs.get("reaction_type", "point")
                for j_annot, j_attrs in annots_below.items():
                    corresponding_annotations[j_annot].setdefault(
                        "correspondents_above", []
                    )
                    j_attrs = annots_below[j_annot]
                    j_page = j_annot.page
                    j_geom = j_attrs["geometry"]
                    j_rank = j_attrs["rank"]
                    j_tag = j_attrs["tag"]
                    j_rxn_type = j_attrs.get("reaction_type", "point")
                    if (
                        i_geom.geom_type == "LineString"
                        and j_geom.geom_type == "LineString"
                    ):
                        continue  # No correspondence between lines
                    correspondence_ratio = geom_ops.check_corresponds(i_geom, j_geom)
                    if (
                        correspondence_ratio and i_rank >= j_rank
                    ):  # Same rank allowed to transfer in correspondents (e.g. column to column)
                        corr_below = Correspondent(
                            correspondence_ratio,
                            j_geom,
                            j_tag,
                            other_rank=j_rank,
                            other_reaction_type=j_rxn_type,
                        )
                        corr_above = Correspondent(
                            correspondence_ratio,
                            i_geom,
                            i_attrs["tag"],
                            other_rank=i_rank,
                            other_reaction_type=i_rxn_type,
                        )
                        if (
                            corr_below
                            not in corresponding_annotations[i_annot][
                                "correspondents_below"
                            ]
                        ):
                            corresponding_annotations[i_annot][
                                "correspondents_below"
                            ].append(corr_below)
                        if (
                            corr_above
                            not in corresponding_annotations[j_annot][
                                "correspondents_above"
                            ]
                        ):
                            corresponding_annotations[j_annot][
                                "correspondents_above"
                            ].append(corr_above)

                    else:
                        # Populate empty fields for annotations with no correspondents
                        corresponding_annotations[i_annot].setdefault(
                            "correspondents_below", []
                        )
                        corresponding_annotations[i_annot].setdefault(
                            "correspondents_above", []
                        )
        else:
            correspondents_above = {}
            annots_last = annots_by_page[page]
            if len(descending_pages) == 1:
                correspondents_above = (
                    {}
                )  # There are no correspondents above or below on a single page
            for i_annot, i_attrs in annots_last.items():
                i_tag = i_attrs["tag"]
                corresponding_annotations[i_annot]["correspondents_above"] = (
                    correspondents_above.get(i_tag, [])
                )
                corresponding_annotations[i_annot]["correspondents_below"] = []
        if prev_page is None:
            prev_page = page
    return corresponding_annotations


def trim_cantilevers(element: Element, abs_tol: Optional[float] = 0.02):
    """
    Mutates the geometry in node elements so that any cantilevers which are
    below self.cantilever_abs_tol are removed from
    the geometry and the geometry spans exactly from support to support.
    """
    new_element = deepcopy(element)
    geometry = joist_prototype = element.geometry
    if geometry.geom_type == "LineString" and element.extent_line is None:
        orig_support_geoms = [ib.other_geometry for ib in element.intersections_below]
        support_geoms = geom_ops.clean_polygon_supports(orig_support_geoms, geometry)
        support_geoms = geom_ops.sort_supports(geometry, support_geoms)
        ordered_geom = LineString(
            geom_ops.order_nodes_positive(
                [Point(geometry.coords[0]), Point(geometry.coords[-1])]
            )
        )
        cantilevers = geom_ops.get_cantilever_segments(
            ordered_geom, support_geoms, abs_tol=abs_tol
        )
        # ordered_geom = LineString(
        #     geom_ops.order_nodes_positive(
        #         [Point(geometry.coords[0]), Point(geometry.coords[-1])]
        #     )
        # )
        # cantilevers = geom_ops.get_cantilever_segments(
        #     ordered_geom, intersection_points, rel_tol=rel_tol, abs_tol=abs_tol
        # )
        start_point, end_point = Point(ordered_geom.coords[0]), Point(
            ordered_geom.coords[-1]
        )

        if (cantilevers["A"] == 0.0) and (cantilevers["A"] != cantilevers["A_orig"]):
            start_point = cantilevers["A_intersection"]
        if (cantilevers["B"] == 0.0) and (cantilevers["B"] != cantilevers["B_orig"]):
            end_point = cantilevers["B_intersection"]
        new_geometry = LineString([start_point, end_point])  # type: ignore
        new_element.geometry = new_geometry
        intersection_checks = [
            new_geometry.intersects(ib.other_geometry)
            for ib in element.intersections_below
        ]
        new_intersections_below = [
            Intersection(
                intersecting_region=ib.intersecting_region,
                other_tag=ib.other_tag,
                other_overlap=ib.other_overlap,
                other_geometry=ib.other_geometry,
                other_reaction_type=ib.other_reaction_type,
                other_extents=ib.other_extents,
                other_index=ib.other_index,
            )
            for ib in element.intersections_below
        ]
        new_element.intersections_below = new_intersections_below
    return new_element


def align_frames_to_centroids(element: Element):
    """
    Mutates the geometry in node elements so that any cantilevers which are
    below self.cantilever_abs_tol are removed from
    the geometry and the geometry spans exactly from support to support.
    """
    new_element = deepcopy(element)
    geometry = element.geometry
    if geometry.geom_type != "LineString":
        return new_element
    start_point, end_point = geometry.coords
    start_point, end_point = Point(start_point), Point(end_point)
    start_support = None
    end_support = None
    new_start_point = None
    new_end_point = None
    if geometry.geom_type == "LineString":
        new_intersections = []
        for ib in new_element.intersections_below:
            ib: Intersection
            support_geom = ib.other_geometry
            support_reaction_type = ib.other_reaction_type
            intersecting_region = ib.intersecting_region
            if support_geom.contains(start_point):
                start_support = support_geom
            if support_geom.contains(end_point):
                end_support = support_geom
            overlap_region = ib.other_overlap
            if support_geom.geom_type == "Polygon" and support_reaction_type == "point":
                # intersecting_region = support_geom.centroid
                intersecting_region = geom_ops.get_projected_support_centroid(
                    geometry, support_geom
                )
                overlap_region = support_geom.intersection(geometry)
            elif (
                support_geom.geom_type == "Polygon"
                and support_reaction_type == "linear"
            ):
                intersecting_region = geom_ops.get_projected_support_centerline(
                    geometry, support_geom
                )
                overlap_region = support_geom.intersection(geometry)

            if support_geom == start_support:
                new_start_point = intersecting_region
            elif support_geom == end_support:
                new_end_point = intersecting_region
            new_intersection = Intersection(
                intersecting_region,
                ib.other_geometry,
                ib.other_tag,
                other_overlap=overlap_region,
                other_index=ib.other_index,
                other_reaction_type=ib.other_reaction_type,
                other_extents=ib.other_extents,
            )
            new_intersections.append(new_intersection)
        if new_start_point is None:
            new_start_point = start_point
        if new_end_point is None:
            new_end_point = end_point
        new_geom = LineString([new_start_point, new_end_point])
        new_element.intersections_below = new_intersections
        new_element.geometry = new_geom
    return new_element


def filter_correspondents(
    tagged_annotations: dict[Annotation, dict],
) -> dict[Annotation, dict]:
    """
    Applies the filter from prioritize_correspondents and back-propagates
    the filtered correspodents to the correspondents_above:
    """
    copied_annotations = tagged_annotations.copy()
    correspondents_above = {}
    for annot, annot_attrs in tagged_annotations.items():
        element_family = annot_attrs["tag"][0]
        filtered_below = prioritize_correspondents(
            annot_attrs["correspondents_below"], element_family
        )
        copied_annotations[annot]["correspondents_below"] = filtered_below
        for corr in filtered_below:
            corr: Correspondent
            overlap_ratio = corr.overlap_ratio
            overlap_area = overlap_ratio * annot_attrs["geometry"].area
            above_overlap_ratio = overlap_area / corr.other_geometry.area
            above_tag = corr.other_tag
            corr_above = Correspondent(
                overlap_ratio=above_overlap_ratio,
                other_geometry=annot_attrs["geometry"],
                other_tag=annot_attrs["tag"],
                other_rank=annot_attrs["rank"],
                other_reaction_type=annot_attrs["reaction_type"],
            )
            correspondents_above.setdefault(above_tag, [])
            correspondents_above[above_tag].append(corr_above)

    for copied_annot, copied_annot_attrs in copied_annotations.items():
        copied_annot_attrs["correspondents_above"] = correspondents_above.get(
            copied_annot_attrs["tag"], []
        )

    return copied_annotations


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
