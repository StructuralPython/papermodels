from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Any, Optional
import warnings
import numpy as np
from shapely import (
    LineString,
    Point,
    MultiLineString,
    Polygon,
    MultiPoint,
    convex_hull,
    GeometryCollection,
    box,
    set_precision,
)
import shapely.ops as ops

from papermodels.datatypes.element import (
    Element,
    Intersection,
    trim_cantilevers,
    align_frames_to_centroids,
)
from papermodels.geometry import geom_ops
import load_distribution as ld

from rich import print


@dataclass
class CollectorTribModel:
    element: Element
    trib_width: float = 1.0
    reaction_type: str = "linear"
    use_subelements: bool = False
    cantilever_tolerance: float = 1e-2

    def __post_init__(self):
        if self.element.extent_polygon is not None:
            self.use_subelements = True

    def __call__(self):
        """
        Generates a representative trib area for the joist prototype.

        Assumptions:
        - The supports are assumed to be orthogonal
        - IF EXTENT LINES ARE USED, then cantilevers are not supported
            (joist prototypes with extent lines using this model will
            have their cantilevers trimmed off...currently).
        - The loading is consistent for all joists within the spread
        - The joist represents a one-of-many similar elements within
            the spread and the spread thus represents a linear reaction
            over the supports.
        """
        e = self.element
        geom = e.geometry
        collector_extents = e.get_collector_extents(relative=False)
        left = []
        right = []
        for extent in collector_extents.values():
            p0_relation = geom_ops.relate_point_to_line(extent[0], geom)
            if p0_relation in (("left", "above"), ("right", "above")):
                left.append(extent[0])
                right.append(extent[1])
            elif p0_relation in (("left", "below"), ("right", "below")):
                left.append(extent[1])
                right.append(extent[0])
        left_dist = [geom.distance(extent) for extent in left]
        right_dist = [geom.distance(extent) for extent in right]
        left_minimum_idx = left_dist.index(min(left_dist))
        right_minimum_idx = right_dist.index(min(right_dist))

        joist_vector = np.abs(geom_ops.get_direction_vector(geom))
        # This is one of the places where orthogonality is assumed
        joist_orientation = None
        if joist_vector[0] > joist_vector[1]:
            joist_orientation = "horizontal"
        elif joist_vector[1] > joist_vector[0]:
            joist_orientation = "vertical"
        else:
            print(f"JOIST ORIENTATION VERIANT: {geom=}")

        if joist_orientation == "vertical":
            minx = left[left_minimum_idx].coords[0][0]
            miny = geom.coords[0][1]
            maxx = right[right_minimum_idx].coords[0][0]
            maxy = geom.coords[1][1]
        elif joist_orientation == "horizontal":
            minx = geom.coords[0][0]
            miny = left[left_minimum_idx].coords[0][1]
            maxx = geom.coords[1][0]
            maxy = right[right_minimum_idx].coords[0][1]
        trib_area = box(minx, miny, maxx, maxy)
        # trib_area = e.geometry.buffer(self.trib_width/2)
        if not self.use_subelements:
            collector_element = Element(
                e.geometry,
                e.tag,
                0,
                e.intersections_above,
                e.intersections_below,
                e.correspondents_above,
                e.correspondents_below,
                e.plane_id,
                e.element_type,
                e.subelements,
                trib_area=trib_area,
                reaction_type="linear",
                kwargs=e.kwargs,
                extent_line=e.extent_line,
            )
        else:
            ext_poly = e.extent_polygon
            joist_prototype = e.geometry
            support_lines = {}
            for ib in e.intersections_below:
                if ib.other_geometry.geom_type == "Polygon":
                    support_line = geom_ops.get_rectangle_centerline(ib.other_geometry)
                    support_lines.update({support_line: ib.other_tag})
                elif ib.other_geometry.geom_type == "LineString":
                    support_lines.update({ib.other_geometry: ib.other_tag})
            support_geoms = [ib.other_geometry for ib in e.intersections_below]

            # 1. Find polygon extent edges that intersect with joist prototype
            # These will be our extent boundaries for the length of the prototype
            poly_edge_points = list(
                zip(ext_poly.exterior.coords, ext_poly.exterior.coords[1:])
            )
            start_edge = None
            end_edge = None
            start_point, end_point = geom_ops.order_nodes_positive(
                [Point(coord) for coord in geom.coords]
            )
            for edge_points in poly_edge_points:
                edge = LineString(edge_points)
                if edge.intersects(geom) and (
                    start_point.distance(edge) < end_point.distance(edge)
                ):
                    start_edge = edge
                elif edge.intersects(geom) and (
                    end_point.distance(edge) < start_point.distance(edge)
                ):
                    end_edge = edge
            # Either the start or end edge should work since
            # orthogonality is assumed.

            # 2. Find support geoms which intersect with start and end edges
            start_supports = []
            end_supports = []
            joist_vector = np.abs(geom_ops.get_direction_vector(geom))
            # This is one of the places where orthogonality is assumed
            joist_orientation = None
            if joist_vector[0] > joist_vector[1]:
                joist_orientation = "horizontal"
            elif joist_vector[1] > joist_vector[0]:
                joist_orientation = "vertical"
            else:
                print(f"JOIST ORIENTATION VERIANT: {geom=}")
            support_centroids = [
                support_geom.centroid for support_geom in support_geoms
            ]
            if joist_orientation == "horizontal":
                start_support = min(support_centroids, key=lambda x: x.coords[0][0])
                end_support = max(support_centroids, key=lambda x: x.coords[0][0])
                start_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        start_support.coords[0][0],
                        support_line.coords[0][0],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]
                end_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        end_support.coords[0][0],
                        support_line.coords[0][0],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]
            elif joist_orientation == "vertical":
                start_support = min(support_centroids, key=lambda x: x.coords[0][1])
                end_support = max(support_centroids, key=lambda x: x.coords[0][1])
                start_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        start_support.coords[0][1],
                        support_line.coords[0][1],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]
                end_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        end_support.coords[0][1],
                        support_line.coords[0][1],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]

            # 2b. Get intermediate supports
            intermediate_support_lines = []
            for support_line in support_lines:
                if support_line not in start_supports + end_supports:
                    intermediate_support_lines.append(support_line)

            # 3. Generate overlap regions

            overlap_polys = []
            for start_support in start_supports:
                for end_support in end_supports:
                    overlap_poly = None
                    pa0, pa1 = start_support.coords
                    pb0, pb1 = end_support.coords
                    if joist_orientation == "vertical":
                        overlap_region = ld.get_overlap_coords(
                            pa0[0], pa1[0], pb0[0], pb1[0]
                        )
                        if overlap_region is not None:
                            overlap_poly = box(
                                overlap_region[0], pa0[1], overlap_region[1], pb1[1]
                            )
                    elif joist_orientation == "horizontal":
                        a_sort = sorted([pa0, pa1], key=lambda x: x[1])
                        b_sort = sorted([pb0, pb1], key=lambda x: x[1])
                        pa0 = a_sort[0]
                        pa1 = a_sort[1]
                        pb0 = b_sort[0]
                        pb1 = b_sort[1]
                        overlap_region = ld.get_overlap_coords(
                            pa0[1], pa1[1], pb0[1], pb1[1]
                        )
                        if overlap_region is not None:
                            overlap_poly = box(
                                pa0[0], overlap_region[0], pb1[0], overlap_region[1]
                            )
                    if overlap_poly is not None:
                        overlap_within_extent = ext_poly.intersection(overlap_poly)
                        overlap_polys.append(overlap_within_extent)

            # 5. Do overlap polys intersect with intermediate supports?
            #    if so, break the support as required.
            revised_poly_overlaps = []
            for overlap_poly in set(overlap_polys):
                split_polys = []
                for intermediate_support in intermediate_support_lines:
                    if intermediate_support.intersects(overlap_poly):
                        inter_coords = intermediate_support.coords
                        poly_splits = geom_ops.split_polygon(
                            overlap_poly, joist_orientation, inter_coords
                        )
                        split_polys += poly_splits
                if not split_polys:
                    revised_poly_overlaps.append(overlap_poly)
                else:
                    revised_poly_overlaps += split_polys
            sorted_poly_overlaps = sorted(
                revised_poly_overlaps,
                key=lambda x: (x.centroid.coords[0][0], x.centroid.coords[0][1]),
            )

            joist_prototype_geometries = []
            for overlap_poly in sorted_poly_overlaps:
                overlap_poly: Polygon
                overlap_edge_points = list(
                    zip(overlap_poly.exterior.coords, overlap_poly.exterior.coords[1:])
                )
                for pi, pj in overlap_edge_points:
                    edge_ls = LineString([pi, pj])
                    if geom_ops.check_2d_linestring_parallel(
                        edge_ls, start_edge, tol=0.01
                    ):
                        # Need to translate the original joist prototype to the new position
                        new_joist = geom_ops.translate_joist_to_point(
                            joist_prototype,
                            joist_orientation,
                            intersection_point=edge_ls.centroid,
                        )
                        trimmed_joist = new_joist.intersection(overlap_poly)
                        # We only need to hit one edge of the overlap so we can break here
                        break
                joist_prototype_geometries.append(trimmed_joist)

            # 7. Create an Element for each new joist prototype geometries
            subelements = []
            sorted_joist_geoms = joist_prototype_geometries
            for idx, joist_geom in enumerate(sorted_joist_geoms):
                intersections = []
                total_new_subs = len(joist_prototype_geometries)
                z_fill_qty = math.floor(math.log10(total_new_subs))
                index = f"{idx}".zfill(z_fill_qty)
                subelement_tag = f"{e.tag}-{index}"
                trib_area = sorted_poly_overlaps[idx]
                assert joist_geom.intersects(trib_area)
                for support_geom in support_geoms:
                    support_overlap = None
                    support_intersection = joist_geom.intersection(
                        support_geom, grid_size=1e-3
                    )
                    if not support_intersection:
                        continue
                    if support_geom.geom_type == "Polygon":
                        support_overlap = support_intersection
                        support_line = geom_ops.get_rectangle_centerline(
                            support_geom
                        )  # geom_ops.clean_polygon_supports([support_geom], joist_geom)[0]
                        intersecting_region = support_geom.exterior.intersection(
                            joist_geom, grid_size=1e-3
                        )
                    elif support_geom.geom_type == "LineString":
                        support_line = support_geom
                        intersecting_region = support_intersection

                    tag = support_lines[support_line]
                    # HERE: Previous behaviour was to intersect with the centerline but that is no longer a requirement
                    # Joist geom needs to be rebuilt to ensure it hits the wall centerline

                    # if intersecting_region.is_empty:
                    #     continue
                    intersection = Intersection(
                        intersecting_region=intersecting_region,
                        other_geometry=support_geom,
                        other_tag=tag,
                        other_overlap=support_overlap,
                        other_reaction_type=(
                            "linear" if support_geom.geom_type == "Polygon" else "point"
                        ),
                    )
                    intersections.append(intersection)
                subelement = Element(
                    geometry=joist_geom,
                    tag=subelement_tag,
                    rank=e.rank,
                    intersections_above=[],
                    intersections_below=intersections,
                    correspondents_above=[],
                    correspondents_below=[],
                    plane_id=e.plane_id,
                    element_type=e.element_type,
                    subelements=None,
                    trib_area=trib_area,
                    reaction_type="linear",
                    kwargs=e.kwargs,
                    # extent_polygon=e.extent_polygon,
                )
                aligned_subelement = align_frames_to_centroids(subelement)
                subelements.append(aligned_subelement)

            # 8. Return subelements
            collector_element = Element(
                e.geometry,
                tag=e.tag,
                rank=e.rank,
                intersections_above=e.intersections_above,
                intersections_below=e.intersections_below,
                correspondents_above=e.correspondents_above,
                correspondents_below=e.correspondents_below,
                plane_id=e.plane_id,
                element_type=e.element_type,
                subelements=subelements,
                trib_area=e.trib_area,
                reaction_type="linear",
                kwargs=e.kwargs,
                extent_line=e.extent_line,
            )
        return collector_element


def collector_trib_model(
    element: Element, trib_width: float, reaction_type: str = "linear"
) -> Element:
    """
    An alias for CollectorTribModel.__call__() for temporary
    backwards compatibility.
    """
    model = CollectorTribModel(element, trib_width, reaction_type)
    return model()


class JoistArrayModel:
    """
    Models a spread of joists over a region where the distance
    between the supports may vary linearly.
    """

    def __init__(
        self,
        element: Optional[Element] = None,
        spacing: float = 1,
        # joist_id: str,
        # joist_prototype: LineString,
        initial_offset: float | int = 0.0,
        joist_at_start: bool = True,
        joist_at_end: bool = False,
        cantilever_tolerance: float = 1e-1,
    ):
        self.joist_prototype = LineString(
            geom_ops.get_start_end_nodes(element.geometry)
        )
        self.element = element
        self.extent_polygon = element.extent_polygon
        self.joist_supports = {}
        for ib in element.intersections_below:
            if (
                not ib.other_geometry.intersects(self.joist_prototype)
                and not self.extent_polygon
            ):
                # This condition can exist when extent lines are used
                continue
            if ib.other_geometry.geom_type == "Polygon":
                support = geom_ops.clean_polygon_supports(
                    [ib.other_geometry], self.joist_prototype, self.extent_polygon
                )
                self.joist_supports.update({support[0]: ib.other_reaction_type})
            else:
                support = ib.other_geometry
                self.joist_supports.update({support: ib.other_reaction_type})

        # TODO: FIX THIS HACK THAT SEEMS TO ALTER THE DIRECTION OF SOME CANTILEVERS
        if self.extent_polygon:
            support_centroids = [geom.centroid for geom in self.joist_supports.keys()]
            ordered_centroids = geom_ops.order_nodes_positive(support_centroids)
            ordered_supports = []
            for point in ordered_centroids:
                for support in self.joist_supports.keys():
                    if point == support.centroid:
                        ordered_supports.append(support)
                        break
            self._supports = ordered_supports
        else:
            self._supports = geom_ops.sort_supports(
                self.joist_prototype, self.joist_supports.keys()
            )

        self.joist_support_tags = [ib.other_tag for ib in element.intersections_below]
        self.id = element.tag
        self.plane_id = element.plane_id
        self.elem_kwargs = element.kwargs
        self.spacing = (
            spacing  # Need to include this in the legend and thus, the Element
        )
        self.initial_offset = float(initial_offset)
        self._joist_prototype = self.joist_prototype
        self._cantilever_tolerance = cantilever_tolerance
        self.use_subelements = True
        try:
            self._extents = geom_ops.get_joist_extents(
                self.joist_prototype,
                self._supports,
                trib_area=None,
                extent_polygon=self.extent_polygon,
            )
            # self._extents = geom_ops.get_joist_extents(
            #     self.joist_prototype, self.joist_supports, trib_area=self.extent_polygon, extent_polygon=self.extent_polygon
            # )
        except AssertionError as e:
            raise AssertionError(
                f"No intersection within joist extents: {element.tag=}"
            )
        self._cantilevers = geom_ops.get_cantilever_segments(
            self.joist_prototype, self._supports, abs_tol=0
        )
        self.vector_parallel = geom_ops.get_direction_vector(self.joist_prototype)

        self.vector_normal = geom_ops.rotate_90_vector(self.vector_parallel, ccw=True)

        self.joist_at_start = joist_at_start
        self.joist_at_end = joist_at_end
        self.joist_locations = geom_ops.get_joist_locations(
            self.get_extent_edge("start"),
            self.get_extent_edge("end"),
            self.spacing,
            self.initial_offset,
            self.joist_at_start,
        )
        self.joist_geoms = [
            self.generate_joist_geom(idx) for idx, _ in enumerate(self.joist_locations)
        ]
        self.joist_trib_widths = [
            self.get_joist_trib_widths(idx) for idx, _ in enumerate(self.joist_geoms)
        ]
        self.joist_trib_areas = [
            self.generate_trib_area(idx) for idx, _ in enumerate(self.joist_geoms)
        ]

    # def __repr__(self):
    #     return class_representation(self)

    @classmethod
    def create_subelements(
        cls,
        element: Element,
        extents: Optional[Polygon] = None,
        spacing: Optional[float] = 1.0,
        initial_offset: float | int = 0.0,
        joist_at_start: bool = True,
        joist_at_end: bool = False,
        cantilever_tolerance: float = 1e-2,
    ) -> JoistArrayModel:
        if element.geometry.geom_type != "LineString":
            return None
        joist_array = cls(
            element,
            spacing,
            initial_offset,
            joist_at_start,
            joist_at_end,
            cantilever_tolerance,
            extents,
        )
        # joist_array.show_svg()
        return joist_array.to_subelements()

    def to_subelements(self):
        """
        An alias for __call__ for temporary backwards compatibility
        """
        self()

    def __call__(self) -> list[Element]:
        """
        Returns the sub-joists in the JoistArray (self) as Element
        """
        e = self.element
        subelements = []
        for idx, joist_geom in enumerate(self.joist_geoms):
            if joist_geom is None:
                continue
            trib_area = self.joist_trib_areas[idx]
            sub_id = f"{self.id}-{idx}"
            intersections_below = []
            for sup_idx, support_geom in enumerate(self.joist_supports):
                other_tag = self.joist_support_tags[sup_idx]
                intersection_attrs = geom_ops.get_intersection(
                    joist_geom, support_geom, other_tag
                )
                if intersection_attrs is None:
                    continue
                other_reaction_type = self.joist_supports[support_geom]
                intersection_below = Intersection(
                    *intersection_attrs, other_reaction_type=other_reaction_type
                )
                intersections_below.append(intersection_below)
            subelement = Element(
                joist_geom,
                sub_id,
                intersections_below=intersections_below,
                intersections_above=[],
                correspondents_below=[],
                correspondents_above=[],
                plane_id=self.plane_id,
                element_type="collector",
                subelements=None,
                trib_area=trib_area,
                kwargs=self.elem_kwargs,
                # extent_polygon=self.extent_polygon,
            )
            subelements.append(subelement)
        new_element = Element(
            e.geometry,
            tag=e.tag,
            rank=e.rank,
            intersections_above=e.intersections_above,
            intersections_below=e.intersections_below,
            correspondents_above=e.correspondents_above,
            correspondents_below=e.correspondents_below,
            plane_id=e.plane_id,
            element_type=e.element_type,
            subelements=subelements,
            trib_area=e.trib_area,
            reaction_type="linear",
            kwargs=e.kwargs,
            extent_line=e.extent_line,
        )
        return new_element

    def generate_joist_geom(self, index: int):
        """
        Returns i, j coordinates of the joist in the JoistArray at the position
        of 'index'. Raises IndexError if 'index' is not within the joist array
        extents given the spacing.

        'index': joists are numbered from 0 (first joist, at joist extent) and
            go to n, the last joist in the array.
        """
        start_centroid = self.get_extent_edge("start").centroid
        try:
            joist_distance = self.joist_locations[index]
        except IndexError as e:
            raise IndexError(
                f"Joist index {index} is beyond the extent of the joist array for {self.id}. "
                f"Last index is {len(self.joist_locations) - 1} @ {self.joist_locations[-1]}"
            ) from None

        if index != 0 and index != len(self.joist_locations) - 1:
            new_centroid = geom_ops.project_node(
                start_centroid, -self.vector_normal, joist_distance  # orig -ve
            )
            system_bounds = geom_ops.get_system_bounds(
                self._joist_prototype, list(self._supports)
            )
            projection_distance = geom_ops.get_magnitude(system_bounds)
            ray_ai = geom_ops.project_node(
                new_centroid, self.vector_parallel, projection_distance  # orig +ve
            )
            ray_aj = geom_ops.project_node(
                new_centroid, -self.vector_parallel, projection_distance  # orig -ve
            )
            ray_a = LineString([ray_ai, ray_aj])

            ray_bj = geom_ops.project_node(
                new_centroid, self.vector_parallel, projection_distance  # orig +ve
            )
            ray_bi = geom_ops.project_node(
                new_centroid, -self.vector_parallel, projection_distance  # orig +ve
            )
            ray_b = LineString([ray_bi, ray_bj])
            intersecting_supports = [
                support
                for support in self._supports
                if support.intersects(ray_a | ray_b)
            ]
            try:
                sorted_supports = geom_ops.sort_supports(
                    ray_a | ray_b, intersecting_supports
                )
            except AssertionError:
                # This condition is hit when a particular joist geometry cannot be
                # supported upon the supplied supports. It is skipped and a gap in
                # the joist array will be result. The solution is for the designer
                # to adjust their support conditions so the joist will not be
                # skipped.
                warnings.warn(
                    f"No geometry generated for a joist with index={index}. Designer to "
                    "review support conditions and make required adjustments."
                )
                return None
            support_a_loc = ray_a.intersection(sorted_supports[0], grid_size=1e-3)
            support_b_loc = ray_b.intersection(sorted_supports[-1], grid_size=1e-3)
            end_a = support_a_loc
            end_b = support_b_loc
        # These clauses req'd to deal with floating point error possible
        # on the end joists (occurs after performing project_node)
        elif index == 0:
            end_a = support_a_loc = self._extents[0][0]
            end_b = support_b_loc = self._extents[-1][0]
        elif index == len(self.joist_locations) - 1:
            end_a = support_a_loc = self._extents[0][1]
            end_b = support_b_loc = self._extents[-1][1]

        cant_a = self._cantilevers["A"]
        cant_b = self._cantilevers["B"]
        if cant_a and cant_a >= self._cantilever_tolerance:
            end_a = geom_ops.project_node(
                support_a_loc, -self.vector_parallel, self._cantilevers["A"]
            )
        if cant_b and cant_b >= self._cantilever_tolerance:
            end_b = geom_ops.project_node(
                support_b_loc, self.vector_parallel, self._cantilevers["B"]
            )
        joist_geom = set_precision(LineString([end_a, end_b]), grid_size=1e-3)
        # With diagonal supports, it is possible the the new geom sliiiightly misses the supports. WTF. Why?
        # This is a hack to fix that...TODO: BETTER SOLUTION FOR THIS HACK?
        if not all([joist_geom.intersects(support) for support in self._supports]):
            end_a = geom_ops.project_node(
                end_a, -self.vector_parallel, self._cantilever_tolerance / 10
            )
            end_b = geom_ops.project_node(
                end_b, self.vector_parallel, self._cantilever_tolerance / 10
            )
            joist_geom = LineString([end_a, end_b])
        if joist_geom.length <= self._cantilever_tolerance:
            return None
        return joist_geom

    def get_extent_edge(self, edge: str = "start"):
        """
        Gets the "joist" that would exist at the edge of the array

        'edge': one of {'start', 'end'}
        """
        if edge == "start":
            node_i = self._extents[0][0]
            node_j = self._extents[1][0]
        elif edge == "end":
            node_i = self._extents[0][1]
            node_j = self._extents[1][1]
        return LineString([node_i, node_j])

    def get_joist_trib_widths(self, index) -> tuple[float, float]:
        """
        Returns the trib widths of the the joist at 'index'. The trib
        widths are a tuple representing the left and right width,
        respectively.
        """
        if index < 0:
            # Convert -ve index lookup to a +ve index lookup
            index = len(self.joist_locations) + index
        if index == 0:  # The first joist
            spacing_right = self.joist_locations[1] - self.joist_locations[0]
            trib_widths = (0.0, spacing_right / 2.0)
        elif index == len(self.joist_locations) - 1:  # The last joist
            spacing_left = self.joist_locations[-1] - self.joist_locations[-2]
            trib_widths = (spacing_left / 2.0, 0.0)
        else:
            spacing_left = self.joist_locations[index] - self.joist_locations[index - 1]
            spacing_right = (
                self.joist_locations[index + 1] - self.joist_locations[index]
            )
            trib_widths = (spacing_left / 2.0, spacing_right / 2.0)
        return trib_widths

    def generate_trib_area(self, index: int) -> Polygon:
        """
        Returns a tuple of Polygon representing the tributary area of the 'joist' based on the
        given 'trib_widths'
        """
        joist = self.joist_geoms[index]
        if joist is None:
            return None
        trib_widths = self.joist_trib_widths[index]
        i_node, j_node = joist.boundary.geoms  # Point, Point
        trib_left, trib_right = trib_widths  # float, float

        # Left - # TODO: Can I not just buffer the joist? I guess that if the joist is on an
        # angle then extents won't capture the angle.
        if trib_left != 0.0:
            i_left = geom_ops.project_node(i_node, self.vector_normal, trib_left)
            j_left = geom_ops.project_node(j_node, self.vector_normal, trib_left)
            trib_area_left = convex_hull(MultiPoint([i_left, j_left, j_node, i_node]))
        else:
            trib_area_left = Polygon()

        # Right
        if trib_right != 0.0:
            i_right = geom_ops.project_node(i_node, -self.vector_normal, trib_right)
            j_right = geom_ops.project_node(j_node, -self.vector_normal, trib_right)
            trib_area_right = convex_hull(
                MultiPoint([i_right, j_right, j_node, i_node])
            )
        else:
            trib_area_right = Polygon()
        return trib_area_left | trib_area_right

    def show_svg(self, use_ipython_display: bool = True):
        """
        Returns a GeometryCollection containing:
            - Joists
            - Joist Trib Areas
            - Joist Supports

        For manual visual review
        """
        from IPython.display import display

        display(
            GeometryCollection(
                self.joist_geoms + self.joist_trib_areas + self.joist_supports
            )
        )
