from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Any, Optional
import pycba as cba
import numpy as np
from shapely import (
    LineString,
    Point,
    MultiLineString,
    Polygon,
    MultiPoint,
    convex_hull,
    GeometryCollection,
)
import shapely.ops as ops

from papermodels.datatypes.element import Element, Intersection
from papermodels.geometry import geom_ops

from rich import print
from IPython.display import display


def collector_trib_model(
        element: Element, 
        trib_width: float = 1.0,
        reaction_type: str = "linear"
    ):
    """
    Generates a representative trib area for the joist prototype.

    Assumptions:
    - The supports are assumed to be orthogonal
    - The loading is consistent for all joists within the spread
    - The joist represents a one-of-many similar elements within
        the spread and the spread thus represents a linear reaction
        over the supports.
    """
    e = element
    geom = e.geometry
    trib_area = geom.buffer(
        distance=trib_width/2.0,
        cap_style="flat",
    )
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
    )
    return collector_element



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
        cantilever_tolerance: float = 1e-2,
    ):
        self.joist_prototype = LineString(geom_ops.get_start_end_nodes(element.geometry))
        try:
            self.joist_supports = geom_ops.clean_polygon_supports([ib.other_geometry for ib in element.intersections_below], self.joist_prototype)
        except AssertionError:
            raise AssertionError(f"No intersection at cleaned_support: {element.tag=}. Is geometry right on the edge of the support?")
        self.joist_support_tags = [ib.other_tag for ib in element.intersections_below]
        self.id = element.tag
        self.plane_id = element.plane_id
        self.spacing = spacing  # Need to include this in the legend and thus, the Element
        self.initial_offset = float(initial_offset)
        self._joist_prototype = self.joist_prototype
        self._cantilever_tolerance = cantilever_tolerance
        try:
            self._extents = geom_ops.get_joist_extents(self.joist_prototype, self.joist_supports)
        except AssertionError as e:
            raise AssertionError(f"No intersection within joist extents: {element.tag=}")

        self._supports = geom_ops.sort_supports(self.joist_prototype, self.joist_supports)
        self._cantilevers = geom_ops.get_cantilever_segments(self.joist_prototype, self._supports)
        self.vector_parallel = geom_ops.get_direction_vector(self.joist_prototype)
        self.vector_normal = geom_ops.rotate_90_vector(self.vector_parallel, ccw=True)
        self.joist_at_start = float(joist_at_start)
        self.joist_at_end = float(joist_at_end)
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
            self.get_joist_trib_widths(idx)
            for idx, _ in enumerate(self.joist_locations)
        ]
        self.joist_trib_areas = [
            self.generate_trib_area(idx) for idx, _ in enumerate(self.joist_locations)
        ]
    # def __repr__(self):
    #     return class_representation(self)

    @classmethod
    def create_subelements(
        cls,
        element: Element,
        spacing: float,
        initial_offset: float | int = 0.0,
        joist_at_start: bool = True,
        joist_at_end: bool = False,
        cantilever_tolerance: float = 1e-2,
    ) -> JoistArrayModel:
        if element.geometry.geom_type != "LineString":
            return None
        joist_array = cls(
            element, spacing, initial_offset, joist_at_start, joist_at_end, cantilever_tolerance
        )
        # joist_array.show_svg()
        return joist_array.to_subelements()
    

    def to_subelements(self) -> list[Element]:
        """
        Returns the sub-joists in the JoistArray (self) as Element
        """
        subelements = []
        for idx, joist_geom in enumerate(self.joist_geoms):
            trib_area = self.joist_trib_areas[idx]
            sub_id = f"{self.id}-{idx}"
            # other_tag = self.joist_support_tags[idx]
            intersections_below = []
            for sup_idx, support_geom in enumerate(self.joist_supports):
                other_tag = self.joist_support_tags[sup_idx]
                intersection_attrs = geom_ops.get_intersection(joist_geom,  support_geom, other_tag)
                intersection_below = Intersection(
                    *intersection_attrs
                )
                intersections_below.append(intersection_below)
            element = Element(
                joist_geom,
                sub_id,
                intersections_below=intersections_below,
                intersections_above=[],
                correspondents_below=[],
                correspondents_above=[],
                plane_id=self.plane_id,
                element_type="collector",
                subelements=None,
                trib_area = trib_area,
            )
            subelements.append(element)
        return subelements
    

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
                start_centroid, -self.vector_normal, joist_distance # orig -ve
            )

            system_bounds = geom_ops.get_system_bounds(
                self._joist_prototype, list(self._supports)
            )
            projection_distance = geom_ops.get_magnitude(system_bounds)
            ray_aj = geom_ops.project_node(
                new_centroid, -self.vector_parallel, projection_distance # orig -ve
            )
            ray_a = LineString([new_centroid, ray_aj])
            ray_bj = geom_ops.project_node(
                new_centroid, self.vector_parallel, projection_distance # orig +ve
            )
            ray_b = LineString([new_centroid, ray_bj])
            # display(GeometryCollection([start_centroid, self.get_extent_edge("start"), new_centroid, self._supports[0], self._supports[-1],]))
            support_a_loc = ray_a.intersection(self._supports[0])
            support_b_loc = ray_b.intersection(self._supports[-1])

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
        
        if self._cantilevers["A"]:
            end_a = geom_ops.project_node(
                support_a_loc, -self.vector_parallel, self._cantilevers["A"]
            )
        if self._cantilevers["B"]:
            end_b = geom_ops.project_node(
                support_b_loc, self.vector_parallel, self._cantilevers["B"]
            )
        return LineString([end_a, end_b])

    def get_extent_edge(self, edge: str = "start"):
        """
        Gets the "joist" that would exist at the edge of the array

        'edge': one of {'start', 'end'}
        """
        if edge == "start":
            node_i = self._extents[0][0]
            node_j = self._extents[-1][0]
        elif edge == "end":
            node_i = self._extents[0][1]
            node_j = self._extents[-1][1]
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
        trib_widths = self.joist_trib_widths[index]
        i_node, j_node = joist.boundary.geoms  # Point, Point
        trib_left, trib_right = trib_widths  # float, float

        # Left - # TODO: Can I not just buffer the joist? I guess that if the joist is on an 
        # angle then extents won't capture the angle.
        if trib_left != 0.0:
            i_left = geom_ops.project_node(i_node, -self.vector_normal, trib_left)
            j_left = geom_ops.project_node(j_node, -self.vector_normal, trib_left)
            trib_area_left = convex_hull(MultiPoint([i_left, j_left, j_node, i_node]))
        else:
            trib_area_left = Polygon()

        # Right
        if trib_right != 0.0:
            i_right = geom_ops.project_node(i_node, self.vector_normal, trib_right)
            j_right = geom_ops.project_node(j_node, self.vector_normal, trib_right)
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
        return display(GeometryCollection(self.joist_geoms + self.joist_trib_areas + self.joist_supports))
        


# @dataclass
# class Joist:
#     """
#     Models a joist with a uniform load of
#     'w' on all spans of the joist that exist.

#                  w
#     ||||||||||||||||||||||||||||
#     ----------------------------
#         ^                  ^
#         R1                 R2
#     < a ><      span      >< b >
#     """

#     span: float | Any
#     a: float | Any = 0.0
#     b: float | Any = 0.0

#     def __post_init__(self):
#         L = [self.a, self.span, self.b]
#         EI = [1e3, 1e3, 1e3]
#         R = [
#             0.0,
#             0.0,
#             -1.0,
#             0.0,
#             -1.0,
#             0.0,
#             0.0,
#             0.0,
#         ]

#         if self.a == 0:
#             L.pop(0)
#             EI.pop(0)
#             R.pop(0)
#             R.pop(0)
#         if self.b == 0:
#             L.pop()
#             EI.pop()
#             R.pop()
#             R.pop()

#         self._pycba_model = cba.BeamAnalysis(
#             L,
#             EI,
#             R,
#         )
#         for idx, _ in enumerate(L):
#             self._pycba_model.add_udl(idx + 1, 1)  # 1-based idx

#     def get_r1(self):
#         self._pycba_model.analyze()
#         total_r1 = self._pycba_model._beam_results.R[0]
#         total_load = self.get_total_load()
#         return round(total_r1 / total_load, 9)

#     def get_r2(self):
#         self._pycba_model.analyze()
#         total_r2 = self._pycba_model._beam_results.R[1]
#         total_load = self.get_total_load()
#         return round(total_r2 / total_load, 9)

#     def get_total_load(self):
#         w = 1
#         total_load_a = w * self.a
#         total_load_span = w * self.span
#         total_load_b = w * self.b
#         total_load = sum([total_load_a, total_load_span, total_load_b])
#         return total_load


