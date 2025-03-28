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
from papermodels.datatypes.utils import class_representation
from papermodels.geometry import geom_ops

from rich import print
from IPython.display import display


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
        self.joist_prototype = LineString(get_start_end_nodes(element.geometry))
        self.joist_supports = self.clean_polygon_supports([ib.other_geometry for ib in element.intersections_below])

        self.joist_support_tags = [ib.other_tag for ib in element.intersections_below]
        self.id = element.tag
        self.plane_id = element.plane_id
        self.spacing = spacing  # Need to include this in the legend and thus, the Element
        self.initial_offset = float(initial_offset)
        self._joist_prototype = self.joist_prototype
        self._cantilever_tolerance = cantilever_tolerance
        self._extents = get_joist_extents(self.joist_prototype, self.joist_supports)

        self._supports = determine_support_order(self.joist_prototype, self.joist_supports)
        self._cantilevers = get_cantilever_segments(self.joist_prototype, self._supports)
        self.vector_parallel = get_direction_vector(self.joist_prototype)
        self.vector_normal = rotate_90(self.vector_parallel, ccw=True)
        self.joist_at_start = float(joist_at_start)
        self.joist_at_end = float(joist_at_end)
        self.joist_locations = get_joist_locations(
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
        joist_array = cls(
            element, spacing, initial_offset, joist_at_start, joist_at_end, cantilever_tolerance
        )
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
                plane_id=self.plane_id,
                element_type="collector",
                subelements=None,
                trib_area = trib_area,
            )
            subelements.append(element)
        return subelements
    

    def clean_polygon_supports(self, support_geoms: list[LineString | Polygon]):
        """
        Converts any Polygon in support_geoms into LineStrings. The LineStrings
        are created depending on where the joist prototype lands within the polygon.

        Assumption: the Polygon represents a single rectangle which represents a 
        wall or something similar.

        The resulting LineString will either be located on the inside face of the
        rectangle support or along the centerline.

        Generating the centerline assumes that the Polygon is a rectangle. Results
        will be unpredictable for Polygons of other shapes.
        """
        cleaned_supports = []
        for support_geom in support_geoms:
            if support_geom.geom_type == "Polygon":
                support_lines = geom_ops.explode_polygon(support_geom)
                support_intersections = self.joist_prototype.intersects(np.array(support_lines))
                if sum(support_intersections) == 1: # Intersects on one edge only
                    intersecting_line_index = int(support_intersections.nonzero()[0][0])
                    support_line = support_lines[intersecting_line_index]
                    assert support_line.intersects(self.joist_prototype)
                elif sum(support_intersections) == 2:
                    support_line = geom_ops.get_rectangle_centerline(support_geom)
                    assert support_line.intersects(self.joist_prototype)
                cleaned_supports.append(support_line)
            else:
                cleaned_supports.append(support_geom)
        return cleaned_supports

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
            new_centroid = project_node(
                start_centroid, self.vector_normal, joist_distance # orig -ve
            )

            system_bounds = get_system_bounds(
                self._joist_prototype, list(self._supports.values())
            )
            projection_distance = get_magnitude(system_bounds)
            ray_aj = project_node(
                new_centroid, -self.vector_parallel, projection_distance # orig -ve
            )
            ray_a = LineString([new_centroid, ray_aj])
            ray_bj = project_node(
                new_centroid, self.vector_parallel, projection_distance # orig +ve
            )
            ray_b = LineString([new_centroid, ray_bj])
            # display(GeometryCollection([start_centroid, self._supports['A'], ray_a, self._supports['B'], self._extents['A'][0], self._extents['B'][0]]))
            support_a_loc = ray_a.intersection(self._supports["A"])
            support_b_loc = ray_b.intersection(self._supports["B"])

            end_a = support_a_loc
            end_b = support_b_loc

        # These clauses req'd to deal with floating point error possible
        # on the end joists (occurs after performing project_node)
        elif index == 0:
            end_a = support_a_loc = self._extents["A"][0]
            end_b = support_b_loc = self._extents["B"][0]
        elif index == len(self.joist_locations) - 1:
            end_a = support_a_loc = self._extents["A"][1]
            end_b = support_b_loc = self._extents["B"][1]

        if self._cantilevers["A"]:
            end_a = project_node(
                support_a_loc, -self.vector_parallel, self._cantilevers["A"]
            )
        if self._cantilevers["B"]:
            end_b = project_node(
                support_b_loc, self.vector_parallel, self._cantilevers["B"]
            )
        return LineString([end_a, end_b])

    def get_extent_edge(self, edge: str = "start"):
        """
        Gets the "joist" that would exist at the edge of the array

        'edge': one of {'start', 'end'}
        """
        if edge == "start":
            node_i = self._extents["A"][0]
            node_j = self._extents["B"][0]
        elif edge == "end":
            node_i = self._extents["A"][1]
            node_j = self._extents["B"][1]
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
            trib_widths = (0.0, self.joist_locations[1] / 2.0)
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
            i_left = project_node(i_node, self.vector_normal, trib_left)
            j_left = project_node(j_node, self.vector_normal, trib_left)
            trib_area_left = convex_hull(MultiPoint([i_left, j_left, j_node, i_node]))
        else:
            trib_area_left = Polygon()

        # Right
        if trib_right != 0.0:
            i_right = project_node(i_node, -self.vector_normal, trib_right)
            j_right = project_node(j_node, -self.vector_normal, trib_right)
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
        return GeometryCollection(self.joist_geoms + self.joist_trib_areas + self.joist_supports)
        


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


def get_joist_extents(
    joist_prototype: LineString, joist_supports: list[LineString]
) -> dict[str, tuple[Point, Point]]:
    """
    Returns the extents for the supports "A" and "B". Each extent is represented by a tuple of
    Point objects which represent the "i" (start) and "j" (end) locations on the supports
    given in 'joist_supports' which support the 'joist_prototype'.

    'joist_supports' is a list of two LineString where each LineString only has one line segment
        (the relevant line segment which provides the support to 'joist_prototype')
    """
    supports_bbox = get_system_bounds(joist_prototype, joist_supports)
    
    magnitude_max = get_magnitude(supports_bbox)
    joist_vector = get_direction_vector(joist_prototype)
    ordered_supports = determine_support_order(joist_prototype, joist_supports)
    a_support, b_support = ordered_supports["A"], ordered_supports["B"]

    ai_node, aj_node = get_start_end_nodes(a_support)

    ai_to_b_jnode = project_node(ai_node, joist_vector, magnitude_max)
    ai_to_b_ray = LineString([ai_node, ai_to_b_jnode])
    aj_to_b_jnode = project_node(aj_node, joist_vector, magnitude_max)
    aj_to_b_ray = LineString([aj_node, aj_to_b_jnode])

    bi_node, bj_node = get_start_end_nodes(b_support)

    bi_to_a_jnode = project_node(bi_node, -joist_vector, magnitude_max)
    bi_to_a_ray = LineString([bi_node, bi_to_a_jnode])
    bj_to_a_jnode = project_node(bj_node, -joist_vector, magnitude_max)
    bj_to_a_ray = LineString([bj_node, bj_to_a_jnode])

    extents_a = [ai_node, aj_node]
    extents_b = [bi_node, bj_node]

    # Project A onto B
    if ai_to_b_ray.intersects(b_support):
        extents_b[0] = ai_to_b_ray & b_support
    if aj_to_b_ray.intersects(b_support):
        extents_b[1] = aj_to_b_ray & b_support

    # Project B onto A
    if bi_to_a_ray.intersects(a_support):
        extents_a[0] = bi_to_a_ray & a_support
    if bj_to_a_ray.intersects(a_support):
        extents_a[1] = bj_to_a_ray & a_support

    return {"A": tuple(extents_a), "B": tuple(extents_b)}


def get_cantilever_segments(
    joist_prototype: LineString,
    ordered_supports: dict[str, LineString],
    tolerance: float = 1e-1,
) -> dict[str, float]:
    """
    Returns a dictionary containing the cantilever lengths over-hanging supports "A" and
    "B", respectively. Returns a length of 0.0 if the length is less than the tolerance.
    """
    joist_interior = convex_hull(MultiLineString([geom for geom in ordered_supports.values()]))
    cantilevers = ops.split(joist_prototype, joist_interior) - joist_interior
    cantilever_segments = {"A": 0.0, "B": 0.0}
    if isinstance(cantilevers, LineString):
        split_a = cantilevers
        split_b = Point() # A geometry of length 0
    elif hasattr(cantilevers, "geoms"):
        split_a, split_b = cantilevers.geoms
    if split_a.distance(ordered_supports['A']) < split_a.distance(ordered_supports['B']):
        cantilever_segments = {"A": split_a.length, "B": split_b.length}
    else:
        cantilever_segments = {"A": split_b.length, "B": split_a.length}
    return cantilever_segments


def get_system_bounds(
    joist_prototype: LineString, joist_supports: list[LineString]
) -> tuple[float, float, float, float]:
    """
    Returns the minx, miny, maxx, maxy bounding box of all the LineStrings in 'joist_supports',
    taken as a group.
    """
    all_lines = MultiLineString(joist_supports + [joist_prototype])
    bbox = all_lines.bounds
    return bbox


def get_magnitude(bounds: tuple[float, float, float, float]) -> float:
    """
    Returns the distance of the "min" and "max" coordinates described in 'bounds'

    'bounds': represents the minx, miny, maxx, maxy values of the "min" and "max" coordinates
    """
    minx, miny, maxx, maxy = bounds
    delta_y = maxy - miny
    delta_x = maxx - minx
    magnitude = (delta_y**2 + delta_x**2) ** 0.5
    return magnitude


def get_joist_locations(
    start_edge: LineString,
    end_edge: LineString,
    spacing: float,
    initial_offset: float,
    joist_at_start: bool,
) -> list[float]:
    """
    Returns a list of location offsets (starting from 0.0)
    """
    distance = start_edge.distance(end_edge)
    distance_remaining = distance
    joist_locs = []
    if joist_at_start:
        joist_locs.append(0.0)
        if initial_offset:
            joist_locs.append(initial_offset)
            distance_remaining -= initial_offset
    else:
        if initial_offset:
            joist_locs.append(initial_offset)
            distance_remaining -= initial_offset
    while distance_remaining > spacing:
        distance_remaining -= spacing
        joist_locs.append(distance - distance_remaining)
    else:
        joist_locs.append(distance)
    return joist_locs


def get_direction_vector(ls: LineString) -> np.ndarray:
    """
    Returns a numpy array representing the normalized +ve direction vector of the LineString 'ls'.

    'ls': A LineString with two or more points. If there are more than two points, it
        is assumed that all points are co-linear.
    """
    i_node, j_node = get_start_end_nodes(ls)
    column_vector = np.array(j_node.xy) - np.array(i_node.xy)
    column_vector_norm = np.linalg.norm(column_vector)
    parallel_vector =  column_vector / column_vector_norm
    return parallel_vector
    # return column_vector.T[0] # Return a flat, 1D vector


def determine_support_order(
    joist_prototype: LineString, supports: list[LineString]
) -> dict[str, LineString]:
    """
    Returns a dict identifying which support is "A" and which is "B".

    The "A" and "B" supports are arranged so that the vector of the joist spanning
    between them is going to be in the +ve direction (positive X bias). See the
    docstring for get_start_end_nodes for more explanation of the +ve vector direction.
    """
    # TODO: THIS FUNCTION CAN ENABLE HAVING JOISTS WITH MORE THAN TWO SUPPORTS
    # This function should still return the {"A": .., "B": ...} dict but should
    # ignore supports in between A and B. The intermediate supports are re-captured
    # in .to_subelements().

    all_supports = MultiLineString(supports)
    joist_a_node, joist_b_node = order_nodes_positive(
        (joist_prototype & all_supports).geoms
    )
    support_a, support_b = supports
    if joist_a_node.buffer(1e-6).intersects(support_a):
        return {"A": support_a, "B": support_b}
    elif joist_b_node.buffer(1e-6).intersects(support_a):
        return {"A": support_b, "B": support_a}


def get_start_end_nodes(ls: LineString) -> tuple[Point, Point]:
    """
    Returns the "i" and "j" nodes for the coordinates comprising 'ls' in such a way that
    it produces a +ve vector when j_node - i_node is performed. See docstring for
    order_nodes_positive for more information about the +ve vector.

    'ls': A LineString whose points are assumed to be co-linear if there are more than two.
    """
    first_coord = Point(ls.coords[0])
    last_coord = Point(ls.coords[-1])
    return order_nodes_positive([first_coord, last_coord])


def order_nodes_positive(nodes: list[Point]) -> list[Point]:
    """
    Returns the 'i_node' and 'j_node' in the order of "A" and "B" node where "A"
    and "B" node generate a +ve vector when B - A.

    A +ve vector is a vector that has an angle, theta, when measured from horizontal, with
    the following range: -pi / 2 < theta <= pi/2. This can also be thought of as a vector
    with a "positive x bias" because such a vector will never point in the -ve x direction.
    """
    # return sorted(nodes, key=lambda x: x.coords[0][0])
    i_node, j_node = nodes
    ix, iy = i_node.coords[0]
    jx, jy = j_node.coords[0]

    delta_y = jy - iy
    delta_x = jx - ix

    if -math.pi / 2 < math.atan2(delta_y, delta_x) <= math.pi / 2:
        return i_node, j_node
    else:
        return j_node, i_node


def project_node(node: Point, vector: np.ndarray, magnitude: float):
    """
    Returns a Point representing 'node' projected along 'vector' for a distance of
    'magnitude'.

    'node': a point in 2D or 3D space
    'vector': a normalized vector in 2D or 3D space
    'magnitude': the distance along 'vector' that 'node' should be projected
    """
    scaled_vector = vector * magnitude
    projected_node = np.array(node.xy) + scaled_vector
    return Point(projected_node)


def rotate_90(v: np.ndarray, precision: int = 6, ccw=True) -> tuple[float, float]:
    """
    Rotate the vector components, 'x1' and 'y1' by 90 degrees.

    'precision': round result to this many decimal places
    'ccw': if True, rotate counter-clockwise (clockwise, otherwise)
    """
    v_angle = np.arctan2(v[1], v[0])

    if ccw:
        if 0 < v_angle <= math.pi / 2: # Positive x-bias
            angle = -math.pi / 2
        else:
            angle = math.pi/2
    else:
        if 0 < v_angle <= math.pi / 2: # Positive x-bias
            angle = math.pi/2
        else:
            angle = -math.pi / 2

    rot = np.array(
        [
            [round(math.cos(angle), precision), -round(math.sin(angle), precision)],
            [round(math.sin(angle), precision), round(math.cos(angle), precision)],
        ]
    )
    return rot @ v
