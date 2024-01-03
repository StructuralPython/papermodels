from dataclasses import dataclass
import math
from typing import Any, Optional
import pycba as cba
import numpy as np
from shapely import LineString, Point, MultiLineString
import shapely.ops as ops

from papermodels.datatypes.utils import class_representation


class JoistArray:
    """
    Models a spread of joists over a region where the distance
    between the supports may vary linearly.
    """

    def __init__(
        self,
        joist_prototype: LineString,
        joist_supports: list[LineString],
        spacing: float | int,
        initial_offset: float | int = 0.0,
        joist_at_start: bool = True,
        joist_at_end: bool = False,
        cantilever_tolerance: float = 1e-2,
    ):
        self.spacing = spacing
        self.initial_offset = initial_offset
        self._joist_prototype = joist_prototype
        self._cantilever_tolerance = cantilever_tolerance
        self._extents = get_joist_extents(joist_prototype, joist_supports)
        self._supports = determine_support_order(joist_prototype, joist_supports)
        self._cantilevers = get_cantilever_segments(joist_prototype, self._supports)
        self.vector_parallel = get_direction_vector(joist_prototype)
        self.vector_normal = rotate_90(self.vector_parallel, ccw=False)
        self.joist_at_start = joist_at_start
        self.joist_at_end = joist_at_end

    # def __repr__(self):
    #     return class_representation(self)

    def generate_joist(self, index: int):
        """
        Returns i, j coordinates of the joist in the JoistArray at the position
        of 'index'. Raises IndexError if 'index' is not within the joist array
        extents given the spacing.

        'index': joists are numbered from 0 (first joist, at joist extent) and
            go to n, the last joist in the array.
        """
        start_centroid = self.get_extent_edge("start").centroid
        joist_distance = self.get_joist_location(index)
        new_centroid = project_node(start_centroid, -self.vector_normal, joist_distance)
        # return new_centroid
        system_bounds = get_system_bounds(
            self._joist_prototype, list(self._supports.values())
        )
        projection_distance = get_magnitude(system_bounds)
        ray_aj = project_node(new_centroid, -self.vector_parallel, projection_distance)
        ray_a = LineString([new_centroid, ray_aj])
        ray_bj = project_node(new_centroid, self.vector_parallel, projection_distance)
        ray_b = LineString([new_centroid, ray_bj])
        support_a_loc = ray_a.intersection(self._supports["A"])
        support_b_loc = ray_b.intersection(self._supports["B"])
        end_a = support_a_loc
        end_b = support_a_loc
        if self._cantilevers["A"]:
            end_a = project_node(
                support_a_loc, -self.vector_parallel, self._cantilevers["A"]
            )
        if self._cantilevers["B"]:
            end_b = project_node(
                support_b_loc, self.vector_parallel, self._cantilevers["B"]
            )
        return LineString([end_a, end_b])

    def index_in_array(self, index: int):
        """
        Returns True if
        """
        pass

    def get_joist_location(self, index: int) -> float:
        """
        Returns the position of the joist centroid projected
        along self.vector_normal
        """
        joist_start = 0.0
        if not self.joist_at_start:
            if not self.initial_offset:
                joist_start = self.spacing
            else:
                joist_start = self.initial_offset

        joist_next = self.spacing
        if self.initial_offset:
            joist_next = self.initial_offset

        start_sequence = [joist_start]
        if index >= 1:
            start_sequence += [joist_next]
        # exclude_start_joist = int(not self.joist_at_start)
        remaining_indexes = max(0, index - (len(start_sequence) - 1))
        print(index, remaining_indexes)
        remaining_sequence = [self.spacing] * remaining_indexes
        total_sequence = start_sequence + remaining_sequence
        return sum(total_sequence)

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


@dataclass
class Joist:
    """
    Models a joist with a uniform load of
    'w' on all spans of the joist that exist.

                 w
    ||||||||||||||||||||||||||||
    ----------------------------
        ^                  ^
        R1                 R2
    < a ><      span      >< b >
    """

    span: float | Any
    a: float | Any = 0.0
    b: float | Any = 0.0

    def __post_init__(self):
        L = [self.a, self.span, self.b]
        EI = [1e3, 1e3, 1e3]
        R = [
            0.0,
            0.0,
            -1.0,
            0.0,
            -1.0,
            0.0,
            0.0,
            0.0,
        ]

        if self.a == 0:
            L.pop(0)
            EI.pop(0)
            R.pop(0)
            R.pop(0)
        if self.b == 0:
            L.pop()
            EI.pop()
            R.pop()
            R.pop()

        self._pycba_model = cba.BeamAnalysis(
            L,
            EI,
            R,
        )
        for idx, _ in enumerate(L):
            self._pycba_model.add_udl(idx + 1, 1)  # 1-based idx

    def get_r1(self):
        self._pycba_model.analyze()
        total_r1 = self._pycba_model._beam_results.R[0]
        total_load = self.get_total_load()
        return round(total_r1 / total_load, 9)

    def get_r2(self):
        self._pycba_model.analyze()
        total_r2 = self._pycba_model._beam_results.R[1]
        total_load = self.get_total_load()
        return round(total_r2 / total_load, 9)

    def get_total_load(self):
        w = 1
        total_load_a = w * self.a
        total_load_span = w * self.span
        total_load_b = w * self.b
        total_load = sum([total_load_a, total_load_span, total_load_b])
        return total_load


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
    splits_a = ops.split(joist_prototype, ordered_supports["A"])
    splits_b = ops.split(joist_prototype, ordered_supports["B"])
    for geom_a in splits_a.geoms:
        for geom_b in splits_b.geoms:
            if not geom_a.intersects(geom_b):
                return {
                    "A": 0.0 if geom_a.length < tolerance else geom_a.length,
                    "B": 0.0 if geom_b.length < tolerance else geom_b.length,
                }


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


def get_direction_vector(ls: LineString) -> np.ndarray:
    """
    Returns a numpy array representing the normalized +ve direction vector of the LineString 'ls'.

    'ls': A LineString with two or more points. If there are more than two points, it
        is assumed that all points are co-linear.
    """
    i_node, j_node = get_start_end_nodes(ls)
    column_vector = np.array(j_node.xy) - np.array(i_node.xy)
    column_vector_norm = np.linalg.norm(column_vector)
    return column_vector / column_vector_norm
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

    all_supports = MultiLineString(supports)
    joist_a_node, joist_b_node = order_nodes_positive(
        *(joist_prototype & all_supports).geoms
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
    return order_nodes_positive(first_coord, last_coord)


def order_nodes_positive(i_node: Point, j_node: Point) -> tuple[Point, Point]:
    """
    Returns the 'i_node' and 'j_node' in the order of "A" and "B" node where "A"
    and "B" node generate a +ve vector when B - A.

    A +ve vector is a vector that has an angle, theta, when measured from horizontal, with
    the following range: -pi / 2 < theta <= pi/2. This can also be thought of as a vector
    with a "positive x bias" because such a vector will never point in the -ve x direction.
    """
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
    # print(projected_node, scaled_vector)
    return Point(projected_node)


def rotate_90(v: np.ndarray, precision: int = 6, ccw=True) -> tuple[float, float]:
    """
    Rotate the vector components, 'x1' and 'y1' by 90 degrees.

    'precision': round result to this many decimal places
    'ccw': if True, rotate counter-clockwise (clockwise, otherwise)
    """
    if ccw:
        angle = math.pi / 2
    else:
        angle = -math.pi / 2

    rot = np.array(
        [
            [round(math.cos(angle), precision), -round(math.sin(angle), precision)],
            [round(math.sin(angle), precision), round(math.cos(angle), precision)],
        ]
    )
    return rot @ v
