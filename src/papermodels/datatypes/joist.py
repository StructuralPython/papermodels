from dataclasses import dataclass
import math
from typing import Any, Optional
import pycba as cba
import numpy as np
from shapely import LineString, Point, MultiLineString


@dataclass
class JoistArray:
    """
    Models a spread of joists over a region where the distance
    between the supports may vary linearly.
    """

    cant_a: float
    cant_b: float
    extents_a: tuple[Point, Point]
    extents_b: tuple[Point, Point]
    spacing: float
    initial_offset: float
    vector_normal: np.ndarray


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

    from IPython.display import display
    from shapely import GeometryCollection

    display(GeometryCollection([a_support, b_support, ai_to_b_ray, aj_to_b_ray]))

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
    return Point(projected_node)
