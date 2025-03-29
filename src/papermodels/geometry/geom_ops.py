import math
from typing import Optional, Union
import numpy as np

from shapely import (
    Point,
    MultiPoint,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
    convex_hull,
)
import shapely.ops as ops
import shapely.affinity as aff

Geometry = Union[LineString, Polygon]
IntersectingGeometry = Union[Point, LineString]

def get_intersection(
    above: Geometry, below: Geometry, j_tag: str
) -> Optional[tuple[str, IntersectingGeometry, Geometry]]:
    """
    Returns the details of the intersection
    """
    # intersecting_region = above.intersection(below)
    i_type = above.geom_type
    j_type = below.geom_type
    if i_type == "LineString" and j_type == "Polygon":
        intersecting_region = above.intersection(below.exterior)
    elif i_type == "Polygon" and j_type == "LineString":
        intersecting_region = below.intersection(above.exterior)
    else:
        intersecting_region = above.intersection(below)

    all_linestrings = i_type == j_type == "LineString"
    if intersecting_region.geom_type == "Point" and all_linestrings:
        return (intersecting_region, below, j_tag)
    elif intersecting_region.geom_type == "MultiPoint": # Line enters and exits a polygon boundary
        if (
            (i_type == "Polygon" and j_type == "LineString")
            or
            (i_type == "LineString" and j_type == "Polygon")
        ):
            point = intersecting_region.centroid
            assert above.contains(point)
            return (point, below, j_tag)
        else:
            raise ValueError(
                "Could not get intersecting region for MultiPoint. Should not see this error.\n"
                f"{above.wkt=} | {below.wkt=}"
            )
    elif intersecting_region.geom_type == "Point": # LineString and Polygon intersection @ boundary
        return (intersecting_region, below, j_tag)
    else:
        return None


def check_corresponds(above: Union[LineString, Polygon], below: Union[LineString, Polygon]) -> float:
    """
    Returns the ratio of overlap between geometry above and the geometry below.

    A return value of 1.0 represents full correspondence with above and below
    A return value of 0.0 indicates no correspondence with above and below
    A return value in between represents an off-set between the two

    If 'above' and 'below' are Polygons: the ratio represents (above & below).area / below.area.
    If 'above' and 'below' are LineString: the ratio represents (above & below).length / below.length.
    If 'above' is a Polygon and 'below' is a LineString: the ratio represents the 1.0 - distance(above.centroid, below)

    In all cases, a return value of 1.0 represents a "full bearing ratio" (100% of the area of the
    element corresponds with the element on the plane below).  This ratio represents the accuracy
    of the alignment of the sketch from plane to plane and does not necessarily represent the
    bearing area at a connection. For example, a ratio of 1.0 between two polygons representing
    columns may indicate that there is no "slope" in the column and that the bottom of the column has
    been sketched so that it is directly under the top of the column. 
    """
    intersecting_region = above.intersection(below)
    a_type = above.geom_type
    b_type = below.geom_type
    c_type = intersecting_region.geom_type
    if intersecting_region is None:
        return 0.0
    elif a_type == b_type == c_type == 'LineString':
        return intersecting_region.length / below.length
    elif a_type == b_type == c_type == "Polygon" :
        return intersecting_region.area / below.area
    else:
        return 0.0
    

def get_local_intersection_ordinates(start_node: Point, intersections: list[Point]) -> list[float]:
    """
    Returns the relative distances of the Points in 'intersections' relative to the 'start_node'.
    """
    return [
        start_node.distance(intersection) for intersection in intersections
    ]


def get_linestring_start_node(ls: LineString) -> Point:
    """
    Returns a Point representing the starting node of the 'ls' LineString
    when the nodes are ordered with a +ve X bias.
    """
    coords_a, coords_b = ls.coords
    ordered_coords = order_nodes_positive(Point(coords_a), Point(coords_b))
    start_coord = ordered_coords[0]
    return start_coord

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
    supports = MultiLineString([ordered_supports["A"], ordered_supports["B"]])
    cantilever_segments = {"A": 0.0, "B": 0.0}
    for geom_a in splits_a.geoms:
        if isinstance(geom_a & supports, Point):
            cantilever_segments["A"] = (
                0.0 if geom_a.length < tolerance else geom_a.length
            )

    for geom_b in splits_b.geoms:
        if isinstance(geom_b & supports, Point):
            cantilever_segments["B"] = (
                0.0 if geom_a.length < tolerance else geom_a.length
            )
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


def rotate_to_horizontal(line: LineString, geoms: list[Geometry]):
    """
    Rotate the line so that it is horizonatla. Bring the geomswiith it
    """
    i_end, j_end = get_start_end_nodes(line)
    ix, iy = i_end.coords[0]
    jx, jy = j_end.coords[0]

    delta_y = jy - iy
    delta_x = jx - ix

    angle = math.atan2(delta_y, delta_x)

    rotated_line = aff.translate(aff.rotate(line, -angle, origin=i_end, use_radians=True), xoff=-ix)
    rotated_geoms = [aff.translate(aff.rotate(geom, -angle, origin=i_end, use_radians=True), xoff=-ix) for geom in geoms]

    return rotated_line, rotated_geoms

def explode_polygon(p: Polygon) -> list[LineString]:
    """
    Explodes the exterior of the polygon in to a list of individual line segments
    """
    ext_ls = LineString(p.exterior)
    exploded = [LineString(tup) for tup in zip(ext_ls.coords, ext_ls.coords[1:])]
    return exploded

def get_rectangle_centerline(p: Polygon) -> LineString:
    """
    Returns the centerline of the Polygon 'p' assuming that 'p' represents
    a regular rectangle with a long dimension and a short dimension
    """
    rectangle_edges = explode_polygon(p)
    sorted_edges = sorted(rectangle_edges, key=lambda x: x.length)
    short_edges = sorted_edges[:2]
    edge1, edge2 = short_edges
    center_line = LineString([edge1.centroid, edge2.centroid])
    return center_line
    

def calculate_trapezoid_area_sums(
    member_loads: list[list[list[tuple]]]
) -> list[float]:
    """
    Returns a list of the sums of the areas of the trapezoids
    in 'traps'
    """
    member_polys = []
    for polygon_load in member_loads:
        polygon_loads = []
        for inner_pair in polygon_load:
            start, end = inner_pair
            start_x, start_y = start
            end_x, end_y = end
            h = end_x - start_x
            b2 = start_y
            b1 = end_y
            trap_area = trapezoid_area(h, b2, b1)
            polygon_loads.append(trap_area)
        member_polys.append(sum(polygon_loads))
    return member_polys
        


def trapezoid_area(h: float, b2: float, b1: float) -> float:
    """
    Returns the area of the trapezoid.
    """
    area = (b1 + b2) / 2 * h
    return area


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


def create_linestring(points: list[tuple]) -> LineString:
    return LineString(points)


def create_multipoint(points: list[tuple]) -> MultiPoint:
    return MultiPoint(points)


def create_polygon(points: list[tuple]) -> Polygon:
    return Polygon(points)


def create_multipolygon(polygons: list[Polygon]) -> MultiPolygon:
    return MultiPolygon(polygons)


def create_convex_hull(points: list[Point]) -> Polygon:
    return convex_hull(points)
