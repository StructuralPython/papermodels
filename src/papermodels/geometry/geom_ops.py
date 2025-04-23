from decimal import Decimal
import math
from typing import Optional, Union
import numpy as np
from numpy.typing import ArrayLike
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
from papermodels.datatypes.exceptions import GeometryError

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
        if intersecting_region.is_empty:
            intersecting_region = below.intersection(above)
    else:
        intersecting_region = above.intersection(below)
    if intersecting_region.is_empty:
        return None
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
            return (point, below, j_tag)
        else:
            raise ValueError(
                "Could not get intersecting region for MultiPoint. Should not see this error.\n"
                f"{above.wkt=} | {below.wkt=}"
            )
    elif intersecting_region.geom_type == "LineString":
        return (intersecting_region, below, j_tag)
    elif intersecting_region.geom_type == "Point": # LineString and Polygon intersection @ boundary
        return (intersecting_region, below, j_tag)
    elif intersecting_region.geom_type == "Polygon": # Polygon point/line load intersecting with another polygon
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
    ordered_coords = order_nodes_positive([Point(coords_a), Point(coords_b)])
    start_coord = ordered_coords[0]
    return start_coord



def clean_polygon_supports(support_geoms: list[LineString | Polygon], joist_prototype: LineString):
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
            support_lines = explode_polygon(support_geom)
            support_intersections = joist_prototype.intersects(np.array(support_lines))
            if sum(support_intersections) == 1: # Intersects on one edge only
                intersecting_line_index = int(support_intersections.nonzero()[0][0])
                support_line = support_lines[intersecting_line_index]
                # Ensure there are no missing intersections on the support line
                assert support_line.intersects(joist_prototype)
            # elif sum(support_intersections) == 2:
            elif sum(support_intersections) == 0:
                assert support_geom.intersects(support_lines)
                raise GeometryError(f"The geometry {support_geom.wkt} does not intersect {joist_prototype.wkt}")
            elif sum(support_intersections) == 2:
                # Ensure there are no missing intersections on the support line
                # Can sometimes be caused by a joist intersecting with a column
                # (joists should not be "supported" by columns)
                support_line = get_rectangle_centerline(support_geom)
                assert support_line.intersects(joist_prototype)
            cleaned_supports.append(support_line)
        else:
            cleaned_supports.append(support_geom)
    return cleaned_supports


def get_joist_extents(
    joist_prototype: LineString, joist_supports: list[LineString], eps: float = 1e-6
) -> dict[str, tuple[Point, Point]]:
    """
    Returns the extents for the supports "A" and "B". Each extent is represented by a tuple of
    Point objects which represent the "i" (start) and "j" (end) locations on the supports
    given in 'joist_supports' which support the 'joist_prototype'.

    'joist_supports' is a list of LineString where each LineString only has one line segment
        (the relevant line segment which provides the support to 'joist_prototype')
    'eps' is a small tolerance amount to deal with floating point error in the extent
        calcualtion.
    """
    supports_bbox = get_system_bounds(joist_prototype, joist_supports)
    magnitude_max = get_magnitude(supports_bbox)
    joist_vector = get_direction_vector(joist_prototype).flatten()
    joist_origin, joist_end = get_start_end_nodes(joist_prototype)
    joist_origin = np.array(joist_origin.coords[0])
    left_coords = []
    right_coords = []
    for joist_support in joist_supports:
        start_coord, end_coord = joist_support.coords
        start_coord_rotation = cross_product_2d(joist_vector, np.array(start_coord) - joist_origin)
        end_coord_rotation = cross_product_2d(joist_vector, np.array(end_coord) - joist_origin)
        if 0.0 < start_coord_rotation:
            left_coords.append(start_coord)
        elif start_coord_rotation < 0.0:
            right_coords.append(start_coord)
        if 0.0 < end_coord_rotation:
            left_coords.append(end_coord)
        elif end_coord_rotation < 0.0:
            right_coords.append(end_coord)
    closest_left_coord = min(left_coords, key=lambda x: Point(x).distance(joist_prototype))
    closest_right_coord = min(right_coords, key=lambda x: Point(x).distance(joist_prototype))
    closest_left_distance = Point(closest_left_coord).distance(joist_prototype)
    closest_right_distance = Point(closest_right_coord).distance(joist_prototype)
    joist_vector_normal = rotate_90_vector(joist_vector, ccw=True)
    joist_left = LineString([
        project_node(Point(joist_origin), joist_vector_normal, magnitude=closest_left_distance - eps),
        project_node(Point(joist_end), joist_vector_normal, magnitude=closest_left_distance - eps)
    ])
    joist_right = LineString([
        project_node(Point(joist_origin), -joist_vector_normal, magnitude=closest_right_distance - eps),
        project_node(Point(joist_end), -joist_vector_normal, magnitude=closest_right_distance - eps)
    ])
    ordered_joist_supports = sort_supports(joist_prototype, joist_supports)
    extents = []
    # from IPython.display import display
    # from shapely import GeometryCollection
    for support_linestring in ordered_joist_supports:
        left_extent = support_linestring.intersection(joist_left)
        right_extent = support_linestring.intersection(joist_right)
        # display(GeometryCollection([joist_left, support_linestring]))

        # Make sure the intersection geometries are not empty before proceeding
        # If one or more is empty, there is a problem that needs investigating
        assert not left_extent.is_empty
        assert not right_extent.is_empty
        extents.append((left_extent, right_extent))
    return extents


def get_cantilever_segments(
    joist_prototype: LineString,
    ordered_supports: list[LineString],
    tolerance: float = 1e-1,
) -> dict[str, float]:
    """
    Returns a dictionary containing the cantilever lengths over-hanging supports "A" and
    "B", respectively. Returns a length of 0.0 if the length is less than the tolerance.
    """
    joist_interior = convex_hull(MultiLineString([geom for geom in ordered_supports]))
    cantilevers = ops.split(joist_prototype, joist_interior) - joist_interior
    cantilever_segments = {"A": 0.0, "B": 0.0}
    if isinstance(cantilevers, LineString):
        split_a = cantilevers
        split_b = Point() # A geometry of length 0
    elif hasattr(cantilevers, "geoms"):
        split_a, split_b = cantilevers.geoms
    if split_a.distance(ordered_supports[0]) < split_a.distance(ordered_supports[-1]):
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
    column_vector = np.array(j_node.coords[0]) - np.array(i_node.coords[0])
    column_vector_norm = np.linalg.norm(column_vector)
    parallel_vector =  column_vector / column_vector_norm
    return parallel_vector
    # return column_vector.T[0] # Return a flat, 1D vector


def sort_supports(
    joist_prototype: LineString, supports: list[LineString]
) -> dict[str, LineString]:
    """
    Returns a list of the supports arranged so that the vector of the joist spanning
    between them is going to be in the +ve direction (positive X bias). See the
    docstring for get_start_end_nodes for more explanation of the +ve vector direction.
    """
    all_supports = MultiLineString(supports)
    ordered_intersections = order_nodes_positive(
        (joist_prototype & all_supports).geoms
    )
    ordered_supports = []
    for point in ordered_intersections:
        for linestring in supports:
            if linestring.buffer(1e-6).intersects(point):
                ordered_supports.append(linestring)
    return ordered_supports


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


def order_nodes_positive(points: list[Point]) -> tuple[Point]:
    """
    Returns the 'i_node' and 'j_node' in the order of "A" and "B" node where "A"
    and "B" node generate a +ve vector when B - A.

    A +ve vector is a vector that has an angle, theta, when measured from horizontal, with
    the following range: -pi / 2 < theta <= pi/2. This can also be thought of as a vector
    with a "positive x bias" because such a vector will never point in the -ve x direction.
    """
    return tuple(sorted(points, key=lambda x: x.coords[0]))


def project_node(node: Point, vector: np.ndarray, magnitude: float):
    """
    Returns a Point representing 'node' projected along 'vector' for a distance of
    'magnitude'.

    'node': a point in 2D or 3D space
    'vector': a normalized vector in 2D or 3D space
    'magnitude': the distance along 'vector' that 'node' should be projected
    """
    scaled_vector = vector * magnitude
    projected_node = np.array(node.coords[0]) + scaled_vector
    return Point(projected_node)


def scale_vertices(
    vertices: list[Decimal],
    scale: Decimal,
    paper_origin: Optional[tuple[Decimal, Decimal]] = None,
    round_precision: int = 4
) -> tuple[Decimal | float]:
    """
    Scale the vertices in relation to the origin or in relation to 'paper_origin'.

    If 'paper_origin' is provided, then the annotation coordinates will have their origin reset
    to 'paper_origin'. Note that 'paper_origin' is the unscaled coordinate space (i.e. in points)
    """
    if paper_origin is not None:
        offset_x = paper_origin[0]
        offset_y = paper_origin[1]
        vertices = _translate_vertices(vertices, offset_x, offset_y)

    scaled_vertices = [round(vertex * scale, round_precision) for vertex in vertices]
    return tuple(scaled_vertices)


def _translate_vertices(
    vertices: list[Decimal], offset_x: float, offset_y: float
) -> tuple[Decimal]:
    """
    Returns a list of float representing 'verticies' translated by 'offset_x' and 'offset_y'.
    """
    vertices_floats = [float(vertex) for vertex in vertices]
    coord_array = np.array(_group_vertices(vertices_floats))
    offset_array = np.array([offset_x, offset_y])
    translated_array = coord_array + offset_array
    flattened_array = flatten_vertex_array(translated_array)
    return flattened_array



def _group_vertices(vertices: list[Decimal | float], close=False) -> list[tuple[Decimal, Decimal]]:
    """
    Returns a list of (x, y) tuples from a list of vertices in the format of:
    'x1 y1 x2 y2 x3 y3 ... xn yn'
    """
    grouped_vertices = []
    coordinates = []
    for idx, ordinate in enumerate(vertices):
        if idx % 2:
            coordinates.append(ordinate)
            grouped_vertices.append(coordinates)
            coordinates = []
        else:
            coordinates.append(ordinate)
    if close:
        grouped_vertices.append(grouped_vertices[0])

    return grouped_vertices


def vertices_to_array(vertices: list[Decimal]) -> ArrayLike:
    """
    Returns a numpy array representing 'vertices' but reshaped to (n, 2)
    """
    return np.array(_group_vertices(vertices), dtype=float)


def flatten_vertex_array(v: ArrayLike, precision=6) -> tuple[Decimal]:
    """
    Returns a flattened version of 'v' in the format of 
    (x1, y1, x2, y2, x3, y3, ..., xn, yn) where 'v' is either
    a row or column-based vector of shape (2, n) or (n, 2) rounded to
    'precision'.

    The returned array is in the format of a tuple of Decimal objects
    for use in pdf annotations
    """
    return tuple([round(Decimal(x), precision) for x in v.flatten()])



def _group_vertices_str(vertices: str, close=False) -> str:
    """
    Returns a list of (x, y) tuples from a list of vertices in the format of:
    'x1 y1 x2 y2 x3 y3 ... xn yn'
    """
    acc = []
    coordinates = []
    for idx, ordinate in enumerate(vertices):
        if idx % 2:
            coordinates.append(f"{ordinate}")
            acc.append(" ".join(coordinates))
            coordinates = []
        else:
            coordinates.append(f"{ordinate}")
    if close:
        acc.append(acc[0])
    return ", ".join(acc)

def rotate_90_vector(v: ArrayLike, precision: int = 6, ccw=True) -> tuple[float, float]:
    """
    Rotate the vector components, 'x1' and 'y1' by 90 degrees.

    'precision': round result to this many decimal places
    'ccw': if True, rotate counter-clockwise (clockwise, otherwise)
    """
    # v_angle = np.arctan2(v[1], v[0])
    if ccw:
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


def rotate_90_coords(v: ArrayLike, precision: int = 6, ccw=True) -> tuple[float, float]:
    """
    Rotate the vector components, 'x1' and 'y1' by 90 degrees.

    'precision': round result to this many decimal places
    'ccw': if True, rotate counter-clockwise (clockwise, otherwise)
    """
    # v_angle = np.arctan2(v[1], v[0])
    if ccw:
            angle = math.pi/2
    else:
            angle = -math.pi / 2
    rot = np.array(
        [
            [round(math.cos(angle), precision), -round(math.sin(angle), precision)],
            [round(math.sin(angle), precision), round(math.cos(angle), precision)],
        ]
    )
    return v @ rot


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
    a regular rectangle with a long dimension and a short dimension.
    The LineString is created with a +ve X-bias.
    """
    rectangle_edges = explode_polygon(p)
    sorted_edges = sorted(rectangle_edges, key=lambda x: x.length)
    short_edges = sorted_edges[:2]
    edge1, edge2 = short_edges
    start, end = order_nodes_positive([edge1.centroid, edge2.centroid])
    center_line = LineString([start, end])
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


def get_vector_angle(v1, v2) -> float:
    """
    Returns the angle between two vectors
    """
    num = np.dot(v1, v2)
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(num / denom)
    return angle

def cross_product_2d(v1, v2):
    return v1[0] * v2[1] - v2[0] * v1[1]

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
