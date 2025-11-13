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
    box,
    convex_hull,
    intersects,
    union,
    GeometryCollection,
    set_precision,
    intersection_all,
)
import shapely.ops as ops
import shapely.affinity as aff
from papermodels.datatypes.exceptions import GeometryError
import load_distribution as ld

Geometry = Union[LineString, Polygon]
IntersectingGeometry = Union[Point, LineString]


def get_intersection(
    above: Geometry,
    below: Geometry,
    below_tag: str,
    above_extent_polygon: Optional[Polygon] = None,
) -> Optional[tuple[IntersectingGeometry, Geometry, str, Optional[Geometry]]]:
    """
    Returns the details of the intersection
    """
    # intersecting_region = above.intersection(below)
    i_type = above.geom_type
    j_type = below.geom_type
    i_extent = above_extent_polygon
    overlap_region = None  # Overlap region is like the "raw" intersecting region
    if i_extent and j_type == "Polygon":
        # Goal: calculate the intersecting region as being along the centerline
        # of the linear polygon support so that, down the line, it becomes easy
        # to calculate the extents from the intersecting region
        intersecting_region = i_extent.intersection(below, grid_size=1e-3)
        overlap_region = intersecting_region
        if not intersecting_region.is_empty:
            inter_centerline = get_rectangle_centerline(intersecting_region)
            support_centerline = get_rectangle_centerline(below)
            # Project the intersecting region centerline onto the support centerline
            _, projected_a = ops.nearest_points(
                Point(inter_centerline.coords[0]), support_centerline
            )
            _, projected_b = ops.nearest_points(
                Point(inter_centerline.coords[1]), support_centerline
            )
            intersecting_region = LineString([projected_a, projected_b])
    elif i_extent and j_type == "LineString":
        intersecting_region = i_extent.intersection(below, grid_size=1e-3)
        overlap_region = intersecting_region
    elif i_type == "LineString" and j_type == "Polygon":
        # intersecting_region = above.intersection(below.exterior)
        start, end = Point(above.coords[0]), Point(above.coords[1])
        if below.contains(start):
            intersecting_region = start
        elif below.contains(end):
            intersecting_region = end
        else:
            intersecting_region = above.intersection(below.exterior, grid_size=1e-3)
        overlap_region = above.intersection(below, grid_size=1e-3)
    elif i_type == "Polygon" and j_type == "LineString":
        start, end = Point(below.coords[0]), Point(below.coords[1])
        if above.contains(start):
            intersecting_region = start
        elif above.contains(end):
            intersecting_region = end
        else:
            intersecting_region = below.intersection(above.exterior, grid_size=1e-3)
        overlap_region = below.intersection(above, grid_size=1e-3)
        if intersecting_region.is_empty:
            intersecting_region = below.intersection(above, grid_size=1e-3)
    else:
        intersecting_region = above.intersection(below, grid_size=1e-3)
        if intersecting_region.length != 0.0:
            overlap_region = intersecting_region  # We do not want a point overlap
    if intersecting_region.is_empty:
        return None
    all_linestrings = i_type == j_type == "LineString"
    if intersecting_region.geom_type == "Point" and all_linestrings:
        return (intersecting_region, below, below_tag, overlap_region)
    elif (
        intersecting_region.geom_type == "MultiPoint"
    ):  # Line enters and exits a polygon boundary
        if (i_type == "Polygon" and j_type == "LineString") or (
            i_type == "LineString" and j_type == "Polygon"
        ):
            point = intersecting_region.centroid
            return (point, below, below_tag, overlap_region)
        else:
            raise ValueError(
                "Could not get intersecting region for MultiPoint. Should not see this error.\n"
                f"{above.wkt=} | {below.wkt=}"
            )
    elif intersecting_region.geom_type == "LineString":
        return (intersecting_region, below, below_tag, overlap_region)
    elif (
        intersecting_region.geom_type == "Point"
    ):  # LineString and Polygon intersection @ boundary
        return (intersecting_region, below, below_tag, overlap_region)
    elif (
        intersecting_region.geom_type == "Polygon"
    ):  # Polygon point/line load intersecting with another polygon
        return (intersecting_region, below, below_tag, overlap_region)
    else:
        return None


def check_corresponds(
    above: Union[LineString, Polygon], below: Union[LineString, Polygon]
) -> float:
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
    elif a_type == b_type == c_type == "LineString":
        return intersecting_region.length / below.length
    elif a_type == b_type == c_type == "Polygon":
        return intersecting_region.area / below.area
    else:
        return 0.0


def get_local_intersection_ordinates(
    start_node: Point, intersections: list[Point]
) -> list[float]:
    """
    Returns the relative distances of the Points in 'intersections' relative to the 'start_node'.
    """
    return [start_node.distance(intersection) for intersection in intersections]


def get_linestring_start_node(ls: LineString) -> Point:
    """
    Returns a Point representing the starting node of the 'ls' LineString
    when the nodes are ordered with a +ve X bias.
    """
    coords_a, coords_b = ls.coords
    ordered_coords = order_nodes_positive([Point(coords_a), Point(coords_b)])
    start_coord = ordered_coords[0]
    return start_coord


def clean_polygon_supports(
    support_geoms: list[LineString | Polygon],
    joist_prototype: LineString,
    extent_polygon: Optional[Polygon] = None,
):
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
            if sum(support_intersections) == 1:  # Intersects on one edge only
                intersecting_line_index = int(support_intersections.nonzero()[0][0])
                support_line = support_lines[intersecting_line_index]
                center_line = get_rectangle_centerline(support_geom)

                # Ensure the center_line is not parallel with the joist
                j_vec = get_direction_vector(joist_prototype)
                s_vec = get_direction_vector(center_line)
                if np.dot(j_vec, s_vec) == 1.0:
                    center_line = get_rectangle_centerline(
                        support_geom, on_long_edge=True
                    )

                # If the joist intersects, then that is the support we will use
                if joist_prototype.intersects(center_line):
                    support_line = center_line
            elif sum(support_intersections) == 0 and not extent_polygon:
                # assert support_geom.intersects(support_lines).any()

                raise GeometryError(
                    f"The geometry {support_geom.wkt} does not intersect {joist_prototype.wkt}"
                )
            elif sum(support_intersections) == 0:
                support_line = get_rectangle_centerline(support_geom)
            elif sum(support_intersections) == 2:
                # Ensure there are no missing intersections on the support line
                # Can sometimes be caused by a joist intersecting with a column
                # (joists should not be "supported" by columns)
                support_line = get_rectangle_centerline(support_geom)
                if not support_line.intersects(joist_prototype):
                    support_line = get_rectangle_centerline(
                        support_geom, on_long_edge=True
                    )
                assert support_line.intersects(joist_prototype)
            cleaned_supports.append(support_line)
        else:
            cleaned_supports.append(support_geom)
    return cleaned_supports


def get_projected_support_centroid(
    frame_element_geometry: LineString,
    polygon_point_support: Polygon,
) -> Point:
    """
    Returns a point that represents the centroid of the polygon point support
    projected onto the vector of the frame_element_geometry. The purpose of this function
    is to effectively "snap" the frame_element_geometry to the centroids of the its polygon
    point supports.

    This effect is desireable when papermodels is used to create "design spans" where
    the frame elements are intended to span from center-of-support to center-of-support.
    """
    fg = frame_element_geometry
    magnitude_max = 3 * fg.length
    direction_vector = get_direction_vector(fg).flatten()
    orig_joist_origin, orig_joist_end = get_start_end_nodes(fg)
    # rfg -> revised_frame_geometry
    rfg_start = project_node(orig_joist_origin, -direction_vector, magnitude_max / 2)
    rfg_end = project_node(orig_joist_end, direction_vector, magnitude_max / 2)
    rfg = LineString([rfg_start, rfg_end])
    centroid = polygon_point_support.centroid
    distance = rfg.project(centroid)
    projected_centroid = rfg.interpolate(distance)
    if not polygon_point_support.contains(projected_centroid):
        raise GeometryError("Projected support centroid is outside of the polygon.")
    return projected_centroid


def get_projected_support_centerline(
    frame_element_geometry: LineString,
    polygon_linear_support: Polygon,
) -> Point:
    """
    Returns a point that represents the either the direct intersection of the frame_element_geometry
    with the centerline or, if no intersection is present, the projection of the nearest end coordinate
    to the centerline onto the center line. The purpose of this function
    is to effectively "snap" the frame_element_geometry to the centroids of the its polygon
    linear supports.

    This effect is desireable when papermodels is used to create "design spans" where
    the frame elements are intended to span from center-of-support to center-of-support.
    """
    fg = frame_element_geometry
    centerline = get_rectangle_centerline(polygon_linear_support)
    if fg.intersects(centerline):
        intersection_point = fg.intersection(centerline)
    else:
        intersection_point, fg_start = ops.nearest_points(centerline, fg)
    if not centerline.intersects(intersection_point):
        raise GeometryError(
            f"Projected support centroid is outside of the polygon: {centerline=} | {intersection_point=}."
        )
    return intersection_point


def get_joist_extents(
    joist_prototype: LineString,
    joist_supports: list[LineString],
    trib_area: Optional[Polygon] = None,
    extent_polygon: Optional[Polygon] = None,
    eps: float = 1e-6,
) -> list[tuple[Point, Point]]:
    """
    Returns the extents for the supports "A" and "B". Each extent is represented by a tuple of
    Point objects which represent the "i" (start) and "j" (end) locations on the supports
    given in 'joist_supports' which support the 'joist_prototype'.

    'joist_supports' is a list of LineString where each LineString only has one line segment
        (the relevant line segment which provides the support to 'joist_prototype')
    'trib_area' if passed, the intersection of the trib area and the support geoms
        will be used to determine the extent locations.
    'extent_polygon', an optional parameter which changes the behaviour of this function
        for situations when the extents of the total extent polygon are required, regardless
        of the support locations within it.
    'eps' is a small tolerance amount to deal with floating point error in the extent
        calcualtion.
    """
    if extent_polygon is not None:
        supports_bbox = get_system_bounds(
            joist_prototype, joist_supports, normal=False, extent_polygon=extent_polygon
        )
        minx, miny, maxx, maxy = supports_bbox
        joist_vector = np.abs(get_direction_vector(joist_prototype))
        # This is one of the places where orthogonality is assumed
        joist_orientation = None
        if joist_vector[0] > joist_vector[1]:
            joist_orientation = "horizontal"
        elif joist_vector[1] > joist_vector[0]:
            joist_orientation = "vertical"
        else:
            print(f"JOIST ORIENTATION VERIANT: {joist_prototype=}")
        if joist_orientation == "horizontal":
            extents = [
                (Point(minx, maxy), Point(minx, miny)),
                (Point(maxx, maxy), Point(maxx, miny)),
            ]
        if joist_orientation == "vertical":
            extents = [
                (Point(minx, miny), Point(maxx, miny)),
                (Point(minx, maxy), Point(maxx, maxy)),
            ]
        return extents

    if trib_area is not None:
        supports_bbox = trib_area.bounds
    else:
        supports_bbox = get_system_bounds(
            joist_prototype, joist_supports, normal=True, extent_polygon=extent_polygon
        )

    magnitude_max = get_magnitude(supports_bbox)
    joist_vector = get_direction_vector(joist_prototype).flatten()
    orig_joist_origin, orig_joist_end = get_start_end_nodes(joist_prototype)
    joist_origin = project_node(orig_joist_origin, -joist_vector, magnitude_max / 2)
    joist_end = project_node(orig_joist_end, joist_vector, magnitude_max / 2)
    extended_joist_prototype = LineString([joist_origin, joist_end])
    joist_origin = np.array(joist_origin.coords[0])
    orig_joist_origin = np.array(orig_joist_origin.coords[0])

    left_coords = []
    right_coords = []
    support_intersection = intersection_all(joist_supports)
    for joist_support in joist_supports:
        joist_support = joist_support.intersection(box(*supports_bbox), grid_size=1e-3)

        start_coord, end_coord = joist_support.coords
        start_coord, end_coord = Point(start_coord), Point(end_coord)

        start_coord_rotation = cross_product_2d(
            joist_vector, np.array(start_coord.coords[0]) - orig_joist_origin
        )
        end_coord_rotation = cross_product_2d(
            joist_vector, np.array(end_coord.coords[0]) - orig_joist_origin
        )

        if start_coord_rotation == 0.0:
            if end_coord_rotation > 0.0:
                right_coords.append(start_coord)
                left_coords.append(end_coord)
            else:
                left_coords.append(start_coord)
                right_coords.append(end_coord)

        elif end_coord_rotation == 0.0:
            if start_coord_rotation > 0.0:
                left_coords.append(start_coord)
                right_coords.append(end_coord)
            else:
                right_coords.append(start_coord)
                left_coords.append(end_coord)
        elif 0.0 < start_coord_rotation:
            left_coords.append(start_coord)
            right_coords.append(end_coord)
        else:
            left_coords.append(end_coord)
            right_coords.append(start_coord)

    left_coords.append(support_intersection)
    right_coords.append(support_intersection)

    closest_left_coord = min_with_none(
        left_coords, key=lambda x: x.distance(extended_joist_prototype)
    )
    closest_right_coord = min_with_none(
        right_coords, key=lambda x: x.distance(extended_joist_prototype)
    )
    closest_left_distance = set_precision(closest_left_coord, grid_size=1e-3).distance(
        extended_joist_prototype
    )
    closest_right_distance = set_precision(
        closest_right_coord, grid_size=1e-3
    ).distance(extended_joist_prototype)
    joist_vector_normal = rotate_90_vector(joist_vector, ccw=True)
    joist_left = set_precision(
        LineString(
            [
                project_node(
                    Point(joist_origin),
                    joist_vector_normal,
                    magnitude=closest_left_distance,
                ),
                project_node(
                    Point(joist_end),
                    joist_vector_normal,
                    magnitude=closest_left_distance,
                ),
            ]
        ),
        grid_size=1e-3,
    )

    joist_right = set_precision(
        LineString(
            [
                project_node(
                    Point(joist_origin),
                    -joist_vector_normal,
                    magnitude=closest_right_distance,
                ),
                project_node(
                    Point(joist_end),
                    -joist_vector_normal,
                    magnitude=closest_right_distance,
                ),
            ]
        ),
        grid_size=1e-3,
    )
    ordered_joist_supports = sort_supports(joist_prototype, joist_supports)
    extents = []
    import shapely.ops as ops

    for support_linestring in ordered_joist_supports:
        support_linestring = set_precision(support_linestring, grid_size=1e-3)

        left_extent = support_linestring.intersection(joist_left)
        right_extent = support_linestring.intersection(joist_right)
        # Make sure the intersection geometries are not empty before proceeding
        # If one or more is empty, there is a problem that needs investigating
        assert not left_extent.is_empty
        assert not right_extent.is_empty
        # Do not "sort_nodes_positive" on these nodes. They are in the correct order.
        extents.append((left_extent, right_extent))
    return extents


def intersecting_support_extents(
    joist_prototype: LineString, supports: list[LineString]
) -> list[tuple[LineString, LineString]]:
    """
    Returns a joist extent taking into account the self-intersection of the supports.
    """
    support_intersections = intersection_all(supports)
    joist_and_supports = supports.copy()
    joist_and_supports = [joist_prototype] + joist_and_supports
    joist_intersections = []
    for support in supports:
        joist_intersections.append(joist_prototype.intersection(support))
    return joist_intersections


def get_cantilever_segments(
    joist_prototype: LineString,
    ordered_supports: list[LineString],
    rel_tol: float = 5e-2,
    abs_tol: Optional[float] = None,
) -> dict[str, float]:
    """
    Returns a dictionary containing the cantilever lengths over-hanging supports "A" and
    "B", respectively. Returns a length of 0.0 if the length is less than the tolerance.

    If 'abs_tol' is given, then 'rel_tol' is ignored
    """
    joist_interior_region = convex_hull(
        GeometryCollection([geom for geom in ordered_supports])
    )
    # NOTE (2025-10-27): An attempt was made to replace 'ordered_supports' with a list of the intersection
    # points of the joist_prototype. The idea being that it is only teh intersection points that defined
    # the support system of cantilevers. However, this proved to carry several unintended consequences to
    # how the geometry interacted such as when a beam intersected with both a wall and a column. It created
    # a strange polygon which was not desired. One solution to this might be to add rules such that the
    # beam cannot intersect a linear reaction element but this creates a problem for when we want to have
    # that rule not apply. It just creates an error, which could be caught, but is not useful.
    # Thus, I went back to the original implementation of using the actual support lines. This allowed
    # all tests to pass again without creating geometry errors.

    # The below commented code is an artifact of the above attempt which is preserved here as a reminder of that
    # implementation in the event it proves to have an advantage over the current one.
    # - CMF

    # if joist_interior_region.geom_type == "LineString":
    #     joist_interior_region = joist_interior_region.buffer(1, cap_style="flat")

    joist_interior = joist_interior_region.intersection(joist_prototype, grid_size=1e-3)
    joist_interior_length = joist_interior.length

    cantilevers = (
        ops.split(joist_prototype, joist_interior_region) - joist_interior_region
    )
    cantilever_segments = {"A": 0.0, "B": 0.0}
    if isinstance(cantilevers, LineString):
        split_a = cantilevers
        split_b = Point()  # A geometry of length 0
    elif hasattr(cantilevers, "geoms"):
        split_a, split_b = cantilevers.geoms
    else:
        split_a, split_b = Point(joist_interior.coords[0]), Point(
            joist_interior.coords[-1]
        )
    a_orig = split_a.length
    b_orig = split_b.length

    if abs_tol is not None:
        split_a = (
            split_a if split_a.length > abs_tol else Point(joist_interior.coords[0])
        )
        split_b = (
            split_b if split_b.length > abs_tol else Point(joist_interior.coords[-1])
        )
    elif rel_tol:
        split_a = (
            split_a
            if (split_a.length / joist_interior_length) > rel_tol
            else Point(joist_interior.coords[0])
        )
        split_b = (
            split_b
            if (split_b.length / joist_interior_length) > rel_tol
            else Point(joist_interior.coords[-1])
        )

    if split_a.distance(ordered_supports[0]) < split_a.distance(ordered_supports[-1]):
        print(joist_prototype.intersection(ordered_supports[0], grid_size=1e-3))
        a_intersection = get_intersection(joist_prototype, ordered_supports[0], "")
        b_intersection = get_intersection(joist_prototype, ordered_supports[-1], "")

        if a_intersection is not None:
            a_intersect = a_intersection[0]
        else:
            a_intersect = None

        if b_intersection is not None:
            b_intersect = b_intersection[0]
        else:
            b_intersect = None

        cantilever_segments = {
            "A": split_a.length,
            "A_intersection": a_intersect,
            "A_orig": a_orig,
            "B": split_b.length,
            "B_intersection": b_intersect,
            "B_orig": b_orig,
        }
    else:
        b_intersection = get_intersection(joist_prototype, ordered_supports[0], "")
        a_intersection = get_intersection(joist_prototype, ordered_supports[-1], "")

        if a_intersection is not None:
            a_intersect = a_intersection[0]
        else:
            a_intersect = None

        if b_intersection is not None:
            b_intersect = b_intersection[0]
        else:
            b_intersect = None

        cantilever_segments = {
            "A": split_b.length,
            "A_intersection": a_intersection,
            "A_orig": b_orig,
            "B": split_a.length,
            "B_intersection": b_intersection,
            "B_orig": a_orig,
        }
    return cantilever_segments


def find_extent_intersections(
    element_geoms: list[LineString], extent_geoms: list[LineString]
) -> list[Optional[LineString]]:
    """
    Returns a list of LineString representing the extent geometries put
    into the order of the element geometries. If an extent geometry intersects
    with an element_geometry, the resulting list will have the extent geometry
    in the same corresponding list position that the element geometry is in.
    If there is no such intersection, then that list position will be None.
    """
    extents_array = np.array(extent_geoms)
    acc = []
    for element_geom in element_geoms:
        mask = intersects(element_geom, extents_array)
        if mask.any():
            extent = extents_array[mask][
                0
            ]  # Assume the first one until a better idea comes
            acc.append(extent)
        else:
            acc.append(None)
    return acc


def create_extent_polygon(
    element_geom: LineString, extent_geom: Optional[LineString] = None
) -> Polygon:
    """
    Returns a Polygon representing the bounding box of the union
    of 'element_geom' and 'extent_geom'
    """
    if extent_geom is None:
        return None
    return box(*(union(element_geom, extent_geom).bounds))


def split_polygon(
    polygon: Polygon, joist_orientation: str, split_points: list[tuple[float, float]]
) -> list[Polygon]:
    split_locations = []
    for split_point in split_points:
        if joist_orientation == "vertical":
            split_location = split_point[0]
        elif joist_orientation == "horizontal":
            split_location = split_point[1]
        if Point(split_point).within(polygon):
            split_locations.append(split_location)
    polygons = polygon_splitter(polygon.bounds, split_locations, joist_orientation)
    return polygons


def polygon_splitter(
    poly_bounds: tuple[float, float, float, float],
    split_locations: list[float],
    joist_orientation: str,
) -> list[tuple[float, float, float, float]]:
    """
    Returns a list of box bounds describing the sub-boxes remaining after splitting
    """
    xmin, ymin, xmax, ymax = poly_bounds
    sub_polys = []
    if joist_orientation == "horizontal":
        origin = ymin
        for sl in split_locations:
            sub_poly = box(xmin, origin, xmax, sl)
            sub_polys.append(sub_poly)
            origin = sl
        else:
            sub_poly = box(xmin, origin, xmax, ymax)
            sub_polys.append(sub_poly)
        if sub_polys:
            return sub_polys
        return [box(xmin, ymin, xmax, ymax)]
    elif joist_orientation == "vertical":
        origin = xmin
        for sl in split_locations:
            sub_poly = box(origin, ymin, sl, ymax)
            sub_polys.append(sub_poly)
            origin = sl
        else:
            sub_poly = box(origin, ymin, xmax, ymax)
            sub_polys.append(sub_poly)
        if sub_polys:
            return sub_polys
        return [box(xmin, ymin, xmax, ymax)]


def translate_joist_to_point(
    joist_geom: LineString, joist_orientation: str, intersection_point: Point
) -> LineString:
    """
    Returns a LineString representing 'joist_geom' translated so that it intersects with 'intersection_point'
    """
    point_i, point_j = joist_geom.coords
    ix, iy = point_i
    jx, jy = point_j
    ipx, ipy = intersection_point.coords[0]
    if joist_orientation == "horizontal":
        return LineString([(ix, ipy), (jx, ipy)])
    elif joist_orientation == "vertical":
        return LineString([(ipx, iy), (ipx, jy)])


def get_system_bounds(
    joist_prototype: LineString,
    joist_supports: list[LineString],
    normal: bool = True,
    extent_polygon: Optional[Polygon] = None,
) -> tuple[float, float, float, float]:
    """
    Returns the minx, miny, maxx, maxy bounding box of all the LineStrings in 'joist_supports',
    taken as a group.
    """
    if normal:
        all_lines = MultiLineString(joist_supports + [joist_prototype])
        bbox = all_lines.bounds
        return bbox
    else:
        joist_vector = np.abs(get_direction_vector(joist_prototype))
        # This is one of the places where orthogonality is assumed
        joist_orientation = None
        if joist_vector[0] > joist_vector[1]:
            joist_orientation = "horizontal"
        elif joist_vector[1] > joist_vector[0]:
            joist_orientation = "vertical"
        else:
            print(f"JOIST ORIENTATION VERIANT: {joist_prototype=}")

        overlap_polys = []
        overlap_poly = None
        for start_support in joist_supports:
            for end_support in joist_supports:
                if start_support == end_support:
                    continue
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
                    overlap_region = ld.get_overlap_coords(
                        pa0[1], pa1[1], pb0[1], pb1[1]
                    )
                    if overlap_region is not None:
                        overlap_poly = box(
                            pa0[0], overlap_region[0], pb1[0], overlap_region[1]
                        )
                if overlap_poly is not None:
                    if extent_polygon is not None:
                        overlap_poly = extent_polygon.intersection(overlap_poly)
                    overlap_polys.append(overlap_poly)
        return MultiPolygon(overlap_polys).bounds


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
    parallel_vector = column_vector / column_vector_norm
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
    joist_intersections = joist_prototype.intersection(all_supports, grid_size=1e-3)
    assert joist_intersections.geom_type != "Point"
    assert not joist_intersections.is_empty
    ordered_intersections = order_nodes_positive(joist_intersections.geoms)
    ordered_supports = []
    for point in ordered_intersections:
        for linestring in supports:
            if linestring.buffer(1e-3).intersects(
                point
            ):  # If using 1e-3 grid size, use 1e-3 buffer
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


def relate_point_to_line(point: Point, line: LineString) -> tuple:
    """
    REturns a tuple of ('left'/'right', 'above'/'below') to describe
    where teh point is in relation to the line
    """
    try:
        slope, intercept = ld.get_slope_and_intercept(*line.coords)
    except ZeroDivisionError:
        slope = 1
        intercept = float("inf")
    xp, yp = point.coords[0]
    yl = slope * xp + intercept
    delta_y = yl - yp
    if delta_y > 0 and slope >= 0:
        return ("right", "below")
    elif delta_y > 0 and slope < 0:
        return ("left", "below")
    elif delta_y < 0 and slope >= 0:
        return ("left", "above")
    elif delta_y < 0 and slope < 0:
        return ("right", "above")


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
    round_precision: int = 4,
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


def _group_vertices(
    vertices: list[Decimal | float], close=False
) -> list[tuple[Decimal, Decimal]]:
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


def rotate_90_coords(v: ArrayLike, precision: int = 6, ccw=True) -> tuple[float, float]:
    """
    Rotate the vector components, 'x1' and 'y1' by 90 degrees.

    'precision': round result to this many decimal places
    'ccw': if True, rotate counter-clockwise (clockwise, otherwise)
    """
    # v_angle = np.arctan2(v[1], v[0])
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
    return v @ rot


def check_2d_linestring_parallel(ls1: LineString, ls2: LineString, tol=1e-6) -> bool:
    """
    Returns True if ls1 and ls2 are parallel within an absolute tolerance
    """
    x1, y1 = get_direction_vector(ls1)
    x2, y2 = get_direction_vector(ls2)
    return math.isclose(abs(x1 * y2 - x2 * y1), 0.0, abs_tol=tol)


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

    rotated_line = aff.translate(
        aff.rotate(line, -angle, origin=i_end, use_radians=True), xoff=-ix
    )
    rotated_geoms = [
        aff.translate(
            aff.rotate(geom, -angle, origin=i_end, use_radians=True), xoff=-ix
        )
        for geom in geoms
    ]

    return rotated_line, rotated_geoms


def explode_polygon(p: Polygon) -> list[LineString]:
    """
    Explodes the exterior of the polygon in to a list of individual line segments
    """
    ext_ls = LineString(p.exterior)
    exploded = [LineString(tup) for tup in zip(ext_ls.coords, ext_ls.coords[1:])]
    return exploded


def get_rectangle_centerline(p: Polygon, on_long_edge: bool = False) -> LineString:
    """
    Returns the centerline of the Polygon 'p' assuming that 'p' represents
    a regular rectangle with a long dimension and a short dimension.
    The LineString is created with a +ve X-bias.
    """
    rectangle_edges = explode_polygon(p)
    sorted_edges = sorted(rectangle_edges, key=lambda x: x.length)
    short_edges = sorted_edges[:2]
    long_edges = sorted_edges[2:]
    if on_long_edge:
        edge1, edge2 = long_edges
    else:
        edge1, edge2 = short_edges
    start, end = order_nodes_positive([edge1.centroid, edge2.centroid])
    center_line = LineString([start, end])
    return center_line


def calculate_trapezoid_area_sums(member_loads: list[list[list[tuple]]]) -> list[float]:
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


def min_with_none(x: list, key=None):
    """
    Returns min but ignoring None
    """
    cleaned = [y for y in x if y is not None]
    return min(cleaned, key=key)


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
