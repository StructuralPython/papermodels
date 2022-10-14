from typing import Optional
from dataclasses import dataclass
from shapely.geometry import (
    Polygon,
    LineString,
    Point,
    MultiPolygon,
    GeometryCollection,
)
import itertools
import more_itertools

## Testing
from IPython.display import display


@dataclass
class Overlap:
    """
    Represents attributes of two line segments, a and b, that are one above the other,
    overlapping.
    Each line segment is defined by the following attributes representing
    each line as the equation y = mx + b
        x0, x1: starting and ending x coordinates
        m, b: slope and y-intercept
    """

    x0: float
    x1: float
    ma: float
    ba: float
    mb: float
    bb: float


# Examples
OL0 = Overlap(x0=-5.0, x1=10.0, ma=-4.0, ba=15.0, mb=2, bb=-2.0)
OL1 = Overlap(x0=12.3, x1=16.3, ma=0.5, ba=6.1, mb=-3.34, bb=2.5)


@dataclass
class Singularity:
    """
    Represents a singularity function with customizable attributes. The
    Singularity can be called like a function to generate a 'y' value
    from a given 'x' value.

    The singularity function is assumed to be linear and defined over a
    certain domain between x0 and x1, a slope of m, and with an initial
    y-value of y0. The resulting y-value for y1 will be rounded to
    'precision'.

    This class is defined as a callable class instead of as a function
    so that the attributes can be altered as required for appropriate
    scaling of the singularity function after it has been defined.
    """

    x0: float
    x1: float
    m: float
    y0: float
    precision: int

    def __call__(self, x: float) -> float:
        """
        Returns a value for 'y' calculated from 'x' and the attributes
        stored in the singularity function class instance.
        """
        y0 = self.y0
        x0 = self.x0
        x1 = self.x1
        m = self.m
        return round((x > x0) * (y0 + m * (x - x0)) * (x <= x1), self.precision)


def get_projected_polygons(p: Polygon) -> list[Polygon]:
    """
    Returns the projected trapezoids corresponding to the polygon,
    'p'.
    """
    overlap_regions = get_overlap_regions(p)
    sorted_regions = sorted(overlap_regions, key=lambda x: x.x0)
    projected_polygons = []
    for overlap_region in sorted_regions:
        projected_polygon = overlap_region_to_polygon(overlap_region)
        projected_polygons.append(projected_polygon)
    return projected_polygons


def overlap_region_to_polygon(ovlp: Overlap) -> Polygon:
    """
    Returns a projected Polygon generated from an Overlap.
    """
    y0, y1 = get_range(ovlp)

    coords = [[ovlp.x0, 0], [ovlp.x0, y0], [ovlp.x1, y1], [ovlp.x1, 0]]
    return Polygon(coords)


def overlap_region_to_singularity(ovlp: Overlap, precision: int = 6) -> callable:
    """
    Returns a singularity function generated from 'ovlp'
    """
    y0, y1 = get_range(ovlp)
    x0 = ovlp.x0
    x1 = ovlp.x1
    m = (y1 - y0) / (x1 - x0)
    singularity_function = Singularity(x0, x1, m, y0, precision)
    return singularity_function


def singularities_to_polygon(los: list[Singularity]) -> Polygon:
    """
    Returns a Polygon in the shape of the singularity function.
    """
    sorted_sings = sorted(los, key=lambda x: x.x1)
    x_acc = []
    for idx, sing in enumerate(sorted_sings):
        if idx == 0:
            x_acc.append(round(sing.x0, sing.precision))
            continue
        x_acc.append(round(sing.x0 + 10**-sing.precision, sing.precision))
        x_acc.append(round(sing.x1, sing.precision))

    # else:
    x_acc.append(round(sing.x1, sing.precision))
    x_acc = sorted(list(set(x_acc)))
    y_vals = [sum([sing(x) for sing in los]) for x in x_acc]

    xy_vals = list(zip(x_acc, y_vals))
    xy_vals = xy_vals + [(x_acc[-1], 0)]
    poly = Polygon(xy_vals)
    return poly


def apply_total_load(sing: Singularity, total_load: float) -> Singularity:
    """
    Returns a Singularity function that has been scaled to represent a trapezoid
    with an area that is equal to the value of 'total_load'.
    """
    delta_x = sing.x1 - sing.x0

    area = total_load
    slope = sing.m
    load_0 = (2 * area - delta_x**2) / (2 * delta_x)

    scaled_sing = Singularity(
        x0=sing.x0, x1=sing.x1, m=slope, y0=load_0, precision=sing.precision
    )
    return scaled_sing


def get_range(ovlp: Overlap) -> tuple[float, float]:
    """
    Returns a tuple representing the values of y0 and y1 calculated
    from the information in 'ovlp' where y0 and y1 correlate with
    the locations of x0 and x1.
    """
    ya0 = ovlp.ma * ovlp.x0 + ovlp.ba
    ya1 = ovlp.ma * ovlp.x1 + ovlp.ba
    yb0 = ovlp.mb * ovlp.x0 + ovlp.bb
    yb1 = ovlp.mb * ovlp.x1 + ovlp.bb

    y0 = max(ya0, yb0) - min(ya0, yb0)
    y1 = max(ya1, yb1) - min(ya1, yb1)
    return y0, y1


def get_y_vals(ovlp: Overlap) -> tuple[float, float]:
    """
    Returns a tuple representing the values of y0 and y1 calculated
    from the information in 'ovlp' where y0 and y1 correlate with
    the locations of x0 and x1.
    """
    ya0 = min(ovlp.ma * ovlp.x0 + ovlp.ba, ovlp.mb * ovlp.x0 + ovlp.bb)
    ya1 = min(ovlp.ma * ovlp.x1 + ovlp.ba, ovlp.mb * ovlp.x1 + ovlp.bb)
    yb0 = max(ovlp.ma * ovlp.x0 + ovlp.ba, ovlp.mb * ovlp.x0 + ovlp.bb)
    yb1 = max(ovlp.ma * ovlp.x1 + ovlp.ba, ovlp.mb * ovlp.x1 + ovlp.bb)
    return ya0, ya1, yb0, yb1


def get_overlap_regions(p: Polygon, display_progress: bool = False) -> list[Overlap]:
    """
    Returns a list of overlapping regions described by the vertices in the polygon, 'p'.
    An overlapping region is an enclosed region in the polygon bounded by an upper y and a lower y.
    """
    vertices = [vertex for vertex in p.exterior.coords]
    edge_indexes = [[k, k + 1] for k in range(len(vertices) - 1)]
    edge_combinations = itertools.combinations(edge_indexes, 2)
    overlapping_regions = []
    is_convex = check_convex_polygon(p)
    for edge_a, edge_b in edge_combinations:
        pa0_idx, pa1_idx = edge_a
        pb0_idx, pb1_idx = edge_b
        pa0, pa1 = vertices[pa0_idx], vertices[pa1_idx]
        pb0, pb1 = vertices[pb0_idx], vertices[pb1_idx]
        if is_convex:
            overlapping_region = get_overlap_region(
                pa0, pa1, pb0, pb1, p, convex=True, show_overlap=display_progress
            )
        else:
            overlapping_region = get_overlap_region(
                pa0, pa1, pb0, pb1, p, convex=False, show_overlap=display_progress
            )
        if overlapping_region is None:
            continue
        overlapping_regions.append(overlapping_region)
    return overlapping_regions


def get_overlap_region(
    pa0: tuple[float, float],
    pa1: tuple[float, float],
    pb0: tuple[float, float],
    pb1: tuple[float, float],
    p: Polygon,
    convex: bool,
    show_overlap: bool = False,
) -> Optional[Overlap]:
    """
    Returns an Overlap object representing the overlapping portion of the line
    segments a and b, represented by their respective points, p0 and p1.
    """
    xa0, _ = pa0
    xa1, _ = pa1
    xb0, _ = pb0
    xb1, _ = pb1
    xa0, xa1 = list(sorted([xa0, xa1]))
    xb0, xb1 = list(sorted([xb0, xb1]))
    xc0, xc1 = None, None
    if not convex:
        xc0, xc1 = get_void_extents(p)
    overlapping_coords = get_overlap_coords(xa0, xa1, xb0, xb1, xc0, xc1)
    if overlapping_coords is None:
        return None
    a_slope, a_intercept = get_slope_and_intercept(pa0, pa1)
    b_slope, b_intercept = get_slope_and_intercept(pb0, pb1)
    x0 = min(overlapping_coords)
    x1 = max(overlapping_coords)
    overlap = Overlap(x0, x1, a_slope, a_intercept, b_slope, b_intercept)
    if convex == False:
        ya0, ya1, yb0, yb1 = get_y_vals(overlap)

        overlap_test_poly = Polygon([[x0, ya0], [x0, yb0], [x1, yb1], [x1, ya1]])
        if show_overlap:
            display(GeometryCollection([overlap_test_poly, p]))
        if not check_poly_contains(overlap_test_poly, p):
            return None
    return overlap


def get_overlap_coords(
    xa0: float,
    xa1: float,
    xb0: float,
    xb1: float,
    xc0: Optional[float],
    xc1: Optional[float],
) -> Optional[tuple[float, float]]:
    """
    Returns a tuple representing the starting and ending x-coordinates of the region of overlap
    between two line segments, a and b, defined solely by each line segment's start and end
    x-coordinates, x0 and x1.

    xc0 and xc1 represent the extents of void space that may incur into the overlap region.
    If xc0 and xc1 are provided, then the returned extents are further bounded by their location.

    If no overlap exists, None is returned.
    """
    # Five possible states captured by four conditions
    # (The fifth state is when the coordinates of xa and xb
    # are the same which will be captured by condition 3 and 4 but by
    # condition 3 first)
    overlap_coords = None
    # 1. Overlapping region with b "ahead" of a
    if xb0 < xa1 < xb1 and xa0 < xb0:
        overlap_coords = (xa1, xb1)
    # 2. Overlapping region with a "ahead" of b
    elif xb0 < xa0 < xb1 and xa1 > xb1:
        overlap_coords = (xa0, xb1)
    # 3. Overlapping region with a "inside" b
    elif xa0 >= xb0 and xa1 <= xb1:
        overlap_coords = (xa0, xa1)
    # 4. Overlapping region with b "inside" a
    elif xb0 >= xa0 and xb1 <= xa1:
        overlap_coords = (xb0, xb1)

    # # Check for void space boundaries
    if overlap_coords is not None and xc0 is not None and xc1 is not None:
        xa, xb = overlap_coords
        if xc0 <= xa and xb <= xc1:
            pass
        elif xa <= xc1 < xb:
            overlap_coords = (xc1, xb)
        elif xc0 <= xb:
            overlap_coords = (xa, xc0)
    return overlap_coords


def get_slope_and_intercept(p0: tuple[float], p1: tuple[float]) -> tuple[float, float]:
    """
    Returns a tuple representing the slope and y-intercept of the line enclosing the
    line segment represented by the points 'p0' and 'p1'.
    """
    if not len(p0) == len(p1) == 2:
        return ValueError("Both p0 and p1 must be tuples of len == 2")
    x0, y0 = p0
    x1, y1 = p1
    slope = (y1 - y0) / (x1 - x0)
    y_intercept = y1 - slope * x1
    return slope, y_intercept


def check_convex_polygon(p: Polygon) -> bool:
    """
    Returns True if the polygon is convex.
    Returns False if the polygon is non-convex.
    """
    convex_hull = p.convex_hull
    if not isinstance(convex_hull - p, (Polygon, MultiPolygon)):
        return True
    return False


def get_void_extents(p: Polygon) -> bool:
    """
    Returns the region of any void space that exists in the concave Polygon, 'p'.
    """
    convex_hull = p.convex_hull
    void_space = convex_hull - p
    minx, miny, maxx, maxy = void_space.bounds
    return minx, maxx


def check_poly_contains(test_p: Polygon, p: Polygon) -> bool:
    """
    Returns True if the Polygon 'test_p' is wholly contained within 'p'.
    """
    test = False
    test_geom = (p | test_p) - p
    if test_geom.is_empty:
        test = True
    elif not isinstance(test_geom, Polygon) and not isinstance(
        test_geom, (GeometryCollection, MultiPolygon)
    ):
        test = True
    elif test_geom.area < 1e-6:
        test = True
    return test
