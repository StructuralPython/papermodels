from __future__ import annotations
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


    def __neg__(self):
        return Singularity(self.x0, self.x1, -self.m, -self.y0, self.precision)

    


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
    eps = 1e-12
    prev_x = 0
    for idx, sing in enumerate(sorted_sings):
        print(f"{sing=}")
        n = sing.precision
        if idx == 0:
            x_acc.append(0.)
        if prev_x != sing.x0:
            x_acc.append(prev_x + eps)
        x_acc.append(sing.x0)
        x_acc.append(sing.x0 + eps)
        if idx == len(sorted_sings) - 1:
            x_acc.append(sing.x1 - eps)
        x_acc.append(sing.x1)
        prev_x = sing.x1

    x_acc.append(sing.x1)
    x_acc = sorted(list(set(x_acc)))
    print(f"{x_acc=}")
    y_acc = [sum([sing(x) for sing in sorted_sings]) for x in x_acc[:-1]]
    y_acc += [0.]
    xy_acc = zip(x_acc, y_acc)
    print(xy_acc)
    return Polygon(xy_acc)



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


def get_singularity_functions(p: Polygon, display_progress: bool = False) -> tuple[list[Singularity]]:
    """
    Returns a 2-tuple, each element a list of Singularity. The first element is the overlapping regions
    of the convex hull of the polygon. The second element is the overlapping regions of the void spaces
    within the convex hull.
    """
    if p == p.convex_hull:
        return (get_overlap_regions(p), [])
    elif p.is_empty:
        return ([], [])
    else:
        void_regions = get_void_regions(p)
        convex_overlaps = get_overlap_regions(p.convex_hull)
        void_overlaps = []
        for void_region in void_regions:
            if void_region == void_region.convex_hull:
                void_overlaps += [-void_overlap for void_overlap in get_overlap_regions(void_region)]
                print("convex_void")
            else:
                convex_voids, negative_voids = get_singularity_functions(void_region)
                print("Hull void + voids")
                void_overlaps += [-void_overlap for void_overlap in convex_voids]
                void_overlaps += [-void_overlap for void_overlap in negative_voids]
    return (convex_overlaps, void_overlaps)



def get_overlap_regions(p: Polygon, display_progress: bool = False) -> list[Singularity]:
    """
    Returns a list of overlapping regions described by the vertices in the CONVEX polygon.
    The overlapping regions are described as a list of Singularity functions.
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
        overlapping_region = get_overlap_region(
            pa0, pa1, pb0, pb1, p, convex=True, show_overlap=display_progress
        )
        if overlapping_region is None:
            continue
        singularity_function = overlap_region_to_singularity(overlapping_region)
        overlapping_regions.append(singularity_function)
    return overlapping_regions


def get_overlap_region(
    pa0: tuple[float, float],
    pa1: tuple[float, float],
    pb0: tuple[float, float],
    pb1: tuple[float, float],
    p: Polygon,
    convex: bool,
    show_overlap: bool = True,
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

    overlapping_coords = get_overlap_coords(xa0, xa1, xb0, xb1)
    if overlapping_coords is None:
        return None
    a_slope, a_intercept = get_slope_and_intercept(pa0, pa1)
    b_slope, b_intercept = get_slope_and_intercept(pb0, pb1)
    x0 = min(overlapping_coords)
    x1 = max(overlapping_coords)
    overlap = Overlap(x0, x1, a_slope, a_intercept, b_slope, b_intercept)
    return overlap


def get_overlap_coords(
    xa0: float,
    xa1: float,
    xb0: float,
    xb1: float,
    void_regions: Optional[list[Polygon]] = None,
    # # xc0: Optional[float],
    # # xc1: Optional[float],
) -> Optional[tuple[float, float]]:
    """
    Returns a tuple representing the starting and ending x-coordinates of the region of overlap
    between two line segments, a and b, defined solely by each line segment's start and end
    x-coordinates, x0 and x1.

    # xc0 and xc1 represent the extents of void space that may incur into the overlap region.
    # If xc0 and xc1 are provided, then the returned extents are further bounded by their location.

    If no overlap exists, None is returned.
    """
    # Five possible states captured by four conditions
    # (The fifth state is when the coordinates of xa and xb
    # are the same which will be captured by condition 3 and 4 but by
    # condition 3 first)
    print(f"a: {xa0}, {xa1}")
    print(f"b: {xb0}, {xb1}")
    overlap_coords = None
    if xa0 == xa1 or xb0 == xb1: # ignore vertical lines
        return overlap_coords
    # 1. Overlapping region with b "ahead" of a
    if xb0 < xa1 < xb1 and xa0 < xb0:
        overlap_coords = (xb0, xa1)
    # 2. Overlapping region with a "ahead" of b
    elif xb0 < xa0 < xb1 and xa1 > xb1:
        overlap_coords = (xa0, xb1)
    # 3. Overlapping region with a "inside" b
    elif xa0 >= xb0 and xa1 <= xb1:
        overlap_coords = (xa0, xa1)
    # 4. Overlapping region with b "inside" a
    elif xb0 >= xa0 and xb1 <= xa1:
        overlap_coords = (xb0, xb1)
    print("overlap: ", overlap_coords)

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


def get_void_regions(p: Polygon) -> list[Polygon]:
    """
    Returns the void regions as a list of Polygon, if any
    void regions exist within 'p'. If no void region is found,
    then an empty list is returned.

    A void region is defined as any region that does not lie
    within the convex hull of 'p'.
    """
    convex_hull = p.convex_hull
    void_regions = convex_hull - p
    if isinstance(void_regions, Polygon):
        return [void_regions]
    elif isinstance(void_regions, MultiPolygon):
        return list(void_regions.geoms)
    else:
        return []


def get_void_extents(void_region: Polygon) -> bool:
    """
    Returns the region of any void space that exists in the concave Polygon, 'p'.
    """
    minx, miny, maxx, maxy = void_region.bounds
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
