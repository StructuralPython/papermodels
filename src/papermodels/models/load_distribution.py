from typing import Optional
from dataclasses import dataclass
from shapely.geometry import Polygon, LineString, Point
import itertools
import more_itertools


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
OL0 = Overlap(x0=-5., x1=10., ma=-4., ba=15., mb=2, bb=-2.)
OL1 = Overlap(x0=12.3, x1=16.3, ma=0.5, ba=6.1, mb=-3.34, bb=2.5)


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
    Returns a Polygon generated from an Overlap.
    """
    ya0 = ovlp.ma * ovlp.x0 + ovlp.ba
    ya1 = ovlp.ma * ovlp.x1 + ovlp.ba
    yb0 = ovlp.mb * ovlp.x0 + ovlp.bb
    yb1 = ovlp.mb * ovlp.x1 + ovlp.bb

    y0 = max(ya0, yb0) - min(ya0, yb0)
    y1 = max(ya1, yb1) - min(ya1, yb1)

    coords = [
        [ovlp.x0, 0],
        [ovlp.x0, y0],
        [ovlp.x1, y1],
        [ovlp.x1, 0]
    ]
    return Polygon(coords)


def get_overlap_regions(p: Polygon) -> list[Overlap]:
    """
    Returns a list of overlapping regions described by the vertices in the polygon, 'p'.
    An overlapping region is an enclosed region in the polygon bounded by an upper y and a lower y.
    """
    vertices = [vertex for vertex in p.exterior.coords]
    edge_indexes = [[k,k+1] for k in range(len(vertices) - 1)]
    edge_combinations = itertools.combinations(edge_indexes, 2)
    overlapping_regions = []
    for edge_a, edge_b in edge_combinations:
        pa0_idx, pa1_idx = edge_a
        pb0_idx, pb1_idx = edge_b
        pa0, pa1 = vertices[pa0_idx], vertices[pa1_idx]
        pb0, pb1 = vertices[pb0_idx], vertices[pb1_idx]
        overlapping_region = get_overlap_region(pa0, pa1, pb0, pb1)
        if overlapping_region is None: continue
        overlapping_regions.append(overlapping_region)
    return overlapping_regions


def get_overlap_region(pa0, pa1, pb0, pb1) -> Optional[Overlap]:
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
    return Overlap(x0, x1, a_slope, a_intercept, b_slope, b_intercept)


def get_overlap_coords(xa0: float, xa1: float, xb0: float, xb1: float) -> Optional[tuple[float, float]]:
    """
    Returns a tuple representing the starting and ending x-coordinates of the region of overlap
    between two line segments, a and b, defined solely by each line segment's start and end 
    x-coordinates, x0 and x1.

    If no overlap exists, None is returned.
    """
    # Five possible states captured by four conditions
    # (The fifth state is when the coordinates of xa and xb
    # are the same which will be captured by condition 3 and 4 but by 
    # condition 3 first)

    # 1. Overlapping region with b "ahead" of a
    if xb0 < xa1 < xb1 and xa0 < xb0:
        return (xa1, xb1)
    # 2. Overlapping region with a "ahead" of b
    elif xb0 < xa0 < xb1 and xa1 > xb1:
        return (xa0, xb1)
    # 3. Overlapping region with a "inside" b
    elif xa0 >= xb0 and xa1 <= xb1:
        return (xa0, xa1)
    # 4. Overlapping region with b "inside" a
    elif xb0 >= xa0 and xb1 <= xa1:
        return (xb0, xb1)
    else:
        return None


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
    

def sort_points_in_polygon(p: Polygon, axis: str = "x", ascending: bool = True) -> list[Point]:
    """
    Returns the list of vertices in 'p' sorted according to 'axis' and 'ascending'.
    """
    sorted_acc = list(set([vertex for vertex in p.exterior.coords]))
    sort_index = 0 if axis == 'x' else 1
    sorted_acc.sort(key=lambda x: x[sort_index])
    return sorted_acc


def min_x_point(p: Polygon) -> int:
    """
    Returns the index of the vertex in 'p' with the minimum value of x coordinate.
    If there are more than one vertex with the same minimum, the vertex
    with lowest y value coordinate will be returned.
    """
    vertices = [vertex for vertex in p.exterior.coords]
    unique_vertices = list(set([vertex for vertex in p.exterior.coords]))
    sorted_vertices = sorted(unique_vertices, key=lambda x: (x[0], x[1]))
    min_vertex = sorted_vertices[0]    
    index = vertices.index(min_vertex)
    return index


def get_poly_edges(p: Polygon, start: tuple[float, float]) -> list[tuple[int, int]]:
    """
    Returns a list of tuples representing the edges of the polygon in clockwise order
    starting from the vertex of 'start'.
    """
    vertices = [vertex for vertex in p.exterior.coords]
    num_vertices = len(vertices)
    start_idx = vertices.index(start)
    vertex_cycler = itertools.cycle(vertices)
    vertex_cycler = more_itertools.consume(vertex_cycler, start_idx)
    
    edges = []
    for idx, vertex in enumerate(vertex_cycler):
        if idx == num_vertices:
            break
        else:
            pass
            
            
        
        
    
    
    
        