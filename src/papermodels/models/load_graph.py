import numpy as np
from shapely import intersects, Geometry

def intersections_by_all_ranks(geoms_by_rank: dict[str, Geometry]):
    """
    Returns a dictionary of intersection arrays by rank.

    The underlying assumption is that all geometries of a lower rank
    can be supported by any of the geometries of a higher rank. 
    Geometries of "Rank:0" can be supported by geometries of "Rank:1",
    "Rank:2", or "Rank:3", etc.

    The returned dictionary will contain keys for n - 1 ranks so if the
    the 'geoms_by_rank' dict goes to "Rank:4", then the return dict will
    go to "Rank:3". This is because the geometries on "Rank:4" are not
    being supported by anything (i.e. they are the end of the load path)
    """
    nodes = []
    for geoms in rank_geoms.values():
        nodes += geoms

    ranks = list(rank_geoms.keys())

    rank_intersections = {}
    for i, rank in enumerate(ranks):
        for sub_rank in ranks[i+1:]:
            a = np.array([rank_geoms[rank]])
            b = np.array([rank_geoms[sub_rank]])
            intersections = intersects(a, b.T)
            rank_intersections.update({rank: intersections})
        
    