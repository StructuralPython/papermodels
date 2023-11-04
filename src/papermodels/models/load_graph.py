import numpy as np
from shapely import intersects, Geometry, Point
from papermodels.datatypes.annotation import Annotation

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


def get_all_ranks(tagged_annotations: dict[Annotation, dict]) -> list[int]:
    """
    Returns a list of all rank values present in 'tagged_annotations'
    """
    acc = []
    for annot_attr in tagged_annotations.values():
        acc.append(annot_attr['rank'])
    return sorted(set(acc))


def get_geometry_intersections(tagged_annotations: dict[Annotation, dict]) -> dict[str, Geometry]:
    """
    Returns a dictionary of 
    """
    annots = list(tagged_annotations.keys())
    intersected_annotations = tagged_annotations.copy()
    for idx, i_annot in enumerate(annots):
        next_annot = min(idx + 1, len(annots) - 1)
        i_attrs = intersected_annotations[i_annot]
        i_rank = i_attrs['rank']
        i_page = i_annot.page
        intersections = []
        for j_annot in annots[next_annot:]:
            j_attrs = intersected_annotations[j_annot]
            j_rank = j_attrs['rank']
            j_page = j_annot.page
            if i_rank < j_rank and i_page == j_page:
                i_geom = i_attrs['geometry']
                j_geom = j_attrs['geometry']
                intersection_point = i_geom & j_geom if j_geom.geom_type != "Polygon" else i_geom & j_geom.exterior
                if intersection_point.is_empty:
                    continue
                elif intersection_point.geom_type == "MultiPoint":
                    intersection_point = Point(
                        np.array(
                            [np.array(geom.coords[0]) for geom in intersection_point.geoms]
                            ).mean(axis=1)
                    )
                else:
                    intersection = (j_attrs['tag'], intersection_point)
                intersections.append(intersection)
        i_attrs['intersections'] = intersections
    return intersected_annotations  
                    



