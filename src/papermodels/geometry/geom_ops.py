from typing import Optional
import numpy as np

from shapely import Point, MultiPoint, LineString, Polygon


def get_intersection(i_geom: LineString, j_geom: LineString | Polygon, j_tag: str) -> Optional[tuple[str, Point, LineString]]:
    """
    Returns the details of the intersection 
    """
    intersection_point = (
        i_geom & j_geom
        if j_geom.geom_type != "Polygon"
        else i_geom & j_geom.exterior
    )
    if intersection_point.is_empty:
        return
    if intersection_point.geom_type == "MultiPoint":
        intersection_point = Point(
            np.array(
                [
                    np.array(geom.coords[0])
                    for geom in intersection_point.geoms
                ]
            ).mean(axis=1)
        )
    intersection = (j_tag, intersection_point, j_geom)
    return intersection