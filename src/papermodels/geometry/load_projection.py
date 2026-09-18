"""
Robust projection of loading areas onto a member.

``load_distribution.get_distributed_loads_from_projected_polygons`` rotates the
member and its loading polygons to horizontal and then derives singularity
functions from the polygon edges. After an arbitrary rotation, an edge that is
geometrically vertical keeps a ~1e-17 x-offset between its endpoints, which the
singularity step turns into a ramp of slope ~1e16; the reconstructed load then
collapses to zero length (e.g. a rectangular trib band about a joist at 30 deg).

This wrapper performs the same rotation itself (identical convention: rotate
about the member's positive-x start node, which lands at x = 0), snaps the
rotated polygons to ``SNAP_GRID`` so vertical edges are exactly vertical, and
hands ``load_distribution`` an exactly-horizontal member so its own rotation is
the identity. Output is unchanged wherever the library was already correct.

This is a local guard; the underlying fix belongs in load_distribution.
"""

from __future__ import annotations

import load_distribution as ld
from shapely import LineString, set_precision

from . import geom_ops

# Far below drawing precision (~1e-3) and far above rotation noise (~1e-16).
SNAP_GRID = 1e-9


def project_loading_areas(member: LineString, applied_loading_areas: list) -> list:
    """
    Drop-in replacement for
    ``ld.get_distributed_loads_from_projected_polygons(member, applied_loading_areas)``.

    'applied_loading_areas': list of (Polygon, LoadingGeometry-or-None) pairs.
    """
    if not applied_loading_areas:
        return []
    load_geoms = [area[0] for area in applied_loading_areas]
    rotated_member, rotated_geoms = geom_ops.rotate_to_horizontal(member, load_geoms)
    (x0, y0), (x1, _) = rotated_member.coords[0], rotated_member.coords[-1]
    horizontal_member = LineString([(0.0, y0), (abs(x1 - x0), y0)])
    snapped = [set_precision(geom, SNAP_GRID) for geom in rotated_geoms]
    snapped_areas = [
        (geom, *area[1:]) for geom, area in zip(snapped, applied_loading_areas)
    ]
    return ld.get_distributed_loads_from_projected_polygons(
        horizontal_member, snapped_areas
    )
