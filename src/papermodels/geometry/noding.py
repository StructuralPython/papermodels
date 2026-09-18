"""
Phase 2 — tolerance extend/trim noding (design §7).

A preprocessing step that repairs "approximately intersecting" geometry so the
arrangement is topologically clean *before* the STRtree / node-registry build
(design §6, step 2).  It is governed by ``noding_abs_tol``, an **absolute**
tolerance in real-world units that is deliberately **separate** from the node
canonicalization tolerance (decided, design §9): endpoint repair may need to
reach wider than crossing-merging.

The operation, per line endpoint:

- Query the spatial index for geometry within ``noding_abs_tol``.
- If the nearest point on a candidate is within ``noding_abs_tol`` (and the
  endpoint is not already touching something), relocate the endpoint onto that
  nearest point — **extending** the line if it fell short, **trimming** it if it
  overshot.  Both are the same move: the endpoint goes to the nearest point on
  the candidate, a correction bounded by ``noding_abs_tol``.

Absolute (not relative) tolerance is deliberate: a relative tolerance would make
endpoint behaviour depend on line length, reintroducing scale-dependence.

**Cascading snaps.** Moving endpoint A onto line L can pull it off line M it
previously touched, so a single sweep is not guaranteed to converge (design §7).
Each pass works from a frozen snapshot and only moves endpoints that are *not*
already coincident with some geometry; the passes iterate to a fixed point.  If
the pass cap is reached first, a warning is emitted and the best-so-far
arrangement is returned — the **cap + warning** decision (§9), not a hard error.

Only ``LineString`` endpoints move.  Polygons (columns) and any ``fixed_ids``
(e.g. wall centerlines, a derived spine we do not want dragged around) are valid
snap *targets* but are never themselves moved.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
from warnings import warn

from shapely import STRtree
from shapely.geometry import LineString, Point
from shapely.geometry.base import BaseGeometry
from shapely.ops import nearest_points

GeomId = str

DEFAULT_MAX_PASSES = 10
# An endpoint within this distance of a candidate is treated as already touching
# it ("anchored") and is left alone.  Kept far below any sensible noding
# tolerance; inputs are grid-snapped (~1e-3) and a fresh snap lands on the
# target to float precision (~1e-12).
COINCIDENCE_EPS = 1e-9


@dataclass
class NodingReport:
    """Outcome of a noding run."""

    passes: int
    converged: bool
    moves: int


def node_geometries(
    geometries: dict[GeomId, BaseGeometry],
    noding_abs_tol: float,
    max_passes: int = DEFAULT_MAX_PASSES,
    fixed_ids: Optional[set[GeomId]] = None,
    suppress_warnings: bool = False,
) -> tuple[dict[GeomId, BaseGeometry], NodingReport]:
    """
    Return ``(repaired_geometries, report)``.

    ``geometries`` is mapped ``GeomId -> geometry`` (walls already reduced to
    their centerline, matching the model's index representation).  Line endpoints
    are extended/trimmed onto nearby geometry within ``noding_abs_tol``.

    ``fixed_ids`` geometries are never moved but remain snap targets.
    """
    if noding_abs_tol <= 0:
        raise ValueError(f"noding_abs_tol must be > 0, got {noding_abs_tol!r}")
    fixed = fixed_ids or set()

    current = dict(geometries)
    total_moves = 0
    passes = 0
    converged = False
    while passes < max_passes:
        passes += 1
        current, moved = _noding_pass(current, noding_abs_tol, fixed)
        total_moves += moved
        if moved == 0:
            converged = True
            break

    if not converged and not suppress_warnings:
        warn(
            f"Tolerance noding did not converge within {max_passes} passes "
            f"(noding_abs_tol={noding_abs_tol}). Proceeding with the best-so-far "
            "arrangement; some near-touches may remain unrepaired."
        )
    return current, NodingReport(passes=passes, converged=converged, moves=total_moves)


def _noding_pass(
    geometries: dict[GeomId, BaseGeometry],
    tol: float,
    fixed: set[GeomId],
) -> tuple[dict[GeomId, BaseGeometry], int]:
    """
    One frozen-snapshot sweep.  Proposals are computed against ``geometries`` and
    applied all at once, so within a pass no endpoint sees another's move.
    """
    ids = list(geometries.keys())
    geoms = [geometries[gid] for gid in ids]
    tree = STRtree(geoms)

    proposals: dict[GeomId, LineString] = {}
    for gid, geom in geometries.items():
        if geom.geom_type != "LineString" or gid in fixed:
            continue
        coords = list(geom.coords)
        new_coords = list(coords)
        changed = False
        for idx in (0, len(coords) - 1):
            snapped = _snap_endpoint(
                Point(coords[idx]), gid, geometries, ids, tree, tol
            )
            if snapped is not None:
                new_coords[idx] = snapped
                changed = True
        if changed:
            proposals[gid] = LineString(new_coords)

    if not proposals:
        return geometries, 0
    result = dict(geometries)
    result.update(proposals)
    return result, len(proposals)


def _snap_endpoint(
    point: Point,
    gid: GeomId,
    geometries: dict[GeomId, BaseGeometry],
    ids: list[GeomId],
    tree: STRtree,
    tol: float,
) -> Optional[tuple[float, float]]:
    """
    Return the ``(x, y)`` the endpoint should move to, or ``None`` if it should
    stay (already touching something, or nothing within tolerance).
    """
    best_xy: Optional[tuple[float, float]] = None
    best_d: Optional[float] = None
    for ci in tree.query(point, predicate="dwithin", distance=tol):
        other_gid = ids[int(ci)]
        if other_gid == gid:
            continue
        _, q = nearest_points(point, geometries[other_gid])
        d = point.distance(q)
        if d <= COINCIDENCE_EPS:
            return None  # already anchored to some geometry; do not move
        if d <= tol and (best_d is None or d < best_d):
            best_d, best_xy = d, (q.x, q.y)
    return best_xy
