"""
Station engine for joist arrays.

A joist array is modelled in the prototype's own frame instead of by projecting
rays and re-testing the result:

- ``u`` is the unit vector along the joist, ``n`` the unit normal (the array
  axis).  A **station** ``s`` is a position along ``n``; ``t`` is a position
  along ``u``.  The joist line at station ``s`` is ``origin + s*n + t*u``.
- The frame is oriented **once**, by a dominant-axis rule on ``n`` (``n`` points
  +x when ``|n_x| > |n_y|``, else +y) with ``u = rotate_cw(n)``.  Inside the
  engine, points are ordered only by their ``t`` / ``s`` projections — never by
  the positive-x bias (``order_nodes_positive``), whose orientation flips at
  vertical where 1e-3 grid snapping produces exact x-ties.  End "A" of a joist
  is always its smaller-``t`` end.
- Every support is a straight line (a beam, or a wall's centerline).  A support
  is present at station ``s`` if it crosses the joist line there *inside the
  array region R*.  Projecting the endpoints of each support's pieces inside R
  onto ``n`` gives **breakpoints**; between two consecutive breakpoints the set
  of supports present is constant, so each interval is classified once.
- Target stations that fall where supports are missing are moved to the
  nearest valid station (a ``bisect`` over the breakpoints, then an outward
  scan limited to the joist's own window) — see :func:`resolve_stations`.
- Crossings are computed once, by interpolation along the support line, and
  snapped through a shared :class:`~papermodels.datatypes.geometry_model.NodeRegistry`;
  the joist's support end-points *are* those canonical nodes.  Nothing is
  re-intersected afterwards.
- Trib areas are bands between midpoints of adjacent final stations, clipped to
  R, so they tile R with no gaps or overlaps.

Only frame-based arithmetic and ``rotate_90_vector`` (characterized in
``tests/test_helper_characterization.py``) are used from the older helpers.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from shapely import LineString, Point, Polygon, MultiPolygon
from shapely.geometry.base import BaseGeometry

from . import geom_ops

# Numerical tolerance for comparing stations/t-values that were derived from the
# same inputs (not a drawing tolerance).
EPS = 1e-9
# Margin (real-world units) used when building "infinite" frame lines/bands.
_FAR = 1e3

MODE_FIXED = "fixed"  # ends from the outer supports + constant cantilevers
MODE_CONTAINER = "container"  # ends from the container boundary


# --------------------------------------------------------------------------
# Frame
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ArrayFrame:
    origin: np.ndarray
    u: np.ndarray  # along the joist
    n: np.ndarray  # array axis

    @classmethod
    def from_prototype(cls, prototype: LineString) -> "ArrayFrame":
        coords = list(prototype.coords)
        p0 = np.asarray(coords[0][:2], dtype=float)
        p1 = np.asarray(coords[-1][:2], dtype=float)
        d = p1 - p0
        length = float(np.linalg.norm(d))
        if length == 0.0:
            raise geom_ops.GeometryError("Joist prototype has zero length.")
        d = d / length
        n = np.asarray(geom_ops.rotate_90_vector(d, ccw=True), dtype=float)
        # Dominant-axis orientation of the array axis (see module docstring).
        if abs(n[0]) > abs(n[1]):
            sign = 1.0 if n[0] > 0 else -1.0
        else:
            sign = 1.0 if n[1] > 0 else -1.0
        n = n * sign
        u = np.asarray(geom_ops.rotate_90_vector(n, ccw=False), dtype=float)
        # Midpoint origin: identical whichever way round the prototype was drawn.
        origin = (p0 + p1) / 2.0
        return cls(origin=origin, u=u, n=n)

    def ts(self, xy) -> tuple[float, float]:
        """(t, s) of a point."""
        v = np.asarray(_xy(xy), dtype=float) - self.origin
        return float(v @ self.u), float(v @ self.n)

    def xy(self, t: float, s: float) -> tuple[float, float]:
        p = self.origin + t * self.u + s * self.n
        return float(p[0]), float(p[1])

    def line(self, s: float, t0: float, t1: float) -> LineString:
        return LineString([self.xy(t0, s), self.xy(t1, s)])

    def band(self, s0: float, s1: float, t0: float, t1: float) -> Polygon:
        return Polygon(
            [self.xy(t0, s0), self.xy(t1, s0), self.xy(t1, s1), self.xy(t0, s1)]
        )

    def s_range(self, geom: BaseGeometry) -> tuple[float, float]:
        ss = [self.ts(c)[1] for c in _all_coords(geom)]
        return min(ss), max(ss)

    def t_range(self, geom: BaseGeometry) -> tuple[float, float]:
        tt = [self.ts(c)[0] for c in _all_coords(geom)]
        return min(tt), max(tt)


def _xy(p) -> tuple[float, float]:
    if isinstance(p, Point):
        return (p.x, p.y)
    return (float(p[0]), float(p[1]))


def _all_coords(geom: BaseGeometry) -> list[tuple[float, float]]:
    if geom.is_empty:
        return []
    gt = geom.geom_type
    if gt == "Polygon":
        return [c[:2] for c in geom.exterior.coords]
    if gt in ("MultiPolygon", "MultiLineString", "MultiPoint", "GeometryCollection"):
        out = []
        for g in geom.geoms:
            out.extend(_all_coords(g))
        return out
    return [c[:2] for c in geom.coords]


# --------------------------------------------------------------------------
# Regions
# --------------------------------------------------------------------------


def extent_region(
    frame: ArrayFrame,
    t_range: tuple[float, float],
    extent_line: LineString,
    s_include: float = 0.0,
) -> Polygon:
    """
    The extent-line array region: the parallelogram (a rectangle in the frame)
    swept by ``t_range`` over the extent line's station range.  Replaces the
    axis-aligned bounding box of ``create_extent_polygon``.

    ``s_include`` (default: the prototype's own station, 0) is always included.
    """
    s0, s1 = frame.s_range(extent_line)
    s0, s1 = min(s0, s_include), max(s1, s_include)
    return frame.band(s0, s1, t_range[0], t_range[1])


def fixed_region(
    frame: ArrayFrame,
    support_a: LineString,
    support_b: LineString,
    s_range: tuple[float, float],
    cant_a: float,
    cant_b: float,
) -> Polygon:
    """
    The plain-prototype array region: the quadrilateral bounded by the outer
    supports ``support_a`` (smaller t) and ``support_b`` over ``s_range``,
    extended outward along ``u`` by the (non-negative) cantilevers.  Handles
    non-parallel supports (e.g. a triangular roof layout) without special cases.
    """
    s0, s1 = s_range
    a0, a1 = _t_on_line(frame, support_a, s0), _t_on_line(frame, support_a, s1)
    b0, b1 = _t_on_line(frame, support_b, s0), _t_on_line(frame, support_b, s1)
    return Polygon(
        [
            frame.xy(a0 - cant_a, s0),
            frame.xy(b0 + cant_b, s0),
            frame.xy(b1 + cant_b, s1),
            frame.xy(a1 - cant_a, s1),
        ]
    )


def clip_support_to_region(
    support_line: LineString, region: BaseGeometry, registry=None
) -> Optional[BaseGeometry]:
    """
    The part of a support line (a beam, or a wall's full-polygon centerline)
    inside ``region``.  Replaces the extent branch of ``geom_ops.get_intersection``,
    which took the oriented-bounding-box centerline of the *clipped wall piece*
    and could turn perpendicular when the clip was narrower than the wall.

    Endpoints are snapped through ``registry`` if given.  Returns None when the
    support does not enter the region.
    """
    piece = support_line.intersection(region)
    if piece.is_empty:
        return None
    if piece.geom_type == "MultiLineString":
        piece = max(piece.geoms, key=lambda g: g.length)
    if piece.geom_type not in ("LineString", "Point"):
        return None
    if registry is None:
        return piece
    if piece.geom_type == "Point":
        return Point(registry.coord[registry.get_or_create((piece.x, piece.y))])
    coords = list(piece.coords)
    first = registry.coord[registry.get_or_create(coords[0][:2])]
    last = registry.coord[registry.get_or_create(coords[-1][:2])]
    return LineString([first] + [c[:2] for c in coords[1:-1]] + [last])


def support_rectangularity(polygon: Polygon) -> float:
    """Area of ``polygon`` over the area of its oriented bounding box (1.0 = rectangle)."""
    from shapely import minimum_rotated_rectangle

    obb = minimum_rotated_rectangle(polygon)
    return polygon.area / obb.area if obb.area else 0.0


# --------------------------------------------------------------------------
# Supports, breakpoints, intervals
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Support:
    gid: str
    line: LineString  # beam line or wall centerline (2 points)


def _t_on_line(frame: ArrayFrame, line: LineString, s: float) -> float:
    """t where the (infinite) support line crosses station s."""
    (t0, s0), (t1, s1) = frame.ts(line.coords[0]), frame.ts(line.coords[-1])
    if s1 == s0:
        raise geom_ops.GeometryError("Support is parallel to the joists.")
    return t0 + (s - s0) * (t1 - t0) / (s1 - s0)


def _crossing_xy(frame: ArrayFrame, line: LineString, s: float) -> tuple[float, float]:
    """
    Where station s crosses the support, interpolated along the support itself
    so the point lies on the support to float precision.
    """
    (x0, y0), (x1, y1) = line.coords[0][:2], line.coords[-1][:2]
    s0, s1 = frame.ts((x0, y0))[1], frame.ts((x1, y1))[1]
    lam = (s - s0) / (s1 - s0)
    return (x0 + lam * (x1 - x0), y0 + lam * (y1 - y0))


def prototype_cantilevers(
    frame: ArrayFrame,
    prototype: LineString,
    supports: list[Support],
    tolerance: float = 0.0,
) -> tuple[float, float]:
    """
    Signed cantilevers (cant_a, cant_b) of the prototype past its outer supports,
    measured in the frame: A at the smaller-t end.  A cantilever below
    ``tolerance`` — including a negative one, where the prototype was drawn to a
    wall face short of the centerline — is 0, so joists always reach their
    outer supports.
    """
    if len(supports) < 2:
        raise geom_ops.GeometryError("A joist prototype needs at least two supports.")
    s_p = frame.ts(prototype.interpolate(0.5, normalized=True))[1]
    ts = sorted(_t_on_line(frame, sup.line, s_p) for sup in supports)
    t_min, t_max = frame.t_range(prototype)
    cant_a, cant_b = ts[0] - t_min, t_max - ts[-1]
    cant_a = cant_a if cant_a >= tolerance and cant_a > 0 else 0.0
    cant_b = cant_b if cant_b >= tolerance and cant_b > 0 else 0.0
    return cant_a, cant_b


def outer_supports(
    frame: ArrayFrame, prototype: LineString, supports: list[Support]
) -> tuple[Support, Support]:
    """The prototype's supports at its smaller-t (A) and larger-t (B) ends."""
    s_p = frame.ts(prototype.interpolate(0.5, normalized=True))[1]
    ordered = sorted(supports, key=lambda sup: _t_on_line(frame, sup.line, s_p))
    return ordered[0], ordered[-1]


def common_s_range(
    frame: ArrayFrame, supports: list[Support]
) -> Optional[tuple[float, float]]:
    """Station range over which every support exists (None if they never overlap)."""
    ranges = [frame.s_range(sup.line) for sup in supports]
    lo = max(r[0] for r in ranges)
    hi = min(r[1] for r in ranges)
    return (lo, hi) if hi > lo + EPS else None


@dataclass
class Interval:
    s0: float
    s1: float
    present: tuple[str, ...]  # support gids present, ordered by t at the midpoint

    @property
    def count(self) -> int:
        return len(self.present)


def usable_supports(frame: ArrayFrame, supports: list[Support]) -> list[Support]:
    """Supports that can bear joists: straight 2-point lines not parallel to u."""
    out = []
    for sup in supports:
        if len(sup.line.coords) != 2:
            raise geom_ops.GeometryError(
                f"Support {sup.gid} must be a single straight segment."
            )
        s0, s1 = frame.s_range(sup.line)
        if s1 - s0 > EPS:
            out.append(sup)
    return out


def support_spans(
    frame: ArrayFrame, supports: list[Support], region: BaseGeometry
) -> dict[str, list[tuple[float, float]]]:
    """Per support, the station ranges over which it lies inside ``region``."""
    spans: dict[str, list[tuple[float, float]]] = {}
    for sup in supports:
        piece = sup.line.intersection(region)
        pieces = (
            list(piece.geoms)
            if piece.geom_type == "MultiLineString"
            else ([piece] if piece.geom_type == "LineString" else [])
        )
        for p in pieces:
            if p.is_empty or p.length <= EPS:
                continue
            spans.setdefault(sup.gid, []).append(frame.s_range(p))
    return spans


def classify_intervals(
    frame: ArrayFrame,
    supports: list[Support],
    region: BaseGeometry,
    s_range: tuple[float, float],
) -> list[Interval]:
    """
    Partition ``s_range`` at every support-span endpoint; each interval carries
    the supports present in it, ordered by t.
    """
    s_min, s_max = s_range
    spans = support_spans(frame, supports, region)
    points = {s_min, s_max}
    for ranges in spans.values():
        for a, b in ranges:
            for v in (a, b):
                if s_min < v < s_max:
                    points.add(v)
    breakpoints = sorted(points)
    by_gid = {sup.gid: sup for sup in supports}
    intervals = []
    for a, b in zip(breakpoints, breakpoints[1:]):
        if b - a <= EPS:
            continue
        mid = (a + b) / 2.0
        present = [
            gid
            for gid, ranges in spans.items()
            if any(lo - EPS <= mid <= hi + EPS for lo, hi in ranges)
        ]
        present.sort(key=lambda gid: _t_on_line(frame, by_gid[gid].line, mid))
        intervals.append(Interval(a, b, tuple(present)))
    if not intervals:  # degenerate zero-width range
        intervals.append(Interval(s_min, s_max, ()))
    return intervals


def interval_at(intervals: list[Interval], s: float) -> int:
    """Index of the interval containing s (boundary stations use the inner side)."""
    starts = [iv.s0 for iv in intervals]
    idx = bisect.bisect_right(starts, s) - 1
    return min(max(idx, 0), len(intervals) - 1)


# --------------------------------------------------------------------------
# Stations
# --------------------------------------------------------------------------


def target_stations(
    s_min: float,
    s_max: float,
    spacing: float,
    initial_offset: float = 0.0,
    joist_at_start: bool = True,
    joist_at_end: bool = True,
) -> list[float]:
    """
    Nominal joist stations over [s_min, s_max].  Joists are laid at ``spacing``
    from the first station; no gap exceeds ``spacing`` (the closing gap to an
    end joist may be shorter).  Replaces ``geom_ops.get_joist_locations``, which
    used a line-to-line distance, ignored ``joist_at_end`` and let the last gap
    reach 1.5x spacing.
    """
    if spacing <= 0:
        raise ValueError(f"spacing must be > 0, got {spacing!r}")
    stations: list[float] = [s_min] if joist_at_start else []
    # Regular joists start one offset (or one spacing) in from the start edge.
    first = s_min + (initial_offset if initial_offset else spacing)
    k = 0
    while (s := first + k * spacing) < s_max - EPS:
        stations.append(s)
        k += 1
    if joist_at_end and (not stations or s_max - stations[-1] > EPS):
        stations.append(s_max)
    return stations


@dataclass
class StationEvent:
    kind: str  # "relocated" | "dropped" | "span_jump"
    target: float
    final: Optional[float] = None
    detail: str = ""


def _valid(count: int, left: Optional[int], right: Optional[int]) -> bool:
    """
    The station-validity rule from the spec: at least two supports, and not a
    support count that differs from both neighbours.  End stations (one
    neighbour) only need two supports — they cannot be told apart from a
    legitimate change in the support pattern.
    """
    if count < 2:
        return False
    if left is None or right is None:
        return True
    return count in (left, right)


def resolve_stations(
    targets: list[float],
    intervals: list[Interval],
    spacing: float,
    min_bearing: float,
) -> tuple[list[tuple[float, float]], list[StationEvent]]:
    """
    Returns ([(target, final_station)], events).  Invalid targets move to the
    nearest valid station within their window (half a spacing either side,
    and strictly between the previous final station and the next target);
    targets with no valid station in the window are dropped.
    """
    counts = [intervals[interval_at(intervals, t)].count for t in targets]
    resolved: list[tuple[float, float]] = []
    events: list[StationEvent] = []
    prev_final = -math.inf
    for i, t in enumerate(targets):
        left = counts[i - 1] if i > 0 else None
        right = counts[i + 1] if i + 1 < len(targets) else None
        idx = interval_at(intervals, t)
        if _valid(intervals[idx].count, left, right) and t > prev_final + EPS:
            resolved.append((t, t))
            prev_final = t
            continue
        lo = max(t - spacing / 2.0, prev_final + EPS)
        hi = t + spacing / 2.0
        if i + 1 < len(targets):
            hi = min(hi, targets[i + 1] - EPS)
        best = _nearest_valid(t, idx, intervals, left, right, lo, hi, min_bearing)
        if best is None:
            events.append(
                StationEvent(
                    "dropped",
                    t,
                    detail=f"no station with >= 2 supports within [{lo:.4f}, {hi:.4f}]",
                )
            )
            continue
        resolved.append((t, best))
        prev_final = best
        events.append(StationEvent("relocated", t, best))
    return resolved, events


def _nearest_valid(t, idx, intervals, left, right, lo, hi, min_bearing):
    """Scan outward from interval ``idx`` (found by bisect) for the closest valid point."""
    best = None
    best_d = math.inf
    for step in range(len(intervals)):
        progressed = False
        for j in {idx - step, idx + step}:
            if not 0 <= j < len(intervals):
                continue
            iv = intervals[j]
            if iv.s1 < lo - EPS or iv.s0 > hi + EPS:
                continue
            progressed = True
            if not _valid(iv.count, left, right):
                continue
            a = max(iv.s0 + min_bearing, lo)
            b = min(iv.s1 - min_bearing, hi)
            if a > b + EPS:
                continue
            cand = min(max(t, a), b)
            d = abs(cand - t)
            if d < best_d - EPS or (abs(d - best_d) <= EPS and cand < best):
                best, best_d = cand, d
        if not progressed and step > 0:
            break
    return best


# --------------------------------------------------------------------------
# Joists
# --------------------------------------------------------------------------


@dataclass
class Crossing:
    gid: str
    node_id: Optional[int]
    xy: tuple[float, float]
    t: float


@dataclass
class GeneratedJoist:
    target: float
    station: float
    geometry: LineString
    crossings: list[Crossing]  # ordered by t (A end first)
    trib_area: Optional[BaseGeometry] = None


def _crossings_at(frame, s, present, by_gid, registry) -> list[Crossing]:
    out = []
    for gid in present:
        xy = _crossing_xy(frame, by_gid[gid].line, s)
        node_id = None
        if registry is not None:
            node_id = registry.get_or_create(xy)
            xy = registry.coord[node_id]
        out.append(Crossing(gid, node_id, xy, frame.ts(xy)[0]))
    out.sort(key=lambda c: c.t)
    return out


def joist_at_station(
    frame: ArrayFrame,
    s: float,
    present: tuple[str, ...],
    supports: list[Support],
    mode: str,
    region: BaseGeometry,
    cant_a: float = 0.0,
    cant_b: float = 0.0,
    registry=None,
) -> Optional[tuple[LineString, list[Crossing]]]:
    """
    The joist at station ``s`` and its canonical support crossings.

    - MODE_FIXED: from end A's support crossing minus ``cant_a`` to end B's plus
      ``cant_b`` (constant cantilevers).
    - MODE_CONTAINER: the joist is the container cut along the station line;
      the ends are the container edges (extended to the outer supports if the
      container stops short of them), so cantilevers and backspans both vary.
    """
    by_gid = {sup.gid: sup for sup in supports}
    crossings = _crossings_at(frame, s, present, by_gid, registry)
    if len(crossings) < 2:
        return None
    t_a, t_b = crossings[0].t, crossings[-1].t
    if mode == MODE_FIXED:
        t_start, t_end = t_a - cant_a, t_b + cant_b
    elif mode == MODE_CONTAINER:
        t_lo, t_hi = frame.t_range(region)
        cut = frame.line(s, t_lo - 1.0, t_hi + 1.0).intersection(region)
        pieces = (
            list(cut.geoms)
            if cut.geom_type == "MultiLineString"
            else ([cut] if cut.geom_type == "LineString" else [])
        )
        spans = sorted(frame.t_range(p) for p in pieces if not p.is_empty)
        if not spans:
            return None
        start_span = next(
            (sp for sp in spans if sp[0] - EPS <= t_a <= sp[1] + EPS), None
        )
        end_span = next((sp for sp in spans if sp[0] - EPS <= t_b <= sp[1] + EPS), None)
        t_start = min(start_span[0] if start_span else t_a, t_a)
        t_end = max(end_span[1] if end_span else t_b, t_b)
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    def end_xy(t, crossing):
        # A joist end that sits exactly on a support is that support's node.
        return crossing.xy if abs(t - crossing.t) <= EPS else frame.xy(t, s)

    geom = LineString([end_xy(t_start, crossings[0]), end_xy(t_end, crossings[-1])])
    return geom, crossings


def trib_bands(
    frame: ArrayFrame,
    stations: list[float],
    region: BaseGeometry,
    s_range: tuple[float, float],
) -> list[BaseGeometry]:
    """
    Trib area per station: the band between the midpoints to its neighbours
    (the array's outer station limits at the ends), clipped to ``region``.
    Neighbouring bands share one boundary value, so the bands tile the region.
    """
    if not stations:
        return []
    s_min, s_max = s_range
    bounds = [s_min] + [(a + b) / 2.0 for a, b in zip(stations, stations[1:])] + [s_max]
    t_lo, t_hi = frame.t_range(region)
    out = []
    for i in range(len(stations)):
        band = frame.band(bounds[i], bounds[i + 1], t_lo - 1.0, t_hi + 1.0)
        out.append(band.intersection(region))
    return out


# --------------------------------------------------------------------------
# Whole array
# --------------------------------------------------------------------------


@dataclass
class ArrayResult:
    frame: ArrayFrame
    region: BaseGeometry
    s_range: tuple[float, float]
    intervals: list[Interval]
    joists: list[GeneratedJoist]
    events: list[StationEvent] = field(default_factory=list)


def build_array(
    frame: ArrayFrame,
    region: BaseGeometry,
    supports: list[Support],
    s_range: tuple[float, float],
    mode: str,
    spacing: float,
    initial_offset: float = 0.0,
    joist_at_start: bool = True,
    joist_at_end: bool = True,
    cant_a: float = 0.0,
    cant_b: float = 0.0,
    min_bearing: float = 1e-3,
    span_jump_tol: Optional[float] = None,
    registry=None,
) -> ArrayResult:
    """Lays out, resolves and builds every joist of one array (see module docstring)."""
    supports = usable_supports(frame, supports)
    intervals = classify_intervals(frame, supports, region, s_range)
    targets = target_stations(
        s_range[0], s_range[1], spacing, initial_offset, joist_at_start, joist_at_end
    )
    resolved, events = resolve_stations(targets, intervals, spacing, min_bearing)
    joists: list[GeneratedJoist] = []
    for target, s in resolved:
        present = intervals[interval_at(intervals, s)].present
        built = joist_at_station(
            frame, s, present, supports, mode, region, cant_a, cant_b, registry
        )
        if built is None:
            events.append(StationEvent("dropped", target, s, "fewer than 2 supports"))
            continue
        geom, crossings = built
        joists.append(GeneratedJoist(target, s, geom, crossings))
    for joist, band in zip(
        joists, trib_bands(frame, [j.station for j in joists], region, s_range)
    ):
        joist.trib_area = band
    tol = spacing / 2.0 if span_jump_tol is None else span_jump_tol
    events.extend(span_jumps(frame, joists, supports, tol))
    return ArrayResult(frame, region, s_range, intervals, joists, events)


def span_jumps(
    frame: ArrayFrame,
    joists: list[GeneratedJoist],
    supports: list[Support],
    tol: float,
) -> list[StationEvent]:
    """
    Flags abrupt changes in the outer supports between consecutive joists.  A
    jump can only happen where an outer support changes identity; the new
    support's crossing is compared with the old support's line extended to the
    same station, so a single angled support never triggers it.
    """
    by_gid = {sup.gid: sup for sup in supports}
    events = []
    for prev, cur in zip(joists, joists[1:]):
        for end, label in ((0, "A"), (-1, "B")):
            old, new = prev.crossings[end], cur.crossings[end]
            if old.gid == new.gid:
                continue
            expected = _t_on_line(frame, by_gid[old.gid].line, cur.station)
            jump = abs(new.t - expected)
            if jump > tol:
                events.append(
                    StationEvent(
                        "span_jump",
                        cur.target,
                        cur.station,
                        f"end {label}: outer support changes {old.gid} -> {new.gid} "
                        f"and moves {jump:.3f} along the joist",
                    )
                )
    return events
