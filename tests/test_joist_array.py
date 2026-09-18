"""
Unit tests for the joist-array station engine (papermodels.geometry.joist_array).

Scenes are built in a local frame (joists vertical, supports horizontal) and then
rotated, so every expectation can be written once and checked at any angle.
"""

import math

import numpy as np
import pytest
from shapely import LineString, Point, Polygon, box, unary_union
from shapely import affinity as aff

from papermodels.datatypes.geometry_model import NodeRegistry
from papermodels.geometry import joist_array as ja

SPACING = 1.0
MB = 1e-3  # min_bearing


def S(gid, coords):
    return ja.Support(gid, LineString(coords))


def assert_tiles(result, region=None):
    region = result.region if region is None else region
    bands = [j.trib_area for j in result.joists]
    assert math.isclose(sum(b.area for b in bands), region.area, rel_tol=1e-9)
    assert unary_union(bands).symmetric_difference(region).area < 1e-9
    for a, b in zip(bands, bands[1:]):
        assert a.intersection(b).area < 1e-9


def fixed_array(proto, supports, s_range=None, region=None, **kw):
    frame = ja.ArrayFrame.from_prototype(proto)
    cant_a, cant_b = ja.prototype_cantilevers(frame, proto, supports)
    a, b = ja.outer_supports(frame, proto, supports)
    s_range = s_range or ja.common_s_range(frame, supports)
    region = region or ja.fixed_region(frame, a.line, b.line, s_range, cant_a, cant_b)
    return ja.build_array(
        frame,
        region,
        supports,
        s_range,
        ja.MODE_FIXED,
        SPACING,
        cant_a=cant_a,
        cant_b=cant_b,
        min_bearing=MB,
        registry=NodeRegistry(1e-6),
        **kw,
    )


# --------------------------------------------------------------------------
# Frame
# --------------------------------------------------------------------------


def test_frame_dominant_axis_rule():
    vert = ja.ArrayFrame.from_prototype(LineString([(5, 0), (5, 10)]))
    np.testing.assert_allclose(vert.n, [1, 0], atol=1e-12)  # array runs left->right
    horiz = ja.ArrayFrame.from_prototype(LineString([(0, 5), (10, 5)]))
    np.testing.assert_allclose(horiz.n, [0, 1], atol=1e-12)  # bottom->top


@pytest.mark.parametrize("dx", [1e-4, -1e-4, 1e-12, -1e-12])
def test_frame_is_stable_either_side_of_vertical(dx):
    """The positive-x bias flips here; the dominant-axis frame must not."""
    f = ja.ArrayFrame.from_prototype(LineString([(5, 0), (5 + dx, 10)]))
    assert f.n[0] > 0.99


def test_frame_ignores_prototype_direction():
    a = ja.ArrayFrame.from_prototype(LineString([(0, 0), (3, 7)]))
    b = ja.ArrayFrame.from_prototype(LineString([(3, 7), (0, 0)]))
    np.testing.assert_allclose(a.n, b.n)
    np.testing.assert_allclose(a.u, b.u)
    np.testing.assert_allclose(a.origin, b.origin)


# --------------------------------------------------------------------------
# Target stations
# --------------------------------------------------------------------------


@pytest.mark.parametrize("length", [10.0, 10.3, 9.95, 0.4, 3.0])
@pytest.mark.parametrize("offset", [0.0, 0.25])
def test_target_stations_never_exceed_spacing(length, offset):
    st = ja.target_stations(0.0, length, 1.0, initial_offset=offset)
    assert st[0] == 0.0 and st[-1] == length
    gaps = np.diff(st)
    assert gaps.max() <= 1.0 + 1e-12
    assert (gaps > 0).all()


def test_target_stations_options():
    assert ja.target_stations(0, 3, 1) == [0, 1, 2, 3]
    assert ja.target_stations(0, 3, 1, joist_at_end=False) == [0, 1, 2]
    assert ja.target_stations(0, 3, 1, joist_at_start=False) == [1, 2, 3]
    assert ja.target_stations(0, 3, 1, initial_offset=0.5) == [0, 0.5, 1.5, 2.5, 3]
    with pytest.raises(ValueError):
        ja.target_stations(0, 3, 0)


# --------------------------------------------------------------------------
# Orthogonal array, canonical crossings, tiling
# --------------------------------------------------------------------------

PROTO = LineString([(5, -0.5), (5, 11.0)])  # cantilevers: 0.5 below, 1.0 above
BOTTOM = S("BOT", [(0, 0), (10, 0)])
TOP = S("TOP", [(0, 10), (10, 10)])


def test_orthogonal_array():
    r = fixed_array(PROTO, [BOTTOM, TOP])
    assert [j.station for j in r.joists] == [float(s) for s in range(-5, 6)]
    for j in r.joists:
        x = j.station + 5
        (x0, y0), (x1, y1) = j.geometry.coords
        assert {round(y0, 12), round(y1, 12)} == {-0.5, 11.0}
        assert math.isclose(x0, x, abs_tol=1e-12) and math.isclose(x1, x, abs_tol=1e-12)
        assert {c.gid for c in j.crossings} == {"BOT", "TOP"}
    assert not r.events
    assert_tiles(r)


def test_crossings_are_registry_nodes_and_joist_ends_are_nodes():
    reg = NodeRegistry(1e-6)
    proto = LineString([(5, 0), (5, 10)])  # no cantilevers: ends are nodes
    frame = ja.ArrayFrame.from_prototype(proto)
    sups = [BOTTOM, TOP]
    a, b = ja.outer_supports(frame, proto, sups)
    s_range = ja.common_s_range(frame, sups)
    region = ja.fixed_region(frame, a.line, b.line, s_range, 0, 0)
    r = ja.build_array(frame, region, sups, s_range, ja.MODE_FIXED, 1.0, registry=reg)
    for j in r.joists:
        for c in j.crossings:
            assert reg.coord[c.node_id] == c.xy
        assert set(j.geometry.coords) == {c.xy for c in j.crossings}


# --------------------------------------------------------------------------
# Gap resolution
# --------------------------------------------------------------------------


def test_gap_under_target_relocates_to_nearest_valid_station():
    # Bottom support has a gap from x=4.8 to 5.3; the target at x=5 falls in it.
    bottom = [S("B1", [(0, 0), (4.8, 0)]), S("B2", [(5.3, 0), (10, 0)])]
    proto = LineString([(2, -0.5), (2, 11.0)])
    frame = ja.ArrayFrame.from_prototype(proto)
    sups = bottom + [TOP]
    s_range = (frame.ts((0, 0))[1], frame.ts((10, 0))[1])
    region = ja.fixed_region(frame, TOP.line, bottom[0].line, s_range, 1.0, 0.5)
    r = ja.build_array(
        frame,
        region,
        sups,
        s_range,
        ja.MODE_FIXED,
        1.0,
        cant_a=1.0,
        cant_b=0.5,
        min_bearing=MB,
    )
    xs = [round(frame.xy(0, j.station)[0], 9) for j in r.joists]
    assert xs == [0, 1, 2, 3, 4, round(4.8 - MB, 9), 6, 7, 8, 9, 10]
    (ev,) = [e for e in r.events if e.kind == "relocated"]
    assert math.isclose(frame.xy(0, ev.target)[0], 5.0, abs_tol=1e-9)
    assert all(len(j.crossings) == 2 for j in r.joists)
    assert_tiles(r)


def test_intermediate_support_that_comes_and_goes_is_not_a_gap():
    mid = S("MID", [(2.5, 5), (7.5, 5)])
    r = fixed_array(PROTO, [BOTTOM, TOP], s_range=(-5, 5))
    r2 = fixed_array(
        PROTO,
        [BOTTOM, TOP, mid],
        s_range=(-5, 5),
        region=r.region,
    )
    assert [j.station for j in r2.joists] == [j.station for j in r.joists]
    assert not [e for e in r2.events if e.kind == "relocated"]
    counts = [len(j.crossings) for j in r2.joists]
    assert counts == [2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2]


def test_support_count_differing_from_both_neighbours_relocates():
    # A short intermediate beam under only the x=5 joist: 3 supports vs 2 and 2.
    short = S("SHORT", [(4.6, 5), (5.4, 5)])
    base = fixed_array(PROTO, [BOTTOM, TOP])
    r = fixed_array(PROTO, [BOTTOM, TOP, short], s_range=(-5, 5), region=base.region)
    (ev,) = [e for e in r.events if e.kind == "relocated"]
    assert math.isclose(ev.target, 0.0, abs_tol=1e-12)
    assert math.isclose(abs(ev.final), 0.4 + MB, abs_tol=1e-9)
    assert all(len(j.crossings) == 2 for j in r.joists)


def test_no_valid_station_in_window_drops_joist():
    bottom = [S("B1", [(0, 0), (3, 0)]), S("B2", [(7, 0), (10, 0)])]
    proto = LineString([(1, -0.5), (1, 11.0)])
    frame = ja.ArrayFrame.from_prototype(proto)
    s_range = (frame.ts((0, 0))[1], frame.ts((10, 0))[1])
    region = ja.fixed_region(frame, TOP.line, bottom[0].line, s_range, 1.0, 0.5)
    r = ja.build_array(
        frame,
        region,
        bottom + [TOP],
        s_range,
        ja.MODE_FIXED,
        1.0,
        cant_a=1.0,
        cant_b=0.5,
        min_bearing=MB,
    )
    xs = [round(frame.xy(0, j.station)[0], 9) for j in r.joists]
    # x=3 is exactly B1's end (zero bearing), so it moves in by min_bearing.
    assert xs == [0, 1, 2, round(3 - MB, 9), 7, 8, 9, 10]
    dropped = [
        round(frame.xy(0, e.target)[0], 9) for e in r.events if e.kind == "dropped"
    ]
    assert dropped == [4, 5, 6]
    assert_tiles(r)


# --------------------------------------------------------------------------
# Non-orthogonal supports, containers, span jumps
# --------------------------------------------------------------------------


def test_diagonal_support_varies_backspan_with_constant_cantilevers():
    diag = S("DIAG", [(0, 2), (10, 12)])  # y = x + 2
    proto = LineString([(5, -0.5), (5, 8.0)])  # 0.5 below BOTTOM, 1.0 above DIAG
    r = fixed_array(proto, [BOTTOM, diag])
    assert len(r.joists) == 11
    for j in r.joists:
        x = frame_x = j.geometry.coords[0][0]
        ys = sorted(c[1] for c in j.geometry.coords)
        assert math.isclose(ys[0], -0.5, abs_tol=1e-9)
        assert math.isclose(ys[1], x + 2 + 1.0, abs_tol=1e-9)
    assert_tiles(r)


def test_container_clips_joists_so_cantilevers_vary():
    # Container bottom edge y=-1; top edge rises from y=11 (x=0) to y=13 (x=10).
    container = Polygon([(0, -1), (10, -1), (10, 13), (0, 11)])
    proto = LineString([(5, -1), (5, 12)])
    frame = ja.ArrayFrame.from_prototype(proto)
    s_range = frame.s_range(container)
    r = ja.build_array(
        frame,
        container,
        [BOTTOM, TOP],
        s_range,
        ja.MODE_CONTAINER,
        1.0,
        registry=NodeRegistry(1e-6),
    )
    assert len(r.joists) == 11
    for j in r.joists:
        x = j.geometry.coords[0][0]
        ys = sorted(c[1] for c in j.geometry.coords)
        assert math.isclose(ys[0], -1.0, abs_tol=1e-9)
        assert math.isclose(ys[1], 11 + 0.2 * x, abs_tol=1e-9)
    assert_tiles(r)


def test_container_short_of_support_extends_joist_to_support_node():
    container = box(0, -1, 10, 9.9)  # stops 0.1 short of TOP's centerline... but
    # TOP must still be inside the region to be present; widen to include it.
    container = Polygon([(0, -1), (10, -1), (10, 10.0), (0, 10.0)])
    proto = LineString([(5, -1), (5, 10)])
    frame = ja.ArrayFrame.from_prototype(proto)
    r = ja.build_array(
        frame,
        container,
        [BOTTOM, TOP],
        frame.s_range(container),
        ja.MODE_CONTAINER,
        1.0,
        registry=NodeRegistry(1e-6),
    )
    for j in r.joists:
        assert math.isclose(max(c[1] for c in j.geometry.coords), 10.0, abs_tol=1e-12)


def test_span_jump_warns_only_for_abrupt_changes():
    step = [S("T1", [(0, 10), (5.5, 10)]), S("T2", [(5.5, 12), (10, 12)])]
    proto = LineString([(2, 0), (2, 10)])
    frame = ja.ArrayFrame.from_prototype(proto)
    s_range = (frame.ts((0, 0))[1], frame.ts((10, 0))[1])
    region = box(0, -1, 10, 13)
    r = ja.build_array(frame, region, [BOTTOM] + step, s_range, ja.MODE_FIXED, 1.0)
    jumps = [e for e in r.events if e.kind == "span_jump"]
    assert len(jumps) == 1 and "T1 -> T2" in jumps[0].detail

    near = [S("T1", [(0, 10), (5.5, 10)]), S("T2", [(5.5, 10.1), (10, 10.1)])]
    r = ja.build_array(frame, region, [BOTTOM] + near, s_range, ja.MODE_FIXED, 1.0)
    assert not [e for e in r.events if e.kind == "span_jump"]


def test_clip_support_to_region_follows_centerline_even_when_narrow():
    wall_centerline = LineString([(0, 0), (10, 0)])
    narrow = box(4.9, -1, 5.1, 1)  # narrower than a typical wall is thick
    reg = NodeRegistry(1e-6)
    piece = ja.clip_support_to_region(wall_centerline, narrow, reg)
    assert piece.geom_type == "LineString"
    (x0, y0), (x1, y1) = piece.coords
    assert y0 == y1 == 0.0 and {round(x0, 9), round(x1, 9)} == {4.9, 5.1}
    assert all(tuple(c) in set(reg.coord.values()) for c in piece.coords)
    assert ja.clip_support_to_region(wall_centerline, box(0, 5, 1, 6)) is None


# --------------------------------------------------------------------------
# Direction: cantilevers never flip, array order never reverses
# --------------------------------------------------------------------------

DIR_ANGLES = list(range(0, 360, 15)) + [89.9, 90.1, 179.9, 180.1, 269.9, 270.1]


def _rotated_scene(angle, reverse=False, jitter=0.0, seed=0):
    rng = np.random.default_rng(seed)

    def rot(g):
        g = aff.rotate(g, angle, origin=(0, 0))
        if jitter:
            g = aff.translate(g, *rng.uniform(-jitter, jitter, 2))
        return g

    proto = rot(PROTO)
    if reverse:
        proto = LineString(list(proto.coords)[::-1])
    sups = [ja.Support(s.gid, rot(s.line)) for s in (BOTTOM, TOP)]
    return proto, sups


@pytest.mark.parametrize("angle", DIR_ANGLES)
def test_cantilevers_keep_their_side_and_length_at_any_angle(angle):
    proto, sups = _rotated_scene(angle)
    by_gid = {s.gid: s.line for s in sups}
    r = fixed_array(proto, sups)
    assert len(r.joists) == 11
    for j in r.joists:
        ends = [Point(c) for c in j.geometry.coords]
        cross = {c.gid: Point(c.xy) for c in j.crossings}
        bot, top = cross["BOT"], cross["TOP"]
        span = np.subtract(top.coords[0], bot.coords[0])
        # The end nearest each support lies OUTSIDE the span, by the drawn amount.
        bot_end = min(ends, key=lambda e: e.distance(bot))
        top_end = min(ends, key=lambda e: e.distance(top))
        assert np.dot(np.subtract(bot_end.coords[0], bot.coords[0]), span) < 0
        assert np.dot(np.subtract(top_end.coords[0], top.coords[0]), span) > 0
        assert math.isclose(bot_end.distance(bot), 0.5, abs_tol=1e-9)
        assert math.isclose(top_end.distance(top), 1.0, abs_tol=1e-9)


@pytest.mark.parametrize("angle", DIR_ANGLES)
def test_array_is_invariant_to_prototype_direction_and_jitter(angle):
    base = fixed_array(*_rotated_scene(angle))
    rev = fixed_array(*_rotated_scene(angle, reverse=True))
    assert [
        j.geometry.equals_exact(k.geometry, 1e-9)
        for j, k in zip(base.joists, rev.joists)
    ] == [True] * 11
    if angle % 90 in (45,):  # the dominant-axis rule's own discontinuity
        return
    for seed in range(3):
        jit = fixed_array(*_rotated_scene(angle, jitter=1e-6, seed=seed))
        assert len(jit.joists) == len(base.joists)
        for j, k in zip(base.joists, jit.joists):
            assert j.geometry.hausdorff_distance(k.geometry) < 1e-5
