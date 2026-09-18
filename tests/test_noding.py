"""
Tests for Phase 2 tolerance extend/trim noding
(``papermodels.geometry.noding``) and its integration into ``GeometryModel``.
"""

import pytest
from shapely.geometry import LineString, box

from papermodels.datatypes.geometry_model import GeometryModel
from papermodels.geometry.noding import node_geometries

TOL = 1e-3
SUPPORT = LineString([(0, 0), (10, 0)])


# --------------------------------------------------------------------------
# node_geometries — module behaviour
# --------------------------------------------------------------------------


def test_extend_endpoint_that_falls_short():
    beam = LineString([(5, 0.0008), (5, 5)])  # 0.0008 gap below the support
    result, report = node_geometries({"S": SUPPORT, "B": beam}, TOL)
    assert result["B"].coords[0] == pytest.approx((5.0, 0.0))
    assert result["B"].distance(SUPPORT) == pytest.approx(0.0, abs=1e-12)
    assert report.converged and report.moves == 1


def test_trim_endpoint_that_overshoots():
    beam = LineString([(5, -0.0007), (5, 5)])  # pokes 0.0007 past the support
    result, report = node_geometries({"S": SUPPORT, "B": beam}, TOL)
    assert result["B"].coords[0] == pytest.approx((5.0, 0.0))
    # No part of the trimmed beam lies below the support.
    assert result["B"].bounds[1] == pytest.approx(0.0, abs=1e-12)
    assert report.converged and report.moves == 1


def test_no_move_when_already_touching():
    beam = LineString([(5, 0), (5, 5)])  # endpoint already on the support
    result, report = node_geometries({"S": SUPPORT, "B": beam}, TOL)
    assert result["B"].coords[0] == (5.0, 0.0)
    assert report.moves == 0 and report.passes == 1 and report.converged


def test_no_move_beyond_tolerance():
    beam = LineString([(5, 0.002), (5, 5)])  # 0.002 gap, tol is 1e-3
    result, report = node_geometries({"S": SUPPORT, "B": beam}, TOL)
    assert result["B"].coords[0] == (5.0, 0.002)
    assert report.moves == 0 and report.converged


def test_snap_onto_polygon_boundary():
    column = box(4, 4, 6, 6)
    beam = LineString([(5, 3.999), (5, 0)])  # tip 0.001 below the column
    result, _ = node_geometries({"C": column, "B": beam}, TOL)
    assert result["B"].coords[0] == pytest.approx((5.0, 4.0))


def test_anchored_endpoint_does_not_move():
    # B's endpoint sits exactly on A but is also within tol of C (a diagonal line
    # passing ~0.0005 from the origin, whose own endpoints are far from
    # everything). Being already anchored to A, B must not be dragged onto C.
    a = LineString([(-10, 0), (10, 0)])
    c = LineString([(-10, -10 + 0.0007), (10, 10 + 0.0007)])
    b = LineString([(0, 0), (0, -20)])
    result, report = node_geometries({"A": a, "C": c, "B": b}, TOL)
    assert list(result["B"].coords) == [(0.0, 0.0), (0.0, -20.0)]
    assert report.moves == 0 and report.converged


def test_fixed_ids_are_targets_but_never_move():
    # A is fixed; B's near endpoint should snap onto A, and A must stay put.
    a = LineString([(0, 0), (0, 10)])
    b = LineString([(0.0008, 5), (5, 5)])
    result, _ = node_geometries({"A": a, "B": b}, TOL, fixed_ids={"A"})
    assert list(result["A"].coords) == [(0.0, 0.0), (0.0, 10.0)]  # unmoved
    assert result["B"].coords[0] == pytest.approx((0.0, 5.0))  # snapped to A


def test_cap_and_warning_on_non_convergence():
    # A move happens in pass 1, but max_passes=1 denies the confirming zero-move
    # pass, so the run reports non-convergence and warns (cap + warning, §9).
    beam = LineString([(5, 0.0008), (5, 5)])
    with pytest.warns(UserWarning, match="did not converge"):
        result, report = node_geometries(
            {"S": SUPPORT, "B": beam}, TOL, max_passes=1
        )
    assert report.passes == 1 and not report.converged


def test_non_convergence_warning_can_be_suppressed():
    beam = LineString([(5, 0.0008), (5, 5)])
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning would raise
        _, report = node_geometries(
            {"S": SUPPORT, "B": beam}, TOL, max_passes=1, suppress_warnings=True
        )
    assert not report.converged


def test_rejects_bad_tol():
    with pytest.raises(ValueError):
        node_geometries({"S": SUPPORT}, 0.0)


# --------------------------------------------------------------------------
# GeometryModel integration
# --------------------------------------------------------------------------


class _Elem:
    def __init__(self, tag, geometry, rank, plane_id=0, reaction_type="point"):
        self.tag = tag
        self.geometry = geometry
        self.rank = rank
        self.plane_id = plane_id
        self.reaction_type = reaction_type


def test_model_noding_bridges_subtolerance_gap():
    # A joist ending 0.0008 short of its support: no crossing without noding,
    # a crossing (and edge) once noding closes the gap.
    els = [
        _Elem("S", LineString([(0, 0), (10, 0)]), rank=1),
        _Elem("J", LineString([(5, 0.0008), (5, 5)]), rank=0),
    ]
    assert GeometryModel.from_elements(els).intersection_edges() == set()
    model = GeometryModel.from_elements(els, noding_abs_tol=TOL)
    assert model.intersection_edges() == {("J", "S")}
    assert model.noding_report.converged


def test_model_noding_holds_walls_fixed():
    # A joist ending just short of a wall (linear polygon). The wall is indexed
    # by its centerline and held fixed; the joist snaps onto the spine.
    wall = _Elem("W", box(4, 0, 6, 10), rank=9, reaction_type="linear")
    joist = _Elem("J", LineString([(5, 10.0008), (5, 20)]), rank=0)
    # centerline of the wall runs x==5, y in [0, 10]; joist tip is 0.0008 above.
    model = GeometryModel.from_elements([wall, joist], noding_abs_tol=TOL)
    assert model.intersection_edges() == {("J", "W")}
    # Wall centerline (the fixed spine) is unchanged by noding.
    assert list(model.geometries["W"].coords) == [(5.0, 0.0), (5.0, 10.0)]


def test_model_without_noding_leaves_report_none():
    els = [_Elem("S", LineString([(0, 0), (10, 0)]), rank=1)]
    assert GeometryModel.from_elements(els).noding_report is None
