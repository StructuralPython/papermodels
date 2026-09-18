"""
Tests for JoistArrayModel (papermodels.datatypes.joist_models) built on the
station engine: synthetic elements outside a graph, and acceptance checks on
the real fixtures inside a GeometryGraph.
"""

import copy
import math
import pathlib
import warnings
from decimal import Decimal

import numpy as np
import pytest
from shapely import LineString, Point, box, unary_union
from shapely import affinity as aff

from papermodels.datatypes.element import Element, Intersection
from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.joist_models import JoistArrayModel

TEST_DATA = pathlib.Path(__file__).parent / "test_data"
QUARTER = Decimal(1) / Decimal(72) * Decimal(4)
EIGHTH = Decimal(1) / Decimal(72) * Decimal(8)

FIXTURES = [
    ("collector_extents.pdf", QUARTER, {}),
    ("sketch_to_scale.pdf", QUARTER, {}),
    ("horiz_extents.pdf", QUARTER, {}),
    ("resi_dormers.pdf", QUARTER, {"cantilever_abs_tol": 0.3}),
    ("frame_collectors_transfers.pdf", EIGHTH, {}),
    ("collector_extents_walls.pdf", QUARTER, {}),
    ("intersections.pdf", QUARTER, {}),
]


# --------------------------------------------------------------------------
# Synthetic element (standalone model)
# --------------------------------------------------------------------------


def _element(proto, supports, tag="J0.0"):
    """A rank-0 prototype with (tag, geometry, reaction_type) supports below."""
    return Element(
        geometry=proto,
        tag=tag,
        rank=0,
        intersections_below=[
            Intersection(Point(0, 0), geom, sup_tag, other_reaction_type=rxn)
            for sup_tag, geom, rxn in supports
        ],
        intersections_above=[],
        correspondents_above=[],
        correspondents_below=[],
        plane_id=0,
    )


WALL = ("WT0.0", box(0, -0.25, 10, 0.25), "linear")  # centerline y=0
BEAM = ("FB0.0", LineString([(0, 10), (10, 10)]), "point")


def test_standalone_array_bears_on_wall_centerline_and_beam():
    e = _element(LineString([(5, -0.1), (5, 10.5)]), [WALL, BEAM])
    model = JoistArrayModel(e, spacing=1.0, cantilever_tolerance=0.2)
    new = model()
    subs = new.subelements
    assert [s.tag for s in subs] == [f"J0.0-{i}" for i in range(11)]
    for sub in subs:
        ys = sorted(c[1] for c in sub.geometry.coords)
        # Prototype stops 0.1 short of the wall centerline (inside the wall):
        # the joist still reaches its support node. 0.5 cantilever past the beam.
        assert math.isclose(ys[0], 0.0, abs_tol=1e-12)
        assert math.isclose(ys[1], 10.5, abs_tol=1e-12)
        by_tag = {ib.other_tag: ib for ib in sub.intersections_below}
        assert set(by_tag) == {"WT0.0", "FB0.0"}
        assert by_tag["WT0.0"].other_geometry.geom_type == "Polygon"
        assert by_tag["WT0.0"].other_reaction_type == "linear"
        for ib in sub.intersections_below:
            assert ib.intersecting_region.geom_type == "Point"
            xy = (ib.intersecting_region.x, ib.intersecting_region.y)
            assert xy in set(model.geometry_model.nodes.coord.values())
            lo, hi = ib.other_extents
            assert 0 <= lo <= hi <= 10
    bands = [s.trib_area for s in subs]
    assert math.isclose(sum(b.area for b in bands), model.result.region.area)
    assert subs[0].trib_area.area == pytest.approx(0.5 * 10.5, abs=1e-5)


@pytest.mark.parametrize("dx", [1e-4, -1e-4, 1e-9, -1e-9])
@pytest.mark.parametrize("reverse", [False, True])
def test_near_vertical_prototype_never_flips_cantilevers(dx, reverse):
    """
    The reported bug: with the positive-x bias, a prototype tilted just past
    vertical projected its cantilevers into the span.
    """
    coords = [(5, -0.5), (5 + dx, 11.0)]  # 0.5 cantilever below, 1.0 above
    if reverse:
        coords = coords[::-1]
    beam_bottom = ("FB0.1", LineString([(0, 0), (10, 0)]), "point")
    e = _element(LineString(coords), [beam_bottom, BEAM])
    subs = JoistArrayModel(e, spacing=1.0, cantilever_tolerance=0.1)().subelements
    assert len(subs) == 11
    xs = [s.geometry.centroid.x for s in subs]
    assert xs == sorted(xs)  # array runs left to right
    for sub in subs:
        ys = sorted(c[1] for c in sub.geometry.coords)
        assert ys[0] == pytest.approx(-0.5, abs=1e-6)
        assert ys[1] == pytest.approx(11.0, abs=1e-6)


def test_prototype_on_one_support_is_a_clear_error():
    e = _element(LineString([(5, -0.5), (5, 5)]), [BEAM])
    with pytest.raises(Exception, match="at least two"):
        JoistArrayModel(e)


# --------------------------------------------------------------------------
# Fixtures (inside a GeometryGraph)
# --------------------------------------------------------------------------


def _assigned(name, scale, kw):
    graph = GeometryGraph.from_pdf_file(TEST_DATA / name, scale=scale, **kw)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
    return graph, [str(w.message) for w in caught]


@pytest.mark.parametrize("name,scale,kw", FIXTURES)
def test_fixture_arrays_are_complete_and_canonical(name, scale, kw):
    graph, messages = _assigned(name, scale, kw)
    assert not messages, messages
    assert not graph.omitted
    coords = set(graph.geometry_model.nodes.coord.values())
    for node in graph.collector_elements:
        element = graph.nodes[node]["element"]
        for sub in element.subelements or []:
            assert len(sub.intersections_below) >= 2, sub.tag
            for ib in sub.intersections_below:
                p = ib.intersecting_region
                assert (p.x, p.y) in coords
            # Generated joists are written back into the model
            assert sub.tag in graph.geometry_model.geometries


def test_extent_line_array_follows_outer_support_lines():
    """
    collector_extents.pdf: SJ0.0 is drawn once with an extent line across a
    run of bottom supports (WT0.1, FB0.0, WT0.0) and top supports (FB0.2,
    WT0.2, FB0.3), with WT0.3 as an intermediate support that comes and goes.
    """
    graph, _ = _assigned("collector_extents.pdf", QUARTER, {})
    subs = graph.nodes["SJ0.0"]["element"].subelements
    bottom, top = {"WT0.1", "FB0.0", "WT0.0"}, {"FB0.2", "WT0.2", "FB0.3"}
    saw_intermediate = False
    for sub in subs:
        tags = {ib.other_tag for ib in sub.intersections_below}
        assert len(tags & bottom) == 1 and len(tags & top) == 1, (sub.tag, tags)
        saw_intermediate |= "WT0.3" in tags
    assert saw_intermediate
    # The gravity-frame processing aligns FB0.2/FB0.3 to column centroids,
    # leaving gaps on the top line at x~10.65-11.06 and x~15.56-16.0. The
    # nominal joists there (x=10.794, 15.794) move to the nearest station
    # with a full support pattern instead of being built on one support.
    xs = [round(s.geometry.centroid.x, 3) for s in subs]
    assert 10.794 not in xs and 15.794 not in xs
    assert 11.064 in xs and 16.001 in xs
    # Trib areas tile the extent region: no gaps next to the moved joists
    bands = [s.trib_area for s in subs]
    union = unary_union(bands)
    assert sum(b.area for b in bands) == pytest.approx(union.area, rel=1e-9)
    assert union.geom_type == "Polygon"  # contiguous


def test_reassigning_behaviour_replaces_generated_joists():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents.pdf", scale=QUARTER
    )
    graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
    first = list(graph.geometry_model.generated["SJ0.0"])
    graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
    assert graph.geometry_model.generated["SJ0.0"] == first
    generated = [g for kids in graph.geometry_model.generated.values() for g in kids]
    assert len(generated) == len(set(generated))


def test_abrupt_span_change_warns_with_tag_and_advice():
    # Top support steps 2.0 further out halfway along the extent: a span jump.
    top_1 = ("FB0.1", LineString([(0, 10), (5.5, 10)]), "point")
    top_2 = ("FB0.2", LineString([(5.5, 12), (10, 12)]), "point")
    e = _element(LineString([(2, -0.5), (2, 10.5)]), [WALL, top_1, top_2])
    e.extent_line = LineString([(0, 5), (10, 5)])
    with pytest.warns(UserWarning, match=r"Joist array J0\.0: abrupt change in span"):
        model = JoistArrayModel(e, spacing=1.0)
    assert [ev.kind for ev in model.events] == ["span_jump"]
    # Generation continues: every joist still bears on the wall and a top beam
    for joist in model.result.joists:
        assert len(joist.crossings) == 2
