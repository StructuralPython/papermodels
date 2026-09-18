"""
Tests for the Phase 1 node-canonicalization prototype
(``papermodels.datatypes.geometry_model``).

Three groups:

1. ``NodeRegistry`` unit behaviour — snapping, bucket-boundary straddling,
   idempotence.
2. ``GeometryModel`` validated against a real fixture (``intersections.pdf``):
   the documented physical crossings must be reproduced.
3. The **jitter test** (design §10) — the acceptance criterion for the whole
   effort: perturbing every input coordinate by 1e-6..1e-9 must leave the
   intersection topology (node set, edge set, incidence structure) invariant.
   Plus idempotence of the build.
"""

import copy
import pathlib
from decimal import Decimal

import numpy as np
import pytest
from shapely import LineString, Point, Polygon, box, transform

from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.geometry_model import (
    GeometryModel,
    NodeRegistry,
    ROLE_BOUNDARY,
    ROLE_ENDPOINT,
    ROLE_INTERIOR,
)

TEST_DATA = pathlib.Path(__file__).parent / "test_data"
QUARTER_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(4)


# --------------------------------------------------------------------------
# NodeRegistry unit tests
# --------------------------------------------------------------------------


def test_node_registry_snaps_within_tolerance():
    reg = NodeRegistry(node_abs_tol=1e-3)
    a = reg.get_or_create((10.0, 20.0))
    # Well within tolerance -> same node, coordinate is "first point wins".
    b = reg.get_or_create((10.0 + 4e-4, 20.0 - 3e-4))
    assert a == b
    assert reg.coord[a] == (10.0, 20.0)
    assert len(reg) == 1


def test_node_registry_distinct_beyond_tolerance():
    reg = NodeRegistry(node_abs_tol=1e-3)
    a = reg.get_or_create((0.0, 0.0))
    b = reg.get_or_create((0.0, 2e-3))  # beyond tol -> new node
    assert a != b
    assert len(reg) == 2


def test_node_registry_snaps_across_bucket_boundary():
    # Two points within tol but on opposite sides of a quantization cell edge.
    # cell size == tol == 0.1: floor(0.099/0.1)=0, floor(0.101/0.1)=1, so the
    # points live in adjacent buckets and only the 9-cell probe reunites them.
    reg = NodeRegistry(node_abs_tol=0.1)
    a = reg.get_or_create((0.099, 0.0))
    b = reg.get_or_create((0.101, 0.0))
    assert a == b
    assert len(reg) == 1


def test_node_registry_snaps_to_nearest():
    reg = NodeRegistry(node_abs_tol=1.0)
    n0 = reg.get_or_create((0.0, 0.0))
    n1 = reg.get_or_create((1.5, 0.0))
    # (0.7, 0) is within tol of BOTH existing nodes; must snap to the nearer.
    assert reg.get_or_create((0.7, 0.0)) == n0
    assert reg.get_or_create((0.9, 0.0)) == n1


def test_node_registry_idempotent():
    reg = NodeRegistry(node_abs_tol=1e-3)
    pts = [(0.0, 0.0), (5.0, 5.0), (5.0 + 2e-4, 5.0), (10.0, 0.0)]
    ids_first = [reg.get_or_create(p) for p in pts]
    coords_first = dict(reg.coord)
    # Re-inserting the same points yields identical ids and coordinates.
    ids_second = [reg.get_or_create(p) for p in pts]
    assert ids_first == ids_second
    assert reg.coord == coords_first


def test_node_registry_accepts_shapely_point():
    reg = NodeRegistry(node_abs_tol=1e-3)
    a = reg.get_or_create(Point(3.0, 4.0))
    b = reg.get_or_create((3.0, 4.0))
    assert a == b


def test_node_registry_rejects_bad_tol():
    with pytest.raises(ValueError):
        NodeRegistry(node_abs_tol=0.0)


# --------------------------------------------------------------------------
# GeometryModel — small synthetic model
# --------------------------------------------------------------------------


class _Elem:
    """Minimal Element stand-in for GeometryModel.from_elements."""

    def __init__(self, tag, geometry, rank, plane_id=0, reaction_type="point"):
        self.tag = tag
        self.geometry = geometry
        self.rank = rank
        self.plane_id = plane_id
        self.reaction_type = reaction_type


def test_model_shares_one_node_for_a_crossing():
    # Two crossing lines and a support line all meeting near one point.
    a = _Elem("A", LineString([(0, 0), (2, 2)]), rank=0)
    b = _Elem("B", LineString([(0, 2), (2, 0)]), rank=1)
    model = GeometryModel.from_elements([a, b])
    # Exactly one canonical node; both geoms incident to it.
    assert len(model.nodes) == 1
    (node_id,) = model.incidence
    incident = {inc.geom_id for inc in model.incident_geoms(node_id)}
    assert incident == {"A", "B"}
    # Direction is by rank: A (0) transfers down to B (1).
    assert model.intersection_edges() == {("A", "B")}


def test_model_respects_plane_separation():
    a = _Elem("A", LineString([(0, 0), (2, 2)]), rank=0, plane_id=0)
    b = _Elem("B", LineString([(0, 2), (2, 0)]), rank=1, plane_id=1)
    model = GeometryModel.from_elements([a, b])
    # Different planes -> no crossing recorded.
    assert model.intersection_edges() == set()
    assert len(model.nodes) == 0


def test_model_wall_indexed_by_centerline():
    # A horizontal joist crossing a vertical wall (linear-reaction polygon).
    joist = _Elem("J", LineString([(-1, 5), (11, 5)]), rank=0)
    wall = _Elem("W", box(4, 0, 6, 10), rank=9, reaction_type="linear")
    model = GeometryModel.from_elements([joist, wall])
    # The wall is indexed by its centerline (x == 5), so the crossing lands on it.
    assert model.geometries["W"].geom_type == "LineString"
    assert model.polygons["W"].geom_type == "Polygon"  # polygon retained
    assert model.intersection_edges() == {("J", "W")}
    (node_id,) = [n for n in model.geom_nodes["J"]]
    (cx, _) = model.nodes.coord[node_id]
    assert cx == pytest.approx(5.0, abs=1e-9)


def test_model_role_classification():
    # Line endpoint landing on a support line's interior.
    support = _Elem("S", LineString([(0, 0), (10, 0)]), rank=1)
    beam = _Elem("B", LineString([(5, 0), (5, 5)]), rank=0)
    model = GeometryModel.from_elements([support, beam])
    (node_id,) = model.incidence
    roles = {inc.geom_id: inc.role for inc in model.incident_geoms(node_id)}
    assert roles["B"] == ROLE_ENDPOINT
    assert roles["S"] == ROLE_INTERIOR


# --------------------------------------------------------------------------
# GeometryModel — validated against the intersections.pdf fixture
# --------------------------------------------------------------------------


def _fixture_elements(name):
    """Raw elements from a fixture PDF (no gravity-frame post-processing)."""
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / name,
        scale=QUARTER_INCH_SCALE,
        process_gravity_frame=False,
    )
    return [graph.nodes[n]["element"] for n in graph.nodes]


# Physical crossings documented by test_geometry_graph.test_intersections_below_above.
_EXPECTED_INTERSECTION_EDGES = {
    ("J0.0", "DB0.0"),
    ("J0.0", "DB0.1"),
    ("J0.1", "DB0.0"),
    ("J0.1", "WT0.0"),
    ("DB0.0", "CT0.0"),
    ("DB0.0", "CT0.1"),
    ("DB0.1", "CT0.2"),
    ("DB0.1", "CT0.3"),
}


def test_model_reproduces_fixture_intersections():
    model = GeometryModel.from_elements(_fixture_elements("intersections.pdf"))
    edges = model.intersection_edges()
    missing = _EXPECTED_INTERSECTION_EDGES - edges
    assert not missing, f"model failed to reproduce crossings: {missing}"


@pytest.mark.parametrize(
    "fixture_name",
    ["intersections.pdf", "collector_extents_walls.pdf"],
)
def test_model_builds_from_fixtures(fixture_name):
    # Every fixture must build without error and produce a non-empty arrangement.
    model = GeometryModel.from_elements(_fixture_elements(fixture_name))
    assert len(model.nodes) > 0
    assert len(model.intersection_edges()) > 0


def test_noding_preserves_fixture_intersections():
    # Enabling Phase 2 noding must not drop the already-clean crossings.
    elements = _fixture_elements("intersections.pdf")
    model = GeometryModel.from_elements(elements, noding_abs_tol=1e-3)
    missing = _EXPECTED_INTERSECTION_EDGES - model.intersection_edges()
    assert not missing, f"noding dropped crossings: {missing}"
    assert model.noding_report is not None and model.noding_report.converged


# --------------------------------------------------------------------------
# Jitter test (acceptance criterion, design §10) + idempotence
# --------------------------------------------------------------------------


def _jitter_geometry(geom, rng, mag):
    """
    Perturb every vertex by up to ``mag`` in each axis.  Identical input
    coordinates (e.g. a polygon ring's shared first/last vertex) get the same
    offset, so rings stay closed and shared vertices move together.
    """
    cache = {}

    def fn(coords):
        out = np.empty_like(coords)
        for i, (x, y) in enumerate(coords):
            key = (round(float(x), 9), round(float(y), 9))
            if key not in cache:
                cache[key] = rng.uniform(-mag, mag, size=2)
            out[i] = np.array([x, y], dtype=float) + cache[key]
        return out

    return transform(geom, fn)


@pytest.mark.parametrize("magnitude", [1e-6, 1e-7, 1e-8, 1e-9])
@pytest.mark.parametrize("noding_abs_tol", [None, 1e-3])
def test_topology_invariant_under_jitter(magnitude, noding_abs_tol):
    elements = _fixture_elements("intersections.pdf")
    base = GeometryModel.from_elements(
        elements, noding_abs_tol=noding_abs_tol
    ).topology_signature()
    rng = np.random.default_rng(1234)
    for _ in range(10):
        jittered = copy.deepcopy(elements)
        for element in jittered:
            element.geometry = _jitter_geometry(element.geometry, rng, magnitude)
        sig = GeometryModel.from_elements(
            jittered, noding_abs_tol=noding_abs_tol
        ).topology_signature()
        assert sig["n_nodes"] == base["n_nodes"]
        assert sig["edges"] == base["edges"]
        assert sig["incidence"] == base["incidence"]


def test_build_is_idempotent():
    elements = _fixture_elements("intersections.pdf")
    m1 = GeometryModel.from_elements(elements)
    m2 = GeometryModel.from_elements(elements)
    # Identical node ids and canonical coordinates both times.
    assert m1.nodes.coord == m2.nodes.coord
    assert m1.topology_signature() == m2.topology_signature()


# --------------------------------------------------------------------------
# Support queries and generated geometry (used by JoistArrayModel)
# --------------------------------------------------------------------------


def _support_scene():
    wall = _Elem(
        "WT", box(0, -0.25, 10, 0.25), rank=2, reaction_type="linear"
    )  # centerline y=0
    beam = _Elem("FB", LineString([(0, 5), (10, 5)]), rank=1)
    column = _Elem("CT", box(4.9, 2.4, 5.1, 2.6), rank=3, reaction_type="point")
    other_plane = _Elem("FB_P1", LineString([(0, 3), (10, 3)]), rank=1, plane_id=1)
    joist = _Elem("J", LineString([(5, -1), (5, 6)]), rank=0)
    return [wall, beam, column, other_plane, joist]


def test_query_supports_filters_plane_rank_and_columns():
    model = GeometryModel.from_elements(_support_scene())
    region = box(4, -1, 6, 6)
    assert model.query_supports(region, plane=0, rank=0) == ["WT", "FB"]
    # Rank filter: nothing above rank 2 bears a joist here (column excluded)
    assert model.query_supports(region, plane=0, rank=2) == []
    assert model.query_supports(region, plane=1, rank=0) == ["FB_P1"]


def test_query_supports_uses_wall_centerline_not_polygon():
    model = GeometryModel.from_elements(_support_scene())
    # Region touches the wall polygon's face (y=0.2) but not its centerline (y=0)
    assert "WT" not in model.query_supports(box(4, 0.1, 6, 0.2), plane=0, rank=0)
    assert model.source_geometry("WT").geom_type == "Polygon"
    assert model.geometries["WT"].geom_type == "LineString"


def test_seed_points_are_reused_first_point_wins():
    seed = (5.0 + 4e-7, 5.0)  # within node_abs_tol of the J/FB crossing
    model = GeometryModel.from_elements(
        _support_scene(), node_abs_tol=1e-6, seed_points=[seed]
    )
    (node,) = model.intersections_below("J")["FB"]
    assert model.nodes.coord[node] == seed


def test_shared_registry_is_used():
    reg = NodeRegistry(1e-6)
    model = GeometryModel.from_elements(_support_scene(), node_abs_tol=1e-6, nodes=reg)
    assert model.nodes is reg and len(reg) > 0


def test_add_geometry_records_incidence_and_reindexes():
    model = GeometryModel.from_elements(_support_scene(), node_abs_tol=1e-6)
    before = model.query_supports(box(7, -1, 8, 6), plane=0, rank=0)
    n_wall = model.nodes.get_or_create((7.5, 0.0))
    n_beam = model.nodes.get_or_create((7.5, 5.0))
    model.add_geometry(
        "J-1",
        LineString([(7.5, 0.0), (7.5, 5.0)]),
        rank=0,
        plane=0,
        crossings={n_wall: "WT", n_beam: "FB"},
    )
    assert model.intersections_below("J-1") == {"WT": {n_wall}, "FB": {n_beam}}
    assert "J-1" in model.intersections_above("WT")
    # Generated joists (rank 0) never become supports; index rebuilt lazily
    assert model.query_supports(box(7, -1, 8, 6), plane=0, rank=0) == before
    assert "J-1" in model._index_geom_ids
    with pytest.raises(ValueError):
        model.add_geometry("J-1", LineString([(0, 0), (1, 1)]), rank=0, plane=0)


def test_graph_geometry_model_holds_every_intersection_node():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "intersections.pdf", scale=QUARTER_INCH_SCALE
    )
    model = graph.geometry_model
    assert model is not None
    assert set(model.geometries) == set(graph.nodes)
    coords = set(model.nodes.coord.values())
    for n in graph.nodes:
        for ib in graph.nodes[n]["element"].intersections_below or []:
            region = ib.intersecting_region
            if region.geom_type == "Point":
                assert (region.x, region.y) in coords


def test_remove_generated_replaces_previous_generation():
    model = GeometryModel.from_elements(_support_scene(), node_abs_tol=1e-6)
    edges_before = model.intersection_edges()
    n_wall = model.nodes.get_or_create((7.5, 0.0))
    n_beam = model.nodes.get_or_create((7.5, 5.0))
    model.add_geometry(
        "J-0",
        LineString([(7.5, 0.0), (7.5, 5.0)]),
        rank=0,
        plane=0,
        crossings={n_wall: "WT", n_beam: "FB"},
        parent="J",
    )
    assert model.generated == {"J": ["J-0"]}
    model.remove_generated("J")
    assert "J-0" not in model.geometries and model.generated == {}
    assert model.intersection_edges() == edges_before
    assert n_wall not in model.geom_nodes["WT"]
    # Regenerating under the same id is allowed after removal
    model.add_geometry("J-0", LineString([(7.5, 0.0), (7.5, 5.0)]), rank=0, plane=0)
