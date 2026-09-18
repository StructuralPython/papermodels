import os
from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt
from papermodels.paper import annotations as an
from papermodels.paper import pdf
from papermodels.datatypes.joist_models import JoistArrayModel
from papermodels.datatypes.element import create_element_filter
import numpy as np
import numpy.testing as npt
from pytest import fixture
import pytest
from pytest_check.context_manager import check
from shapely import Polygon, box, Point, transform
import copy
import pathlib
import fixtures
from decimal import Decimal

EIGHTTH_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(8)
QUARTER_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(4)

TEST_DATA = pathlib.Path(__file__).parent / "test_data"


@fixture()
def load_frame_collectors_transfers():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "frame_collectors_transfers.pdf",
        scale=EIGHTTH_INCH_SCALE,
    )
    return graph


@fixture()
def load_resi_dormers():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "resi_dormers.pdf", scale=QUARTER_INCH_SCALE, cantilever_abs_tol=0.3
    )
    return graph


@fixture()
def load_many_correspondents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "many_correspondents.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def load_collector_extents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents.pdf",
        scale=EIGHTTH_INCH_SCALE,
    )
    return graph


@fixture()
def load_collector_extents_walls():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents_walls.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def load_intersections():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "intersections.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


def test_load_frame_collectors_transfers(load_frame_collectors_transfers):
    assert load_frame_collectors_transfers


def test_collector_assignment_frame_collectors_transfers(
    load_frame_collectors_transfers,
):
    graph = load_frame_collectors_transfers
    steel_joist_arrays = create_element_filter(element_types=["SJ"])
    user_designated_joists = create_element_filter(
        user_defined={"collector behaviour": "array"}
    )
    graph.assign_collector_behaviour(
        JoistArrayModel, spacing=4.0
    )  # Assign all collectors a coarse array
    coarse = len(graph.nodes["SJ0.0"]["element"].subelements)
    graph.assign_collector_behaviour(
        JoistArrayModel, steel_joist_arrays, spacing=1.0
    )  # Re-assign steel joists a finer array
    graph.assign_collector_behaviour(
        JoistArrayModel, user_designated_joists, spacing=1.0
    )
    les = graph.create_loaded_elements()
    # The later, filtered assignments replace the earlier one...
    assert "SJ0.0-9" in les
    assert len(graph.nodes["SJ0.0"]["element"].subelements) > coarse
    assert "WJ0.1-9" in les
    # ...and elements outside the filters keep the first assignment
    assert "WJ0.0-0" in les and "WJ0.0-9" not in les


def test_load_collector_extents(load_collector_extents):
    assert load_collector_extents


def test_load_resi_dormers(load_resi_dormers):
    assert load_resi_dormers


def test_load_collector_extent_walls(load_collector_extents_walls):
    assert load_collector_extents_walls


def test_load_intersections(load_intersections):
    assert load_intersections


def test_resi_dormers_array(load_resi_dormers):
    graph = load_resi_dormers
    # Assign collector behaviour using filter functions
    roof_joist_filter = create_element_filter(element_types=["RJ"])
    all_other_joists_filter = create_element_filter(exclude_element_types=["RJ"])

    # First, assign the default behaviour to everything else
    graph.assign_collector_behaviour(
        JoistArrayModel, filter_function=all_other_joists_filter, spacing=2.0
    )

    # Then assign the special cases. These will overwrite the previously
    # set behaviours for elements that pass the filter
    graph.assign_collector_behaviour(
        JoistArrayModel, filter_function=roof_joist_filter, spacing=1.0
    )

    # Check that joist arrays are created and that their length varies
    les = graph.create_loaded_elements()
    assert les["FB0.3"].model()["loads"]["point_loads"]
    fb03_pl = les["FB0.3"].model()["loads"]["point_loads"]
    # Every joist in the triangular array reaches FB0.3, including the long one
    # at the open end; the zero-backspan joist at the apex is not generated.
    assert len(fb03_pl) == 6
    assert not graph.omitted
    rj_lengths = [
        les[f"RJ0.0-{idx}"].model()["element_attributes"]["length"] for idx in range(6)
    ]
    assert rj_lengths[0] == 5.707
    assert rj_lengths[-1] == 0.864
    # Joists shorten linearly toward the apex (FB0.3 is at 45 degrees; 1.0 spacing)
    steps = [a - b for a, b in zip(rj_lengths, rj_lengths[1:])]
    assert max(steps) - min(steps) < 2e-3


def test_many_correspondents(load_many_correspondents):
    graph = load_many_correspondents
    assert len(graph.nodes["WT2.0"]["element"].correspondents_below) == 1
    assert len(graph.nodes["WB1.1"]["element"].correspondents_above) == 1
    assert graph.nodes["WT2.0"]["element"].correspondents_below[0].other_tag == "WB1.1"
    assert graph.nodes["WB1.1"]["element"].correspondents_above[0].other_tag == "WT2.0"


def test_plot_connectivity(load_collector_extents, capsys):
    # Not currently testing for correct SVG output because
    if "CI" in os.environ:
        assert True
    else:
        graph = load_collector_extents
        graph.plot_connectivity()  # Should write bytes to stdout
        captured = capsys.readouterr()
        assert captured.out is not None


def test_intersections_below_above(load_collector_extents_walls, load_intersections):
    graph_inters = load_intersections
    graph_walls = load_collector_extents_walls

    j0 = graph_walls.nodes["J0.0"]["element"]
    j0_below_tags = [ib.other_tag for ib in j0.intersections_below]
    wt0 = graph_walls.nodes["WT0.0"]["element"]
    wt0_above_tags = [ib.other_tag for ib in wt0.intersections_above]
    db0 = graph_walls.nodes["DB0.0"]["element"]
    db0_above_tags = [ib.other_tag for ib in db0.intersections_above]

    with check:
        assert "WT0.0" in j0_below_tags
    with check:
        assert "DB0.0" in j0_below_tags
    with check:
        assert "J0.0" in wt0_above_tags
    with check:
        assert "J0.0" in db0_above_tags

    j0 = graph_inters.nodes["J0.0"]["element"]
    j0_below_tags = [ib.other_tag for ib in j0.intersections_below]

    j1 = graph_inters.nodes["J0.1"]["element"]
    j1_below_tags = [ib.other_tag for ib in j1.intersections_below]

    wt0 = graph_inters.nodes["WT0.0"]["element"]
    wt0_below_tags = [ib.other_tag for ib in wt0.intersections_below]
    wt0_above_tags = [ib.other_tag for ib in wt0.intersections_above]

    db0 = graph_inters.nodes["DB0.0"]["element"]
    db0_below_tags = [ib.other_tag for ib in db0.intersections_below]
    db0_above_tags = [ib.other_tag for ib in db0.intersections_above]

    db1 = graph_inters.nodes["DB0.1"]["element"]
    db1_below_tags = [ib.other_tag for ib in db1.intersections_below]
    db1_above_tags = [ib.other_tag for ib in db1.intersections_above]

    ct0 = graph_inters.nodes["CT0.0"]["element"]
    ct0_above_tags = [ib.other_tag for ib in ct0.intersections_above]

    ct1 = graph_inters.nodes["CT0.1"]["element"]
    ct1_above_tags = [ib.other_tag for ib in ct1.intersections_above]

    ct2 = graph_inters.nodes["CT0.2"]["element"]
    ct2_above_tags = [ib.other_tag for ib in ct2.intersections_above]

    with check:
        assert "DB0.0" in j0_below_tags
    with check:
        assert "DB0.1" in j0_below_tags
    with check:
        assert "DB0.0" in j1_below_tags
    with check:
        assert "WT0.0" in j1_below_tags
    with check:
        assert "CT0.2" in db1_below_tags
    with check:
        assert "CT0.3" in db1_below_tags
    with check:
        assert "CT0.0" in db0_below_tags
    with check:
        assert "CT0.1" in db0_below_tags
    with check:
        assert "J0.1" in wt0_above_tags
    with check:
        assert "J0.0" in db1_above_tags
    with check:
        assert "J0.1" in db0_above_tags
    with check:
        assert "J0.0" in db0_above_tags
    with check:
        assert "DB0.0" in ct0_above_tags
    with check:
        assert "DB0.0" in ct1_above_tags
    with check:
        assert "DB0.1" in ct2_above_tags


def _jitter_geometry(geom, rng, magnitude):
    """Perturb every vertex by up to 'magnitude'; identical vertices move together."""
    cache = {}

    def fn(coords):
        out = np.empty_like(coords)
        for i, (x, y) in enumerate(coords):
            key = (round(float(x), 9), round(float(y), 9))
            if key not in cache:
                cache[key] = rng.uniform(-magnitude, magnitude, size=2)
            out[i] = np.array([x, y], dtype=float) + cache[key]
        return out

    return transform(geom, fn)


def _graph_edge_set(graph):
    return set((u, v, graph.edges[u, v]["edge_type"]) for u, v in graph.edges)


@pytest.mark.parametrize("magnitude", [1e-6, 1e-8])
@pytest.mark.parametrize("fixture_name", ["intersections.pdf", "sketch_to_scale.pdf"])
def test_production_graph_topology_invariant_under_jitter(fixture_name, magnitude):
    """
    The acceptance criterion for the node-canonicalization effort (design §10),
    applied to the real production pipeline: perturbing every input coordinate
    must leave the load-path graph (node set and typed edge set) unchanged.
    """
    raw = GeometryGraph.from_pdf_file(
        TEST_DATA / fixture_name,
        scale=QUARTER_INCH_SCALE,
        process_gravity_frame=False,
    )
    elements = [raw.nodes[n]["element"] for n in raw.nodes]
    base = GeometryGraph.from_elements(copy.deepcopy(elements))
    base_nodes, base_edges = set(base.nodes), _graph_edge_set(base)
    rng = np.random.default_rng(20260717)
    for _ in range(5):
        jittered = copy.deepcopy(elements)
        for element in jittered:
            element.geometry = _jitter_geometry(element.geometry, rng, magnitude)
        graph = GeometryGraph.from_elements(jittered)
        with check:
            assert set(graph.nodes) == base_nodes
        with check:
            assert _graph_edge_set(graph) == base_edges


def test_intersection_region_shared_below_and_above():
    """
    Each physical crossing is computed once and shared: as built, the 'below'
    view on the lower-rank element and the 'above' view on the higher-rank
    element carry the byte-identical intersecting region (design §6, node
    canonicalization). This is the invariant produced by
    get_geometry_intersections; it is asserted on the raw graph (before the
    gravity-frame post-processing that recomputes regions — a downstream site
    still to be migrated in Phase 3).
    """
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "intersections.pdf",
        scale=QUARTER_INCH_SCALE,
        process_gravity_frame=False,
    )
    checked = 0
    for node_name in graph.nodes:
        element = graph.nodes[node_name]["element"]
        if not element.intersections_below:
            continue
        for below in element.intersections_below:
            other = graph.nodes[below.other_tag]["element"]
            matching_above = [
                above
                for above in (other.intersections_above or [])
                if above.other_tag == element.tag
            ]
            with check:
                assert (
                    matching_above
                ), f"{below.other_tag} has no 'above' view of {element.tag}"
            for above in matching_above:
                with check:
                    assert below.intersecting_region.equals_exact(
                        above.intersecting_region, 0.0
                    )
                checked += 1
    assert checked > 0


def _array_signature(graph):
    """Per collector: its subelement tags and each subelement's support tags."""
    sig = {}
    for node in graph.collector_elements:
        element = graph.nodes[node]["element"]
        sig[node] = [
            (sub.tag, tuple(sorted(ib.other_tag for ib in sub.intersections_below)))
            for sub in element.subelements or []
        ]
    return sig


@pytest.mark.parametrize("magnitude", [1e-6, 1e-8])
@pytest.mark.parametrize(
    "fixture_name",
    ["intersections.pdf", "sketch_to_scale.pdf", "collector_extents.pdf"],
)
def test_joist_arrays_invariant_under_jitter(fixture_name, magnitude):
    """
    Jitter acceptance (design §10) extended to JoistArrayModel: perturbing every
    input coordinate must not change which joists are generated or what each
    one bears on.
    """
    raw = GeometryGraph.from_pdf_file(
        TEST_DATA / fixture_name,
        scale=QUARTER_INCH_SCALE,
        process_gravity_frame=False,
    )
    elements = [raw.nodes[n]["element"] for n in raw.nodes]
    base = GeometryGraph.from_elements(copy.deepcopy(elements))
    base.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
    base_sig = _array_signature(base)
    assert any(base_sig.values())
    rng = np.random.default_rng(20260918)
    for _ in range(3):
        jittered = copy.deepcopy(elements)
        for element in jittered:
            element.geometry = _jitter_geometry(element.geometry, rng, magnitude)
        graph = GeometryGraph.from_elements(jittered)
        graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
        with check:
            assert _array_signature(graph) == base_sig
