from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt
from papermodels.paper import annotations as an
from papermodels.paper import pdf
from papermodels.datatypes.joist_models import JoistArrayModel
import numpy as np
import numpy.testing as npt
from pytest import fixture
from pytest_check import check
from shapely import Polygon, box, Point
import pathlib
import fixtures
import math
from decimal import Decimal

QUARTER_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(4)

TEST_DATA = pathlib.Path(__file__).parent / "test_data"


@fixture()
def load_sketch_to_scale():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "sketch_to_scale.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def load_trib_areas_basic():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "trib-areas-basic.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def load_collector_extents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def load_horiz_extents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "horiz_extents.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def sketch_to_scale_to_array_loaded_elements(load_sketch_to_scale):
    graph = load_sketch_to_scale
    graph.assign_collector_behaviour(JoistArrayModel, spacing=1)
    les = graph.create_loaded_elements()
    return les


@fixture()
def collector_extents_to_array_loaded_elements(load_collector_extents):
    graph = load_collector_extents
    graph.assign_collector_behaviour(JoistArrayModel)
    les = graph.create_loaded_elements()
    return les


@fixture()
def horiz_extents_to_array_loaded_elements(load_horiz_extents):
    graph = load_horiz_extents
    graph.assign_collector_behaviour(JoistArrayModel)
    les = graph.create_loaded_elements()
    return les


@fixture()
def trib_areas_basic_loaded_elements(load_trib_areas_basic):
    graph = load_trib_areas_basic
    les = graph.create_loaded_elements()
    return les


def test_sketch_to_scale_loads(load_sketch_to_scale):
    assert load_sketch_to_scale


def test_sketch_to_scale_creates_array_loaded_elements(
    sketch_to_scale_to_array_loaded_elements,
):
    les = sketch_to_scale_to_array_loaded_elements
    assert les


def test_horiz_extents_loads(horiz_extents_to_array_loaded_elements):
    les = horiz_extents_to_array_loaded_elements
    assert les


def test_kwargs_pass_thru_sketch_to_scale_array(
    sketch_to_scale_to_array_loaded_elements,
):
    les = sketch_to_scale_to_array_loaded_elements
    assert les["J4.0-1"].model()["element_attributes"]["user_defined"] == {
        "slope": "4/12",
        "slope_down": "right",
    }
    assert les["FB2.0"].model()["element_attributes"]["user_defined"] == {
        "user_defined": "data"
    }


def test_joists_loaded_sketch_to_scale(sketch_to_scale_to_array_loaded_elements):
    les = sketch_to_scale_to_array_loaded_elements
    j40_subs = [tag for tag in les if tag.startswith("J4.0-")]
    j40_loads = [
        dl for tag in j40_subs for dl in les[tag].model()["loads"]["distributed_loads"]
    ]
    with check:
        assert {dl["occupancy"] for dl in j40_loads} == {"roof"}
        # All of the roof area over the array is carried by its joists
        assert math.isclose(
            sum(dl["applied_area"] for dl in j40_loads), 187.941, abs_tol=2e-3
        )
        wt40_loads = les["WT4.0"].model()["loads"]["point_loads"]
        assert {pl["transfer_source"].rsplit("-", 1)[0] for pl in wt40_loads} == {
            "J4.0"
        }
        wt40_locations = [pl["location"] for pl in wt40_loads]
        assert min(wt40_locations) == 0.573
        assert max(wt40_locations) == 11.687

        fb1_3 = les["FB1.3"].model()["loads"]["point_loads"]
        fb1_3_locations = [
            pl["location"] for pl in fb1_3 if pl["transfer_source"].startswith("J1.1")
        ]
        assert min(fb1_3_locations) == 0.433
        assert max(fb1_3_locations) == 10.906

        j1_1 = les["J1.1-0"].model()
        assert j1_1["element_geometry"]["supports"][1]["overlap_length"] == 0.392


def test_collector_extent_loads(load_collector_extents):
    assert load_collector_extents


def test_collector_extent_creates_loaded_elements(
    collector_extents_to_array_loaded_elements,
):
    """
    An extent-line joist array (SJ0.0) spreading over supports that come and go:
    WT0.3 is an intermediate support under part of the array only, and the
    outer supports change from WT0.1/FB0.2 to FB0.0/WT0.2 to WT0.0/FB0.3.
    """
    les = collector_extents_to_array_loaded_elements
    structural = {
        "WT0.0",
        "WT0.1",
        "WT0.2",
        "WT0.3",
        "FB0.0",
        "FB0.1",
        "FB0.2",
        "FB0.3",
        "FB0.4",
        "CT0.0",
        "CT0.1",
        "CT0.2",
        "CT0.3",
        "CT0.4",
        "CT0.5",
        "CT0.6",
        "CT0.7",
        "CT0.8",
    }
    joists = {f"SJ0.0-{i}" for i in range(15)} | {f"SJ0.1-{i}" for i in range(6)}
    assert set(les.keys()) == structural | joists

    def point_loads(tag):
        return les[tag].model()["loads"]["point_loads"]

    with check:
        # There are loads present on the intermediate support
        assert point_loads("WT0.3")
        assert {pl["transfer_source"] for pl in point_loads("WT0.3")} == {
            f"SJ0.0-{i}" for i in range(2, 8)
        }
    with check:
        wt01 = point_loads("WT0.1")
        assert wt01[0]["transfer_source"] == "SJ0.0-0"
        assert wt01[0]["location"] == 2.195
        assert [pl["location"] for pl in wt01] == [2.195, 3.195, 4.195, 5.195, 6.195]
    with check:
        fb00 = point_loads("FB0.0")
        assert fb00[0]["transfer_source"] == "SJ0.0-5"
        assert fb00[0]["location"] == 0.379
        assert {pl["transfer_source"] for pl in fb00} == {
            f"SJ0.0-{i}" for i in range(5, 10)
        }


def test_collector_extents_creates_array_loaded_elements(
    collector_extents_to_array_loaded_elements,
):
    les = collector_extents_to_array_loaded_elements
    assert les


def test_wall_point_load_locations(sketch_to_scale_to_array_loaded_elements):
    les = sketch_to_scale_to_array_loaded_elements
    wt4_pt = les["WT4.0"].model()["loads"]["point_loads"]
    acc = []
    for load in wt4_pt:
        acc.append(load["location"])

    with check:
        # Test that all intervals are about the same size by not exceeding the prescribed spacing
        # First, for walls (order-independent: the array direction is a convention)
        acc = sorted(acc)
        joist_intervals = [x[1] - x[0] for x in zip(acc[:-1], acc[1:])]
        assert math.isclose(max(joist_intervals), 1)

        # Test that all intervals are about the same size by not exceeding the prescribed spacing
        # Next, for beams
        fb4_pt = les["FB4.0"].model()["loads"]["point_loads"]
        acc = sorted(load["location"] for load in fb4_pt)
        joist_intervals = [x[1] - x[0] for x in zip(acc[:-1], acc[1:])]
        assert math.isclose(max(joist_intervals), 1)


def test_horiz_extents_joist_extents(horiz_extents_to_array_loaded_elements):
    les = horiz_extents_to_array_loaded_elements
    db0 = les["DB0.0"].model()
    db1 = les["DB0.1"].model()
    locations = sorted(pl["location"] for pl in db0["loads"]["point_loads"])
    assert locations[-1] == 7.83
    assert locations[0] == 3.134
    assert db0["element_geometry"]["supports"][0]["location"] == 0.0
    assert db0["element_geometry"]["supports"][1]["location"] == 11.610


def test_trib_areas_basic_loads(trib_areas_basic_loaded_elements):
    les = trib_areas_basic_loaded_elements
    dl = les["FB0.0"].model()["loads"]["distributed_loads"]
    assert dl[0]["start_loc"] == -0.144
    assert dl[0]["end_loc"] == 20.358
    assert dl[1]["start_loc"] == 8.143
    assert dl[1]["end_loc"] == 10.123
    assert dl[2]["start_loc"] == 10.123
    assert dl[2]["end_loc"] == 10.132
    assert dl[3]["start_loc"] == 10.132
    assert dl[3]["end_loc"] == 12.113
