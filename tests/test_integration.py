from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt
from papermodels.paper import annotations as an
from papermodels.paper import pdf
from papermodels.datatypes.joist_models import JoistArrayModel, CollectorTribModel
import numpy as np
import numpy.testing as npt
from pytest import fixture
from shapely import Polygon, box, Point
import pathlib
import fixtures
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
def load_collector_extents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents.pdf",
        scale=QUARTER_INCH_SCALE,
    )
    return graph


@fixture()
def sketch_to_scale_to_trib_loaded_elements(load_sketch_to_scale):
    graph = load_sketch_to_scale
    graph.assign_collector_behaviour(CollectorTribModel)
    les = graph.create_loaded_elements()
    return les


@fixture()
def sketch_to_scale_to_array_loaded_elements(load_sketch_to_scale):
    graph = load_sketch_to_scale
    graph.assign_collector_behaviour(JoistArrayModel, spacing=1)
    les = graph.create_loaded_elements()
    return les


@fixture()
def collector_extents_to_trib_loaded_elements(load_collector_extents):
    graph = load_collector_extents
    graph.assign_collector_behaviour(CollectorTribModel)
    les = graph.create_loaded_elements()
    return les


def test_sketch_to_scale_loads(load_sketch_to_scale):
    assert load_sketch_to_scale


def test_sketch_to_scale_creates_trib_loaded_elements(
    sketch_to_scale_to_trib_loaded_elements,
):
    assert sketch_to_scale_to_trib_loaded_elements


def test_sketch_to_scale_creates_array_loaded_elements(
    sketch_to_scale_to_array_loaded_elements,
):
    les = sketch_to_scale_to_array_loaded_elements
    assert les


def test_kwargs_pass_thru_sketch_to_scale_trib(sketch_to_scale_to_trib_loaded_elements):
    les = sketch_to_scale_to_trib_loaded_elements
    assert les["J4.0"].model()["element_attributes"]["user_defined"] == {
        "slope": "4/12",
        "slope_down": "right",
    }
    assert les["FB2.0"].model()["element_attributes"]["user_defined"] == {
        "user_defined": "data"
    }


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


def test_joists_loaded_sketch_to_scale(sketch_to_scale_to_trib_loaded_elements):
    les = sketch_to_scale_to_trib_loaded_elements
    assert les["J4.0"].model()["loads"]["distributed_loads"]


def test_collector_extent_loads(load_collector_extents):
    assert load_collector_extents


def test_collector_extent_creates_loaded_elements(
    collector_extents_to_trib_loaded_elements,
):
    les = collector_extents_to_trib_loaded_elements
    assert set(les.keys()) == set(
        [
            "SJ0.0-0",
            "SJ0.0-1",
            "SJ0.0-2",
            "SJ0.0-3",
            "SJ0.0-4",
            "SJ0.0-5",
            "SJ0.0-6",
            "SJ0.1",
            "WT0.0",
            "WT0.1",
            "WT0.2",
            "WT0.3",
            "FB0.0",
            "FB0.2",
            "FB0.1",
            "FB0.3",
            "CT0.4",
            "CT0.5",
            "CT0.8",
            "FB0.4",
            "WT0.0",
            "CT0.0",
            "CT0.6",
            "WT0.2",
            "CT0.2",
            "CT0.7",
            "WT0.1",
            "CT0.1",
            "CT0.3",
        ]
    )
    assert les["WT0.3"].model()["loads"][
        "distributed_loads"
    ]  # There are loads present on the intermediate support

    assert (
        les["WT0.1"].model()["loads"]["distributed_loads"][0]["transfer_source"]
        == "SJ0.0-0"
    )
    assert les["WT0.1"].model()["loads"]["distributed_loads"][0]["start_loc"] == 2.521
    assert les["WT0.1"].model()["loads"]["distributed_loads"][0]["end_loc"] == 3.858

    # This member experiences a splitting that occurs from an intermediate support that is found
    # within its overlap region.
    assert (
        les["FB0.0"].model()["loads"]["distributed_loads"][0]["transfer_source"]
        == "SJ0.0-2"
    )
    assert les["FB0.0"].model()["loads"]["distributed_loads"][0]["start_loc"] == 0.326
    assert les["FB0.0"].model()["loads"]["distributed_loads"][0]["end_loc"] == 3.01
    assert (
        les["FB0.0"].model()["loads"]["distributed_loads"][1]["transfer_source"]
        == "SJ0.0-3"
    )
    assert les["FB0.0"].model()["loads"]["distributed_loads"][1]["start_loc"] == 3.01
    assert les["FB0.0"].model()["loads"]["distributed_loads"][1]["end_loc"] == 3.912
