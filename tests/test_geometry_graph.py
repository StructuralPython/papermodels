from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt
from papermodels.paper import annotations as an
from papermodels.paper import pdf
from papermodels.datatypes.joist_models import JoistArrayModel, CollectorTribModel
from papermodels.datatypes.element import create_element_filter
import numpy as np
import numpy.testing as npt
from pytest import fixture
from shapely import Polygon, box, Point
import pathlib
import fixtures
from decimal import Decimal

EIGHTTH_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(8)

TEST_DATA = pathlib.Path(__file__).parent / "test_data"


@fixture()
def load_frame_collectors_transfers():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "frame_collectors_transfers.pdf",
        scale=EIGHTTH_INCH_SCALE,
    )
    return graph


@fixture()
def load_collector_extents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents.pdf",
        scale=EIGHTTH_INCH_SCALE,
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
        CollectorTribModel
    )  # Assign all collectors trib model
    graph.assign_collector_behaviour(
        JoistArrayModel, steel_joist_arrays, spacing=1.0
    )  # Assign steel joists the array
    graph.assign_collector_behaviour(
        JoistArrayModel, user_designated_joists, spacing=1.0
    )
    les = graph.create_loaded_elements()
    assert (
        "SJ0.0-9" in les
    )  # Confirms that JoistArray behaviour created for steel joists
    assert "WJ0.0" in les  # Confirms that collector_trib behaviour created for WJ0.0
    assert "WJ0.1-9" in les  # Confirms that JoistArray behaviour created for WJ0.1


def test_load_collector_extents(load_collector_extents):
    assert load_collector_extents
