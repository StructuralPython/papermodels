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
        TEST_DATA / "resi_dormers.pdf",
        scale=QUARTER_INCH_SCALE,
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


def test_load_resi_dormers(load_resi_dormers):
    assert load_resi_dormers


def test_resi_dormers_array(load_resi_dormers):
    graph = load_resi_dormers
    # Assign collector behaviour using filter functions
    roof_joist_filter = create_element_filter(element_types=["RJ"])
    all_other_joists_filter = create_element_filter(exclude_element_types=["RJ"])

    # First, assign the default behaviour to everything
    graph.assign_collector_behaviour(
        CollectorTribModel, filter_function=all_other_joists_filter
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
    assert len(fb03_pl) == 7
    rj001 = les["RJ0.0-1"].model()
    rj006 = les["RJ0.0-6"].model()
    assert rj001["element_attributes"]["length"] == 1.094
    assert rj006["element_attributes"]["length"] == 5.869


def test_many_correspondents(load_many_correspondents):
    graph = load_many_correspondents
    assert len(graph.nodes["WT2.0"]["element"].correspondents_below) == 1
    assert len(graph.nodes["WB1.1"]["element"].correspondents_above) == 1
    assert graph.nodes["WT2.0"]["element"].correspondents_below[0].other_tag == "WB1.1"
    assert graph.nodes["WB1.1"]["element"].correspondents_above[0].other_tag == "WT2.0"
