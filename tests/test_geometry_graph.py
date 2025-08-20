from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt
from papermodels.paper import annotations as an
from papermodels.paper import pdf
from papermodels.datatypes.joist_models import JoistArrayModel, collector_trib_model
import numpy as np
import numpy.testing as npt
from pytest import fixture
from shapely import Polygon, box, Point
import pathlib
import fixtures
from decimal import Decimal

EIGHTTH_INCH_SCALE = Decimal(1) / Decimal(72) / Decimal(8)

TEST_DATA = pathlib.Path(__file__).parent / "test_data"

@fixture()
def load_frame_collectors_transfers():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "frame_collectors_transfers.pdf",
        scale=EIGHTTH_INCH_SCALE,
    )
    return graph

def test_load_frame_collectors_transfers(load_frame_collectors_transfers):
    assert load_frame_collectors_transfers
