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
def sketch_to_scale_to_trib_loaded_elements(load_sketch_to_scale):
    graph = load_sketch_to_scale
    graph.assign_collector_behaviour(collector_trib_model, as_subelements=False)
    les = graph.create_loaded_elements()
    return les

@fixture()
def sketch_to_scale_to_array_loaded_elements(load_sketch_to_scale):
    graph = load_sketch_to_scale
    graph.assign_collector_behaviour(JoistArrayModel.create_subelements, as_subelements=True, spacing=1)
    les = graph.create_loaded_elements()
    return les

def test_sketch_to_scale_loads(load_sketch_to_scale):
    assert load_sketch_to_scale

def test_sketch_to_scale_creates_trib_loaded_elements(sketch_to_scale_to_trib_loaded_elements):
    assert sketch_to_scale_to_trib_loaded_elements

def test_sketch_to_scale_creates_array_loaded_elements(sketch_to_scale_to_array_loaded_elements):
    les = sketch_to_scale_to_array_loaded_elements
    assert les

def test_kwargs_pass_thru_sketch_to_scale_trib(sketch_to_scale_to_trib_loaded_elements):
    les = sketch_to_scale_to_trib_loaded_elements
    assert les['J4.0'].model()['element_attributes']['user_defined'] == {"slope": "4/12", "slope_down": "right"}
    assert les['FB2.0'].model()['element_attributes']['user_defined'] == {"user_defined": "data"}


def test_kwargs_pass_thru_sketch_to_scale_array(sketch_to_scale_to_array_loaded_elements):
    les = sketch_to_scale_to_array_loaded_elements
    assert les['J4.0-1'].model()['element_attributes']['user_defined'] == {"slope": "4/12", "slope_down": "right"}
    assert les['FB2.0'].model()['element_attributes']['user_defined'] == {"user_defined": "data"}
