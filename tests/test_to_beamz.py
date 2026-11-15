"""LoadedElement.to_beamz: a member's symbolic load transfer as a beamz9000 beam."""

import math
import pathlib
from decimal import Decimal

import pytest

from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.joist_models import JoistArrayModel

QUARTER_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(4)
TEST_DATA = pathlib.Path(__file__).parent / "test_data"


def _loaded(name: str, joists: bool = False):
    graph = GeometryGraph.from_pdf_file(TEST_DATA / name, scale=QUARTER_INCH_SCALE)
    if joists:
        graph.assign_collector_behaviour(JoistArrayModel, spacing=1)
    return graph.create_loaded_elements()


@pytest.fixture(scope="module")
def sketch():
    return _loaded("sketch_to_scale.pdf")


@pytest.fixture(scope="module")
def sketch_joists():
    return _loaded("sketch_to_scale.pdf", joists=True)


@pytest.fixture(scope="module")
def collector_extents():
    return _loaded("collector_extents.pdf")


def _loads(beam: dict) -> list[tuple]:
    return [(kind, load) for item in beam["loads"] for kind, load in item.items()]


def test_transfers_are_labelled_with_their_source(sketch):
    beam = sketch["FB3.1"].to_beamz()
    assert beam["spans"] == [13.436, 4.875]
    assert beam["supports"] == ["pin", "roller", "roller"]
    assert _loads(beam) == [
        ("point", {"p": 1, "x": 4.98, "label": "CB3.0"}),
        ("dist", {"w": 1, "x1": 0.025, "x2": 18.282, "label": "J3.0"}),
    ]


def test_each_joist_of_an_array_is_a_labelled_point_load(sketch_joists):
    loads = _loads(sketch_joists["FB3.1"].to_beamz())
    labels = [load["label"] for kind, load in loads if kind == "point"]
    assert "CB3.0" in labels
    assert [l for l in labels if l.startswith("J3.0-")] == [
        f"J3.0-{i}" for i in range(20)
    ]
    assert all(load["p"] == 1 for _, load in loads)


def test_loading_areas_are_one_pressure_load_per_stretch(sketch_joists):
    # "big stuff" reaches J3.0-5 in several trapezoid pieces; they become one load.
    loads = _loads(sketch_joists["J3.0-5"].to_beamz())
    assert loads == [
        (
            "pressure",
            {"q": 1, "trib": 1, "x1": 0.0, "x2": 11.062, "label": "residential"},
        ),
        (
            "pressure",
            {"q": 1, "trib": 1, "x1": 2.391, "x2": 9.627, "label": "big stuff"},
        ),
    ]


def test_loads_stay_on_the_beam(sketch, sketch_joists):
    for loaded in (sketch, sketch_joists):
        for elem in loaded.values():
            if elem.geometry.geom_type != "LineString":
                continue
            beam = elem.to_beamz()
            length = sum(beam["spans"])
            for kind, load in _loads(beam):
                xs = [load["x"]] if kind == "point" else [load["x1"], load["x2"]]
                assert all(0 <= x <= length + 1e-9 for x in xs)


def test_terminus_supports_start_the_beam_and_overhangs_are_cantilevers(
    collector_extents,
):
    # FB0.3 bears on CT0.2 at its start and runs past CT0.7 at its end.
    beam = collector_extents["FB0.3"].to_beamz()
    assert beam["supports"] == ["pin", "roller", "free"]
    assert beam["spans"] == [10.704, 1.362]
    # FB0.2 overhangs CT0.8 at its start and bears on FB0.4 at its end.
    beam = collector_extents["FB0.2"].to_beamz()
    assert beam["supports"] == ["free", "pin", "roller"]
    assert math.isclose(sum(beam["spans"]), 7.526)


def test_units_and_beam_options(sketch):
    beam = sketch["FB3.1"].to_beamz(length_unit="ft", node_labels="auto")
    assert beam["spans"] == ["13.436 ft", "4.875 ft"]
    assert _loads(beam)[0] == ("point", {"p": "1 kN", "x": "4.98 ft", "label": "CB3.0"})
    assert beam["node_labels"] == "auto"


def test_only_horizontal_elements(collector_extents):
    wall = next(
        e for e in collector_extents.values() if e.geometry.geom_type == "Polygon"
    )
    with pytest.raises(ValueError, match="only horizontal"):
        wall.to_beamz()


@pytest.mark.parametrize("length_unit", [None, "ft"])
def test_every_member_draws_in_beamz9000(sketch_joists, length_unit):
    bz = pytest.importorskip("beamz9000")
    for tag, elem in sketch_joists.items():
        if elem.geometry.geom_type != "LineString":
            continue
        scheme = {tag: elem.to_beamz(length_unit=length_unit)}
        if length_unit:
            scheme["_environment"] = "structural"
        spec = bz.from_scheme(bz.parse_scheme(scheme)[tag])
        assert [l.label1 for l in spec.loads] == [
            load["label"] for _, load in _loads(scheme[tag])
        ]
        assert not spec.keys  # Symbolic: no load cases, no legend.
