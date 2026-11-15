"""LoadedElement.to_beamz: a loaded element as a beamz9000 scheme beam."""

import math
import pathlib
from decimal import Decimal

import pytest

from papermodels.datatypes.geometry_graph import GeometryGraph

QUARTER_INCH_SCALE = Decimal(1) / Decimal(72) * Decimal(4)
TEST_DATA = pathlib.Path(__file__).parent / "test_data"
OCCUPANCIES = {"f": {"D": 0.015, "L": 0.040}, "a": {"D": 0.015, "L": 0.100}}


@pytest.fixture(scope="module")
def trib_areas_basic():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "trib-areas-basic.pdf", scale=QUARTER_INCH_SCALE
    )
    return graph.create_loaded_elements()


@pytest.fixture(scope="module")
def collector_extents():
    graph = GeometryGraph.from_pdf_file(
        TEST_DATA / "collector_extents.pdf", scale=QUARTER_INCH_SCALE
    )
    return graph.create_loaded_elements()


def _area(load: dict) -> float:
    p = load["pressure"]
    t1 = p.get("trib", p.get("trib_width1"))
    t2 = p.get("trib", p.get("trib_width2"))
    return (t1 + t2) / 2 * (p["x2"] - p["x1"])


def test_beam_with_area_loads(trib_areas_basic):
    beam = trib_areas_basic["FB0.0"].to_beamz(OCCUPANCIES)
    assert beam["spans"] == [20.274]
    assert beam["supports"] == ["pin", "roller"]
    assert beam["occupancies"] == OCCUPANCIES
    first = beam["loads"][0]["pressure"]
    assert first == {"q": "f", "trib": 4.028, "x1": 0.0, "x2": 20.274}
    # The loads stay on the beam.
    for load in beam["loads"]:
        assert 0 <= load["pressure"]["x1"] < load["pressure"]["x2"] <= 20.274


def test_trib_widths_carry_each_loading_areas_share(trib_areas_basic):
    elem = trib_areas_basic["FB0.0"]
    beam = elem.to_beamz(OCCUPANCIES)
    # Loading area "a" lies within the beam, so its pressure loads cover its
    # whole applied area. "f" is trimmed to the beam's ends.
    a_area = sum(_area(l) for l in beam["loads"] if l["pressure"]["q"] == "a")
    applied = next(
        d["applied_area"]
        for d in elem._get_distributed_loads(12)
        if d["occupancy"] == "a"
    )
    assert a_area == pytest.approx(applied, rel=2e-3)
    f_area = sum(_area(l) for l in beam["loads"] if l["pressure"]["q"] == "f")
    assert f_area == pytest.approx(4.028 * 20.274, rel=1e-6)


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


def test_units_and_beam_options(trib_areas_basic):
    beam = trib_areas_basic["FB0.0"].to_beamz(
        {"f": {"D": "15 psf"}, "a": {"D": "15 psf"}},
        length_unit="ft",
        results=["forces", "reactions"],
        node_labels="auto",
    )
    assert beam["spans"] == ["20.274 ft"]
    assert beam["loads"][0]["pressure"]["trib"] == "4.028 ft"
    assert beam["results"] == ["forces", "reactions"] and beam["node_labels"] == "auto"


def test_missing_occupancy_is_named(trib_areas_basic):
    with pytest.raises(ValueError, match="occupancy 'a'"):
        trib_areas_basic["FB0.0"].to_beamz({"f": {"D": 1.0}})


def test_only_horizontal_elements(collector_extents):
    wall = next(
        e for e in collector_extents.values() if e.geometry.geom_type == "Polygon"
    )
    with pytest.raises(ValueError, match="only horizontal"):
        wall.to_beamz()


def test_the_beam_draws_and_analyses_in_beamz9000(trib_areas_basic):
    bz = pytest.importorskip("beamz9000")
    analysis = pytest.importorskip("beamz9000.analysis")
    pytest.importorskip("pycba_tools")
    beam = trib_areas_basic["FB0.0"].to_beamz(OCCUPANCIES, results=["reactions"])
    scheme = bz.parse_scheme({"FB0.0": beam})
    result = analysis.analyze(scheme["FB0.0"])
    total_d = sum(r.fy for r in result.reactions if r.case == "D")
    expected = sum(0.015 * _area(l) for l in beam["loads"])
    assert total_d == pytest.approx(expected, rel=1e-4)
