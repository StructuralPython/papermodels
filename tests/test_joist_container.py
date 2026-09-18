"""
End-to-end tests for the joist_container markup: legend parsing, correlation
to a joist prototype, graph intersections and the container-mode joist array.

Built from synthetic annotations (see fixtures.shapely_to_annotation) through
the same GeometryGraph.from_annotations path used for PDF files.
"""

import dataclasses
import warnings
from decimal import Decimal

import pytest
from shapely import LineString, Polygon, box, unary_union

import fixtures
from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.joist_models import JoistArrayModel

# The synthetic beam has nothing below it; the graph rightly calls it orphaned.
pytestmark = pytest.mark.filterwarnings("ignore:Orphaned element")

RED = (Decimal("1"), Decimal("0"), Decimal("0"))
GREEN = (Decimal("0"), Decimal("1"), Decimal("0"))
BLUE = (Decimal("0"), Decimal("0"), Decimal("1"))
GREY = (Decimal("0.5"), Decimal("0.5"), Decimal("0.5"))

LEGEND = {
    "joist": ("Legend\nType: Joist\nRank: 0", RED),
    "wall": ("Legend\nType: Wall\nRank: 2\nReaction Type: linear", GREEN),
    "beam": ("Legend\nType: Beam\nRank: 1", BLUE),
    "container": ("Legend\nType: Joist Container", GREY),
}


def shapely_to_annotation(geom, **kw):
    """fixtures.shapely_to_annotation with hashable (tuple) vertices/matrix."""
    annot = fixtures.shapely_to_annotation(geom, **kw)
    return dataclasses.replace(
        annot, vertices=tuple(annot.vertices), matrix=tuple(annot.matrix)
    )


WALL = box(0, -0.25, 20, 0.25)  # centerline y = 0
BEAM = LineString([(0, 10), (20, 14)])  # y = 10 + 0.2x
# Top edge rises faster than the beam: the cantilever grows from 1.0 to 2.0.
CONTAINER = Polygon([(0, -1), (20, -1), (20, 16), (0, 11)])  # top y = 11 + 0.25x
PROTOTYPE = LineString([(10, -1), (10, 13.5)])


def _annots(extra=()):
    def a(geom, kind, page=0):
        return shapely_to_annotation(
            geom, page=page, text="", line_color=LEGEND[kind][1]
        )

    legend = [
        shapely_to_annotation(
            (
                box(100 + 5 * i, 100, 102 + 5 * i, 102)
                if kind in ("wall", "container")
                else LineString([(100 + 5 * i, 100), (102 + 5 * i, 100)])
            ),
            text=text,
            line_color=color,
        )
        for i, (kind, (text, color)) in enumerate(LEGEND.items())
    ]
    items = [
        a(WALL, "wall"),
        a(BEAM, "beam"),
        a(CONTAINER, "container"),
        a(PROTOTYPE, "joist"),
    ] + [a(geom, kind) for geom, kind in extra]
    return legend + items


def _graph(extra=(), **kw):
    return GeometryGraph.from_annotations(_annots(extra), scale=Decimal(1), **kw)


def _prototype_node(graph):
    (node,) = [n for n in graph.collector_elements]
    return node


def test_container_is_parsed_and_correlated_to_the_prototype():
    graph = _graph()
    element = graph.nodes[_prototype_node(graph)]["element"]
    assert element.joist_container is not None
    assert element.joist_container.equals(CONTAINER)
    assert element.array_region.equals(CONTAINER)
    # The container is markup, not a structural element
    assert len(graph.nodes) == 3


def test_container_array_varies_cantilever_and_backspan():
    graph = _graph()
    node = _prototype_node(graph)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
    # (The synthetic beam has nothing below it, so an orphan warning is expected.)
    assert not [w for w in caught if "Joist array" in str(w.message)]
    subs = graph.nodes[node]["element"].subelements
    assert len(subs) == 21
    (wall_tag,) = [
        n for n in graph.nodes if graph.nodes[n]["element"].geometry.geom_type == "Polygon"
    ]
    for sub in subs:
        x = sub.geometry.coords[0][0]
        ys = sorted(c[1] for c in sub.geometry.coords)
        assert ys[0] == pytest.approx(-1.0, abs=1e-9)  # container bottom edge
        assert ys[1] == pytest.approx(11 + 0.25 * x, abs=1e-9)  # container top edge
        regions = {
            ib.other_tag: ib.intersecting_region for ib in sub.intersections_below
        }
        assert len(regions) == 2
        assert regions[wall_tag].y == pytest.approx(0.0, abs=1e-12)  # wall centerline
    # Cantilever past the beam grows from 1.0 to 2.0 along the array
    cantilevers = [
        max(c[1] for c in s.geometry.coords) - (10 + 0.2 * s.geometry.coords[0][0])
        for s in subs
    ]
    assert cantilevers[0] == pytest.approx(1.0) and cantilevers[-1] == pytest.approx(
        2.0
    )
    bands = [s.trib_area for s in subs]
    assert sum(b.area for b in bands) == pytest.approx(CONTAINER.area, rel=1e-9)
    assert unary_union(bands).symmetric_difference(CONTAINER).area < 1e-9


def test_container_drawn_to_wall_face_still_bears_on_wall_centerline():
    # Container stops at the wall's inner face (y=0.25); the centerline (y=0)
    # is outside it, but the wall polygon is inside.
    face_container = Polygon([(0, 0.25), (20, 0.25), (20, 16), (0, 11)])
    annots = [a for a in _annots() if not (a.line_color == GREY and a.text == "")] + [
        shapely_to_annotation(face_container, line_color=GREY)
    ]
    graph = GeometryGraph.from_annotations(annots, scale=Decimal(1))
    graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
    subs = graph.nodes[_prototype_node(graph)]["element"].subelements
    assert len(subs) == 21
    for sub in subs:
        assert min(c[1] for c in sub.geometry.coords) == pytest.approx(0.0, abs=1e-12)


def test_container_warnings():
    empty = box(40, 0, 50, 10)
    with pytest.warns(UserWarning, match="holds 0 joist prototypes"):
        _graph(extra=[(empty, "container")])
