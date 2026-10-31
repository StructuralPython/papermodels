from shapely import LineString, Point, box

from papermodels.datatypes.element import Correspondent, Element, Intersection


def test_supporting_geometries():
    wall_a = box(0, -0.1, 0.2, 0.1)
    wall_b = box(4.8, -0.1, 5.0, 0.1)
    column_below = box(2.4, -0.1, 2.6, 0.1)
    beam = Element(
        geometry=LineString([(0, 0), (5, 0)]),
        tag="FB0.0",
        intersections_below=[
            Intersection(Point(0.1, 0), wall_a, "W0.0"),
            Intersection(Point(4.9, 0), wall_b, "W0.1"),
        ],
        correspondents_below=[Correspondent(1.0, column_below, "C0.0", 50)],
    )
    assert beam.supporting_geometries == [wall_a, wall_b, column_below]


def test_supporting_geometries_without_supports():
    element = Element(geometry=LineString([(0, 0), (1, 0)]), tag="J0.0")
    assert element.supporting_geometries == []
