from shapely import Point, LineString, Polygon, Geometry, GeometryCollection, MultiPoint
from shapely.affinity import translate
from shapely import wkt
from math import isclose

from papermodels.geometry import geom_ops

def test_check_corresponds():

    ls1 = LineString([[0, 0], [1, 0]])
    ls2 = LineString([[0, 0], [2, 0]])
    ls3 = LineString([[0, 1], [2, 1]])
    poly1 = Polygon([[0.95, -0.05], [0.95, 0.05], [1.05, 0.05], [1.05, -0.05]])
    poly2 = translate(poly1, yoff=0.01)

    assert geom_ops.check_corresponds(ls1, ls2) == 0.5
    assert geom_ops.check_corresponds(ls2, ls1) == 1.0
    assert geom_ops.check_corresponds(poly1, poly1) == 1.00
    assert geom_ops.check_corresponds(poly1, poly2) == 0.9
    assert geom_ops.check_corresponds(poly2, ls1) == 0.0
    assert geom_ops.check_corresponds(poly1, ls1) == 0.0
    assert geom_ops.check_corresponds(ls1, ls3) == 0.0
    assert geom_ops.check_corresponds(ls3, ls1) == 0.0


def test_get_joist_extents():
    ls1 = LineString([[50, 4], [300,56]])
    ls2 = LineString([[-23, 300], [350, 335]])
    j1 = LineString([[140.0, -23.4], [100.0, 390.3]])
    extents = geom_ops.get_joist_extents(j1, [ls1, ls2])
    assert (
        wkt.dumps(MultiPoint(extents[0] + extents[1]), trim=True, rounding_precision=3) 
        == 'MULTIPOINT (273.716 327.842, 20.981 304.127, 300 56, 50 4)'
    )
    ls1 = LineString([[0, 0], [0, 100]])
    ls2 = LineString([[50, -20], [50, 80]])
    j1 = LineString([[-20, 40], [60, 40]])
    extents = geom_ops.get_joist_extents(j1, [ls1, ls2])
    assert (
        wkt.dumps(MultiPoint(extents[0] + extents[1]), trim=True, rounding_precision=3)
        == 'MULTIPOINT (0 80, 0 0, 50 80, 50 0)'
    )


def test_order_nodes_positive():
    p1 = Point([0, 0])
    p2 = Point([0, 10])
    p3 = Point([10, 10])
    p4 = Point([12, 0])
    p5 = Point([5, -10])
    p6 = Point([0.1, -12])

    assert geom_ops.order_nodes_positive([p6, p1]) == (p1, p6)
    assert geom_ops.order_nodes_positive([p4, p3]) == (p3, p4)
    assert geom_ops.order_nodes_positive([p6, p5]) == (p6, p5)
    assert geom_ops.order_nodes_positive([p2, p1]) == (p1, p2)