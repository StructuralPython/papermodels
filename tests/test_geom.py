from shapely import LineString, Polygon, Geometry, GeometryCollection, MultiPoint
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
    assert wkt.dumps(MultiPoint(extents['A']), rounding_precision=3) == 'MULTIPOINT (20.981 304.127, 273.716 327.842)'
    assert wkt.dumps(MultiPoint(extents['B']), rounding_precision=3) == 'MULTIPOINT (50.000 4.000, 300.000 56.000)'