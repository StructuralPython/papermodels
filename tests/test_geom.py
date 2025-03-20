from shapely import LineString, Polygon, Geometry, GeometryCollection
from shapely.affinity import translate
from math import isclose

from papermodels.geometry.geom_ops import check_corresponds

def test_check_corresponds():

    ls1 = LineString([[0, 0], [1, 0]])
    ls2 = LineString([[0, 0], [2, 0]])
    ls3 = LineString([[0, 1], [2, 1]])
    poly1 = Polygon([[0.95, -0.05], [0.95, 0.05], [1.05, 0.05], [1.05, -0.05]])
    poly2 = translate(poly1, yoff=0.01)

    assert check_corresponds(ls1, ls2) == 0.5
    assert check_corresponds(ls2, ls1) == 1.0
    assert check_corresponds(poly1, poly1) == 1.00
    assert check_corresponds(poly1, poly2) == 0.9
    assert check_corresponds(poly2, ls1) == 0.99
    assert check_corresponds(poly1, ls1) == 1.0
    assert check_corresponds(ls1, ls3) == 0.0
    assert check_corresponds(ls3, ls1) == 0.0