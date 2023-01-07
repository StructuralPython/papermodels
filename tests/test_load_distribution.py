from papermodels.models import load_distribution as ld
from shapely import Polygon, GeometryCollection

import pytest

square_45 = Polygon([
    [0., 5.],
    [5., 10.],
    [10., 5.],
    [5., 0.,],
])

square_pac_man = Polygon([
    [0., 0.],
    [0., 10.],
    [10., 10.],
    [5., 5.,],
    [10., 0.,],
])

square_bow_tie = Polygon([
    [0., 0.],
    [4., 5.],
    [0., 10.],
    [10., 10.],
    [6., 5.,],
    [10., 0.,],
])

square_with_L_hole = Polygon([
    [0., 0.,],
    [0., 10.],
    [10., 10.,],
    [10., 0.],
], 
    holes = [[
    [2., 2.],
    [2., 8.],
    [8., 8.],
    [8., 6.],
    [4., 6.],
    [4., 2.],
]]
)

p0 = (0., 0.)
p1 = (10., 10.)
p2 = (6., 0.)
p3 = (6., 10.)
p4 = (15., 20.)
p5 = (0., 10.)
p6 = (10., 0.)


def test_get_void_regions():
    assert str(ld.get_void_regions(square_with_L_hole)) == '[<POLYGON ((2 8, 8 8, 8 6, 4 6, 4 2, 2 2, 2 8))>]'
    assert str(ld.get_void_regions(square_bow_tie)) == '[<POLYGON ((0 10, 4 5, 0 0, 0 10))>, <POLYGON ((10 0, 6 5, 10 10, 10 0))>]'
    assert str(ld.get_void_regions(square_pac_man)) == '[<POLYGON ((10 0, 5 5, 10 10, 10 0))>]'
    assert ld.get_void_regions(square_45) == []


def test_get_slope_and_intercept():
    assert ld.get_slope_and_intercept(p0, p1) == (1.0, 0.)
    assert ld.get_slope_and_intercept(p2, p1) == (10./4., -15)
    with pytest.raises(ZeroDivisionError):
        ld.get_slope_and_intercept(p2, p3)

def test_get_overlap_coords():
    assert ld.get_overlap_coords(0., 5., 3., 6.) == (3., 5.)
    assert ld.get_overlap_coords(2., 6., 0., 5.) == (2., 5.)
    assert ld.get_overlap_coords(0., 10., 1., 6.) == (1., 6.)
    assert ld.get_overlap_coords(1., 5., 0., 10.) == (1., 5.)
    assert ld.get_overlap_coords(5., 5., 1., 10.) == None

def test_get_overlap_region():
    assert ld.get_overlap_region(p0, p6, p3, p1) == ld.Overlap(6., 10., 0., 0., 0., 10.)
    assert ld.get_overlap_region(p0, p5, p3, p1) == None

def test_get_overlap_regions():
    assert ld.get_overlap_regions(square_45) == [
        ld.Singularity(0., 5., 2.0, 0., 6),
        ld.Singularity(5., 10., -2.0, 10., 6),
    ]

def test_get_singularity_functions():
    assert ld.get_singularity_functions(square_45) == (
        [
            ld.Singularity(0., 5., 2.0, 0., 6),
            ld.Singularity(5., 10., -2.0, 10., 6),
        ],
        []
    )
    assert ld.get_singularity_functions(square_pac_man) == (
        [ld.Singularity(0., 10., 0., 10., 6)],
        [ld.Singularity(5., 10., -2., -0., 6)]
    )
    assert ld.get_singularity_functions(square_bow_tie) == (
        [ld.Singularity(0., 10., 0., 10., 6)],
        [
            ld.Singularity(0., 4., 2.5, -10., 6),
            ld.Singularity(6., 10., -2.5, -0., 6)
        ]
    )
    assert ld.get_singularity_functions(square_with_L_hole) == (
        [ld.Singularity(x0=0.0, x1=10.0, m=0.0, y0=10.0, precision=6)],
        [
            ld.Singularity(x0=4.0, x1=8.0, m=1.0, y0=-6.0, precision=6),
            ld.Singularity(x0=2.0, x1=4.0, m=-0.0, y0=-6.0, precision=6),
            ld.Singularity(x0=4.0, x1=8.0, m=-1.0, y0=4.0, precision=6)
        ]
    )


def test_get_range():
    assert ld.get_range(ld.Overlap(6., 10., 0., 0., 0., 10.)) == (10., 10.)
    assert ld.get_range(ld.Overlap(-10., -2., -1., -11., 1., 11.)) == (2., 18.)


def test_apply_total_load():
    assert ld.apply_total_load(ld.Singularity(x0=4.0, x1=8.0, m=2.0, y0=1, precision=6), 10) == ld.Singularity(x0=4.0, x1=8.0, m=1.0, y0=0.5, precision=6)


def test_singularities_to_polygon():
    square_45_sings = ld.get_singularity_functions(square_45)
    square_pac_man_sings = ld.get_singularity_functions(square_pac_man)
    square_bow_tie_sings = ld.get_singularity_functions(square_bow_tie)
    square_with_L_hole_sings = ld.get_singularity_functions(square_with_L_hole)

    assert ld.singularities_to_polygon(square_45_sings[0] + square_45_sings[1]).wkt == 'POLYGON ((0 0, 0 0, 5 10, 5 10, 10 0, 10 0, 0 0))'
    assert ld.singularities_to_polygon(square_pac_man_sings[0] + square_pac_man_sings[1]).wkt == 'POLYGON ((0 0, 0 10, 5 10, 5 10, 10 0, 10 0, 10 0, 0 0))'
    assert ld.singularities_to_polygon(square_bow_tie_sings[0] + square_bow_tie_sings[1]).wkt == 'POLYGON ((0 0, 0 0, 4 10, 4 10, 6 10, 6 10, 10 0, 10 0, 10 0, 0 0))'
    assert ld.singularities_to_polygon(square_with_L_hole_sings[0] + square_with_L_hole_sings[1]).wkt == 'POLYGON ((0 0, 0 10, 2 10, 2 4, 4 4, 4 8, 8 8, 8 10, 10 10, 10 0, 0 0))'


def test_overlap_region_to_singularity():
    assert ld.overlap_region_to_singularity(ld.Overlap(6., 10., 0., 0., 0., 10.), 6) == ld.Singularity(6., 10., 0., 10., 6)