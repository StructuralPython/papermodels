"""
Characterization tests for the geometry helpers that the JoistArrayModel
refactor keeps relying on (see the helper audit in the refactor plan).

Nothing is reused on trust: each helper the new model depends on is pinned
here across orientations so that a limiting assumption (axis alignment,
positive-x ordering, floating-point collapse) shows up as a failing test
rather than as a silently wrong joist array.
"""

import math

import numpy as np
import pytest
from shapely import LineString, Polygon, Point
from shapely import affinity as aff

from papermodels.geometry import geom_ops
from papermodels.geometry.load_projection import project_loading_areas

ANGLES = [0, 15, 30, 45, 60, 89.9, 90, 90.1, 135, 179.9, 180, 225, 270, 315]


# --------------------------------------------------------------------------
# get_wall_centerline / get_rectangle_centerline (oriented bounding box)
# --------------------------------------------------------------------------


def _wall(length=10.0, thickness=0.5):
    return Polygon(
        [
            (0, -thickness / 2),
            (length, -thickness / 2),
            (length, thickness / 2),
            (0, thickness / 2),
        ]
    )


@pytest.mark.parametrize("angle", ANGLES)
def test_wall_centerline_is_long_axis_at_any_angle(angle):
    wall = aff.rotate(_wall(), angle, origin=(0, 0))
    spine = geom_ops.get_wall_centerline(wall)
    assert math.isclose(spine.length, 10.0, abs_tol=1e-9)
    # Parallel to the wall's long axis
    expected = np.array([math.cos(math.radians(angle)), math.sin(math.radians(angle))])
    (x0, y0), (x1, y1) = spine.coords
    got = np.array([x1 - x0, y1 - y0]) / spine.length
    assert math.isclose(abs(float(np.dot(expected, got))), 1.0, abs_tol=1e-9)
    # Passes through the wall's centroid
    assert spine.distance(wall.centroid) < 1e-9


def test_wall_centerline_absorbs_vertex_noise():
    noisy = Polygon([(0, -0.25), (5, -0.2501), (10, -0.25), (10, 0.25), (0, 0.2499)])
    spine = geom_ops.get_wall_centerline(noisy)
    assert math.isclose(spine.length, 10.0, abs_tol=1e-3)
    assert abs(spine.coords[0][1]) < 1e-3 and abs(spine.coords[1][1]) < 1e-3


def test_wall_centerline_near_square_is_an_axis_not_a_diagonal():
    sq = Polygon([(0, 0), (1.0, 0), (1.0, 0.98), (0, 0.98)])
    spine = geom_ops.get_wall_centerline(sq)
    (x0, y0), (x1, y1) = spine.coords
    assert math.isclose(x0, x1, abs_tol=1e-9) or math.isclose(y0, y1, abs_tol=1e-9)


# --------------------------------------------------------------------------
# rotate_90_vector
# --------------------------------------------------------------------------


@pytest.mark.parametrize("angle", ANGLES)
def test_rotate_90_vector_is_exact_quarter_turn(angle):
    v = np.array([math.cos(math.radians(angle)), math.sin(math.radians(angle))])
    ccw = np.asarray(geom_ops.rotate_90_vector(v, ccw=True))
    cw = np.asarray(geom_ops.rotate_90_vector(v, ccw=False))
    assert math.isclose(float(np.dot(v, ccw)), 0.0, abs_tol=1e-12)
    assert math.isclose(float(np.linalg.norm(ccw)), 1.0, abs_tol=1e-12)
    # ccw is a left turn, cw a right turn
    assert v[0] * ccw[1] - v[1] * ccw[0] > 0
    assert v[0] * cw[1] - v[1] * cw[0] < 0
    np.testing.assert_allclose(ccw, -cw, atol=1e-12)


# --------------------------------------------------------------------------
# Load projection of trib bands (load_distribution, via project_loading_areas)
# --------------------------------------------------------------------------

JOIST = LineString([(0, 0), (10, 0)])
BANDS = {
    # name: (polygon about JOIST, magnitude at x=0, magnitude at x=10) for unit total
    "rect": (Polygon([(0, -0.5), (10, -0.5), (10, 0.5), (0, 0.5)]), 0.1, 0.1),
    "trapezoid": (
        Polygon([(0, -0.5), (10, -0.5), (10, 1.5), (0, 0.5)]),
        1 / 15,
        2 / 15,
    ),
    "triangle": (Polygon([(0, 0), (10, -1), (10, 1)]), 0.0, 0.2),
}


def _endpoints(dist_loads):
    """First and last (x, magnitude) of a single-polygon projection."""
    pairs = dist_loads[0]
    first = pairs[0][0]
    last = pairs[-1][-1]
    return first, last


@pytest.mark.parametrize("angle", ANGLES)
@pytest.mark.parametrize("band", list(BANDS))
def test_projected_band_distribution_is_rotation_invariant(band, angle):
    poly, w0, w10 = BANDS[band]
    joist = aff.rotate(JOIST, angle, origin=(0, 0))
    rotated = aff.rotate(poly, angle, origin=(0, 0))
    dist = project_loading_areas(joist, [(rotated, None)])
    total = geom_ops.calculate_trapezoid_area_sums(dist)[0]
    assert math.isclose(total, 1.0, abs_tol=1e-6)

    (x_first, m_first), (x_last, m_last) = _endpoints(dist)
    assert math.isclose(x_first, 0.0, abs_tol=1e-6)
    assert math.isclose(x_last, 10.0, abs_tol=1e-6)
    # load_distribution measures along the member from its positive-x start
    # node, so the magnitudes are mirrored when the member points in -x.
    start_node, _ = geom_ops.get_start_end_nodes(joist)
    if start_node.distance(Point(0, 0)) < 1e-9:
        expected_first, expected_last = w0, w10
    else:
        expected_first, expected_last = w10, w0
    assert math.isclose(m_first, expected_first, abs_tol=1e-6)
    assert math.isclose(m_last, expected_last, abs_tol=1e-6)
