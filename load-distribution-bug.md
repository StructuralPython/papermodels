# Bug: `load-distribution` projects rotated trib areas to zero load

**Library:** `load-distribution` 0.1.7 (`load_distribution/load_distribution.py`)
**Found:** 2026-09-18, during the JoistArrayModel refactor (helper audit)
**Severity:** High. Load is lost silently: there is no error, and the total comes out as 0.
**Status in papermodels:** worked around locally by
`papermodels.geometry.load_projection.project_loading_areas`. The fix belongs
upstream (see "Fixing it" below).

## Symptom

`get_distributed_loads_from_projected_polygons(member, [(polygon, ...)])`
returns a load of essentially zero length when the member is at some angles.
This affects rectangular and trapezoidal loading polygons. Triangles are fine.

A 10-long member with a 1-wide rectangular band centred on it should give a
uniform load along 0–10 (magnitude 0.1 for a unit total):

```python
import load_distribution as ld
from shapely import LineString, Polygon
from shapely import affinity as aff

member = LineString([(0, 0), (10, 0)])
band = Polygon([(0, -0.5), (10, -0.5), (10, 0.5), (0, 0.5)])

ld.get_distributed_loads_from_projected_polygons(member, [(band, None)])
# 0 deg:  [[[(1e-13, 0.1), (9.9999999999999, 0.1)]]]   correct

m30 = aff.rotate(member, 30, origin=(0, 0))
b30 = aff.rotate(band, 30, origin=(0, 0))
ld.get_distributed_loads_from_projected_polygons(m30, [(b30, None)])
# 30 deg: [[[(0.0, 0.05), (5.55e-17, 0.1)],
#           [(9.9994e-13, 0.1), (1.00005e-12, 0.1)]]]         WRONG
```

At 30° the result is two slivers about 1e-12 long, and the real 0–10 segment
has disappeared. `calculate_trapezoid_area_sums` on it gives **0.0**, so a joist
at that angle carries no load from that band.

Rectangles, 5° steps, failing angles: **20, 30, 105, 205, 255, 275, 280, 285,
305, 325**. The set looks random because it depends on floating-point noise in
the rotation, not on geometry. Trapezoids fail the same way; triangles pass at
every angle tested.

## Root cause

There are two defects. The first creates a spurious ordinate, and the second
turns that into load loss.

### 1. A vertical edge becomes a near-infinite ramp (`get_singularity_functions` / `get_overlap_regions`)

`get_distributed_loads_from_projected_polygons` first calls
`geom_ops.rotate_to_horizontal`, which rotates the member and polygon so the
member lies along +x. After that rotation, an edge that is geometrically
vertical keeps a tiny x-offset between its endpoints:

```
rotated polygon (30 deg):
POLYGON ((5.55e-17 -0.5, 10 -0.4999999999999996, 10 0.5000000000000009,
          -5.55e-17 0.5, 5.55e-17 -0.5))
```

`get_overlap_regions` (line ~316) treats that edge as a real overlap region
spanning `x ∈ [-5.55e-17, 5.55e-17]`, and `overlap_region_to_singularity`
produces a singularity with slope `m ≈ 9.0e15`:

```
Singularity(x0=-5.55e-17, x1=5.55e-17, m=9007199254740992.0, y0=0.0)   <- spurious
Singularity(x0= 5.55e-17, x1=10.0,     m=4.44e-17,           y0=1.0)
```

At 0° there is only one singularity: `Singularity(x0=0.0, x1=10.0, m=0.0, y0=1.0)`.

This step alone doesn't lose any load. `singularities_to_polygon` still gives
the right shape, with area 10.0 and scale ratio 0.1. But the spurious
singularity adds extra x-ordinates near 0:

```
x: [-9.9994e-13, 0.0, 5.55e-17, 9.9994e-13, 1.00006e-12, 9.999999999999, 10.000000000001]
y: [0.0,         0.5, 1.0,      1.0,        1.0,         1.0,            0.0]
```

That is 7 points instead of the usual 4.

### 2. Pairing drops the trailing segment (`get_distributed_loads_from_projected_polygons`, line ~104)

The projected coordinates are turned into `(start, end)` pairs by walking
`projected_poly_coords[1:-1]` two at a time:

```python
for idx, coord in enumerate(projected_poly_coords[1:-1]):
    if idx % 2 == 1:
        inner_pair.append(coord)
        polygon_dist_loads.append(inner_pair)
        inner_pair = []
    else:
        inner_pair.append(coord)
```

This assumes the interior coordinates come in an even number of consecutive
start/end pairs. With the 7 points above, there are 5 interior points. The
pairs become `(0.0 → 5.55e-17)` and `(9.9994e-13 → 1.00006e-12)`, and the last
point, which carries the real segment out to x ≈ 10, is **silently discarded**
as a leftover half-pair. The function never checks that the pairs still cover
the polygon's projected length.

Triangles escape because their sloped edges never produce a sub-1e-12 overlap
region.

## Fixing it

Both should be fixed, because either one alone can lose load:

1. **Ignore degenerate overlap regions.** In `get_overlap_regions` /
   `get_overlap_region`, skip an overlap region whose x-extent is below an
   absolute tolerance (for example `eps`, or 1e-9 in drawing units). That's an
   edge that is vertical after rotation and adds no area. Alternatively, snap
   the rotated geometry (`shapely.set_precision(geom, 1e-9)`) at the end of
   `rotate_to_horizontal`, which is what the papermodels workaround does from
   the outside.
2. **Make the pairing robust.** Build segments from consecutive distinct
   x-ordinates rather than fixed index parity, or at least raise if the number
   of interior points is odd or the pairs don't span the projected polygon's
   x-range. Losing load should never be silent.

**Regression test to add upstream:** for rectangle, trapezoid and triangle
bands about a member rotated through 0–360° (including 89.9°, 90.1°, 179.9°
and 180.1°), the projected load's total area must equal the unit total, and its
start and end must equal 0 and the member length.

## The papermodels workaround

`src/papermodels/geometry/load_projection.py::project_loading_areas` is a
drop-in replacement for `ld.get_distributed_loads_from_projected_polygons`,
used by `LoadedElement._get_distributed_loads` (`datatypes/element.py`). It:

1. does the same rotation (`geom_ops.rotate_to_horizontal`: rotate about the
   member's positive-x start node, which lands at x = 0);
2. snaps the rotated polygons to a 1e-9 grid (`SNAP_GRID`), so vertical edges
   are exactly vertical;
3. passes load-distribution an **exactly** horizontal member, so the library's
   own rotation is the identity and doesn't reintroduce the noise.

Output is identical wherever the library was already correct. The whole
existing papermodels suite passed unchanged when it was switched in.

`tests/test_helper_characterization.py::test_projected_band_distribution_is_rotation_invariant`
pins the behaviour: 3 band shapes × 14 angles, through the workaround.

**Once load-distribution is fixed:**

1. Bump the `load-distribution` requirement in `pyproject.toml`.
2. Point the characterization test at `ld.get_distributed_loads_from_projected_polygons`
   directly and confirm it passes without the workaround.
3. Delete `geometry/load_projection.py`, and call the library directly again in
   `element.py`.
