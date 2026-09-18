# Design: JoistArrayModel on stations and canonical nodes

**Status:** Implemented
**Branch:** `refactor/joist_array_model`
**Author:** Connor Ferster
**Date:** 2026-09-18

---

## 1. Problem

A joist array turns one drawn **joist prototype** into the joists of a spread.
The previous `JoistArrayModel` built each joist by projecting rays from a
translated centroid, intersecting them with copies of the supports, adding a
fixed cantilever, and then testing the result against the supports again. On the
test fixtures this produced:

| Symptom | Fixture evidence |
|---|---|
| A joist lands in a gap between supports and gets the wrong supports, or is silently dropped | `collector_extents.pdf`: `SJ0.0-5` was truncated at WT0.3; `SJ0.0-10` was missing |
| Joists at the array edge fail to build | `resi_dormers.pdf`: `RJ0.0-0` missing and `RJ0.0-6` a GeometryError; `sketch_to_scale.pdf`: `J1.1-10` dropped as "only one support" |
| The last gap exceeds the spacing (up to 1.5x) | `test_wall_point_load_locations` failed with a 1.114 gap at 1.0 spacing |
| Cantilevers projected into the span; arrays generated right to left | Prototypes tilted past vertical (positive-x bias; §3.1) |
| Cantilevers could not vary; supports assumed orthogonal | `"horizontal"`/`"vertical"` branches throughout |
| Not part of node canonicalization | Supports were per-element copies; crossings were recomputed with `grid_size=1e-3` |

The `CollectorTribModel` gave one whole-spread trib area per prototype, so it
could not represent any of this either. It has been removed.

## 2. Goals

1. A joist that lands where supports are missing moves to the nearest valid
   station, and trib areas leave no gaps between adjacent joists.
2. A `joist_container` polygon markup lets cantilevers **and** backspans vary
   along the array.
3. The model reads supports from, and writes joists into, the shared
   `GeometryModel` (STRtree plus `NodeRegistry`), so no crossing is computed twice.

## 3. Design

The engine lives in `src/papermodels/geometry/joist_array.py`, and
`JoistArrayModel` (`datatypes/joist_models.py`) orchestrates it.

### 3.1 One frame, oriented once

Work happens in the prototype's own frame. `u` runs along the joist. `n` is the
normal, which is the **array axis**. A **station** `s` is a position along `n`
and `t` is a position along `u`, so the joist line at station `s` is
`origin + s·n + t·u`. Nothing depends on axis alignment.

The rest of the library uses the **positive-x bias** (`order_nodes_positive`),
which flips its orientation at vertical. Near vertical, the 1e-3 grid snapping
produces exact x-ties, and the ascending-y tiebreak then disagrees with the
joist's direction. That is how cantilevers ended up projected into the span.
The engine avoids it as follows:

- `n` is set by a **dominant-axis rule**: `+x` when `|n_x| > |n_y|`, else `+y`.
  Vertical joists array left to right; horizontal joists array bottom to top.
  `u = rotate_cw(n)`. Every orientation rule has a discontinuity somewhere; this
  one puts it at 45°, where drawings rarely sit and grid ties don't occur.
- The origin is the prototype's midpoint, so the frame is the same whichever way
  the prototype was drawn.
- Inside the engine, points are ordered **only** by their `t`/`s` projections.
  End A is always the smaller-`t` end.
- The positive-x bias is applied only at the output boundary, where each
  subelement's LineString is written in positive-x order for downstream code.

### 3.2 Array region and supports

| Mode | Triggered by | Region R | Station range |
|---|---|---|---|
| container | a `joist_container` polygon around the prototype | the container | the container's |
| extent | an extent line crossing the prototype | the band swept by the prototype and its supports along the extent line (`prototype_extent_region`, not an axis-aligned bbox) | the extent line's |
| plain | neither | the quadrilateral between the two outer supports, plus the cantilevers | where **all** the prototype's supports exist; stops where the outer supports converge (a triangle apex) |

The supports are the parent element's graph supports (`intersections_below`),
read from the shared model: wall centerlines and beam lines, never copies. A
support is **present** at a station if its drawn geometry (the wall polygon)
is inside R there, bounded by its centerline's own extent. So a container drawn
to a wall's face still catches the wall.

The graph build uses the same function (`region_support_crossing`) for the
parent element's intersections, so the graph and the array cannot disagree
about which supports an array bears on. Graph edges were unchanged on all 9
fixtures when this was switched over.

### 3.3 Breakpoints, stations and gap resolution

Projecting every support's endpoints inside R onto `n` gives **breakpoints**.
Between two adjacent breakpoints the set of supports is constant, so each
interval is classified once. Nothing is re-tested per joist.

Target stations start at the array's start edge, at `spacing` apart. No gap
exceeds the spacing, and `joist_at_start`, `joist_at_end` and `initial_offset`
are all honoured.

A target station is **invalid** if its interval has fewer than two supports, or
if its support count differs from **both** neighbours' counts. An invalid
station moves to the nearest valid station, found by a `bisect` over the
breakpoints followed by an outward scan. The scan stays within half a spacing
of the target and never passes a neighbour, and it insets from the edge of the
valid interval by `min_bearing`. If nothing valid is in the window, the joist is
dropped with a warning.

Intermediate supports that come and go (WT0.3 in `collector_extents`) don't
trigger relocation, because their neighbours have the same count.

### 3.4 Joists and canonical crossings

- **extent and plain modes:** each joist runs between the outer support crossings
  at its station, plus the prototype's constant signed cantilevers. A
  cantilever below `cantilever_tolerance`, or a negative one (the prototype
  drawn to a wall face), becomes 0, so every joist reaches its supports.
- **container mode:** each joist is the container cut along its station line,
  extended to the outer supports if the container stops short of them.

Each crossing is interpolated **along the support**, so it lies on the support,
and is snapped through the shared `NodeRegistry`. A joist end that sits on a
support *is* that node. The subelement's `Intersection` records are built
directly from these nodes and are never re-intersected. Each record also
carries its precomputed extent along the support, so the graph doesn't call the
axis-aligned `get_collector_extents` path for subelements.

`GeometryGraph.build_geometry_model()` runs after gravity-frame processing,
which moves beam ends onto column centroids. It seeds the registry with every
existing intersection coordinate. Generated joists are written back with
`add_geometry(..., parent=...)`, and re-assigning collector behaviour replaces
them.

### 3.5 Trib areas

With final stations `s_0 … s_k`, joist *i*'s trib area is the band between the
midpoints to its neighbours (the array's station limits at the two ends),
clipped to R. Neighbouring bands share one boundary value, so the bands **tile R
exactly**, including around relocated joists.

### 3.6 Span jumps

A span can only jump where an outer support changes identity. The new support's
crossing is compared with the old support's line extended to the same station.
If it moves by more than `span_jump_tol` (default: half the spacing), the model
warns with the tag and station and suggests a separate prototype, then
continues.

## 4. Joist container markup

Add a legend entry whose `Type` contains `container`, for example:

```
Legend
Type: Joist Container
```

Then draw a polygon around a joist prototype. The prototype is matched to the
container holding most of its length on the same page. The model warns if a
container holds zero or several prototypes, and if a prototype has both an
extent line and a container (the container wins).

## 5. Helper audit

No existing helper was reused on trust. The audit is pinned by
`tests/test_helper_characterization.py` across orientations.

| Helper | Finding | Outcome |
|---|---|---|
| `create_extent_polygon` | axis-aligned bbox | replaced by `prototype_extent_region`; removed |
| `get_joist_locations` | line-to-line distance; ignored `joist_at_end`; gaps up to 1.5 × spacing | replaced by `target_stations`; removed |
| `get_intersection`, extent branch | OBB centerline of the *clipped* wall piece turns perpendicular for narrow clips | replaced by `region_support_crossing` for arrays |
| `get_wall_centerline` (OBB) | correct at any angle; tolerates noise; near-square gives an axis | kept |
| `rotate_90_vector` | exact quarter turn | kept |
| `load_distribution.get_distributed_loads_from_projected_polygons` | **defect:** after rotation a vertical edge keeps a ~1e-17 x-offset, which becomes a ~1e16 slope, so rectangular and trapezoidal bands about a 30° member project to zero load | guarded locally by `geometry.load_projection.project_loading_areas` (same rotation, snap to 1e-9, horizontal member); **fix upstream** |

## 6. Behaviour changes

- Horizontal joists now array bottom to top, so load lists on the supports come
  out in the opposite order. Order-dependent tests were made order-independent.
- Plain arrays use only the prototype's own supports, and extent arrays use the
  graph's supports inside the region. Both warn about in-region supports the
  graph didn't link.
- Point-reaction polygons linked to a prototype are borne along their long axis,
  with a warning. Drawings made before the `Reaction Type` legend field draw
  walls this way.
- Subelement tags are sequential (`TAG-0 … TAG-n`). Dropped joists no longer
  leave holes in the numbering.

## 7. Open questions

- **The literal neighbour rule** also moves a joist off a short intermediate
  beam that sits under that one joist only (3 supports against 2 and 2). This
  follows the spec, and is pinned by
  `test_support_count_differing_from_both_neighbours_relocates`.
- **Array ends** only require two supports, because the neighbour rule needs two
  neighbours. A wrong-pair end joist is still caught by the span-jump warning.
- **DXF** gained extent and container routing but has no test fixture.
- **The Example notebooks** have unrelated API drift (`LoadedElement.dump_toml`,
  element tags, and a KeyError in Example 5). Those failures also occur on `main`.
