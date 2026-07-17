# Design: Node Canonicalization for the Geometry Model

**Status:** Draft
**Branch:** `node_canon/design_document`
**Author:** Connor Ferster
**Date:** 2026-07-17

---

## 1. Problem

`papermodels` converts marked-up drawings into a directed "geometry graph" where
lower-rank elements (joists, beams) transfer load down to higher-rank elements
(walls, columns, footings). The load path is derived from **where geometries
intersect** and **which geometries correspond** between adjacent planes.

The library works but is fragile under floating-point error. The fragility is not
caused by GEOS/shapely predicates (`intersects`, `contains`) — those are robust.
It is caused by **construction followed by re-testing**: every `intersection()`,
`nearest_points()`, `interpolate()`, `centroid`, and affine rotation produces
coordinates rounded to the nearest double that almost never lie exactly on the
input geometry, and the result is then fed back into an exact predicate or an
exact float comparison, where a branch flips or an `assert` fires.

The deeper structural cause is that **the same physical crossing is computed many
times, independently, and the results are never reconciled.** Concretely, in
`element.py::get_geometry_intersections` (the O(n²) build loop):

- For each ordered pair with `j_rank > i_rank`, `geom_ops.get_intersection` is
  called **twice** — once to store the crossing on `i` as an
  `intersection_below`, once to store it on `j` as an `intersection_above`
  (`element.py:1152-1164`). These two evaluations of one physical point can
  differ in their low bits.
- Each `Intersection` NamedTuple carries a **copy** of `other_geometry`
  (`element.py:43`). The same wall is copied into every element that touches it.
- Downstream methods (`get_transfer_extents`, `align_frames_to_centroids`,
  `_get_support_locations`, `get_collector_extents`) recompute
  `self.geometry.intersection(other_geometry)` **again**, producing further
  independent evaluations of the same crossing.

A single physical node therefore exists as several slightly-different floats
scattered across the model, with nothing forcing them to agree. The scattered
tolerance workarounds in `geom_ops.py` (`grid_size=1e-3`, `eps=1e-6`,
`buffer(1e-3)`, `rel_tol=5e-2`, exact `== 0.0` / `== 1.0` branch tests) patch
individual symptoms but do not remove the cause; the instability simply moves.

## 2. Goal

Establish a **single canonical identity for every node** (intersection point) in
the model, computed **once**, shared by every element incident to it, and
accelerated by spatial indexing so the model can be built and queried
efficiently.

"Intersecting geometries share the same node location" must become true **by
construction**, not by hoping two independent `.intersection()` calls round the
same way.

### Non-goals

- Replacing shapely or moving to exact arithmetic (CGAL, rational coordinates).
  The domain (drawing markup, ~1 mm meaningful precision, mostly rectangles and
  straight lines) does not justify the dependency or performance cost.
- Changing the structural/analysis semantics (rank, reaction types, load
  transfer). This work changes *how geometry identity is managed*, not what a
  beam or column means.
- Changing the `networkx` DiGraph as the element→element structural graph. It
  stays; its edges become backed by shared node ids.

## 3. Key Insight

The spatial index is **not** the fix — it is only acceleration. The fix is
**node canonicalization**: a single registry into which every computed
intersection point is inserted, snapping to an existing node within an absolute
tolerance if one exists, otherwise registering a new node. Every element meeting
at that crossing then references the **same node id**.

This is equivalent to building a *noded arrangement* (a planar graph with shared
vertices) over the input geometry.

## 4. Data-structure choice (and one correction)

A **KD-tree indexes points** and answers nearest-neighbor / radius queries. It
does **not** index extended geometry (line segments, polygons) for intersection
queries. The two halves of this problem need two different structures:

| Need | Structure | Notes |
|------|-----------|-------|
| "Which lines/polygons cross this one?" (broad phase) | **`shapely.STRtree`** (R-tree) | Already available in shapely ≥ 2.0, GEOS-backed. Replaces the O(n²) loop. **No new dependency.** |
| "Is there already a node within `abs_tol` of this point?" | **Quantized-coordinate hash** (preferred) or KD-tree | Nodes are points. A dict keyed by coordinates quantized to `abs_tol` gives O(1) incremental insert/lookup; probe the 9-cell neighborhood to catch points straddling a bucket boundary. A `scipy.spatial.cKDTree` is only worth adding for a one-shot batch clustering pass. |

The "lines/polygons KD-tree" from the original idea is really an `STRtree`. The
"nodes KD-tree" is a real point index, but a quantized hash is a better fit
because node insertion is incremental.

## 5. Proposed Architecture

### 5.1 `GeometryModel` — the single source of truth

A new class owns everything geometry-related in one place:

```
GeometryModel
├── geometries : dict[GeomId, shapely geometry]     # single source of truth, NO copies
├── ranks      : dict[GeomId, int]
├── planes     : dict[GeomId, PlaneId]              # page/plane id
├── index      : shapely.STRtree                    # broad-phase over geometries
├── nodes      : NodeRegistry                        # canonical node identities
└── incidence  : dict[NodeId, list[Incidence]]      # which geoms meet at each node
```

- `Element` stops owning copies of other elements' geometry. It becomes a
  **view**: it holds its own `GeomId` and references `NodeId`s, and answers
  `intersections_below()` / `corresponding_above()` by querying the model.
- The `networkx` DiGraph is retained as the structural graph; its edges are
  backed by shared `NodeId`s.

### 5.2 `NodeRegistry` — canonical node identities

```
NodeRegistry
├── abs_tol : float                                  # absolute snapping tolerance (real-world units)
├── coord   : dict[NodeId, tuple[float, float]]      # canonical coordinate per node
└── _buckets: dict[tuple[int, int], list[NodeId]]    # quantized spatial hash

  get_or_create(point) -> NodeId
      q = quantize(point, abs_tol)
      for candidate in neighbor_buckets(q):          # 9-cell neighborhood
          if distance(point, coord[candidate]) <= abs_tol:
              return candidate                        # snap to existing node
      return register_new(point)                      # otherwise new node
```

Canonicalization can also be done as a **batch clustering** pass (collect all raw
intersection points, cluster within `abs_tol` via union-find over radius-neighbor
pairs, assign one canonical coordinate — e.g. the cluster centroid — per cluster).
The batch form is where a `cKDTree` would earn its place; the incremental
`get_or_create` form needs no new dependency.

### 5.3 `Incidence` — what meets at a node

```
Incidence(geom_id: GeomId, rank: int, role: str)   # role ∈ {"endpoint", "interior", "boundary"}
```

`intersections_above` / `intersections_below` collapse into a single query:
"elements incident to this node, filtered by rank." One structure replaces the
two independently-computed lists maintained today.

## 6. Build Pipeline

```
1. Ingest + snap to grid            (existing: annotations.py:45, set_precision grid=1e-3)
2. Tolerance noding (Phase 2)       (extend/trim endpoints within abs_tol — §7)
3. Build STRtree over geometries
4. Broad phase: per geom, query index for candidates; filter by rank + plane
5. Compute each crossing ONCE; register its point via NodeRegistry.get_or_create
6. Build incidence graph keyed by NodeId; derive above/below by incident ranks
7. Build networkx DiGraph edges from incidence (backed by shared NodeIds)
```

Step 5 is the crux: each physical crossing is evaluated exactly once and its
identity is shared thereafter.

## 7. Phase 2 — Tolerance Extend/Trim (Noding)

A preprocessing step that repairs "approximately intersecting" geometry, governed
by a single **absolute** tolerance `abs_tol` in real-world units:

- For each line endpoint, query the `STRtree` for geometry within `abs_tol`.
- If the nearest point on a candidate is within `abs_tol`:
  - **extend** the line to that point if it falls short, or
  - **trim** the line to that point if it overshoots.

Absolute (not relative) tolerance is deliberate: a relative tolerance makes
endpoint behavior depend on line length, reintroducing the scale-dependence we
are trying to remove.

**Caveat — cascading snaps.** Moving endpoint A onto line L can pull it off line
M it previously touched. A single sweep over a sparse, mostly-orthogonal
structural drawing is usually sufficient, but the operation is not guaranteed to
converge in one pass. The implementation must iterate to a fixed point (or a
capped number of passes, warning if the cap is hit).

**Ordering.** Phase 2 runs *before* the STRtree/node-registry build (steps 3-6),
so the arrangement is already topologically clean when crossings are
canonicalized.

## 8. Migration / Impact

This is a rewrite of the intersection-building layer, not an add-on. Call sites
that today "compute a crossing" must change to "look up a node":

- `element.py::get_geometry_intersections` — the O(n²) double loop → STRtree
  broad phase + single-evaluation canonicalization.
- `Intersection` / `Correspondent` construction — stop copying `other_geometry`;
  reference `GeomId` / `NodeId`.
- Downstream recomputation sites that call `self.geometry.intersection(other)` —
  `get_transfer_extents`, `align_frames_to_centroids`, `_get_support_locations`,
  `get_collector_extents` — read canonical nodes from the model instead.
- Many `geom_ops.py` snapping workarounds become unnecessary once every
  coordinate is a canonical node, and can be removed in a follow-up cleanup.

The existing test fixtures and the geometry-graph tests are the migration safety
net; behavior (graph topology: nodes, edges, intersection counts) must be
preserved.

## 9. Risks & Open Questions

- **Batch vs incremental canonicalization.** Batch clustering gives the most
  stable canonical coordinate (centroid of a cluster) but needs all points up
  front; incremental `get_or_create` is simpler and dependency-free but the
  canonical coordinate is "first point wins." Decide per performance/stability
  needs. *(Leaning incremental + quantized hash to avoid a scipy dependency.)*
- **Bucket-boundary straddling.** Two points within `abs_tol` can fall in
  adjacent buckets; the 9-cell neighborhood probe handles the 2D case and must
  be covered by tests.
- **Single `abs_tol` for both noding and canonicalization?** Likely yes (one
  precision model), but they could be separated if noding needs to be more
  aggressive than crossing-merging.
- **Cascade convergence bound** for Phase 2 on pathological inputs — cap +
  warning vs. hard error.
- **Correspondents** (cross-plane overlap) currently use area/length ratios, not
  point crossings. They benefit from the shared STRtree but not directly from the
  node registry; confirm whether they stay ratio-based.

## 10. Validation Strategy

- **Preserve existing tests.** All current geometry-graph and integration tests
  must pass unchanged.
- **Jitter test.** Perturb every input coordinate in the fixtures by
  ±1e-6…±1e-9 and assert the resulting graph topology (node set, edge set,
  intersection counts) is invariant. This converts "works on my drawings" into a
  checkable robustness property and is the acceptance criterion for the whole
  effort.
- **Idempotence.** Running canonicalization twice yields identical node ids and
  coordinates.

## 11. Proposed Phasing

1. **Phase 0 (this doc).** Design + agreement.
2. **Phase 1.** `GeometryModel` + `NodeRegistry` prototype behind the existing
   API, validated against current fixtures; add the jitter test.
3. **Phase 2.** Tolerance extend/trim noding preprocessing.
4. **Phase 3.** Migrate downstream `Element` methods to read from the model;
   remove obsolete `geom_ops.py` workarounds.
```
