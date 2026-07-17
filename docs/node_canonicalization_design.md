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
- **Walls are stored in the `STRtree` as centerlines, not polygons** (see §5.4).
  The wall's polygon is retained separately for correspondence (§5.5).

### 5.2 `NodeRegistry` — canonical node identities

```
NodeRegistry
├── node_abs_tol : float                             # absolute snapping tolerance (real-world units); separate from noding_abs_tol (§7)
├── coord   : dict[NodeId, tuple[float, float]]      # canonical coordinate per node
└── _buckets: dict[tuple[int, int], list[NodeId]]    # quantized spatial hash

  get_or_create(point) -> NodeId
      q = quantize(point, node_abs_tol)
      for candidate in neighbor_buckets(q):          # 9-cell neighborhood
          if distance(point, coord[candidate]) <= node_abs_tol:
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

### 5.4 Wall representation: polygon in, centerline in the index

Walls are drawn and stored as **polygons**. This is deliberate — a polygon gives
a forgiving, area-based way to determine *correspondence* between levels
(adjacent PDF pages), which is robust to the wall being drawn slightly
differently from one plane to the next (§5.5).

But a polygon is the wrong thing to intersect for load-path purposes: what a
joist or beam actually bears on is the wall's **centerline**. The historical
sticking point has been computing that centerline. The existing
`geom_ops.get_rectangle_centerline` explodes the polygon's own edges and assumes
a clean, effectively axis-aligned 4-vertex rectangle; it degrades on rotated
walls, walls with vertex noise, or walls with extra/clipped vertices, and the
surrounding orientation logic then falls back to an "is it horizontal or
vertical?" test that fails for anything off-axis.

**Generalized approach (implemented): oriented bounding box.**
`geom_ops.get_wall_centerline` fits the **minimum-area rotated bounding
rectangle** (OBB, `shapely.minimum_rotated_rectangle`) around the wall polygon,
then connects the midpoints of the OBB's two short (opposite) end faces:

- The OBB is a true rectangle regardless of the input's orientation, so the
  result is correct at **any angle** — no axis-alignment assumption.
- Because the OBB is fitted rather than read off the raw edges, it **absorbs
  vertex noise, extra vertices, and clipped corners**; an approximately
  rectangular wall yields the same spine as the clean rectangle would.
- Selecting an **opposite** edge pair (not merely "the two shortest edges")
  means a near-square wall degenerates to a valid long-axis spine instead of a
  diagonal.

Validated on axis-aligned, 30°-rotated, noisy (extra vertex + jitter),
clipped-corner, and near-square walls: all produce a correct long-axis
centerline of the expected length.

In the build pipeline, each wall element contributes its **centerline** to the
`STRtree` and the node registry (so crossings are computed against the spine),
while its polygon is retained on the element for correspondence.

### 5.5 Correspondents stay ratio-based (decided)

Cross-plane *correspondence* remains computed from **polygon area/length
overlap ratios**, not from shared nodes. This is a deliberate decision, not a
limitation:

- Walls (and columns) keep their polygon precisely so correspondence can be
  forgiving to plane-to-plane drawing differences.
- The overlap **ratio** is load-bearing information: when multiple elements
  converge at one location between levels, the ratio provides the basis for
  **choosing a specific load path** (see `prioritize_correspondents` in
  `element.py`). Collapsing correspondence to a shared node would discard the
  signal used to disambiguate.

Correspondents therefore benefit from the shared `STRtree` broad phase (fewer
pairwise checks) but not from the node registry. The node registry governs
*intersections* (same-plane crossings); correspondence governs *between-plane*
transfer and stays ratio-based.

## 6. Build Pipeline

```
1. Ingest + snap to grid            (existing: annotations.py:45, set_precision grid=1e-3)
1a. Wall polygons -> centerlines    (geom_ops.get_wall_centerline; polygon kept for correspondence — §5.4)
2. Tolerance noding (Phase 2)       (extend/trim endpoints within abs_tol — §7)
3. Build STRtree over geometries    (walls contribute their centerline)
4. Broad phase: per geom, query index for candidates; filter by rank + plane
5. Compute each crossing ONCE; register its point via NodeRegistry.get_or_create
6. Build incidence graph keyed by NodeId; derive above/below by incident ranks
7. Build networkx DiGraph edges from incidence (backed by shared NodeIds)
```

Step 5 is the crux: each physical crossing is evaluated exactly once and its
identity is shared thereafter.

## 7. Phase 2 — Tolerance Extend/Trim (Noding)

A preprocessing step that repairs "approximately intersecting" geometry, governed
by an **absolute** tolerance `noding_abs_tol` in real-world units. This is
**separate** from the canonicalization tolerance `node_abs_tol` (§5.2), so
endpoint repair can reach wider than crossing-merging (decided, §9):

- For each line endpoint, query the `STRtree` for geometry within `noding_abs_tol`.
- If the nearest point on a candidate is within `noding_abs_tol`:
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
- Wall-handling call sites that use `get_rectangle_centerline` with the
  axis-aligned assumption (and the `"horizontal"`/`"vertical"` orientation
  fallbacks around it) migrate to `get_wall_centerline`. `get_rectangle_centerline`
  may be retained for genuinely axis-aligned internal box geometry, or removed
  once no longer relied upon.

The existing test fixtures and the geometry-graph tests are the migration safety
net; behavior (graph topology: nodes, edges, intersection counts) must be
preserved.

## 9. Risks & Open Questions

- **Batch vs incremental canonicalization — decided: incremental now, batch
  later.** Phase 1 ships the incremental `get_or_create` + quantized hash form
  (dependency-free; canonical coordinate is "first point wins"). Batch clustering
  (union-find + cluster centroid, requiring a scipy `cKDTree`) is kept as a future
  option to be introduced only if stability testing shows "first point wins"
  produces unstable canonical coordinates.
- **Bucket-boundary straddling.** Two points within `abs_tol` can fall in
  adjacent buckets; the 9-cell neighborhood probe handles the 2D case and must
  be covered by tests.
- **Tolerances — decided: separate.** Noding (Phase 2 extend/trim) and node
  canonicalization use **distinct** tolerances rather than a single shared one,
  so endpoint repair can reach wider than crossing-merging without forcing the
  two to move together. (E.g. `node_abs_tol` for canonicalization,
  `noding_abs_tol` for Phase 2.) Both remain absolute, real-world units.
- **Cascade convergence bound — decided: cap + warning.** Phase 2 iterates to a
  fixed point up to a capped number of passes; if the cap is hit it emits a
  warning and proceeds with the best-so-far arrangement rather than raising.
  Pathological inputs degrade gracefully instead of aborting the build.
- **Correspondents** (cross-plane overlap) stay **ratio-based** — decided; see
  §5.5. The overlap ratio is retained deliberately because it is the basis for
  choosing a load path when multiple elements converge. They benefit from the
  shared STRtree broad phase but not from the node registry.
- **Wall centerline endpoints from the OBB — decided: accept overhang.** OBB
  endpoints extend to the bounding-box extent, which for a clipped/irregular wall
  can slightly overhang the true polygon end. This is accepted as desirable for a
  bearing spine; it should still be confirmed (not clamped) against the
  extent/bearing calculations that consume it.

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

1. **Phase 0 (this doc).** Design + agreement. ✅
2. **Phase 1 (done).** `GeometryModel` + `NodeRegistry` prototype behind the
   existing API, validated against current fixtures; jitter test added.
   Implemented in `src/papermodels/datatypes/geometry_model.py`, tested in
   `tests/test_geometry_model.py`:
   - `NodeRegistry` — incremental `get_or_create`, quantized-hash buckets,
     9-cell neighborhood probe, nearest-wins snapping, "first point wins"
     canonical coordinate (the "incremental now, batch later" decision, §9).
   - `GeometryModel.from_elements` — STRtree broad phase, walls indexed by
     centerline (polygon retained), single-evaluation crossings canonicalized
     through the registry, incidence keyed by `NodeId`, rank-derived
     above/below queries. Not yet wired into `GeometryGraph` (parallel
     prototype).
   - Validation: reproduces all documented crossings of the
     `intersections.pdf` fixture; **jitter test** (perturb every coordinate by
     1e-6…1e-9) shows the intersection topology — node set, edge set, incidence
     multiset — is invariant; build is idempotent.
3. **Phase 2.** Tolerance extend/trim noding preprocessing (`noding_abs_tol`).
4. **Phase 3.** Migrate downstream `Element` methods to read from the model;
   remove obsolete `geom_ops.py` workarounds.
