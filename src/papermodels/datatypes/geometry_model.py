"""
Node canonicalization for the geometry model (Phase 1 prototype).

This module implements the two structures described in
``docs/node_canonicalization_design.md``:

- ``NodeRegistry`` — canonical point identities.  Every computed intersection
  point is inserted via :meth:`NodeRegistry.get_or_create`, which snaps to an
  existing node within ``node_abs_tol`` (an absolute, real-world tolerance) if
  one exists, otherwise registers a new node.  The registry uses an incremental
  quantized-coordinate hash (``get_or_create`` is O(1) amortized) and probes the
  9-cell neighborhood so points straddling a bucket boundary still snap.  This
  is the "incremental now, batch later" decision (design §9): no scipy
  dependency; the canonical coordinate is "first point wins".

- ``GeometryModel`` — the single source of truth over the geometry.  It owns the
  geometries (no copies), their ranks and planes, a ``shapely.STRtree`` for the
  broad phase, the ``NodeRegistry``, and an incidence map (which geometries meet
  at each node).  "Intersecting geometries share the same node location" becomes
  true *by construction*: each physical crossing is evaluated exactly once and
  its identity is shared thereafter.

Scope (Phase 1): this is a prototype that lives *behind* the existing API — it
is not yet wired into ``GeometryGraph``.  It governs same-plane *intersections*
(the node registry's job).  Cross-plane *correspondence* stays ratio-based in
the existing layer (design §5.5) and is intentionally not modelled here.

Walls (linear-reaction polygons) contribute their **centerline** to the index
and node registry (design §5.4) while their polygon is retained for future
correspondence use.
"""

from __future__ import annotations

import math
from typing import Iterable, NamedTuple, Optional

from shapely import STRtree, Point
from shapely.geometry.base import BaseGeometry

from ..geometry import geom_ops
from ..geometry.noding import NodingReport, node_geometries

# A geometry's identity in the model.  In practice this is the element ``tag``.
GeomId = str
# A canonical node's identity.  Small, dense, assigned incrementally.
NodeId = int
# A plane (page) identity; whatever ``Element.plane_id`` holds.
PlaneId = object

# Role of a geometry at a node it is incident to.
ROLE_ENDPOINT = "endpoint"
ROLE_INTERIOR = "interior"
ROLE_BOUNDARY = "boundary"

DEFAULT_NODE_ABS_TOL = 1e-3


def _as_xy(point) -> tuple[float, float]:
    """Return ``(x, y)`` from a shapely Point or a coordinate tuple."""
    if isinstance(point, Point):
        return (point.x, point.y)
    x, y = point[0], point[1]
    return (float(x), float(y))


class NodeRegistry:
    """
    Canonical point identities under an absolute snapping tolerance.

    Points inserted within ``node_abs_tol`` of an existing node snap to it (the
    nearest one); otherwise a new node id is created.  The canonical coordinate
    is that of the first point registered for the node ("first point wins").
    """

    def __init__(self, node_abs_tol: float = DEFAULT_NODE_ABS_TOL):
        if node_abs_tol <= 0:
            raise ValueError(f"node_abs_tol must be > 0, got {node_abs_tol!r}")
        self.node_abs_tol = float(node_abs_tol)
        # Canonical coordinate per node.
        self.coord: dict[NodeId, tuple[float, float]] = {}
        # Quantized spatial hash: cell -> node ids whose coord falls in that cell.
        self._buckets: dict[tuple[int, int], list[NodeId]] = {}
        self._next_id = 0

    def _cell(self, x: float, y: float) -> tuple[int, int]:
        t = self.node_abs_tol
        return (math.floor(x / t), math.floor(y / t))

    def get_or_create(self, point) -> NodeId:
        """
        Return the id of the canonical node within ``node_abs_tol`` of ``point``,
        creating a new node (at ``point``) if none exists.

        The 9-cell neighborhood is probed: two points within ``node_abs_tol`` can
        differ by at most ``node_abs_tol`` in each axis, so their quantized cells
        differ by at most one cell in each axis.  This is what catches points
        straddling a bucket boundary (design §9).
        """
        x, y = _as_xy(point)
        cx, cy = self._cell(x, y)
        tol2 = self.node_abs_tol * self.node_abs_tol
        best_id: Optional[NodeId] = None
        best_d2: Optional[float] = None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for cand in self._buckets.get((cx + dx, cy + dy), ()):
                    ex, ey = self.coord[cand]
                    d2 = (ex - x) ** 2 + (ey - y) ** 2
                    if d2 <= tol2 and (best_d2 is None or d2 < best_d2):
                        best_id, best_d2 = cand, d2
        if best_id is not None:
            return best_id
        node_id = self._next_id
        self._next_id += 1
        self.coord[node_id] = (x, y)
        self._buckets.setdefault((cx, cy), []).append(node_id)
        return node_id

    def __len__(self) -> int:
        return len(self.coord)

    def __contains__(self, node_id: NodeId) -> bool:
        return node_id in self.coord


class Incidence(NamedTuple):
    """A geometry that meets at a node, tagged by its rank and its role there."""

    geom_id: GeomId
    rank: int
    role: str  # one of ROLE_ENDPOINT / ROLE_INTERIOR / ROLE_BOUNDARY


def _role(geom: BaseGeometry, x: float, y: float, tol: float) -> str:
    """Classify how ``geom`` meets the point ``(x, y)``."""
    p = Point(x, y)
    if geom.geom_type == "LineString":
        coords = geom.coords
        for c in (coords[0], coords[-1]):
            if p.distance(Point(c)) <= tol:
                return ROLE_ENDPOINT
        return ROLE_INTERIOR
    if geom.geom_type == "Polygon":
        if p.distance(geom.exterior) <= tol:
            return ROLE_BOUNDARY
        return ROLE_INTERIOR
    return ROLE_INTERIOR


def _crossing_points(geom: BaseGeometry) -> list[tuple[float, float]]:
    """
    Reduce an intersection result to the point(s) that become canonical nodes.

    A crossing that is a point is one node; a segment (collinear overlap or a
    line passing through a polygon) contributes its endpoints; a region
    contributes its representative point.
    """
    gt = geom.geom_type
    if gt == "Point":
        return [(geom.x, geom.y)]
    if gt == "MultiPoint":
        return [(p.x, p.y) for p in geom.geoms]
    if gt == "LineString":
        coords = list(geom.coords)
        return [tuple(coords[0]), tuple(coords[-1])]
    if gt == "MultiLineString":
        out: list[tuple[float, float]] = []
        for ls in geom.geoms:
            coords = list(ls.coords)
            out.append(tuple(coords[0]))
            out.append(tuple(coords[-1]))
        return out
    if gt in ("Polygon", "MultiPolygon"):
        rp = geom.representative_point()
        return [(rp.x, rp.y)]
    if gt == "GeometryCollection":
        out = []
        for g in geom.geoms:
            out.extend(_crossing_points(g))
        return out
    rp = geom.representative_point()
    return [(rp.x, rp.y)]


class GeometryModel:
    """
    Single source of truth over the geometry (design §5.1).

    Owns geometries (no copies), ranks, planes, an ``STRtree`` broad-phase index,
    a ``NodeRegistry`` of canonical node identities, and an incidence map keyed by
    ``NodeId``.  Same-plane crossings are computed once, at build time, and their
    node identities are shared by every incident geometry.
    """

    def __init__(
        self,
        node_abs_tol: float = DEFAULT_NODE_ABS_TOL,
        nodes: Optional[NodeRegistry] = None,
    ):
        self.node_abs_tol = float(node_abs_tol)
        # GeomId -> geometry used for the index (walls store their centerline).
        self.geometries: dict[GeomId, BaseGeometry] = {}
        # GeomId -> retained polygon (walls, point-reaction polygons).
        self.polygons: dict[GeomId, BaseGeometry] = {}
        self.ranks: dict[GeomId, int] = {}
        self.planes: dict[GeomId, PlaneId] = {}
        self.reaction_types: dict[GeomId, str] = {}
        # A caller may share its registry so that nodes created here and nodes
        # created elsewhere in the build (e.g. element intersections) are one set.
        self.nodes = nodes if nodes is not None else NodeRegistry(node_abs_tol)
        # NodeId -> list[Incidence] (one entry per distinct incident geom).
        self.incidence: dict[NodeId, list[Incidence]] = {}
        # GeomId -> set of NodeIds it is incident to.
        self.geom_nodes: dict[GeomId, set[NodeId]] = {}
        self.index: Optional[STRtree] = None
        self._index_geom_ids: list[GeomId] = []
        # Set when geometry is added after the index was built; the STRtree is
        # immutable, so it is rebuilt lazily on the next query.
        self._index_stale = False
        # Set when built with noding enabled (Phase 2); None otherwise.
        self.noding_report: Optional[NodingReport] = None

    # -- construction -------------------------------------------------------

    @classmethod
    def from_elements(
        cls,
        elements: Iterable,
        node_abs_tol: float = DEFAULT_NODE_ABS_TOL,
        noding_abs_tol: Optional[float] = None,
        noding_max_passes: int = 10,
        suppress_warnings: bool = False,
        nodes: Optional[NodeRegistry] = None,
        seed_points: Iterable = (),
    ) -> "GeometryModel":
        """
        Build a model from ``Element``-like objects.  Each must expose
        ``tag``, ``geometry``, ``rank``, ``plane_id`` and ``reaction_type``.

        If ``noding_abs_tol`` is given, Phase 2 tolerance extend/trim noding runs
        over the index geometries *before* the crossing build (design §6/§7): line
        endpoints within ``noding_abs_tol`` of nearby geometry are snapped onto it
        so the arrangement is topologically clean when crossings are
        canonicalized.  Wall centerlines are held fixed (valid snap targets but
        never moved).  ``noding_abs_tol`` is separate from ``node_abs_tol`` (§9).

        ``nodes`` shares an existing registry.  ``seed_points`` are registered
        before any crossing is built; because the registry is "first point
        wins", crossings computed here then reuse those exact coordinates (e.g.
        the canonical intersection regions already stored on the elements).
        """
        self = cls(node_abs_tol, nodes=nodes)
        for point in seed_points:
            self.nodes.get_or_create(point)
        for element in elements:
            self._add_element(element)
        if noding_abs_tol is not None:
            self._apply_noding(noding_abs_tol, noding_max_passes, suppress_warnings)
        self._build_index()
        self._build_incidence()
        return self

    def _apply_noding(
        self, noding_abs_tol: float, max_passes: int, suppress_warnings: bool
    ) -> None:
        # Walls contribute a derived centerline spine; keep those fixed so noding
        # snaps other lines onto them without dragging the spine around (§7).
        fixed_ids = {
            gid
            for gid in self.geometries
            if self.reaction_types.get(gid) == "linear" and gid in self.polygons
        }
        self.geometries, self.noding_report = node_geometries(
            self.geometries,
            noding_abs_tol,
            max_passes=max_passes,
            fixed_ids=fixed_ids,
            suppress_warnings=suppress_warnings,
        )

    def _add_element(self, element) -> None:
        gid = element.tag
        geom = element.geometry
        self.ranks[gid] = element.rank
        self.planes[gid] = element.plane_id
        self.reaction_types[gid] = element.reaction_type
        if geom.geom_type == "Polygon":
            self.polygons[gid] = geom
            # Walls (linear reaction) are indexed by their centerline; a
            # point-reaction polygon (e.g. a column) is indexed as-is (§5.4).
            if element.reaction_type == "linear":
                index_geom = geom_ops.get_wall_centerline(geom)
            else:
                index_geom = geom
        else:
            index_geom = geom
        self.geometries[gid] = index_geom
        self.geom_nodes[gid] = set()

    def _build_index(self) -> None:
        self._index_geom_ids = list(self.geometries.keys())
        self.index = STRtree([self.geometries[gid] for gid in self._index_geom_ids])
        self._index_stale = False

    def _ensure_index(self) -> None:
        if self.index is None or self._index_stale:
            self._build_index()

    def _build_incidence(self) -> None:
        ids = self._index_geom_ids
        geoms = [self.geometries[gid] for gid in ids]
        seen_pairs: set[tuple[GeomId, GeomId]] = set()
        for i, gid_i in enumerate(ids):
            geom_i = geoms[i]
            for j in self.index.query(geom_i):
                gid_j = ids[int(j)]
                if gid_i == gid_j:
                    continue
                pair = (gid_i, gid_j) if gid_i <= gid_j else (gid_j, gid_i)
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                if self.planes[gid_i] != self.planes[gid_j]:
                    continue
                self._register_crossing(gid_i, geom_i, gid_j, geoms[int(j)])

    def _register_crossing(self, gid_i: GeomId, geom_i, gid_j: GeomId, geom_j) -> None:
        crossing = geom_i.intersection(geom_j)
        if crossing.is_empty:
            return
        tol = self.node_abs_tol
        for x, y in _crossing_points(crossing):
            node_id = self.nodes.get_or_create((x, y))
            cx, cy = self.nodes.coord[node_id]
            self.geom_nodes[gid_i].add(node_id)
            self.geom_nodes[gid_j].add(node_id)
            incs = self.incidence.setdefault(node_id, [])
            self._touch(incs, gid_i, geom_i, cx, cy, tol)
            self._touch(incs, gid_j, geom_j, cx, cy, tol)

    def _touch(self, incs, gid, geom, x, y, tol) -> None:
        if any(inc.geom_id == gid for inc in incs):
            return
        incs.append(Incidence(gid, self.ranks[gid], _role(geom, x, y, tol)))

    # -- writing generated geometry -----------------------------------------

    def add_geometry(
        self,
        gid: GeomId,
        geom: BaseGeometry,
        rank: int,
        plane: PlaneId,
        reaction_type: str = "point",
        crossings: Optional[dict[NodeId, GeomId]] = None,
    ) -> None:
        """
        Write a generated geometry (e.g. a joist produced by a joist array) into
        the model.

        ``crossings`` maps each canonical node the new geometry meets (already
        obtained from ``self.nodes``) to the geometry it meets there; incidence
        is recorded for both sides.  Crossings are *not* recomputed: the caller
        constructed the geometry from those nodes.
        """
        if gid in self.geometries:
            raise ValueError(f"Geometry id {gid!r} already exists in the model.")
        self.ranks[gid] = rank
        self.planes[gid] = plane
        self.reaction_types[gid] = reaction_type
        self.geometries[gid] = geom
        self.geom_nodes[gid] = set()
        self._index_stale = True
        tol = self.node_abs_tol
        for node_id, other_gid in (crossings or {}).items():
            if node_id not in self.nodes:
                raise KeyError(f"Node {node_id} is not in the registry.")
            x, y = self.nodes.coord[node_id]
            self.geom_nodes[gid].add(node_id)
            self.geom_nodes[other_gid].add(node_id)
            incs = self.incidence.setdefault(node_id, [])
            self._touch(incs, gid, geom, x, y, tol)
            self._touch(incs, other_gid, self.geometries[other_gid], x, y, tol)

    # -- queries ------------------------------------------------------------

    def source_geometry(self, gid: GeomId) -> BaseGeometry:
        """The element's own geometry: the polygon for walls, else the index geometry."""
        return self.polygons.get(gid, self.geometries[gid])

    def is_bearing_support(self, gid: GeomId) -> bool:
        """
        True if a joist can bear on ``gid``: a line (beam) or a linear-reaction
        polygon (wall, indexed by its centerline).  Point-reaction polygons
        (columns/posts) are not joist supports.
        """
        if gid in self.polygons:
            return self.reaction_types.get(gid) == "linear"
        return self.geometries[gid].geom_type == "LineString"

    def query_supports(
        self, region: BaseGeometry, plane: PlaneId, rank: int
    ) -> list[GeomId]:
        """
        Geometry ids that can support a member of ``rank`` on ``plane`` and whose
        index geometry (wall centerline / beam line) intersects ``region``.

        Returned in model insertion order so results are deterministic.
        """
        self._ensure_index()
        hits = set(int(i) for i in self.index.query(region, predicate="intersects"))
        out = []
        for pos, gid in enumerate(self._index_geom_ids):
            if pos not in hits:
                continue
            if self.planes[gid] != plane or not self.ranks[gid] > rank:
                continue
            if self.is_bearing_support(gid):
                out.append(gid)
        return out

    def incident_geoms(self, node_id: NodeId) -> list[Incidence]:
        return self.incidence.get(node_id, [])

    def intersections_below(self, gid: GeomId) -> dict[GeomId, set[NodeId]]:
        """Higher-rank geoms sharing a node with ``gid`` (gid transfers down)."""
        return self._incident_by_rank(gid, above=False)

    def intersections_above(self, gid: GeomId) -> dict[GeomId, set[NodeId]]:
        """Lower-rank geoms sharing a node with ``gid`` (they transfer to gid)."""
        return self._incident_by_rank(gid, above=True)

    def _incident_by_rank(self, gid: GeomId, above: bool) -> dict[GeomId, set[NodeId]]:
        rank = self.ranks[gid]
        out: dict[GeomId, set[NodeId]] = {}
        for node_id in self.geom_nodes[gid]:
            for inc in self.incidence[node_id]:
                if inc.geom_id == gid:
                    continue
                if (inc.rank < rank) if above else (inc.rank > rank):
                    out.setdefault(inc.geom_id, set()).add(node_id)
        return out

    def intersection_edges(self) -> set[tuple[GeomId, GeomId]]:
        """Directed edges ``(lower_rank, higher_rank)`` sharing a node."""
        edges: set[tuple[GeomId, GeomId]] = set()
        for gid in self.geometries:
            for other in self.intersections_below(gid):
                edges.add((gid, other))
        return edges

    def topology_signature(self) -> dict:
        """
        A jitter-invariant summary of the built arrangement (design §10).

        Comparing two signatures answers "same node set, same edge set, same
        incidence structure?" independent of the exact floating-point
        coordinates, which is the acceptance criterion for the whole effort.
        """
        incidence_multiset = sorted(
            tuple(sorted((inc.geom_id, inc.rank) for inc in incs))
            for incs in self.incidence.values()
        )
        return {
            "n_nodes": len(self.nodes),
            "edges": self.intersection_edges(),
            "incidence": incidence_multiset,
        }
