"""
Joist models: turn one drawn joist prototype into the joists of an array.

``JoistArrayModel`` is the joist model; it replaces the former
``CollectorTribModel`` (removed), whose whole-spread trib area could not
represent support gaps, variable spans or joist containers.
"""

from __future__ import annotations
from typing import Optional
import warnings

from shapely import LineString, Point, GeometryCollection

from papermodels.datatypes.element import (
    Element,
    Intersection,
    NODE_ABS_TOL,
)
from papermodels.datatypes.geometry_model import GeometryModel
from papermodels.geometry import geom_ops
from papermodels.geometry import joist_array as ja


class JoistArrayModel:
    """
    Models a spread of joists from one drawn joist prototype.

    The array is laid out by the station engine in
    ``papermodels.geometry.joist_array`` (see its module docstring), in one of
    three modes:

    - **container** — the prototype lies in a ``joist_container`` polygon; each
      joist is the container cut along its station line, so cantilevers and
      backspans can both vary.
    - **extent** — the prototype has an extent line; joists run between the
      outer supports at each station with the prototype's constant cantilevers.
    - **plain** — neither; the array runs over the station range where all of
      the prototype's supports exist, with constant cantilevers.

    Joists that would land where supports are missing are moved to the nearest
    valid station (or dropped, with a warning); trib areas tile the array region.

    Supports are the parent element's graph supports (``intersections_below``),
    read from the shared ``GeometryModel`` (wall centerlines, beam lines).  Every
    joist/support crossing is a canonical node of the model's ``NodeRegistry``,
    and the generated joists are written back into the model.
    """

    use_subelements = True

    def __init__(
        self,
        element: Element,
        spacing: float = 1,
        initial_offset: float | int = 0.0,
        joist_at_start: bool = True,
        joist_at_end: bool = True,
        cantilever_tolerance: float = 1e-1,
        geometry_model: Optional[GeometryModel] = None,
        min_bearing: float = 1e-3,
        span_jump_tol: Optional[float] = None,
        suppress_warnings: bool = False,
    ):
        if element.geometry.geom_type != "LineString":
            raise geom_ops.GeometryError(
                f"Joist prototype {element.tag} must be a LineString."
            )
        self.element = element
        self.id = element.tag
        self.plane_id = element.plane_id
        self.elem_kwargs = element.kwargs
        self.spacing = spacing
        self.initial_offset = float(initial_offset)
        self.joist_at_start = joist_at_start
        self.joist_at_end = joist_at_end
        self._cantilever_tolerance = cantilever_tolerance
        self.min_bearing = min_bearing
        self.suppress_warnings = suppress_warnings
        self.geometry_model = geometry_model or _standalone_model(element)

        self.joist_prototype = element.geometry
        self.frame = ja.ArrayFrame.from_prototype(self.joist_prototype)
        self.supports = self._graph_supports()
        self.prototype_supports = self._supports_at_prototype()
        self.cantilevers = ja.prototype_cantilevers(
            self.frame,
            self.joist_prototype,
            self.prototype_supports,
            tolerance=cantilever_tolerance,
        )
        self.mode, region, s_range = self._region()
        self._warn_on_unlinked_supports(region)
        cant_a, cant_b = self.cantilevers
        self.result = ja.build_array(
            self.frame,
            region,
            self.supports,
            s_range,
            ja.MODE_CONTAINER if self.mode == "container" else ja.MODE_FIXED,
            spacing,
            initial_offset=self.initial_offset,
            joist_at_start=self.joist_at_start,
            joist_at_end=self.joist_at_end,
            cant_a=cant_a,
            cant_b=cant_b,
            min_bearing=min_bearing,
            span_jump_tol=span_jump_tol,
            registry=self.geometry_model.nodes,
        )
        self.events = self.result.events
        self._warn_events()

    # -- construction helpers ------------------------------------------------

    def _graph_supports(self) -> list[ja.Support]:
        model = self.geometry_model
        supports = []
        for ib in self.element.intersections_below or []:
            gid = ib.other_tag
            if gid not in model.geometries or any(sup.gid == gid for sup in supports):
                continue
            if model.is_bearing_support(gid):
                line = model.geometries[gid]
            elif gid in model.polygons:
                # A point-reaction polygon the graph connects to the prototype.
                # Drawings made before the 'Reaction Type' legend field draw
                # walls this way; like the previous JoistArrayModel, bear on it
                # along its long axis (a post parallel to the joists is dropped
                # below as non-bearing).
                line = geom_ops.get_wall_centerline(model.polygons[gid])
                self._warn(
                    f"Joist array {self.id}: support {gid} is a polygon with a point "
                    "reaction type; the joists bear on it along its long axis. If it "
                    "is a wall, add 'Reaction Type: Linear' to its legend entry."
                )
            else:
                continue
            supports.append(ja.Support(gid, line, model.source_geometry(gid)))
        return ja.usable_supports(self.frame, supports)

    def _supports_at_prototype(self) -> list[ja.Support]:
        """The supports the prototype itself bears on (present at its station)."""
        s_p = self.frame.ts(self.joist_prototype.interpolate(0.5, normalized=True))[1]
        present = [
            sup
            for sup in self.supports
            if (lambda r: r[0] - ja.EPS <= s_p <= r[1] + ja.EPS)(
                self.frame.s_range(sup.line)
            )
        ]
        if len(present) < 2:
            raise geom_ops.GeometryError(
                f"The joist prototype {self.id} bears on {len(present)} support(s) at "
                "its own location; it needs at least two. Check that it extends onto "
                "both of its supports in the source sketch."
            )
        return present

    def _region(self):
        frame = self.frame
        container = getattr(self.element, "joist_container", None)
        if container is not None:
            return "container", container, frame.s_range(container)
        cant_a, cant_b = self.cantilevers
        if self.element.extent_line is not None:
            s0, s1 = frame.s_range(self.element.extent_line)
            s_range = (min(s0, 0.0), max(s1, 0.0))
            t_range = self._extent_t_range(s_range, cant_a, cant_b)
            region = frame.band(s_range[0], s_range[1], *t_range)
            return "extent", region, s_range
        s_range = ja.common_s_range(frame, self.prototype_supports)
        if s_range is None:
            raise geom_ops.GeometryError(
                f"The supports of joist prototype {self.id} do not overlap along the array."
            )
        sup_a, sup_b = ja.outer_supports(
            frame, self.joist_prototype, self.prototype_supports
        )
        s_p = frame.ts(self.joist_prototype.interpolate(0.5, normalized=True))[1]
        s_range, at_start, at_end = ja.converging_s_range(
            frame, sup_a.line, sup_b.line, s_range, s_p
        )
        # No zero-backspan joist where the outer supports meet.
        self.joist_at_start = self.joist_at_start and not at_start
        self.joist_at_end = self.joist_at_end and not at_end
        region = ja.fixed_region(frame, sup_a.line, sup_b.line, s_range, cant_a, cant_b)
        return "plain", region, s_range

    def _extent_t_range(self, s_range, cant_a, cant_b) -> tuple[float, float]:
        """
        t-range of the extent region: the prototype plus every support crossing
        within the extent's station range, widened by the cantilevers so the
        joists (and their trib) are fully inside.
        """
        frame = self.frame
        band = frame.band(s_range[0], s_range[1], -_FAR_T, _FAR_T)
        ts = list(frame.t_range(self.joist_prototype))
        for sup in self.supports:
            piece = sup.line.intersection(band)
            if not piece.is_empty:
                ts.extend(frame.t_range(piece))
        return min(ts) - cant_a - ja.REGION_PAD, max(ts) + cant_b + ja.REGION_PAD

    def _warn_on_unlinked_supports(self, region) -> None:
        """
        Extent/container arrays: bearing supports crossing the region that the
        graph did not link to this prototype (their load path would be lost).
        Plain arrays use only the prototype's own supports by definition.
        """
        if self.mode == "plain":
            return
        linked = {sup.gid for sup in self.supports}
        model = self.geometry_model
        found = model.query_supports(region, self.plane_id, self.element.rank or 0)
        candidates = [
            ja.Support(gid, model.geometries[gid], model.source_geometry(gid))
            for gid in found
            if gid not in linked and gid not in model.generated.get(self.id, [])
        ]
        spans = ja.support_spans(
            self.frame, ja.usable_supports(self.frame, candidates), region
        )
        missing = [
            gid
            for gid, ranges in spans.items()
            if any(hi - lo > self.min_bearing for lo, hi in ranges)
        ]
        if missing:
            self._warn(
                f"Joist array {self.id}: supports {missing} lie inside the array region but "
                "are not connected to the prototype in the graph; joists will not bear on them."
            )

    def _warn(self, msg: str) -> None:
        if not self.suppress_warnings:
            warnings.warn(msg)

    def _warn_events(self) -> None:
        for ev in self.events:
            if ev.kind == "dropped":
                self._warn(
                    f"Joist array {self.id}: no joist could be placed near station "
                    f"{ev.target:.3f} ({ev.detail})."
                )
            elif ev.kind == "span_jump":
                self._warn(
                    f"Joist array {self.id}: abrupt change in span at station "
                    f"{ev.final:.3f} ({ev.detail}). Consider using a separate joist "
                    "prototype for each span."
                )

    # -- results -------------------------------------------------------------

    @property
    def joist_geoms(self) -> list[LineString]:
        return [j.geometry for j in self.result.joists]

    @property
    def joist_trib_areas(self) -> list:
        return [j.trib_area for j in self.result.joists]

    @property
    def joist_locations(self) -> list[float]:
        return [j.station for j in self.result.joists]

    def __call__(self) -> Element:
        """
        Returns a copy of the prototype element with one subelement per joist.
        """
        e = self.element
        model = self.geometry_model
        model.remove_generated(self.id)
        # The lines the engine bore on (wall centerlines, beam lines)
        support_lines = {sup.gid: sup.line for sup in self.supports}
        subelements = []
        for idx, joist in enumerate(self.result.joists):
            sub_id = f"{self.id}-{idx}"
            geometry = LineString(
                geom_ops.order_nodes_positive([Point(c) for c in joist.geometry.coords])
            )
            intersections_below = []
            for crossing in joist.crossings:
                support_line = support_lines[crossing.gid]
                source = model.source_geometry(crossing.gid)
                intersections_below.append(
                    Intersection(
                        intersecting_region=Point(crossing.xy),
                        other_geometry=source,
                        other_tag=crossing.gid,
                        other_overlap=_bearing_overlap(geometry, source),
                        other_reaction_type=model.reaction_types[crossing.gid],
                        other_extents=_trib_extent_on_support(
                            support_line, joist.trib_area, crossing.xy
                        ),
                    )
                )
            subelements.append(
                Element(
                    geometry,
                    sub_id,
                    rank=e.rank,
                    intersections_below=intersections_below,
                    intersections_above=[],
                    correspondents_below=[],
                    correspondents_above=[],
                    plane_id=self.plane_id,
                    element_type="collector",
                    subelements=None,
                    trib_area=joist.trib_area,
                    kwargs=self.elem_kwargs,
                )
            )
            model.add_geometry(
                sub_id,
                geometry,
                rank=e.rank if e.rank is not None else 0,
                plane=self.plane_id,
                crossings={c.node_id: c.gid for c in joist.crossings},
                parent=self.id,
            )
        return Element(
            e.geometry,
            tag=e.tag,
            rank=e.rank,
            intersections_above=e.intersections_above,
            intersections_below=e.intersections_below,
            correspondents_above=e.correspondents_above,
            correspondents_below=e.correspondents_below,
            plane_id=e.plane_id,
            element_type=e.element_type,
            subelements=subelements,
            trib_area=e.trib_area,
            reaction_type="linear",
            kwargs=e.kwargs,
            extent_line=e.extent_line,
        )

    def show_svg(self, use_ipython_display: bool = True):
        """
        Displays (or returns) a GeometryCollection of the joists, their trib
        areas and their supports, for manual visual review.
        """
        collection = GeometryCollection(
            self.joist_geoms
            + [t for t in self.joist_trib_areas if t is not None]
            + [sup.line for sup in self.supports]
        )
        if not use_ipython_display:
            return collection
        from IPython.display import display

        display(collection)


# Far t-extent used to build a band before clipping (real-world units).
_FAR_T = 1e4


def _bearing_overlap(joist: LineString, support) -> Optional[LineString]:
    """
    The joist's bearing length on a wall: its segment inside the wall polygon
    (reported as the support's overlap_length). None for line supports (beams),
    which a joist crosses at a point.
    """
    if support.geom_type != "Polygon":
        return None
    overlap = joist.intersection(support)
    if overlap.is_empty or overlap.length == 0.0:
        return None
    if overlap.geom_type == "MultiLineString":
        overlap = max(overlap.geoms, key=lambda g: g.length)
    return overlap


def _trib_extent_on_support(
    support_line: LineString, trib_area, crossing_xy
) -> tuple[float, float]:
    """
    (start, end) of the joist's trib band along the support, measured from the
    support's positive-x start node (the convention of Element.get_collector_extents).
    Falls back to the crossing itself when the band does not overlap the support.
    """
    start, _ = geom_ops.get_start_end_nodes(support_line)
    piece = support_line.intersection(trib_area) if trib_area is not None else None
    if piece is None or piece.is_empty or piece.length == 0.0:
        d = start.distance(Point(crossing_xy))
        return (d, d)
    ds = [start.distance(Point(c)) for c in ja._all_coords(piece)]
    return (min(ds), max(ds))


def _standalone_model(element: Element) -> GeometryModel:
    """
    A GeometryModel over just this element and its supports, for using the
    joist array outside of a GeometryGraph.
    """

    class _Proxy:
        def __init__(self, tag, geometry, rank, plane_id, reaction_type):
            self.tag = tag
            self.geometry = geometry
            self.rank = rank
            self.plane_id = plane_id
            self.reaction_type = reaction_type

    rank = element.rank or 0
    proxies = [
        _Proxy(element.tag, element.geometry, rank, element.plane_id, "point")
    ] + [
        _Proxy(
            ib.other_tag,
            ib.other_geometry,
            rank + 1,
            element.plane_id,
            ib.other_reaction_type,
        )
        for ib in element.intersections_below or []
    ]
    return GeometryModel.from_elements(proxies, node_abs_tol=NODE_ABS_TOL)
