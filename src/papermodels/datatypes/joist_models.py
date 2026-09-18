from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Any, Optional
import warnings
import numpy as np
from shapely import (
    LineString,
    Point,
    MultiLineString,
    Polygon,
    MultiPoint,
    convex_hull,
    GeometryCollection,
    box,
    set_precision,
)
import shapely.ops as ops

from papermodels.datatypes.element import (
    Element,
    Intersection,
    NODE_ABS_TOL,
    trim_cantilevers,
    align_frames_to_centroids,
)
from papermodels.datatypes.geometry_model import GeometryModel
from papermodels.geometry import geom_ops
from papermodels.geometry import joist_array as ja
import load_distribution as ld

from rich import print


@dataclass
class CollectorTribModel:
    element: Element
    trib_width: float = 1.0
    reaction_type: str = "linear"
    use_subelements: bool = False
    cantilever_tolerance: float = 1e-2

    def __post_init__(self):
        if self.element.extent_polygon is not None:
            self.use_subelements = True

    def __call__(self):
        """
        Generates a representative trib area for the joist prototype.

        Assumptions:
        - The supports are assumed to be orthogonal
        - IF EXTENT LINES ARE USED, then cantilevers are not supported
            (joist prototypes with extent lines using this model will
            have their cantilevers trimmed off...currently).
        - The loading is consistent for all joists within the spread
        - The joist represents a one-of-many similar elements within
            the spread and the spread thus represents a linear reaction
            over the supports.
        """
        e = self.element
        geom = e.geometry
        collector_extents = e.get_collector_extents(relative=False)
        left = []
        right = []
        for extent in collector_extents.values():
            p0_relation = geom_ops.relate_point_to_line(extent[0], geom)
            if p0_relation in (("left", "above"), ("right", "above")):
                left.append(extent[0])
                right.append(extent[1])
            elif p0_relation in (("left", "below"), ("right", "below")):
                left.append(extent[1])
                right.append(extent[0])
        left_dist = [geom.distance(extent) for extent in left]
        right_dist = [geom.distance(extent) for extent in right]
        left_minimum_idx = left_dist.index(min(left_dist))
        right_minimum_idx = right_dist.index(min(right_dist))

        joist_vector = np.abs(geom_ops.get_direction_vector(geom))
        # This is one of the places where orthogonality is assumed
        joist_orientation = None
        if joist_vector[0] > joist_vector[1]:
            joist_orientation = "horizontal"
        elif joist_vector[1] > joist_vector[0]:
            joist_orientation = "vertical"
        else:
            print(f"JOIST ORIENTATION VERIANT: {geom=}")

        if joist_orientation == "vertical":
            minx = left[left_minimum_idx].coords[0][0]
            miny = geom.coords[0][1]
            maxx = right[right_minimum_idx].coords[0][0]
            maxy = geom.coords[1][1]
        elif joist_orientation == "horizontal":
            minx = geom.coords[0][0]
            miny = left[left_minimum_idx].coords[0][1]
            maxx = geom.coords[1][0]
            maxy = right[right_minimum_idx].coords[0][1]
        trib_area = box(minx, miny, maxx, maxy)
        # trib_area = e.geometry.buffer(self.trib_width/2)
        if not self.use_subelements:
            collector_element = Element(
                e.geometry,
                e.tag,
                0,
                e.intersections_above,
                e.intersections_below,
                e.correspondents_above,
                e.correspondents_below,
                e.plane_id,
                e.element_type,
                e.subelements,
                trib_area=trib_area,
                reaction_type="linear",
                kwargs=e.kwargs,
                extent_line=e.extent_line,
            )
        else:
            ext_poly = e.extent_polygon
            joist_prototype = e.geometry
            support_lines = {}
            for ib in e.intersections_below:
                if ib.other_geometry.geom_type == "Polygon":
                    support_line = geom_ops.get_rectangle_centerline(ib.other_geometry)
                    support_lines.update({support_line: ib.other_tag})
                elif ib.other_geometry.geom_type == "LineString":
                    support_lines.update({ib.other_geometry: ib.other_tag})
            support_geoms = [ib.other_geometry for ib in e.intersections_below]

            # 1. Find polygon extent edges that intersect with joist prototype
            # These will be our extent boundaries for the length of the prototype
            poly_edge_points = list(
                zip(ext_poly.exterior.coords, ext_poly.exterior.coords[1:])
            )
            start_edge = None
            end_edge = None
            start_point, end_point = geom_ops.order_nodes_positive(
                [Point(coord) for coord in geom.coords]
            )
            for edge_points in poly_edge_points:
                edge = LineString(edge_points)
                if edge.intersects(geom) and (
                    start_point.distance(edge) < end_point.distance(edge)
                ):
                    start_edge = edge
                elif edge.intersects(geom) and (
                    end_point.distance(edge) < start_point.distance(edge)
                ):
                    end_edge = edge
            # Either the start or end edge should work since
            # orthogonality is assumed.

            # 2. Find support geoms which intersect with start and end edges
            start_supports = []
            end_supports = []
            joist_vector = np.abs(geom_ops.get_direction_vector(geom))
            # This is one of the places where orthogonality is assumed
            joist_orientation = None
            if joist_vector[0] > joist_vector[1]:
                joist_orientation = "horizontal"
            elif joist_vector[1] > joist_vector[0]:
                joist_orientation = "vertical"
            else:
                print(f"JOIST ORIENTATION VERIANT: {geom=}")
            support_centroids = [
                support_geom.centroid for support_geom in support_geoms
            ]
            if joist_orientation == "horizontal":
                start_support = min(support_centroids, key=lambda x: x.coords[0][0])
                end_support = max(support_centroids, key=lambda x: x.coords[0][0])
                start_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        start_support.coords[0][0],
                        support_line.coords[0][0],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]
                end_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        end_support.coords[0][0],
                        support_line.coords[0][0],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]
            elif joist_orientation == "vertical":
                start_support = min(support_centroids, key=lambda x: x.coords[0][1])
                end_support = max(support_centroids, key=lambda x: x.coords[0][1])
                start_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        start_support.coords[0][1],
                        support_line.coords[0][1],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]
                end_supports = [
                    support_line
                    for support_line in support_lines
                    if math.isclose(
                        end_support.coords[0][1],
                        support_line.coords[0][1],
                        abs_tol=self.cantilever_tolerance,
                    )
                ]

            # 2b. Get intermediate supports
            intermediate_support_lines = []
            for support_line in support_lines:
                if support_line not in start_supports + end_supports:
                    intermediate_support_lines.append(support_line)

            # 3. Generate overlap regions

            overlap_polys = []
            for start_support in start_supports:
                for end_support in end_supports:
                    overlap_poly = None
                    pa0, pa1 = start_support.coords
                    pb0, pb1 = end_support.coords
                    if joist_orientation == "vertical":
                        overlap_region = ld.get_overlap_coords(
                            pa0[0], pa1[0], pb0[0], pb1[0]
                        )
                        if overlap_region is not None:
                            overlap_poly = box(
                                overlap_region[0], pa0[1], overlap_region[1], pb1[1]
                            )
                    elif joist_orientation == "horizontal":
                        a_sort = sorted([pa0, pa1], key=lambda x: x[1])
                        b_sort = sorted([pb0, pb1], key=lambda x: x[1])
                        pa0 = a_sort[0]
                        pa1 = a_sort[1]
                        pb0 = b_sort[0]
                        pb1 = b_sort[1]
                        overlap_region = ld.get_overlap_coords(
                            pa0[1], pa1[1], pb0[1], pb1[1]
                        )
                        if overlap_region is not None:
                            overlap_poly = box(
                                pa0[0], overlap_region[0], pb1[0], overlap_region[1]
                            )
                    if overlap_poly is not None:
                        overlap_within_extent = ext_poly.intersection(overlap_poly)
                        overlap_polys.append(overlap_within_extent)

            # 5. Do overlap polys intersect with intermediate supports?
            #    if so, break the support as required.
            revised_poly_overlaps = []
            for overlap_poly in set(overlap_polys):
                split_polys = []
                for intermediate_support in intermediate_support_lines:
                    if intermediate_support.intersects(overlap_poly):
                        inter_coords = intermediate_support.coords
                        poly_splits = geom_ops.split_polygon(
                            overlap_poly, joist_orientation, inter_coords
                        )
                        split_polys += poly_splits
                if not split_polys:
                    revised_poly_overlaps.append(overlap_poly)
                else:
                    revised_poly_overlaps += split_polys
            sorted_poly_overlaps = sorted(
                revised_poly_overlaps,
                key=lambda x: (x.centroid.coords[0][0], x.centroid.coords[0][1]),
            )

            joist_prototype_geometries = []
            for overlap_poly in sorted_poly_overlaps:
                overlap_poly: Polygon
                overlap_edge_points = list(
                    zip(overlap_poly.exterior.coords, overlap_poly.exterior.coords[1:])
                )
                for pi, pj in overlap_edge_points:
                    edge_ls = LineString([pi, pj])
                    if geom_ops.check_2d_linestring_parallel(
                        edge_ls, start_edge, tol=0.01
                    ):
                        # Need to translate the original joist prototype to the new position
                        new_joist = geom_ops.translate_joist_to_point(
                            joist_prototype,
                            joist_orientation,
                            intersection_point=edge_ls.centroid,
                        )
                        trimmed_joist = new_joist.intersection(overlap_poly)
                        # We only need to hit one edge of the overlap so we can break here
                        break
                joist_prototype_geometries.append(trimmed_joist)

            # 7. Create an Element for each new joist prototype geometries
            subelements = []
            sorted_joist_geoms = joist_prototype_geometries
            for idx, joist_geom in enumerate(sorted_joist_geoms):
                intersections = []
                total_new_subs = len(joist_prototype_geometries)
                z_fill_qty = math.floor(math.log10(total_new_subs))
                index = f"{idx}".zfill(z_fill_qty)
                subelement_tag = f"{e.tag}-{index}"
                trib_area = sorted_poly_overlaps[idx]
                assert joist_geom.intersects(trib_area)
                for support_geom in support_geoms:
                    support_overlap = None
                    support_intersection = joist_geom.intersection(
                        support_geom, grid_size=1e-3
                    )
                    if not support_intersection:
                        continue
                    if support_geom.geom_type == "Polygon":
                        support_overlap = support_intersection
                        support_line = geom_ops.get_rectangle_centerline(
                            support_geom
                        )  # geom_ops.clean_polygon_supports([support_geom], joist_geom)[0]
                        intersecting_region = support_geom.exterior.intersection(
                            joist_geom, grid_size=1e-3
                        )
                    elif support_geom.geom_type == "LineString":
                        support_line = support_geom
                        intersecting_region = support_intersection

                    tag = support_lines[support_line]
                    # HERE: Previous behaviour was to intersect with the centerline but that is no longer a requirement
                    # Joist geom needs to be rebuilt to ensure it hits the wall centerline

                    # if intersecting_region.is_empty:
                    #     continue
                    intersection = Intersection(
                        intersecting_region=intersecting_region,
                        other_geometry=support_geom,
                        other_tag=tag,
                        other_overlap=support_overlap,
                        other_reaction_type=(
                            "linear" if support_geom.geom_type == "Polygon" else "point"
                        ),
                    )
                    intersections.append(intersection)
                subelement = Element(
                    geometry=joist_geom,
                    tag=subelement_tag,
                    rank=e.rank,
                    intersections_above=[],
                    intersections_below=intersections,
                    correspondents_above=[],
                    correspondents_below=[],
                    plane_id=e.plane_id,
                    element_type=e.element_type,
                    subelements=None,
                    trib_area=trib_area,
                    reaction_type="linear",
                    kwargs=e.kwargs,
                    # extent_polygon=e.extent_polygon,
                )
                aligned_subelement = align_frames_to_centroids(subelement)
                subelements.append(aligned_subelement)

            # 8. Return subelements
            collector_element = Element(
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
        return collector_element


def collector_trib_model(
    element: Element, trib_width: float, reaction_type: str = "linear"
) -> Element:
    """
    An alias for CollectorTribModel.__call__() for temporary
    backwards compatibility.
    """
    model = CollectorTribModel(element, trib_width, reaction_type)
    return model()


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
            if gid not in model.geometries or not model.is_bearing_support(gid):
                continue
            if any(sup.gid == gid for sup in supports):
                continue
            supports.append(ja.Support(gid, model.geometries[gid]))
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
            ja.Support(gid, model.geometries[gid])
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
        subelements = []
        for idx, joist in enumerate(self.result.joists):
            sub_id = f"{self.id}-{idx}"
            geometry = LineString(
                geom_ops.order_nodes_positive([Point(c) for c in joist.geometry.coords])
            )
            intersections_below = []
            for crossing in joist.crossings:
                support_line = model.geometries[crossing.gid]
                intersections_below.append(
                    Intersection(
                        intersecting_region=Point(crossing.xy),
                        other_geometry=model.source_geometry(crossing.gid),
                        other_tag=crossing.gid,
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
