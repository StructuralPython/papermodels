from __future__ import annotations
from typing import Optional, TypeAlias, Callable
from copy import deepcopy
from decimal import Decimal
import pathlib
import networkx as nx
import hashlib
import json

from papermodels.datatypes.element import (
    Element,
    LoadedElement,
    trim_cantilevers,
    align_frames_to_centroids,
)
from shapely import Point, LineString, Polygon
from ..geometry import geom_ops as geom
from ..datatypes.element import (
    Correspondent,
    Intersection,
)
from ..paper.annotations import (
    Annotation,
    scale_annotations,
    parse_annotations,
    parsed_annotations_to_loading_geometry,
    filter_annotations,
    tag_parsed_annotations,
    assign_page_id_to_annotations,
    annotation_to_shapely,
)
from ..paper.plot import plot_annotations, plot_elements
from ..paper import pdf
from ..paper import dxf
from ..datatypes.exceptions import AnnotationError
from rich.progress import track
from rich import print
import numpy.typing as npt

Rule: TypeAlias = Callable


def TRANSFER_LINES_CANNOT_INTERSECT_WITH_LINEAR_POLYGONS(
    e: Element, inter: Intersection
):
    return not (
        (
            e.rank > 0
            and e.geometry.geom_type == "LineString"
            and inter.other_geometry.geom_type == "Polygon"
            and inter.other_reaction_type == "linear"
        )
        or (
            e.rank > 0
            and e.geometry.geom_type == "Polygon"
            and e.reaction_type == "linear"
            and inter.other_geometry.geom_type == "LineString"
        )
    )


class GeometryGraph(nx.DiGraph):
    """
    A class to represent a connected geometry system in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

    Can be used to generate a GeometryGraph.

    The node_hash is how changes to the graph nodes can be tracked.
    """

    def __init__(
        self,
        do_not_process: bool = False,
        cantilever_abs_tol: Optional[float] = 0.2,
    ):
        super().__init__()
        self.do_not_process = do_not_process
        self.node_hash = None
        self.loading_geometries = None
        self.parsed_annotations = None
        self.raw_annotations = None
        self.legend_entries = None
        self.pdf_path = None
        self.cantilever_abs_tol: Optional[float] = cantilever_abs_tol

    @property
    def collector_elements(self):
        return [
            node_name
            for node_name in self.nodes
            if not list(self.predecessors(node_name))
            if self.nodes[node_name]["element"].rank == 0
        ]

    @property
    def transfer_elements(self):
        return [
            node_name for node_name in self.nodes if list(self.predecessors(node_name))
        ]

    @property
    def orphaned_elements(self):
        return [
            node_name
            for node_name in self.nodes
            if (
                (len(list(self.successors(node_name))) < 2)
                and self.nodes[node_name]["element"].geometry.geom_type == "LineString"
            )
            or (
                (len(list(self.successors(node_name))) < 1)
                and self.nodes[node_name]["element"].geometry.geom_type == "Polygon"
            )
        ]

    @property
    def contiguous_elements(self):
        orphans = self.orphaned_elements
        return [node_name for node_name in self.nodes if node_name not in orphans]

    @classmethod
    def from_elements(
        cls,
        elements: list[Element],
        do_not_process: bool = False,
        cantilever_abs_tol: Optional[float] = 0.2,
        intersection_rules: Optional[list[Rule | callable]] = [
            TRANSFER_LINES_CANNOT_INTERSECT_WITH_LINEAR_POLYGONS
        ],
    ) -> GeometryGraph:
        """
        Returns a GeometryGraph (networkx.DiGraph) based upon the intersections and correspondents
        of the 'elements'.
        """
        g = cls(cantilever_abs_tol=cantilever_abs_tol)
        elements_copy = deepcopy(elements)
        if intersection_rules is None:
            intersection_rules = []
        for element in elements_copy:
            hash = hashlib.sha256(str(element).encode()).hexdigest()
            start_coord = None
            if element.geometry.geom_type == "LineString":
                coords_a, coords_b = element.geometry.coords
                ordered_coords = geom.order_nodes_positive(
                    [Point(coords_a), Point(coords_b)]
                )
                start_coord = ordered_coords[0]
            g.add_node(
                element.tag,
                element=element,
                sha256=hash,
                start_coord=start_coord,
                # b_coord=ordered_coords[1]
            )
            if element.correspondents_below is not None:
                for correspondent in element.correspondents_below:
                    j_tag = correspondent.other_tag
                    g.add_edge(element.tag, j_tag, edge_type="correspondent")
            if element.intersections_below is not None:
                filtered_intersections = []
                for intersection in element.intersections_below:
                    passes_intersection_rules = []
                    for intersection_rule in intersection_rules:
                        passes_intersection_rules.append(
                            intersection_rule(element, intersection)
                        )
                    if all(passes_intersection_rules):
                        j_tag = intersection.other_tag
                        g.add_edge(element.tag, j_tag, edge_type="intersection")
                        filtered_intersections.append(intersection)
                element.intersections_below = filtered_intersections
            if element.tag in g.collector_elements:
                for correspondent in element.correspondents_above:
                    j_tag = correspondent.other_tag
                    g.add_edge(j_tag, element.tag, edge_type="correspondent")

        for node in g.collector_elements:
            g.nodes[node]["element"].element_type = "collector"

        for node in g.transfer_elements:
            g.nodes[node]["element"].element_type = "transfer"

        if do_not_process:
            return g

        g.align_frames_to_centroids()
        g.trim_cantilevers()
        g.remove_excess_correspondent_load_paths()
        g.add_intersection_indexes_below()
        g.add_intersection_indexes_above()

        return g

    def align_frames_to_centroids(self):
        """
        Aligns the ends of frame elements so that they start and end on the centroids
        of posts and walls (centerlines).
        """
        # Only execute on transfer elements because collector elements will be modified
        # when collector behaviour is assigned to them.
        transfer_nodes = self.transfer_elements
        contiguous_nodes = self.contiguous_elements
        nodes_to_align = set(transfer_nodes) & set(contiguous_nodes)
        sorted_nodes = nx.topological_sort(self)

        for node_name in sorted_nodes:
            if node_name not in nodes_to_align:
                continue
            node = self.nodes[node_name]
            element = node["element"]
            new_element = align_frames_to_centroids(element)
            node["element"] = new_element

    def trim_cantilevers(self):
        """
        Trims cantilevers if they are within the tolerance
        """
        # Only execute on transfer elements because collector elements will be modified
        # when collector behaviour is assigned to them.
        transfer_nodes = self.transfer_elements
        contiguous_nodes = self.contiguous_elements

        for node_name in set(transfer_nodes) & set(contiguous_nodes):
            node = self.nodes[node_name]
            element = node["element"]
            new_element = trim_cantilevers(element, abs_tol=self.cantilever_abs_tol)
            node["element"] = new_element

    def remove_excess_correspondent_load_paths(self):
        """
        Removes edges from the graph for the following conditions:

        1. A Polygon node that has more than one "correspondent" edge. The "correspondent"
            edge leading to a node with rank 0 is prioritized. If no node with rank 0
            is present, the edge leading to the node with the larger overlap ratio is
            prioritized. This prevents a correspondent load above from transferring
            to both its "end load" and to a correspondent below at the same time. The load
            path from above should terminate at the "end load" and the end load should transfer
            to whatever it intersects with.
        2. A Polygon node, with a "point" reaction type, that has an "intersection" edge
            and one or more "correspondent" edges. If an "intersection" edge is present,
            then the "correspondent" edges will be removed. This represents the condition
            of platform-framing where a post will land on the floor framing, and transfer
            through it to the supporting post below and the floor framing should receive
            a "crushing" load from the posts above and below without the posts directly
            transferring to each other.
        3. A LineString node that has intersections_below with other LineStrings and the
            intersecting_region (Point) is within a Polygon that is also an intersection_below.
            This represents the condition where a frame element connects to another frame element,
            fully transfering to the other frame element, at the same location where the
            supporting frame element is transferring to a column. This rule is intended to prevent
            both frame elements transferring their load to a column in addition to the supported
            frame element transferring load to the supporting frame element. The correct load
            path should be |FB0.1 -> FB0.2 -> column| instead of |FB0.1 -> column| with
            |FB0.1 -> FB0.2 -> column| also.

        Modifications to the implementation of this function can adjust how load paths are
        conceptually created. For example, to implement baloon framing, the second rule
        can be omitted.
        """
        sorted_nodes = nx.topological_sort(self)
        for node in sorted_nodes:
            element = self.nodes[node]["element"]
            dependents = list(self.successors(node))
            dependent_edges = [(node, dep) for dep in dependents]
            edge_properties = [
                self.edges[edge]["edge_type"] for edge in dependent_edges
            ]
            intersection_below_points = [
                inter.intersecting_region
                for inter in element.intersections_below
                if inter.other_geometry.geom_type == "LineString"
            ]
            intersection_below_polygons = [
                (inter.other_tag, inter.other_geometry)
                for inter in element.intersections_below
                if inter.other_geometry.geom_type == "Polygon"
            ]
            intersection_points_in_polygon_below = hits = []
            for pt in intersection_below_points:
                for inter_id, inter_poly in intersection_below_polygons:
                    if inter_poly.contains(pt):
                        hits.append((node, inter_id))
            # Rule 1
            if (
                element.geometry.geom_type == "Polygon"
                and edge_properties.count("correspondent") > 1
                and element.reaction_type == "point"
            ):
                dep_to_keep = None
                max_overlap = 0.0
                for idx, dep in enumerate(dependents):
                    # Keep the rank 0
                    if self.nodes[dep]["element"].rank == 0:
                        dep_to_keep = idx
                        break
                    else:
                        # Or find the correspondent with the largest overlap ratio
                        try:
                            dep_overlap_ratio = next(
                                (
                                    corr.overlap_ratio
                                    for corr in element.correspondents_below
                                    if corr.other_tag == dep
                                )
                            )
                        except StopIteration:
                            raise ValueError(
                                f"Number of dependents does not match number of correspondents_below.\n"
                                "This can happen if a polygon element is corresponding with more than one "
                                "polygons on the page below and the ranks of the polygons below are not "
                                "quite right. If you are intending to transfer out this polygon to a frame member "
                                "then check to make sure that the transfer element below has a rank of 0.\n"
                                f"{node=}\n{dependents=}\n{element.correspondents_below=}"
                            )
                        if dep_overlap_ratio > max_overlap:
                            dep_to_keep = idx
                            max_overlap = dep_overlap_ratio

                # Remove the edges
                if dep_to_keep is not None:
                    for idx, edge in enumerate(dependent_edges):
                        if idx != dep_to_keep:
                            self.remove_edge(*edge)

            # Rule 2
            if (
                element.geometry.geom_type == "Polygon"
                and element.reaction_type == "point"
                and "intersection" in edge_properties
                and element.rank == 0
            ):
                dep_to_keep = None
                secondary_dep_to_keep = None
                # We only want to keep one intersection
                if edge_properties.count("intersection") > 1:
                    # Go through each intersection
                    dep_overlap_length = 0.0
                    dep_overlap_area = 0.0
                    for idx, dep in enumerate(dependents):
                        if edge_properties[idx] == "intersection":
                            dependent_geometry = self.nodes[dep]["element"].geometry
                            element_geometry = element.geometry
                            if dependent_geometry.geom_type == "LineString":
                                overlap_length = element_geometry.intersection(
                                    dependent_geometry
                                ).length
                                if overlap_length > dep_overlap_length:
                                    dep_to_keep = idx
                                    dep_overlap_length = overlap_length
                            elif dependent_geometry.geom_type == "Polygon":
                                dep_element = self.nodes[dep]["element"]
                                if (
                                    dep_element.reaction_type == "point"
                                ):  # Points should transfer to points
                                    overlap_area = element_geometry.intersection(
                                        dep_element.geometry
                                    ).area
                                    if overlap_area > dep_overlap_area:
                                        secondary_dep_to_keep = idx
                                        dep_overlap_area = overlap_area
                else:
                    dep_to_keep = edge_properties.index("intersection")

                if dep_to_keep is None and secondary_dep_to_keep is not None:
                    dep_to_keep = secondary_dep_to_keep

                for idx, edge in enumerate(dependent_edges):
                    if idx != dep_to_keep:
                        self.remove_edge(*edge)

            # Rule 3
            if (
                element.geometry.geom_type == "LineString"
                and intersection_points_in_polygon_below
            ):
                for edge in intersection_points_in_polygon_below:
                    self.remove_edge(*edge)

    def add_intersection_indexes_below(self):
        sorted_nodes = nx.topological_sort(self)
        orphaned_nodes = self.orphaned_elements
        for node in sorted_nodes:
            if node in orphaned_nodes:
                continue
            node_attrs = self.nodes[node]
            element: Element = node_attrs["element"]
            dependents = list(self.successors(node))
            dependent_intersections = get_dependent_intersections(element, dependents)
            dependent_correspondents = get_dependent_correspondents(element, dependents)
            if not dependent_intersections and not dependent_correspondents:
                continue
            if element.geometry.geom_type == "Polygon":  # node geometry is polygon
                updated_intersections_below = []
                all_extents = {}
                if element.reaction_type == "linear":
                    all_extents = element.get_transfer_extents()
                for intersection in dependent_intersections:
                    extents = all_extents.get(intersection.other_tag)
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        self.nodes[intersection.other_tag]["element"].geometry,
                        intersection.other_tag,
                        intersection.other_overlap,
                        0,
                        intersection.other_reaction_type,
                        other_extents=extents,
                    )
                    updated_intersections_below.append(new_intersection)

                updated_correspondents_below = []
                for correspondent in dependent_correspondents:
                    extents = all_extents.get(correspondent.other_tag)
                    new_correspondent = Correspondent(
                        correspondent.overlap_ratio,
                        correspondent.other_geometry,
                        correspondent.other_tag,
                        correspondent.other_rank,
                        correspondent.other_reaction_type,
                        extents,
                    )
                    updated_correspondents_below.append(new_correspondent)
                    element.correspondents_below = updated_correspondents_below
            else:  # For LineStrings
                intersection_below_local_coords = get_local_coords(
                    node_attrs["start_coord"], dependent_intersections
                )
                sorted_below_ints = sorted(
                    intersection_below_local_coords, key=lambda x: x[0]
                )
                _, other_tags_below = zip(*sorted_below_ints)
                updated_intersections_below = []
                extents = {}
                if node in self.collector_elements:
                    if element.subelements is not None:
                        for subelem in element.subelements:
                            start_coord, _ = geom.order_nodes_positive(
                                subelem.geometry.boundary.geoms
                            )
                            sub_dependent_intersections = get_dependent_intersections(
                                subelem, dependents
                            )
                            sub_local_coords = get_local_coords(
                                start_coord, sub_dependent_intersections
                            )
                            sub_id = subelem.tag
                            subextents = subelem.get_collector_extents()
                            sub_sorted_below_ints = sorted(
                                sub_local_coords, key=lambda x: x[0]
                            )
                            if len(sub_sorted_below_ints) < 2:
                                raise ValueError(
                                    f"It seems that this element only has one support: {node}"
                                )
                            _, sub_other_tags_below = zip(*sub_sorted_below_ints)
                            sub_updated_intersections_below = []
                            for sub_intersection in subelem.intersections_below:
                                sub_other_tag = sub_intersection.other_tag
                                sub_local_index = sub_other_tags_below.index(
                                    sub_other_tag
                                )
                                new_sub_intersection = Intersection(
                                    sub_intersection.intersecting_region,
                                    self.nodes[sub_other_tag]["element"].geometry,
                                    sub_intersection.other_tag,
                                    sub_intersection.other_overlap,
                                    sub_local_index,
                                    other_reaction_type=element.reaction_type,
                                    other_extents=subextents[
                                        sub_intersection.other_tag
                                    ],
                                )
                                sub_updated_intersections_below.append(
                                    new_sub_intersection
                                )
                            subelem.intersections_below = (
                                sub_updated_intersections_below
                            )
                        updated_intersections_below = element.intersections_below
                    else:
                        if element.reaction_type == "linear":
                            extents = element.get_collector_extents()
                        for intersection in dependent_intersections:
                            other_tag = intersection.other_tag
                            local_index = other_tags_below.index(other_tag)
                            new_intersection = Intersection(
                                intersection.intersecting_region,
                                self.nodes[intersection.other_tag]["element"].geometry,
                                intersection.other_tag,
                                intersection.other_overlap,
                                local_index,
                                other_reaction_type=intersection.other_reaction_type,
                                other_extents=extents.get(other_tag, None),
                            )
                            updated_intersections_below.append(new_intersection)

                else:
                    if len(sorted_below_ints) < 2:
                        raise ValueError(
                            f"It seems that this element only has one support: {node}"
                        )
                    for intersection in dependent_intersections:
                        other_tag = intersection.other_tag
                        local_index = other_tags_below.index(other_tag)
                        new_intersection = Intersection(
                            intersection.intersecting_region,
                            self.nodes[intersection.other_tag]["element"].geometry,
                            intersection.other_tag,
                            intersection.other_overlap,
                            local_index,
                            other_reaction_type=intersection.other_reaction_type,
                            other_extents=extents.get(other_tag, None),
                        )
                        updated_intersections_below.append(new_intersection)
                element.intersections_below = updated_intersections_below
            element.intersections_below = updated_intersections_below
            self.nodes[node]["element"] = element

    def add_intersection_indexes_above(self):
        sorted_nodes = nx.topological_sort(self)
        transfer_elements = self.transfer_elements
        for node in transfer_elements:
            element = self.nodes[node]["element"]
            element_tag = element.tag
            predecessors = list(self.predecessors(node))
            predecessor_intersections = [
                intersection
                for intersection in element.intersections_above
                if intersection.other_tag in predecessors
            ]
            predecessor_correspondents = [
                correspondent
                for correspondent in element.correspondents_above
                if correspondent.other_tag in predecessors
            ]
            if not predecessor_intersections and not predecessor_correspondents:
                continue
            indexed_intersections_above = []
            for intersection in predecessor_intersections:
                other_tag = intersection.other_tag
                element_above: Element = self.nodes[other_tag]["element"]
                above_dependents = list(self.successors(other_tag))
                element_above_dependent_intersections = [
                    intersection
                    for intersection in element_above.intersections_below
                    if node in above_dependents
                ]
                above_intersections_below = {
                    above_intersection_below.other_tag: (
                        above_intersection_below.other_index,
                        above_intersection_below.other_extents,
                    )
                    for above_intersection_below in element_above_dependent_intersections
                }
                if element_above.subelements is None:
                    local_index = above_intersections_below[element_tag][0]
                    other_extents = above_intersections_below[element_tag][1]
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        element.geometry,
                        intersection.other_tag,
                        intersection.other_overlap,
                        local_index,
                        element_above.reaction_type,
                        other_extents=other_extents,
                    )
                    indexed_intersections_above.append(new_intersection)
                else:
                    for subelem_above in element_above.subelements:
                        # 1. Find subelements above that actually intersect with this (below) element
                        intersections_this_element = []
                        for above_inter_below in subelem_above.intersections_below:
                            if above_inter_below.other_tag == element_tag:
                                other_extents = above_inter_below.other_extents
                                intersections_this_element.append(
                                    Intersection(
                                        above_inter_below.intersecting_region,
                                        element.geometry,
                                        subelem_above.tag,
                                        above_inter_below.other_overlap,
                                        above_inter_below.other_index,
                                        subelem_above.reaction_type,
                                        other_extents=other_extents,
                                    )
                                )
                        indexed_intersections_above += intersections_this_element

            element.intersections_above = indexed_intersections_above

            indexed_correspondents_above = []
            for correspondent in predecessor_correspondents:
                other_tag = correspondent.other_tag
                element_above: Element = self.nodes[other_tag]["element"]
                above_dependents = list(self.successors(other_tag))
                element_above_dependent_correspondents = [
                    correspondent
                    for correspondent in element_above.correspondents_below
                    if element.tag in above_dependents
                ]
                above_correspondents_below = {
                    above_correspondent_below.other_tag: above_correspondent_below.other_extents
                    for above_correspondent_below in element_above_dependent_correspondents
                }
                other_extents = above_correspondents_below[element_tag]
                if element_above.subelements is None:
                    new_correspondent = Correspondent(
                        correspondent.overlap_ratio,
                        correspondent.other_geometry,
                        correspondent.other_tag,
                        correspondent.other_rank,
                        correspondent.other_reaction_type,
                        other_extents,
                    )
                    indexed_correspondents_above.append(new_correspondent)
            element.correspondents_above = indexed_correspondents_above

            self.nodes[node]["element"] = element

    def assign_collector_behaviour(
        self,
        element_constructor: callable,
        filter_function: Optional[callable] = None,
        *args,
        **kwargs,
    ) -> list[Element]:
        """
        Returns a list of Element to be assigned to element.subelements for elements
        that have element_type == "collector".

        'element_constructor': A callable class that will create a new Element populated with
            a trib_area attribute and possibly other modified attributes. The new Element
            will replace the existing element at that node.

            The class must have a .use_subelements attribute that must be set at time of
            initialization so that this function can query it and determine whether the
            resulting elements should be created as subelements on this element.

            The signature of the __call__ method of the class should be one of either:

            def __call__(element: Element, [*args, **kwargs]) -> Element
            -or
            def __call__(element: Element, [*args, **kwargs]) -> list[Element]

            A constructor function that returns a list[Element] should be used with
            as_subelements = True so that the returned list of Element will be assigned
            to the .subelements attribute of the Element

        'filter_function': A callable with the following function signature:

            def filter_function(element: Element) -> bool

            All elements that return True from the filter_function will have the
            element_constructor function called on them. Any functions that return
            False will have no collector behaviour assigned.

            You can use the papermodels.datatypes.element.create_element_filter
            to readily create such a filter function.

            If None, then all elements will be assigned the collector behaviour.

        '*args' and '**kwargs': These are passed through to the 'element_constructor'
            function
        """
        collectors = self.collector_elements

        for node in collectors:
            node_attrs = self.nodes[node]
            node_element: Element = node_attrs["element"]
            if filter_function is not None:
                # try:
                filter_passes = filter_function(node_element)
                # except:
                #     raise Exception("There was an exception generated during element filtering.")
            else:
                filter_passes = True

            if filter_passes:
                # If an incorrect geometry type makes its way into the element_constructor
                # e.g. a polygon is being entered as a joist in a joist-based element_constructor
                # then the element_constructor should return None
                # This prevents an error from being thrown if, for example, unconnected elements
                # are drawn. Unconnected elements have no precedents therefore they are (currently)
                # being categorized as collectors. However, I think incompatible geometries
                # should simply be ignored and not included as part of the processing.
                callable_instance = element_constructor(
                    node_element,
                    cantilever_tolerance=self.cantilever_abs_tol,
                    *args,
                    **kwargs,
                )
                new_elem = callable_instance()
                node_attrs["element"] = new_elem
        self.add_intersection_indexes_below()
        self.add_intersection_indexes_above()

    @classmethod
    def from_dxf_file(
        cls,
        dxf_filepath: pathlib.Path | str,
        legend_table: dict | pathlib.Path | str,
        page_idx: int = 0,
        scale: Decimal = Decimal(1.0),
        polygonize_layers: Optional[list[str]] = None,
        debug: bool = False,
        progress: bool = False,
        do_not_process: bool = False,
        show_skipped: bool = False,
    ):
        """
        Returns a GeometryGraph built from the geometric entities (LINE, LWPOLYLINE, INSERT)
        contained within 'dxf_filepath'.


        'legend_table': A dict (or a path to a JSON file) that maps layer names to
            annotation text properties
        'scale': An optional scale to be applied to the annotations. If not provided,
            the units of the annotations will be in PDF points where 1 point == 1 /72 inch
        'debug':  When True, will provide verbose documentation of the annotation parsing
            process to assist in reviewing errors and geometry inconsistencies.
        'progress': When True, a progress bar will be displayed
        'do_not_process': Reads the file and adds annotations to the graph but does not
            process the connectivity. Useful for debugging and plotting prior to processing.
        'show_skipped': Shows the skipped annotations that occured during pdf.load_pdf_annotations
        """
        if isinstance(legend_table, (str, pathlib.Path)):
            legend_table_path = pathlib.Path(legend_table)
            with open(legend_table_path, "r") as file:
                legend_table = json.load(file)
        annotations = dxf.load_dxf_annotations(
            dxf_filepath, page_idx, polygonize_layers=polygonize_layers
        )
        scaled_annotations = scale_annotations(
            annotations, scale=scale, paper_origin=(0, 0)
        )

        # parsed_annotations = parse_annotations(
        #     scaled_annots_in_page, legend_entries, legend_identifier
        # )

        load_entries = {}
        trib_area_entries = {}
        structural_element_entries = {}
        parsed_annotations = {}
        raw_annotations = {}
        for idx, scaled_annot in enumerate(scaled_annotations):
            raw_annot = annotations[idx]
            layer_name = scaled_annot.text
            annot_attrs = legend_table.get(layer_name)
            if annot_attrs is None:
                continue
            annot_attrs = {k.lower(): v for k, v in annot_attrs.items()}

            annot_attrs.update({"extent_polygon": None})
            existing_annot_tag = annot_attrs.get("tag", None)
            annot_geom = annotation_to_shapely(scaled_annot)
            annot_attrs["geometry"] = annot_geom
            annot_attrs["page_label"] = scaled_annot.page
            annot_attrs["tag"] = existing_annot_tag

            if "extent" in annot_attrs["type"]:
                parsed_annotations.update({scaled_annot: annot_attrs})
            else:
                annot_attrs.setdefault("reaction_type", "point")
                annot_attrs["reaction_type"] = annot_attrs["reaction_type"].lower()
                if (
                    annot_geom.geom_type == "Polygon"
                    and annot_attrs["reaction_type"] == "linear"
                ):
                    annot_attrs["length"] = geom.get_rectangle_centerline(
                        annot_geom
                    ).length
                parsed_annotations.update({scaled_annot: annot_attrs})

            parsed_annotations[scaled_annot] = annot_attrs
            raw_annotations[raw_annot] = annot_attrs

            if "occupancy" in annot_attrs:
                load_entries.update({scaled_annot: annot_attrs})
            elif "type" in annot_attrs and "hole" in annot_attrs["type"].lower():
                load_entries.update({scaled_annot: annot_attrs})
            elif "type" in annot_attrs and "trib area" in annot_attrs["type"].lower():
                trib_area_entries.update({scaled_annot: annot_attrs})
            elif "type" not in annot_attrs:
                continue
            else:
                structural_element_entries.update({scaled_annot: annot_attrs})

        elements = Element.from_parsed_annotations(
            structural_element_entries, trib_area_entries
        )
        graph = cls.from_elements(elements, do_not_process=do_not_process)
        graph.parsed_annotations = tag_parsed_annotations(parsed_annotations)
        graph.raw_annotations = tag_parsed_annotations(raw_annotations)
        graph.legend_entries = {}
        graph.loading_geometries = parsed_annotations_to_loading_geometry(load_entries)
        return graph

    def parse_annotations():
        parsed_annotations = {}
        for legend_item in legend:
            legend_properties = {
                prop: getattr(legend_item, prop) for prop in properties_to_match
            }
            matching_annots = filter_annotations(annots, legend_properties)
            annot_attributes = parse_legend(legend_item.text, legend_identifier)
            for annot in matching_annots:
                if annot in legend:
                    continue
                annot_kwargs = parse_annot_kwargs(annot.text)
                existing_annot_tag = annot_kwargs.get("tag", None)
                annot_geom = annotation_to_shapely(annot)
                annot_attrs = {}
                annot_attrs["geometry"] = annot_geom
                annot_attrs["page_label"] = annot.page
                annot_attrs["tag"] = existing_annot_tag
                for annot_key, annot_attr in annot_attributes.items():
                    annot_attrs[annot_key] = str_to_int(
                        annot_attr.split("<")[0]
                    )  # .split() to remove trailing HTML tags

                # Run tests for this first
                # annot_attrs["rank"] = int(annot_attributes["rank"])
                if "extent" in annot_attrs["type"]:
                    parsed_annotations.update({annot: annot_attrs})
                else:
                    annot_attrs.setdefault("reaction_type", "point")
                    annot_attrs["reaction_type"] = annot_attrs["reaction_type"].lower()
                    if (
                        annot_geom.geom_type == "Polygon"
                        and annot_attrs["reaction_type"] == "linear"
                    ):
                        annot_attrs["length"] = geom_ops.get_rectangle_centerline(
                            annot_geom
                        ).length
                    parsed_annotations.update({annot: annot_attrs | annot_kwargs})

    def unassigned_collectors(self) -> list[str]:
        """
        Returns a list of str that represents collector nodes who do not currently
        have collector behaviour assigned meaning that they do not have a populated
        .trib_area attribute.
        """
        unassigned_acc = []
        for collector_node in self.collector_elements:
            node_element = self.nodes[collector_node]["element"]
            if not node_element.trib_area:
                unassigned_acc.append(collector_node)
        return unassigned_acc

    @classmethod
    def from_pdf_file(
        cls,
        pdf_filepath: pathlib.Path | str,
        legend_identifier: str = "legend",
        scale: Optional[Decimal] = None,
        cantilever_abs_tol: Optional[float] = 0.2,
        debug: bool = False,
        progress: bool = False,
        do_not_process: bool = False,
        save_tagged_pdf_file: bool = False,
        tag_pdf_file_mode: str = "append",
        show_skipped: bool = False,
    ):
        """
        Returns a GeometryGraph built from that annotations in the provided PDF file
        at 'pdf_filepath'.

        The provided annotations are parsed into four different categories:
            0. Legend entries - All legend entries must contain the 'legend_identifier'
                as the first piece of text in their text property. Legend entries need
                only appear on ONE page of the PDF document.
            1. Structural elements - All legend entries for structural elements must have
                the legend identifier, a "Type" field (e.g. "Type: <value>"), and a "Rank"
                field (e.g. "Rank: <integer>")
            2. Area load elements - All legend entries for area load elements must have
                an "Occupancy" field (e.g. "Occupancy: <value>")
            3. Trib area elements - All legend entries for trib area elements must have
                a "Type" field with a value of "trib" (e.g. "Type: trib")
            4. Origin elements - All origin elements (max. 1 per page) must have the
                word "origin" as their text element. THE WORD "origin" CANNOT BE USED
                AS PART OF ANY OTHER LEGEND ENTRY (e.g. in the Type, Rank, or Occupancy
                fields)

        'annotations': the list of Annotations
        'legend_identifier': the str used in the text attribute of the PDF annotation to
            indicates a given geometry is part of the legend.
        'scale': An optional scale to be applied to the annotations. If not provided,
            the units of the annotations will be in PDF points where 1 point == 1 /72 inch
        'debug':  When True, will provide verbose documentation of the annotation parsing
            process to assist in reviewing errors and geometry inconsistencies.
        'progress': When True, a progress bar will be displayed
        'do_not_process': Reads the file and adds annotations to the graph but does not
            process the connectivity. Useful for debugging and plotting prior to processing.
        'show_skipped': Shows the skipped annotations that occured during pdf.load_pdf_annotations
        """
        annotations = pdf.load_pdf_annotations(pdf_filepath, show_skipped)
        graph = cls.from_annotations(
            annotations,
            legend_identifier,
            scale=scale,
            do_not_process=do_not_process,
            cantilever_abs_tol=cantilever_abs_tol,
        )
        graph.pdf_path = pathlib.Path(pdf_filepath).resolve()
        return graph

    @classmethod
    def from_annotations(
        cls,
        annotations: list[Annotation],
        legend_identifier: str = "legend",
        scale: Optional[Decimal] = None,
        cantilever_abs_tol: Optional[float] = 0.02,
        # area_load_properties: Optional[dict] = None,
        # trib_area_properties: Optional[dict] = None,
        debug: bool = False,
        progress: bool = False,
        do_not_process: bool = False,
    ):
        """
        Returns a GeometryGraph built from the provided annotations.

        The provided annotations are parsed into four different categories:
            0. Legend entries - All legend entries must contain the 'legend_identifier'
                as the first piece of text in their text property. Legend entries need
                only appear on ONE page of the PDF document.
            1. Structural elements - All legend entries for structural elements must have
                the legend identifier, a "Type" field (e.g. "Type: <value>"), and a "Rank"
                field (e.g. "Rank: <integer>")
            2. Area load elements - All legend entries for area load elements must have
                an "Occupancy" field (e.g. "Occupancy: <value>")
            3. Trib area elements - All legend entries for trib area elements must have
                a "Type" field with a value of "trib" (e.g. "Type: trib")
            4. Origin elements - All origin elements (max. 1 per page) must have the
                word "origin" as their text element. THE WORD "origin" CANNOT BE USED
                AS PART OF ANY OTHER LEGEND ENTRY (e.g. in the Type, Rank, or Occupancy
                fields)

        'annotations': the list of Annotations
        'legend_identifier': the str used in the text attribute of the PDF annotation to
            indicates a given geometry is part of the legend.
        'scale': An optional scale to be applied to the annotations. If not provided,
            the units of the annotations will be in PDF points where 1 point == 1 /72 inch
        'debug':  When True, will provide verbose documentation of the annotation parsing
            process to assist in reviewing errors and geometry inconsistencies.
        'progress': When True, a progress bar will be displayed
        'do_not_process': Reads teh fille and adds annotations to the graph but does not
            process connectivity. Useful for debugging.
        """
        annots = annotations
        page_ids = sorted(set([annot.page for annot in annots]), reverse=True)
        legend_entries = [
            annot for annot in annotations if legend_identifier in annot.text.lower()
        ]
        non_legend_entries = [
            annot
            for annot in annotations
            if legend_identifier not in annot.text.lower()
        ]
        page_entries = [annot for annot in annotations if "page" in annot.text.lower()]
        origin_entries = [
            annot for annot in annotations if "origin" in annot.text.lower()
        ]
        if page_entries:
            if len(page_entries) != len(origin_entries):
                raise AnnotationError(
                    "An 'origin' annotation must be present for each 'page' annotation. "
                    f"{len(page_entries)=} | {len(origin_entries)=}"
                )
            other_annots = [annot for annot in annotations if annot not in page_entries]
            annots_by_page = assign_page_id_to_annotations(other_annots, page_entries)
        else:
            annots_by_page = [
                [annot for annot in non_legend_entries if annot.page == page_id]
                for page_id in page_ids
            ]
        load_entries = {}
        trib_area_entries = {}
        extent_entries = {}
        structural_element_entries = {}
        parsed_annotations_acc = {}
        raw_annotations_acc = {}
        for annots_in_page in annots_by_page:
            if scale is not None:
                scaled_annots_in_page = scale_annotations(annots_in_page, scale)

            # Separate annotation types
            parsed_annotations = parse_annotations(
                scaled_annots_in_page, legend_entries, legend_identifier
            )
            raw_annotations = parse_annotations(
                annots_in_page, legend_entries, legend_identifier
            )
            parsed_annotations_acc = parsed_annotations | parsed_annotations_acc
            raw_annotations_acc = raw_annotations | raw_annotations_acc
            for annot, annot_attrs in parsed_annotations.items():
                if "occupancy" in annot_attrs:
                    load_entries.update({annot: annot_attrs})
                elif "trib area" in annot_attrs.get("type", "").lower():
                    trib_area_entries.update({annot: annot_attrs})
                elif "extent" in annot_attrs.get("type", "").lower():
                    extent_entries.update({annot: annot_attrs})
                else:
                    structural_element_entries.update({annot: annot_attrs})
            structural_element_entries = correlate_extents(
                structural_element_entries, extent_entries
            )

        if not structural_element_entries and not legend_entries:
            raise ValueError(
                "No structural element entities were found.\n"
                "Did you forget to create some legend entries?"
            )

        elif not structural_element_entries:
            raise ValueError(
                "No structural element entities were found.\n"
                "Do the annotation properties match the legend entry properties?"
            )

        elements = Element.from_parsed_annotations(
            structural_element_entries, trib_area_entries
        )
        graph = cls.from_elements(
            elements,
            cantilever_abs_tol=cantilever_abs_tol,
            do_not_process=do_not_process,
        )
        graph.parsed_annotations = tag_parsed_annotations(parsed_annotations_acc)
        graph.raw_annotations = tag_parsed_annotations(raw_annotations_acc)
        graph.legend_entries = legend_entries
        graph.loading_geometries = parsed_annotations_to_loading_geometry(load_entries)
        return graph

    def plot_connectivity(self, filepath: Optional[pathlib.Path | str] = None) -> None:
        """
        Using GraphViz, this function plots the connectivity, a representation of
        a load path, between members where each member is represented as a node in
        the graph and each connection is represented by a directional edge (arrow).

        If 'filepath' is provided, then the resulting SVG image of the graph will
        be saved to disk.

        If 'filepath' is None, then display of the SVG will be attempted through
        IPython.display.

        If the designer is not within a Jupyter-like environment, then the SVG
        string will be returned.

        Requires GraphViz to be independently installed. Installation instructions
        here: https://graphviz.org/download/
        """
        try:
            plotting = nx.drawing.nx_agraph.to_agraph(self)
        except:
            raise ImportError(
                "The GraphViz application is missing. Install at https://graphviz.org/download/"
            )
        plotting.layout(prog="dot")
        svg_string = plotting.draw(format="svg")
        if filepath:
            with open(filepath, "wb") as file:
                file.write(svg_string)
                return
        else:
            try:
                from IPython.display import SVG, display
            except ImportError:
                print(svg_string)
                return
            try:
                get_ipython
            except NameError:
                display(svg_string)
                return
            display(SVG(svg_string))

    def plot_annotations(
        self, page_idx: int, figsize: tuple[float, float] = (8, 8), dpi: int = 150
    ):
        """
        Plots all annotations in self.parsed_annotations that are on the 'page_idx'
        """
        annots = {
            annot: attrs
            for annot, attrs in self.parsed_annotations.items()
            if annot.page == page_idx and annot not in self.legend_entries
        }
        return plot_annotations(annots, figsize, dpi, plot_tags=True)

    def plot_elements(
        self,
        plane_id: int,
        figsize: tuple[float, float] = (8, 8),
        dpi: int = 150,
        plot_trib_areas: bool = False,
        plot_extent_polygons: bool = False,
        plot_tags: bool = False,
    ):
        """
        Plots all elements in the graph that are on 'page_idx'
        """
        elements = [self.nodes[node_name]["element"] for node_name in self.nodes.keys()]
        elements_on_page = [e for e in elements if e.plane_id == plane_id]
        return plot_elements(
            elements_on_page,
            figsize,
            dpi,
            plot_trib_areas=plot_trib_areas,
            plot_extent_polygons=plot_extent_polygons,
            plot_tags=plot_tags,
        )

    def create_loaded_elements(self) -> dict[str, LoadedElement]:
        """
        Returns a list of LoadedElement, each with 'loading_areas' applied.

        # TODO: Is there a way to include trib areas? A dict of trib areas where the key is the Element.geometry
        # and the value is the trib area Polygon? Perhaps a way to specify a buffer value (or a left/right) value
        # to generate one from thg Element.geometry and some integers?

        # HERE: Need to find a way to add raw load annotations to the graph so that they cann
        # automatically sort themselves by plane_id so that the right loads go to the right Elements
        """
        collector_elements = self.collector_elements
        loading_geoms = self.loading_geometries
        loading_geoms_by_plane = {}
        for loading_geom in loading_geoms:
            lg_plane = loading_geom.plane_id
            loading_geoms_by_plane.setdefault(lg_plane, [])
            loading_geoms_by_plane[lg_plane].append(loading_geom)

        loaded_elements = {}
        for node in nx.topological_sort(self):
            node_attrs = self.nodes[node]
            element = node_attrs["element"]
            element.element_type = (
                "collector" if node in collector_elements else "transfer"
            )
            element_plane_id = node_attrs["element"].plane_id
            loading_geoms_on_plane = loading_geoms_by_plane.get(element_plane_id, [])
            # Using predecessors and successors allows us to easily remove incorrect edges
            # that main be contained with individual elements. Specifically, the correspondents
            # above do not get included. The graph is defined entirely from edges pointing "downward".
            predecessors = list(self.predecessors(node))
            successors = list(self.successors(node))
            if element.element_type == "collector" and element.subelements is not None:
                for sub_elem in element.subelements:
                    le = LoadedElement.from_element_with_loads(
                        sub_elem,
                        loading_geoms=loading_geoms_on_plane,
                        predecessors=predecessors,
                        successors=successors,
                    )
                    loaded_elements.update({sub_elem.tag: le})
            else:
                le = LoadedElement.from_element_with_loads(
                    node_attrs["element"],
                    loading_geoms=loading_geoms_on_plane,
                    predecessors=predecessors,
                    successors=successors,
                )
                loaded_elements.update({node: le})
        return loaded_elements

    def export_tagged_pdf(
        self, export_path: Optional[pathlib.Path | str] = None, mode: str = "append"
    ) -> None:
        """
        Returns None. Generates a copy of the PDF file at self.pdf_path
        with the element tags added to the text field of each annotation
        which represents a structural element (e.g. "tag: FB0.1")
        'export_path' - If not provided, the export is stored in the same
            directory as self.pdf_path with "-tagged" appended to the filename.
        'mode' - One of {"append", "replace"}. If "append", the tag field is
        appended to the end of the existing annotation using a new line character
        as a separator.
        """
        if mode.lower() == "append":
            append = True
        elif mode.lower() == "replace":
            append = False
        else:
            raise ValueError(
                f'tag_pdf_file_mode must be one of {"append", "replace"}, not {mode=}'
            )
        if export_path is None:
            new_filename = f"{self.pdf_path.stem}-tagged{self.pdf_path.suffix}"
            export_path = pathlib.Path(self.pdf_path).with_name(new_filename)
        pdf.update_pdf_annotations(
            self.pdf_path, self.raw_annotations, export_path, append
        )

    def hash_nodes(self):
        """
        Returns None. Sets the value of self.node_hash based on the hashed values of
        the nodes.
        """
        nodes_from_top = nx.topological_sort(self)
        hashes = []
        for node_name in nodes_from_top:
            element_hash = self.nodes[node_name]["sha256"]
            hashes.append(element_hash)
        graph_hash = hashlib.sha256(str(tuple(hashes)).encode()).hexdigest()
        self.node_hash = graph_hash


def get_local_coords(
    start_coord: Point | tuple[float, float],
    dependent_intersections: list[Intersection],
) -> list[tuple[float, str]]:
    start_coord = Point(start_coord)
    intersection_below_local_coords = []
    for intersection in dependent_intersections:
        below_local_coord = start_coord.distance(intersection.intersecting_region)
        intersection_below_local_coords.append(
            (below_local_coord, intersection.other_tag)
        )
    return intersection_below_local_coords


def get_dependent_intersections(
    element: Element, graph_dependents: list[str]
) -> list[Intersection]:
    acc = []
    for intersection in element.intersections_below:
        if intersection.other_tag in graph_dependents:
            acc.append(intersection)
    return acc


def get_dependent_correspondents(
    element: Element, graph_dependents: list[str]
) -> list[Correspondent]:
    dependent_correspondents = [
        correspondent
        for correspondent in element.correspondents_below
        if correspondent.other_tag in graph_dependents
    ]
    return dependent_correspondents


def correlate_extents(
    element_annots: dict[Annotation, dict], extent_annots: dict[Annotation, dict]
) -> dict:
    """
    Returns a copy of 'element_annots' with value dicts that have been updated
    to include an 'extent_polygon' for any element_annots that have been drawn
    with an extent line or polygon.
    """
    element_geoms = [annot_attrs["geometry"] for annot_attrs in element_annots.values()]
    element_annot_keys = [annot for annot in element_annots.keys()]
    extent_geoms = [annot_attrs["geometry"] for annot_attrs in extent_annots.values()]

    matched_extents = geom.find_extent_intersections(element_geoms, extent_geoms)
    element_annots_copy = deepcopy(element_annots)
    for idx, matched_extent in enumerate(matched_extents):
        annot = element_annot_keys[idx]
        element_geom = element_geoms[idx]
        element_annots_copy[annot]["extent_line"] = matched_extent
    return element_annots_copy
