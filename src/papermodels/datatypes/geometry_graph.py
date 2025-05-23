from __future__ import annotations
from typing import Optional
from copy import deepcopy
from decimal import Decimal
import pathlib
import networkx as nx
import hashlib

from papermodels.datatypes.element import Element, LoadedElement
from shapely import Point, LineString, Polygon
from ..geometry import geom_ops as geom
from ..datatypes.element import Correspondent, Intersection, get_collector_extents, get_transfer_extents
from ..paper.annotations import (
    Annotation, 
    scale_annotations, 
    parse_annotations, 
    parsed_annotations_to_loading_geometry,
    filter_annotations,
    tag_parsed_annotations,
    assign_page_id_to_annotations
)
from ..paper.plot import plot_annotations
from ..paper import pdf
from ..datatypes.exceptions import AnnotationError
from rich.progress import track
from rich import print
import numpy.typing as npt


class GeometryGraph(nx.DiGraph):
    """
    A class to represent a connected geometry system in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

    Can be used to generate a GeometryGraph.

    The node_hash is how changes to the graph nodes can be tracked.
    """

    def __init__(self, do_not_process: bool = False):
        super().__init__()
        self.do_not_process = do_not_process
        self.node_hash = None
        self.loading_geometries = None
        self.parsed_annotations = None
        self.raw_annotations = None
        self.legend_entries = None
        self.pdf_path = None

    @property
    def collector_elements(self):
        return [
            node for node in self.nodes if not list(self.predecessors(node))
        ]
    
    @property
    def transfer_elements(self):
        return [
            node for node in self.nodes if list(self.predecessors(node))
        ]

    @classmethod
    def from_elements(
        cls, elements: list[Element], do_not_process: bool = False
    ) -> GeometryGraph:
        """
        Returns a GeometryGraph (networkx.DiGraph) based upon the intersections and correspondents
        of the 'elements'.
        """
        g = cls()
        elements_copy = deepcopy(elements)
        for element in elements_copy:
            hash = hashlib.sha256(str(element).encode()).hexdigest()
            start_coord = None
            if element.geometry.geom_type == "LineString":
                coords_a, coords_b = element.geometry.coords
                ordered_coords = geom.order_nodes_positive([Point(coords_a), Point(coords_b)])
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
                for intersection in element.intersections_below:
                    j_tag = intersection.other_tag
                    g.add_edge(element.tag, j_tag, edge_type="intersection")
            if element.tag in g.collector_elements:
                for correspondent in element.correspondents_above:
                    j_tag = correspondent.other_tag
                    g.add_edge(j_tag, element.tag, edge_type="correspondent")

        for node in g.collector_elements:
            g.nodes[node]['element'].element_type = "collector"
        
        for node in g.transfer_elements:
            g.nodes[node]['element'].element_type = "transfer"
        
        g.add_intersection_indexes_below()
        g.add_intersection_indexes_above()
        
        if do_not_process:
            return g

        g.remove_excess_correspondent_load_paths()
        return g
    

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
            then the "correspondent" edges will be remove. This represents the condition
            of platform-framing where a post will land on the floor framing, and transfer
            through it to the supporting post below.

        Modifications to the implementation of this function can adjust how load paths are
        conceptually created. For example, to implement baloon framing, the second rule
        can be omitted.
        """
        sorted_nodes = nx.topological_sort(self)
        for node in sorted_nodes:
            element = self.nodes[node]['element']
            dependents = list(self.successors(node))
            dependent_edges = [(node, dep) for dep in dependents]
            edge_properties = [self.edges[edge]['edge_type'] for edge in dependent_edges]
            # Rule 1
            if element.geometry.geom_type == "Polygon" and edge_properties.count("correspondent") > 1:
                dep_to_keep = None
                max_overlap = 0.0
                for idx, dep in enumerate(dependents):
                    # Keep the rank 0
                    if self.nodes[dep]['element'].rank == 0:
                        dep_to_keep = idx
                        break
                    else:
                        # Or find the correspondent with the largest overlap ratio
                        try:
                            dep_overlap_ratio = next((corr.overlap_ratio for corr in element.correspondents_below if corr.other_tag==dep))
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
                            dependent_geometry = self.nodes[dep]['element'].geometry
                            element_geometry = element.geometry
                            if dependent_geometry.geom_type == "LineString":
                                overlap_length = element_geometry.intersection(dependent_geometry).length
                                if overlap_length > dep_overlap_length:
                                    dep_to_keep = idx
                                    dep_overlap_length = overlap_length
                            elif dependent_geometry.geom_type == "Polygon":
                                dep_element = self.nodes[dep]['element']
                                if dep_element.reaction_type == "point": # Points should transfer to points
                                    overlap_area = element_geometry.intersection(dep_element.geometry).area
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
                

    def add_intersection_indexes_below(self):
        sorted_nodes = nx.topological_sort(self)
        for node in sorted_nodes:
            node_attrs = self.nodes[node]
            element: Element = node_attrs['element']
            if node_attrs['start_coord'] is None: # node geometry is polygon
                updated_intersections_below = []
                all_extents = {}
                if element.reaction_type == "linear":
                    all_extents = get_transfer_extents(element)
                for intersection in element.intersections_below:
                    extents = all_extents.get(intersection.other_tag)
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        intersection.other_geometry,
                        intersection.other_tag,
                        0,
                        intersection.other_reaction_type,
                        other_extents=extents
                    )
                    updated_intersections_below.append(new_intersection)
            else:
                start_coord = Point(node_attrs['start_coord'])
                intersection_below_local_coords = []
                for intersection in node_attrs['element'].intersections_below:
                    below_local_coord = start_coord.distance(intersection.intersecting_region)
                    intersection_below_local_coords.append((below_local_coord, intersection.other_tag))
                sorted_below_ints = sorted(intersection_below_local_coords, key=lambda x: x[0])
                if len(sorted_below_ints) < 2:
                    raise ValueError(f"It seems that this element only has one support: {node}")
                _, other_tags_below = zip(*sorted_below_ints)
                updated_intersections_below = []
                extents = {}
                if element.element_type == "collector" and element.reaction_type == "linear":
                    extents = get_collector_extents(element)
                for intersection in element.intersections_below:
                    other_tag = intersection.other_tag

                    local_index = other_tags_below.index(other_tag)
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        intersection.other_geometry,
                        intersection.other_tag,
                        local_index,
                        other_reaction_type=intersection.other_reaction_type,
                        other_extents=extents.get(other_tag, None)
                    )
                    updated_intersections_below.append(new_intersection)
                if node in self.collector_elements and element.subelements is not None:

                    for subelem in element.subelements:
                        sub_updated_intersections_below = []
                        for sub_intersection in subelem.intersections_below:
                            sub_other_tag = sub_intersection.other_tag
                            sub_local_index = other_tags_below.index(sub_other_tag)
                            new_sub_intersection = Intersection(
                                sub_intersection.intersecting_region,
                                sub_intersection.other_geometry,
                                sub_intersection.other_tag,
                                sub_local_index
                            )
                            sub_updated_intersections_below.append(new_sub_intersection)
                        subelem.intersections_below = sub_updated_intersections_below

            element.intersections_below = updated_intersections_below
            self.nodes[node]['element'] = element


    def add_intersection_indexes_above(self):
        sorted_nodes = nx.topological_sort(self)
        transfer_elements = [node for node in sorted_nodes if list(self.predecessors(node))]

        for node in transfer_elements:
            indexed_intersections_above = []
            element = self.nodes[node]['element']
            element_tag = element.tag
            for intersection in self.nodes[node]['element'].intersections_above:
                other_tag = intersection.other_tag
                element_above: Element = self.nodes[other_tag]['element']
                above_intersections_below = {
                    above_intersection_below.other_tag: (
                        above_intersection_below.other_index,
                        above_intersection_below.other_extents
                    )
                    for above_intersection_below in element_above.intersections_below
                }
                local_index = above_intersections_below[element_tag][0]
                other_extents = above_intersections_below[element_tag][1]
                if element_above.subelements is None:
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        intersection.other_geometry,
                        intersection.other_tag,
                        local_index,
                        element_above.reaction_type,
                        other_extents=other_extents
                    )
                    indexed_intersections_above.append(new_intersection)
                else:
                    for subelem_above in element_above.subelements:
                        sub_intersection = [inter for inter in subelem_above.intersections_below if inter.other_tag == element_tag][0]
                        new_sub_intersection = Intersection(
                            sub_intersection.intersecting_region,
                            subelem_above.geometry,
                            subelem_above.tag,
                            local_index,
                            element_above.reaction_type,
                            # No collector_extents because this is an element with subelements (all points)
                        )
                        indexed_intersections_above.append(new_sub_intersection)


            element.intersections_above = indexed_intersections_above
            self.nodes[node]['element'] = element


    def assign_collector_behaviour(
            self, 
            element_constructor: callable, 
            as_subelements: bool,
            *args,
            **kwargs
    ) -> list[Element]:
        """
        Returns a list of Element to be assigned to element.subelements for elements
        that have element_type == "collector".

        'subelement_class': This should be a callable with the following signature:
            def subelement_constructor(element: Element, [*args, **kwargs]) -> list[Element]

            Where *args, and **kwargs can be any additional parameters that are defined
            for the callable.
        '*args' and '**kwargs': These are passed through to 'subelement_constructor'
        """
        # TODO: Add the ability to filter collector elements to assign different
        # behaviours based on rules
        collectors = self.collector_elements
        for node in collectors:
            node_attrs = self.nodes[node]
            node_element = node_attrs['element']
            # If an incorrect geometry type makes its way into the element_constructor
            # e.g. a polygon is being entered as a joist in a joist-based element_constructor
            # then the element_constructor should return None
            # This prevents an error from being thrown if, for example, unconnected elements
            # are drawn. Unconnected elements have no precedents therefore they are (currently)
            # being categorized as collectors. However, I think incompatible geometries
            # should simply be ignored and not included as part of the processing.
            new_elem = element_constructor(node_element, *args, **kwargs)
            if new_elem is not None:
                if as_subelements:
                    node_element.subelements = new_elem
                else:
                    # assert isinstance(new_elem, Element) # Not an iterable of multiple elements
                    node_attrs['element'] = new_elem
        self.add_intersection_indexes_below()
        self.add_intersection_indexes_above()


    @classmethod
    def from_pdf_file(
        cls,
        pdf_filepath: pathlib.path | str,
        legend_identifier: str = "legend",
        scale: Optional[Decimal] = None,
        debug: bool = False,
        progress: bool = False,
        do_not_process: bool = False,
        save_tagged_pdf_file: bool = False,
        tag_pdf_file_mode: str = "append",
        show_skipped: bool = False,
    ):
        """
        Returns a GeometryGraph built from that annotations in the provided PDF file
        at 'filepath'.

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
        graph = cls.from_annotations(annotations, legend_identifier, scale=scale, do_not_process=do_not_process)
        graph.pdf_path = pathlib.Path(pdf_filepath).resolve()
        return graph



    @classmethod
    def from_annotations(
        cls,
        annotations: list[Annotation],
        legend_identifier: str = "legend",
        scale: Optional[Decimal] = None,
        # area_load_properties: Optional[dict] = None,
        # trib_area_properties: Optional[dict] = None,
        debug: bool = False,
        progress: bool = False,
        do_not_process: bool = False
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
        legend_entries = [annot for annot in annotations if legend_identifier in annot.text.lower()]
        non_legend_entries = [annot for annot in annotations if legend_identifier not in annot.text.lower()]
        page_entries = [annot for annot in annotations if "page" in annot.text.lower()]
        origin_entries = [annot for annot in annotations if "origin" in annot.text.lower()]
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
                [annot for annot in non_legend_entries if annot.page == page_id] for page_id in page_ids
            ]
        load_entries = {}
        trib_area_entries = {}
        structural_element_entries = {}
        parsed_annotations_acc = {}
        raw_annotations_acc = {}
        for annots_in_page in annots_by_page:
            if scale is not None:
                scaled_annots_in_page = scale_annotations(annots_in_page, scale)

            # Separate annotation types
            parsed_annotations = parse_annotations(scaled_annots_in_page, legend_entries, legend_identifier)
            raw_annotations = parse_annotations(annots_in_page, legend_entries, legend_identifier)
            parsed_annotations_acc = parsed_annotations | parsed_annotations_acc
            raw_annotations_acc = raw_annotations | raw_annotations_acc
            for annot, annot_attrs in parsed_annotations.items():
                if "occupancy" in annot_attrs:
                    load_entries.update({annot: annot_attrs})
                elif "type" in annot_attrs and "trib area" in annot_attrs['type'].lower():
                    trib_area_entries.update({annot: annot_attrs})
                else:
                    structural_element_entries.update({annot: annot_attrs})

        elements = Element.from_parsed_annotations(structural_element_entries, trib_area_entries)
        graph = cls.from_elements(elements, do_not_process=do_not_process)
        graph.parsed_annotations = tag_parsed_annotations(parsed_annotations_acc)
        graph.raw_annotations = tag_parsed_annotations(raw_annotations_acc)
        graph.legend_entries = legend_entries
        graph.loading_geometries = parsed_annotations_to_loading_geometry(load_entries)
        return graph


    def plot_connectivity(self):
        return nx.draw_spectral(self, with_labels=True)
    

    def plot_annotations(self, page_idx: int, figsize: tuple[float, float]=(8, 8), dpi: int = 150):
        """
        Plots all annotations in self.parsed_annotations that are on the 'page_idx'
        """
        annots = {
            annot: attrs 
            for annot, attrs in self.parsed_annotations.items() 
            if annot.page == page_idx and annot not in self.legend_entries
        }
        return plot_annotations(annots, figsize, dpi, plot_tags=True)


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
            element = node_attrs['element']
            element.element_type = "collector" if node in collector_elements else "transfer"
            element_plane_id = node_attrs['element'].plane_id
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
                        successors=successors
                    )
                    loaded_elements.update({sub_elem.tag: le})
            else:
                le = LoadedElement.from_element_with_loads(
                    node_attrs['element'], 
                    loading_geoms=loading_geoms_on_plane, 
                    predecessors=predecessors, 
                    successors=successors
                )
                loaded_elements.update({node: le})
        return loaded_elements
            

    def export_tagged_pdf(
        self,
        export_path: Optional[pathlib.Path | str] = None,
        mode: str = "append"
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
        elif mode.lower() == 'replace':
            append=False
        else:
            raise ValueError(f'tag_pdf_file_mode must be one of {"append", "replace"}, not {mode=}')
        if export_path is None:
            new_filename = f"{self.pdf_path.stem}-tagged{self.pdf_path.suffix}"
            export_path = pathlib.Path(self.pdf_path).with_name(new_filename)
        pdf.update_pdf_annotations(self.pdf_path, self.raw_annotations, export_path, append)


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