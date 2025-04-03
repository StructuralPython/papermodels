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
from ..datatypes.element import Correspondent, Intersection
from ..paper.annotations import (
    Annotation, 
    scale_annotations, 
    parse_annotations, 
    parsed_annotations_to_loading_geometry
)
from ..paper.pdf import load_pdf_annotations
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

    def __init__(self):
        super().__init__()
        self.node_hash = None
        self.loading_geometries = None

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
        cls, elements: list[Element]
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
                ordered_coords = geom.order_nodes_positive(Point(coords_a), Point(coords_b))
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
                    # print("correspondent: ", j_tag)
                    g.add_edge(element.tag, j_tag)
            if element.intersections_below is not None:
                for intersection in element.intersections_below:
                    j_tag = intersection.other_tag
                    g.add_edge(element.tag, j_tag)
                
        g.add_intersection_indexes_below()
        g.add_intersection_indexes_above()

        return g
    
    def add_intersection_indexes_below(self):
        sorted_nodes = nx.topological_sort(self)
        for node in sorted_nodes:
            node_attrs = self.nodes[node]
            element = node_attrs['element']
            if node_attrs['start_coord'] is None: # node geometry is polygon
                updated_intersections_below = []
                for intersection in element.intersections_below:
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        intersection.other_geometry,
                        intersection.other_tag,
                        0
                    )
                    updated_intersections_below.append(new_intersection)
            else:
                start_coord = Point(node_attrs['start_coord'])
                intersection_below_local_coords = []
                for intersection in node_attrs['element'].intersections_below:
                    below_local_coord = start_coord.distance(intersection.intersecting_region)
                    intersection_below_local_coords.append((below_local_coord, intersection.other_tag))
                sorted_below_ints = sorted(intersection_below_local_coords, key=lambda x: x[0])
                _, other_tags_below = zip(*sorted_below_ints)
                updated_intersections_below = []
                for intersection in element.intersections_below:
                    other_tag = intersection.other_tag
                    local_index = other_tags_below.index(other_tag)
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        intersection.other_geometry,
                        intersection.other_tag,
                        local_index
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
                element_above = self.nodes[other_tag]['element']
                above_intersections_below = {
                    above_intersection_below.other_tag: above_intersection_below.other_index
                    for above_intersection_below in element_above.intersections_below
                }
                local_index = above_intersections_below[element_tag]
                if element_above.subelements is None:
                    new_intersection = Intersection(
                        intersection.intersecting_region,
                        intersection.other_geometry,
                        intersection.other_tag,
                        local_index
                    )
                    indexed_intersections_above.append(new_intersection)
                else:
                    for subelem_above in element_above.subelements:
                        sub_intersection = [inter for inter in subelem_above.intersections_below if inter.other_tag == element_tag][0]
                        new_sub_intersection = Intersection(
                            sub_intersection.intersecting_region,
                            subelem_above.geometry,
                            subelem_above.tag,
                            local_index
                        )
                        indexed_intersections_above.append(new_sub_intersection)


            element.intersections_above = indexed_intersections_above
            self.nodes[node]['element'] = element


    def generate_subelements(self, subelement_constructor: callable, *args, **kwargs) -> list[Element]:
        """
        Returns a list of Element to be assigned to element.subelements for elements
        that have element_type == "collector".

        'subelement_class': This should be a callable with the following signature:
            def subelement_constructor(element: Element, [*args, **kwargs]) -> list[Element]

            Where *args, and **kwargs can be any additional parameters that are defined
            for the callable.
        '*args' and '**kwargs': These are passed through to 'subelement_constructor'
        """
        collectors = self.collector_elements
        for node in collectors:
            node_attrs = self.nodes[node]
            node_element = node_attrs['element']
            subs = subelement_constructor(node_element, *args, **kwargs)
            node_element.subelements = subs
        self.add_intersection_indexes_below()
        self.add_intersection_indexes_above()

    @classmethod
    def from_pdf_file(
        cls,
        filepath: pathlib.path | str,
        legend_identifier: str = "legend",
        scale: Optional[Decimal] = None,
        debug: bool = False,
        progress: bool = False
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
        """
        annotations = load_pdf_annotations(filepath)
        return cls.from_annotations(annotations, legend_identifier, scale=scale)
        


    @classmethod
    def from_annotations(
        cls,
        annotations: list[Annotation],
        legend_identifier: str = "legend",
        scale: Optional[Decimal] = None,
        # area_load_properties: Optional[dict] = None,
        # trib_area_properties: Optional[dict] = None,
        debug: bool = False,
        progress: bool = False
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
        """
        annots = annotations
        page_ids = sorted(set([annot.page for annot in annots]), reverse=True)
        legend_entries = [annot for annot in annotations if legend_identifier in annot.text.lower()]
        non_legend_entries = [annot for annot in annotations if legend_identifier not in annot.text.lower()]
        annots_by_page = [
            [annot for annot in non_legend_entries if annot.page == page_id] for page_id in page_ids
        ]
        load_entries = {}
        trib_area_entries = {}
        structural_element_entries = {}
        for annots_in_page in annots_by_page:

            # Scale taking origin into account
            origin_annots = [annot for annot in annots_in_page if "origin" in annot.text.lower()]
            origin_on_page = None
            if len(origin_annots) > 1:
                raise ValueError(f"There should only be one origin per page. Found: {origin_annots=}")
            elif origin_annots == 1:
                origin_annot = origin_annots[0]
                origin_on_page = parse_origin_annotation(origin_annot)

            if scale is not None:
                annots_in_page = scale_annotations(annots_in_page, scale, origin_on_page)

            # Separate annotation types
            parsed_annotations = parse_annotations(annots_in_page, legend_entries, legend_identifier)
            for annot, annot_attrs in parsed_annotations.items():
                if "occupancy" in annot_attrs:
                    load_entries.update({annot: annot_attrs})
                elif "type" in annot_attrs and annot_attrs['type'].lower() == "trib":
                    trib_area_entries.update({annot: annot_attrs})
                else:
                    structural_element_entries.update({annot: annot_attrs})

        elements = Element.from_parsed_annotations(structural_element_entries)
        graph = cls.from_elements(elements)
        graph.loading_geometries = parsed_annotations_to_loading_geometry(load_entries)
        return graph


    def plot(self):
        return nx.draw_spectral(self, with_labels=True)


    def create_loaded_elements(self) -> list[LoadedElement]:
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

        loaded_elements = []
        for node in self.nodes:
            node_attrs = self.nodes[node]
            element = node_attrs['element']
            element.element_type = "collector" if node in collector_elements else "transfer"
            element_plane_id = node_attrs['element'].plane_id
            loading_geoms_on_plane = loading_geoms_by_plane.get(element_plane_id, [])
            if element.element_type == "collector" and element.subelements is not None:
                for sub_elem in element.subelements:
                    le = LoadedElement.from_element_with_loads(sub_elem, loading_geoms=loading_geoms_on_plane)
                    loaded_elements.append(le)
            else:
                le = LoadedElement.from_element_with_loads(node_attrs['element'], loading_geoms=loading_geoms_on_plane)
                loaded_elements.append(le)
        return loaded_elements
            


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