from __future__ import annotations
from typing import Optional
from copy import deepcopy
import networkx as nx
import hashlib

from papermodels.datatypes.element import Element, LoadedElement
from shapely import Point, LineString, Polygon
from ..geometry import geom_ops as geom
from ..datatypes.element import Correspondent, Intersection
from rich.progress import track
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
                    intersection_below.other_tag: intersection_below.other_index
                    for intersection_below in element_above.intersections_below
                }
                local_index = above_intersections_below[element_tag]
                new_intersection = Intersection(
                    intersection.intersecting_region,
                    intersection.other_geometry,
                    intersection.other_tag,
                    local_index
                )
                indexed_intersections_above.append(new_intersection)
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


    def create_loaded_elements(self, loading_geoms: list[tuple[Polygon, npt.ArrayLike]]) -> list[LoadedElement]:
        """
        Returns a list of LoadedElement, each with 'loading_areas' applied.
        
        # TODO: Is there a way to include trib areas? A dict of trib areas where the key is the Element.geometry
        # and the value is the trib area Polygon? Perhaps a way to specify a buffer value (or a left/right) value
        # to generate one from thg Element.geometry and some integers?

        # HERE: Need to find a way to add raw load annotations to the graph so that they cann
        # automatically sort themselves by plane_id so that the right loads go to the right Elements
        """
        loading_geoms_by_plane = {}
        for loading_geom in loading_geoms:
            lg_plane = loading_geom.plane_id
            loading_geoms_by_plane.setdefault(lg_plane, [])
            loading_geoms_by_plane[lg_plane].append(loading_geom)

        loaded_elements = []
        for node in self.nodes:
            node_attrs = self.nodes[node]
            element_plane_id = node_attrs['element'].plane_id
            loading_geoms_on_plane = loading_geoms_by_plane.get(element_plane_id, [])
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



    def apply_loading(self, loads: list[LoadElement] | dict[str, LoadElement]) -> None:
        """
        Mutates the graph to include the loads applied.
        """
