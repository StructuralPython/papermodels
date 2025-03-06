from __future__ import annotations
from typing import Optional
from copy import deepcopy
import networkx as nx
import hashlib

from papermodels.datatypes.element import Element
from .loads import LoadElement, LoadArray
from rich.progress import track


class LoadGraph(nx.DiGraph):
    """
    A class to represent a connected geometry system in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

    Can be used to generate a GeometryGraph.

    The node_hash is how changes to the graph nodes can be tracked.
    """

    def __init__(self):
        super().__init__()
        self.node_hash = None

    @classmethod
    def from_elements(
        cls, elements: list[Element], floor_elevations: Optional[dict] = None
    ) -> LoadGraph:
        """
        Returns a GeometryGraph (networkx.DiGraph) based upon the intersections and correspondents
        of the 'elements'.
        """
        g = cls()
        for element in elements:
            hash = hashlib.sha256(str(element).encode()).hexdigest()
            g.add_node(element.tag, element=element, sha256=hash)
            if element.correspondents is not None:
                for correspondent in element.correspondents:
                    g.add_edge(element.tag, correspondent)
            if element.intersections is not None:
                for intersecting_element_tag in element.intersections.keys():
                    g.add_edge(element.tag, intersecting_element_tag)
        return g

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



# def 
