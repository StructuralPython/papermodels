from __future__ import annotations
from typing import Optional
from copy import deepcopy
import networkx as nx
import hashlib
from papermodels.datatypes.element import (
    Element,
    element_to_beam_model,
    element_to_joist_model,
    get_tag_type,
    get_normalized_coordinate,
)
from rich.progress import track


class GeometryGraph(nx.DiGraph):
    """
    A class to represent a connected geometry system in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

    Can be used to generate a LoadGraph.

    The node_hash is how changes to the graph nodes can be tracked.
    """

    def __init__(self):
        super().__init__()
        self.node_hash = None


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