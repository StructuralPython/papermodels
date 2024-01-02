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


class LoadGraph(nx.DiGraph):
    """
    A class to represent a load path in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

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
        Returns a LoadGraph (networkx.DiGraph) based upon the intersections and correspondents
        of the 'elements'.
        """
        top_down_elements = sorted(elements, key=lambda x: x.page, reverse=True)
        g = cls()
        for element in top_down_elements:
            hash = hashlib.sha256(str(element).encode()).hexdigest()
            g.add_node(element.tag, element=element, sha256=hash)
            for correspondent in element.correspondents:
                g.add_edge(element.tag, correspondent)
            for intersection in element.intersections:
                g.add_edge(element.tag, intersection[0])
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

    def compile_load_distribution_model(self):
        """
        Returns None. Adds an empty 'load_distribution' dict to each node.
        """
        for node in self.nodes:
            dist = {}
            for pred in self.pred[node]:
                dist.update({pred: {}})
                dist[pred].update(dict(self.succ[node]))
            self.nodes[node]["load_distribution"] = dist

    def compile_distribution_functions(
        self,
        prefix_dict: dict = {
            "FB": "beam",
            "DB": "beam",
            "W": "wall",
            "CPL": "point_load",
            "WLL": "line_load",
            "J": "joist",
            "C": "column",
        },
        track_progress=False,
    ):
        """
        Returns None. Populates the empty 'load_distribution' dict in
        each node.
        """
        if not track_progress:
            progress_tracker = lambda x: x
        else:
            progress_tracker = track
        for node in progress_tracker(self.nodes):
            element = self.nodes[node]["element"]
            tag_prefix = get_tag_type(element.tag)
            if prefix_dict[tag_prefix] == "beam":
                model = element_to_beam_model(element)
                self.nodes[node]["model"] = model
                for load_source, support_tags in self.nodes[node][
                    "load_distribution"
                ].items():
                    # Factor this out to a new function
                    # Note: this function needs to sort out the different
                    # source element types (e.g. joist, column point load, wall, etc.)
                    # so that they apply an appropriate load to the beam.
                    model_copy = deepcopy(self.nodes[node]["model"])
                    source_element = self.nodes[load_source]["element"]
                    intersecting_point = source_element.get_intersection(element.tag)
                    source_loc = get_normalized_coordinate(element, intersecting_point)
                    model_copy.add_member_pt_load(
                        node, "Fy", -1, x=source_loc, case="pass"
                    )
                    model_copy.analyze_linear(check_statics=False)
                    for support_tag in support_tags:
                        reaction = model_copy.Nodes[support_tag].RxnFY["Pass"]
                        support_tags[support_tag] = reaction
            elif prefix_dict[tag_prefix] == "joist":
                model = element_to_joist_model(element)
                self.nodes[node]["model"] = model
            elif prefix_dict[tag_prefix] == "column":
                pass
            elif prefix_dict[tag_prefix] == "wall":
                pass
            elif prefix_dict[tag_prefix] == "point_load":
                pass
            elif prefix_dict[tag_prefix] == "line_load":
                pass
