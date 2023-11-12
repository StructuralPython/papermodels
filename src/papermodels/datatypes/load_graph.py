from copy import deepcopy
import networkx as nx
from hashlib import sha256
from papermodels.datatypes.element import (
    element_to_beam_model,
    element_to_joist_model,
    get_tag_type
)

class LoadGraph(nx.DiGraph):
    """
    A class to represent a load path in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

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
            element_hash = self.nodes[node_name]['sha256']
            hashes.append(element_hash)
        graph_hash = sha256(str(tuple(hashes)).encode()).hexdigest()
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
            self.nodes[node]['load_distribution'] = dist


    def compile_distribution_functions(self, prefix_dict: dict = {
        "FB": "beam",
        "DB": "beam",
        "W": "wall",
        "CPL": "point_load",
        "WLL": "line_load",
        "J": "joist",
    }):
        """
        Returns None. Populates the empty 'load_distribution' dict in 
        each node.
        """
        for node in self.nodes:
            element = self.nodes['element']
            tag_prefix = get_tag_type(element.tag)
            if prefix_dict[tag_prefix] == "beam":
                model = element_to_beam_model(element)
                for load_source in self.nodes['load_distribution'].keys():
                    model_copy = deepcopy(model)
                    source_element = self.nodes[load_source]['element']
                    # TODO: Find the location of the source element on the model
                    model.add_member_pt_load(node, "Fy", 1, x=1, )
                self.nodes[node]['model'] = model
            elif prefix_dict[tag_prefix] == "joist":
                model = element_to_joist_model(element)
                self.nodes[node]['model'] = model
            elif prefix_dict[tag_prefix] == "column":
                pass
            elif prefix_dict[tag_prefix] == "wall":
                pass
            elif prefix_dict[tag_prefix] == "point_load":
                pass
            elif prefix_dict[tag_prefix] == "line_load":
                pass
