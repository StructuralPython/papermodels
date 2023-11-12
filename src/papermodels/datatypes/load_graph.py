from dataclasses import dataclass
import networkx as nx
from hashlib import sha256

@dataclass
class LoadGraph(nx.DiGraph):
    """
    A class to represent a load path in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

    The node_hash is how changes to the graph nodes can be tracked.
    """
    node_hash: str


    def hash_nodes(self):
        """
        Returns None. Sets the value of self.node_hash based on the hashed values of
        the nodes.
        """
