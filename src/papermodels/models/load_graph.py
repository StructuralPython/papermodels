from typing import Optional
import hashlib
import numpy as np
from shapely import intersects, Geometry, Point
from papermodels.datatypes.annotation import Annotation
from papermodels.datatypes.element import Element
from papermodels.models.load_factors import nbcc_2020_combos
from papermodels.datatypes.load_graph import LoadGraph
import parse
import networkx as nx

import math

def get_graph_model_from_elements(elements: list[Element], floor_elevations: Optional[dict] = None) -> nx.DiGraph:
    """
    Returns a networkx DiGraph based upon the intersections and correspondents
    of the 'elements'.
    """
    top_down_elements = sorted(elements, key=lambda x: x.page, reverse=True)
    g = LoadGraph()
    for element in top_down_elements:
        hash = hashlib.sha256(str(element).encode()).hexdigest()
        g.add_node(element.tag, element=element, sha256=hash)
        for correspondent in element.correspondents:
            g.add_edge(element.tag, correspondent)
        for intersection in element.intersections:
            g.add_edge(element.tag, intersection[0])
    return g




            

                    



