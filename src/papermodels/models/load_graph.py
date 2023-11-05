from typing import Optional
import numpy as np
from shapely import intersects, Geometry, Point
from papermodels.datatypes.annotation import Annotation
from papermodels.datatypes.element import Element
import parse
import networkx as nx

def get_graph_model_from_elements(elements: list[Element], floor_elevations: Optional[dict] = None) -> nx.DiGraph:
    """
    Returns a networkx DiGraph based upon the intersections and correspondents
    of the 'elements'.
    """
    top_down_elements = sorted(elements, key=lambda x: x.page, reverse=True)
    g = nx.DiGraph()
    for element in top_down_elements:
        g.add_node(element.tag, element=element)
        for intersection in element.intersections:
            g.add_edge(element.tag, intersection[0])
        for correspondent in element.correspondents:
            g.add_edge(element.tag, correspondent)
    return g


def get_new_correspondent_tag(this_element_tag: str, corresponding_element_tag) -> str:
    """
    Returns a new correspondent tag for an element that is connected below
    """
    format = "{type_tag}{page_tag:d}.{enum_tag:d}"
    result = parse.parse(format, this_element_tag)


def get_elements_by_page(elements: list[Element]) -> dict[int, list[Element]]:
    """
    Returns 'elements' sorted by page
    """
    elements_by_page = {}
    for element in elements:
        elements_on_page = elements_by_page.get(element.page, [])
        elements_on_page.append(element)
        elements_by_page[element.page] = elements_on_page
    return elements_by_page


            

                    



