from typing import Optional
import numpy as np
from shapely import intersects, Geometry, Point
from papermodels.datatypes.annotation import Annotation
from papermodels.datatypes.element import Element
from papermodels.models.load_factors import nbcc_2020_combos
from papermodels.datatypes.joist import Joist
from PyNite import FEModel3D
import parse
import networkx as nx
import math

def get_graph_model_from_elements(elements: list[Element], floor_elevations: Optional[dict] = None) -> nx.DiGraph:
    """
    Returns a networkx DiGraph based upon the intersections and correspondents
    of the 'elements'.
    """
    top_down_elements = sorted(elements, key=lambda x: x.page, reverse=True)
    g = nx.DiGraph()
    for element in top_down_elements:
        for correspondent in element.correspondents:
            g.add_edge(element.tag, correspondent)
        for intersection in element.intersections:
            g.add_edge(element.tag, intersection[0])
    return g


def element_to_beam_model(element: Element) -> FEModel3D:
    """
    Returns an FEModel3D for the data in 'element'
    """
    model = FEModel3D()
    model.add_material("default", 1, 1, 1, 1)
    element_length = element.geometry.length
    start_point = Point(element.geometry.coords[0])
    add_start_node = True
    add_end_node = True
    start_node = "N0"
    model.add_node(start_node, 0, 0, 0)
    for node_name, point in element.intersections:
        x_coord = start_point.distance(point)
        if math.isclose(x_coord, element_length, abs_tol=1e-3):
            last_node = node_name
            add_end_node = False
        if math.isclose(x_coord, 0, abs_tol=1e-3):
            model.Nodes.pop(start_node)
            start_node = node_name
        model.add_node(node_name, x_coord, 0, 0)
        model.def_support(node_name, support_DY=True)
    model.def_support(node_name, 1, 1, 1, 1, 1, 0)

    if add_end_node:
        last_node = "Nn"
        model.add_node(last_node, element_length, 0, 0)
    model.add_member(element.tag, "N0", last_node, "default", 1, 1, 1, 1)

    combos = nbcc_2020_combos()
    for combo_name, load_combo in combos.items():
        model.add_load_combo(combo_name, factors=load_combo)
    return model


def element_to_joist_model(element: Element, w: float = 0.) -> Joist:
    """
    Returns a Joist object based on the data in 'element'
    """
    try:
        r1, r2 = element.intersections
    except ValueError:
        raise ValueError(f"Joists currently need to have two supports. {element.tag=} | {element.intersections=}")
    elem_geometry = element.geometry
    i_end = Point(elem_geometry.coords[0])
    j_end = Point(elem_geometry.coords[1])
    length = elem_geometry.length
    
    # R1 should be closest to I-end
    r1_geom = r1[1]
    r2_geom = r2[1]
    if i_end.distance(r1_geom) > i_end.distance(r2_geom):
        r1, r2 = r2, r1
    
    span = r1_geom.distance(r2_geom)
    a_cantilever = round(abs(i_end.distance(r1_geom)), 6)
    b_cantilever = round(abs(r2_geom.distance(j_end)), 6)
    return Joist(span, a_cantilever, b_cantilever)


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


            

                    



