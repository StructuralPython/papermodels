from __future__ import annotations
from typing import Optional
import math
from copy import deepcopy
import networkx as nx
import hashlib
import pathlib

from shapely import Point
from PyNite import FEModel3D

from papermodels.datatypes.element import (
    Element,
    get_tag_type,
    get_normalized_coordinate,
)
from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.element_models import BeamModel
from papermodels.datatypes.element_models import ColumnModel
from papermodels.datatypes.element_models import WallModel
from papermodels.datatypes.joist_models import JoistArrayModel
from rich.progress import track
from rich import print


class LoadGraph(nx.DiGraph):
    """
    A class to represent a load path in a graph. Inherits from networkx.DiGraph
    and adds a .node_hash attribute for storing a hash of all the nodes.

    The node_hash is how changes to the graph nodes can be tracked.
    """

    def __init__(self):
        super().__init__()
        self.geometry_hash = None

    @classmethod
    def from_geometry_graph(cls, graph: GeometryGraph):
        g = cls()
        g.geometry_hash = graph.node_hash
        for node, node_data in graph.nodes.items():
            element = node_data["element"]
            if element.type.lower() == "joist":
                model = JoistArrayModel.from_element(element)
            elif "beam" in element.type.lower():
                model = BeamModel.from_element(element)
            elif "column" in element.type.lower():
                model = ColumnModel.from_element(element)
            elif "wall" in element.type.lower():
                model = WallModel.from_element(element)
            else:
                print(element.type.lower())
            g.add_node(node, model=model)

        for u, v in graph.edges:
            edge_data = graph.edges[(u, v)]
            g.add_edge(u, v, edge_data=edge_data)
        return g

    @classmethod
    def from_beam_files(cls, beam_file_dir: pathlib.Path):
        g = cls()
        g.geometry_hash = graph.node_hash
        for node, node_data in graph.nodes.items():
            element = node_data["element"]
            if element.type.lower() == "joist":
                model = JoistArrayModel(element)
            elif "beam" in element.type.lower():
                model = BeamModel(element)
            elif "column" in element.type.lower():
                model = ColumnModel(element)
            elif "wall" in element.type.lower():
                model = WallModel(element)
            else:
                print(element.type.lower())
            g.add_node(node, model=model)

        for u, v in graph.edges:
            edge_data = graph.edges[(u, v)]
            g.add_edge(u, v, edge_data=edge_data)
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
                model = element_to_joist_array(element)
                self.nodes[node]["model"] = model
            elif prefix_dict[tag_prefix] == "column":
                pass
            elif prefix_dict[tag_prefix] == "wall":
                pass
            elif prefix_dict[tag_prefix] == "point_load":
                pass
            elif prefix_dict[tag_prefix] == "line_load":
                pass


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
    model.add_load_combo("Pass", {"pass": 1.0})
    return model


# def element_to_joist_model(element: Element, w: float = 0.) -> Joist:
#     """
#     Returns a Joist object based on the data in 'element'
#     """
#     try:
#         r1, r2 = element.intersections
#     except ValueError:
#         raise ValueError(f"Joists currently need to have two supports. {element.tag=} | {element.intersections=}")
#     elem_geometry = element.geometry
#     i_end = Point(elem_geometry.coords[0])
#     j_end = Point(elem_geometry.coords[1])
#     length = elem_geometry.length

#     # R1 should be closest to I-end
#     r1_geom = r1[1]
#     r2_geom = r2[1]
#     if i_end.distance(r1_geom) > i_end.distance(r2_geom):
#         r1, r2 = r2, r1

#     span = r1_geom.distance(r2_geom)
#     a_cantilever = round(abs(i_end.distance(r1_geom)), 6)
#     b_cantilever = round(abs(r2_geom.distance(j_end)), 6)
#     return Joist(span, a_cantilever, b_cantilever)


def element_to_joist_array(
    joist_element: Element,
    initial_offset: int | float = 0,
    joist_at_start: bool = True,
    joist_at_end: bool = False,
    cantilever_tolerance: float = 0.01,
) -> JoistArrayModel:
    """
    Returns a Joist object based on the data in 'element'
    """
    supports = [inter[2] for inter in joist_element.intersections]
    joist_array = JoistArrayModel(
        joist_element.tag,
        joist_element.geometry,
        supports,
        joist_element.spacing,
        initial_offset,
        joist_at_start,
        joist_at_end,
        cantilever_tolerance,
    )
    # try:
    #     r1, r2 = element.intersections
    # except ValueError:
    #     raise ValueError(f"Joists currently need to have two supports. {element.tag=} | {element.intersections=}")
    # elem_geometry = element.geometry
    # i_end = Point(elem_geometry.coords[0])
    # j_end = Point(elem_geometry.coords[1])
    # length = elem_geometry.length

    # # R1 should be closest to I-end
    # r1_geom = r1[1]
    # r2_geom = r2[1]
    # if i_end.distance(r1_geom) > i_end.distance(r2_geom):
    #     r1, r2 = r2, r1

    # span = r1_geom.distance(r2_geom)
    # a_cantilever = round(abs(i_end.distance(r1_geom)), 6)
    # b_cantilever = round(abs(r2_geom.distance(j_end)), 6)
    # return Joist(span, a_cantilever, b_cantilever)
