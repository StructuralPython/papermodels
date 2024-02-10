from dataclasses import dataclass
from enum import IntEnum
import pathlib
from typing import Optional, Any
from papermodels.datatypes.element import Element

from PyNite import FEModel3D



@dataclass
class PyNiteFEModel:
    structured_element_data: dict

    def create_model(self):
        """
        Returns a FEModel3D of a beam that has the attributes described in 'element_data' and 'A', 'J', 'nu', and 'rho'.
        """
        element_data = self.structured_element_data
        element_data["Nodes"] = get_node_locations(
            list(element_data["Supports"].keys()), element_data["L"]
        )
        beam_model = FEModel3D()
        for node_name, x_coord in element_data["Nodes"].items():
            beam_model.add_node(node_name, x_coord, 0, 0)
            support_type = element_data["Supports"].get(x_coord, None)
            if support_type == "P":
                beam_model.def_support(node_name, True, True, True, True, True, False)
            elif support_type == "R":
                beam_model.def_support(node_name, False, True, False, False, False, False)
            elif support_type == "F":
                beam_model.def_support(node_name, True, True, True, True, True, True)

        shear_modulus = calc_shear_modulus(element_data["E"], element_data["nu"])
        beam_model.add_material(
            "Beam Material",
            element_data["E"],
            shear_modulus,
            element_data["nu"],
            element_data["rho"],
        )

        # The variable 'node_name' retains the last value from the for loop
        # Which gives us the j-node
        beam_model.add_member(
            element_data["Name"],
            "N0",
            node_name,
            material="Beam Material",
            Iy=element_data["Iy"],
            Iz=element_data["Iz"],
            J=element_data["J"],
            A=element_data["A"],
        )
        load_cases = []
        for load in element_data["Loads"]:
            if load["Type"] == "Point":
                beam_model.add_member_pt_load(
                    element_data["Name"],
                    load["Direction"],
                    load["Magnitude"],
                    load["Location"],
                    case=load["Case"],
                )
                if load["Case"] not in load_cases:
                    load_cases.append(load["Case"])
            elif load["Type"] == "Dist":
                beam_model.add_member_dist_load(
                    element_data["Name"],
                    load["Direction"],
                    load["Start Magnitude"],
                    load["End Magnitude"],
                    load["Start Location"],
                    load["End Location"],
                    case=load["Case"],
                )
                if load["Case"] not in load_cases:
                    load_cases.append(load["Case"])

        for load_case in load_cases:
            beam_model.add_load_combo(load_case, {load_case: 1.0})
        return beam_model



def get_node_locations(
    support_locations: list[float], beam_length: float
) -> dict[str, float]:
    """
    Returns a dict representing the node names and node coordinates required for
    a beam with support locations at 'support_locations' and a beam of length of
    'beam_length'.

    Each node name and node location is unique. Nodes are numbered sequentially from
    left to right starting at "N0".
    """
    nodes_to_create = support_locations[:]  # make a copy
    if 0.0 not in support_locations:
        nodes_to_create.append(0.0)
    if beam_length not in support_locations:
        nodes_to_create.append(beam_length)

    node_locations = {}
    for idx, sup_loc in enumerate(sorted(nodes_to_create)):
        node_locations.update({f"N{idx}": sup_loc})
    return node_locations


def calc_shear_modulus(E: float, nu: float) -> float:
    """
    Returns the shear modulus based on Poisson's ratio and
    the elastic modulus

    'E': Elastic modulus
    'nu': Poisson's ratio
    """
    G = E / (2 * (1 + nu))
    return G