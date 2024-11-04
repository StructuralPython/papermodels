from dataclasses import dataclass
from enum import IntEnum
import pathlib
from typing import Optional, Any
from papermodels.datatypes.element import Element

from PyNite import FEModel3D


@dataclass
class AnalysisModel:
    structured_element_data: dict
    analysis_model: Optional[Any] = None
    analyzed: bool = False
    reactions: Optional[dict] = None

    def create_model(self):
        raise NotImplemented

    def add_load(self):
        raise NotImplemented

    def analyze(self):
        raise NotImplemented

    def get_reactions(self):
        raise NotImplemented

    def get_forces(self, force_type: str):
        raise NotImplemented

    def get_reaction_type(self):
        raise NotImplemented


@dataclass
class PyNiteFEModel:
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
                beam_model.def_support(
                    node_name, False, True, False, False, False, False
                )
            elif support_type == "F":
                beam_model.def_support(node_name, True, True, True, True, True, True)

        shear_modulus = calc_shear_modulus(element_data["E"], element_data["nu"])
        beam_model.add_material(
            element_data["material_name"],
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
            material=element_data["material_name"],
            Iy=element_data["Iy"],
            Iz=element_data["Iz"],
            J=element_data["J"],
            A=element_data["A"],
        )

        self.analysis_model = beam_model

    def add_load(self, load_data: dict) -> None:
        """
        Adds one load to the model in self.analysis_model.
        Sets self.analyzed to False to force a re-analysis.
        """
        load_cases = list(self.analysis_model.LoadCombos.keys())
        if load_data["Type"] == "Point":
            self.analysis_model.add_member_pt_load(
                self.structured_element_data["Name"],
                load_data["Direction"],
                load_data["Magnitude"],
                load_data["Location"],
                case=load_data["Case"],
            )
            if load_data["Case"] not in load_cases:
                load_cases.append(load_data["Case"])

        elif load_data["Type"] == "Dist":
            self.analysis_model.add_member_dist_load(
                self.structured_element_data["Name"],
                load_data["Direction"],
                load_data["Start Magnitude"],
                load_data["End Magnitude"],
                load_data["Start Location"],
                load_data["End Location"],
                case=load_data["Case"],
            )
            if load_data["Case"] not in load_cases:
                load_cases.append(load_data["Case"])
        self.analyzed = False
        self.reactions = None

    def _add_loads(self) -> None:
        element_data = self.structured_element_data
        load_cases = []
        for load in element_data["Loads"]:
            if load["Type"] == "Point":
                self.analysis_model.add_member_pt_load(
                    element_data["Name"],
                    load["Direction"],
                    load["Magnitude"],
                    load["Location"],
                    case=load["Case"],
                )
                if load["Case"] not in load_cases:
                    load_cases.append(load["Case"])
            elif load["Type"] == "Dist":
                self.analysis_model.add_member_dist_load(
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

        # Hack for lack of implementation of load case reporting in PyNite
        for load_case in load_cases:
            self.analysis_model.add_load_combo(load_case, {load_case: 1.0})

    def analyze(
        self,
        analyze_linear: bool = True,
        check_stability: bool = True,
        check_statics=False,
    ):
        if analyze_linear:
            self.analysis_model.analyze_linear(
                check_stability=check_stability, check_statics=check_statics
            )
        else:
            self.analysis_model.analyze(
                check_stability=check_stability, check_statics=check_statics
            )
        self.analyzed = True
        self.reactions = self.get_reactions()

    def get_reactions(self):
        if not self.analyzed:
            raise UserWarning(
                f"{self.analyzed=}. Re-run .analyze() to ensure results are current."
            )
        model_nodes = self.analysis_model.Nodes

        reactions = {}
        for node_name, node_obj in model_nodes:
            reactions[node_name] = {}
            for reaction_dir in ("Fx", "Fy", "Fz", "Mx", "My", "Mz"):
                reaction_key = f"Rxn{reaction_dir.upper()}"
                reaction_combos = getattr(node_obj, reaction_key)
                reactions[reaction_dir] = reaction_combos
        return reactions

    def get_forces(self, n_points: int = 200):
        if not self.analyzed:
            raise UserWarning(
                f"{self.analyzed=}. Re-run .analyze() to ensure results are current."
            )
        member_name = self.structured_element_data["Name"]
        member_obj = self.analysis_model.Members[member_name]
        load_cases = list(self.analysis_model.LoadCombos.keys())

        forces = {}
        force_actions = ("N", "Vz", "Vy", "Mz", "My", "T")
        for force_action in force_action:
            if force_action == "N":
                forces[force_action] = member_obj.axial_array(
                    n_points=n_points,
                )


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
