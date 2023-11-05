import csv
from typing import Optional
import math
from PyNite import FEModel3D
from papermodels.models.utils import str_to_int, str_to_float, read_csv_file
from papermodels.models.load_factors import nbcc_2020_combos, nbcc_2015_combos, nbcc_2010_combos
    

def read_beam_file(filename: str) -> str:
    """
    Returns data contained in the file, 'filename' as a list of lists.
    It is assumed that the data in the file is "csv-ish", meaning with 
    comma-separated values.
    """
    return read_csv_file(filename)


def separate_lines(file_data: str) -> list[str]:
    """
    Returns the data in 'file_data' separated out into separate lines
    as a list.
    """
    return file_data.split("\n")


def separate_data(file_data: list[str]) -> list[list[str]]:
    """
    Returns the data in 'file_data' separated out into separate lines
    as a list.
    """
    acc = []
    for line in file_data:
        split_line = line.split(", ")
        acc.append(split_line)
    return acc


def parse_beam_attributes(beam_attributes: list[float]) -> dict[str, float]:
    """
    Returns a dictionary of beam attributes according to which attributes are present
    in 'beam_attributes' in the order according to the beam file format BEAM_FORMAT.md,
    Workbook_04 edition. The order of attributes are as follows: Length,E,Iz,[Iy,A,J,nu,rho]
    """
    ATTR_ORDER = ["L", "E", "Iz", "Iy", "A", "J", "nu", "rho"]

    #Implementation 1
    parsed = {}
    for idx, attr in enumerate(ATTR_ORDER):
        try:
            attr_present = beam_attributes[idx]
            parsed.update({attr: attr_present})
        except IndexError:
            parsed.update({attr: 1.0})
    return parsed

    # Implementation 2 - Fancy!
    # defaults = [1.0] * len(ATTR_ORDER)
    # default_attrs = dict(zip(ATTR_ORDER, defaults))
    # attrs_present = dict(zip(ATTR_ORDER, ex2))
    # return default_attrs | attrs_present


def parse_supports(supports: list[str]) -> dict[float, str]:
    """
    Returns a dictionary representing the data in 'supports' separated
    out into a dictionary with support locations as keys and support types
    as values. 

    Assumes 'supports' is in a format that looks like this:
    ['support_loc:support_type', 'support_loc:support_type', etc...]
    e.g. ['1000:P', '3800:R', '6000:R']
    Where the valid support types are one of: P (pinned), F (fixed), or R (roller)
    """
    parsed = {}
    for support in supports:
        sup_loc, sup_type = support.split(":")
        parsed.update({str_to_float(sup_loc): sup_type})
    return parsed


def convert_to_numeric(raw_data: list[list[str]]) -> list[list[float]]:
    """
    Returns a nested list of floats representing all of the numeric string data in
    'raw_data' being converted into a floats. 

    If the data cannot be converted, a ValueError will be raised.
    """
    outer_acc = []
    for line in raw_data:
        inner_acc = []
        for element in line:
            inner_acc.append(str_to_float(element))
        outer_acc.append(inner_acc)
    return outer_acc


def parse_loads(loads_data: list[str|float]) -> list[dict]:
    """
    Returns a the loads in 'loads_data' structured into a list of dicts.
    This implementation is just one implementation among many. I do not 
    make claims about it being "the best" implementation :)
    """
    parsed = []
    for load in loads_data:
        load_type, load_dir = load[0].split(":")
        load_case = load[-1].split(":")[-1]
        if load_type == "POINT":
            magnitude = load[1]
            location = load[2]
            parsed.append(
                {
                    "Type": load_type.title(),
                    "Direction": load_dir.title(),
                    "Magnitude": magnitude,
                    "Location": location,
                    "Case": load_case
                }
            )
                
        elif load_type == "DIST":
            start_magnitude = load[1]
            end_magnitude = load[2]
            start_loc = load[3]
            end_loc = load[4]
            parsed.append(
                {
                    "Type": load_type.title(),
                    "Direction": load_dir.title(),
                    "Start Magnitude": start_magnitude,
                    "End Magnitude": end_magnitude,
                    "Start Location": start_loc,
                    "End Location": end_loc,
                    "Case": load_case
                }
            )
    return parsed


def get_structured_beam_data(raw_data: list[list[str]]) -> dict:
    """
    Returns a dictionary that has string keys describing the attributes of a beam for analysis.
    """
    numeric_beam_data = convert_to_numeric(raw_data)
    beam_name = raw_data[0][0]
    beam_attributes = parse_beam_attributes(numeric_beam_data[1])
    supports = numeric_beam_data[2]
    loads = numeric_beam_data[3:]
    structured_data = {}
    structured_data['Name'] = beam_name
    structured_data.update(beam_attributes)
    structured_data['Supports'] = parse_supports(supports)
    structured_data['Loads'] = parse_loads(loads)
    return structured_data


def get_node_locations(support_locations: list[float], beam_length: float) -> dict[str, float]:
    """
    Returns a dict representing the node names and node coordinates required for
    a beam with support locations at 'support_locations' and a beam of length of 
    'beam_length'.

    Each node name and node location is unique. Nodes are numbered sequentially from
    left to right starting at "N0".
    """
    nodes_to_create = support_locations[:] # make a copy
    if 0.0 not in support_locations:
        nodes_to_create.append(0.0)
    if beam_length not in support_locations:
        nodes_to_create.append(beam_length)
        
    node_locations = {}
    for idx, sup_loc in enumerate(sorted(nodes_to_create)):
        node_locations.update({f"N{idx}": sup_loc})
    return node_locations


def extract_data(data_line: list[str], index: int) -> list[str]:
    """
    Returns the data in 'data_line' located at 'index'.
    """
    return data_line[index].split(", ")


def get_spans(beam_length: float, cant_support_loc: float) -> tuple[float, float]:
    """
    Returns the length of the backspan ("b") and the length of the cantilever ("a")
    given a total beam length and the location of the cantilever support.

    This function assumes that the backspan support is located at 0.0
    """
    b = cant_support_loc
    a = beam_length - cant_support_loc
    return b, a


def build_beam(beam_data: dict) -> FEModel3D:
    """
    Returns a FEModel3D of a beam that has the attributes described in 'beam_data' and 'A', 'J', 'nu', and 'rho'.
    """
    beam_data['Nodes'] = get_node_locations(list(beam_data['Supports'].keys()), beam_data['L'])
    beam_model = FEModel3D()
    for node_name, x_coord in beam_data['Nodes'].items():
        beam_model.add_node(node_name, x_coord, 0, 0)
        support_type = beam_data['Supports'].get(x_coord, None)
        if support_type == "P":
            beam_model.def_support(node_name, True, True, True, True, True, False)
        elif support_type == "R":
            beam_model.def_support(node_name, False, True, False, False, False, False)
        elif support_type == "F":
            beam_model.def_support(node_name, True, True, True, True, True, True)

    shear_modulus = calc_shear_modulus(beam_data['E'], beam_data['nu'])
    beam_model.add_material("Beam Material", beam_data['E'], shear_modulus, beam_data['nu'], beam_data['rho'])
    
    # The variable 'node_name' retains the last value from the for loop
    # Which gives us the j-node
    beam_model.add_member(
        beam_data['Name'], 
        "N0", 
        node_name, 
        material="Beam Material", 
        Iy=beam_data['Iy'], 
        Iz=beam_data['Iz'],
        J=beam_data['J'],
        A=beam_data['A'],
    )
    load_cases = []
    for load in beam_data['Loads']:
        if load['Type'] == "Point":
            beam_model.add_member_pt_load(
                beam_data['Name'],
                load['Direction'],
                load['Magnitude'],
                load['Location'],
                case=load["Case"],
            )
            if load['Case'] not in load_cases:
                load_cases.append(load['Case'])
        elif load['Type'] == "Dist":
            beam_model.add_member_dist_load(
                beam_data['Name'],
                load['Direction'],
                load['Start Magnitude'],
                load['End Magnitude'],
                load['Start Location'],
                load['End Location'],
                case=load['Case']
            )
            if load['Case'] not in load_cases:
                load_cases.append(load['Case'])

    for load_case in load_cases:
        beam_model.add_load_combo(load_case, {load_case: 1.0})
    return beam_model


def load_beam_model(filename: str, add_combos: Optional[str] = None) -> FEModel3D: 
    """
    Returns an FEModel3D beam model representing the beam described in 'filename'

    Current combo collections: 'nbcc_2020', 'nbcc_2015', 'nbcc_2010'
    """
    beam_data_raw = read_beam_file(filename)
    beam_data_structured = get_structured_beam_data(beam_data_raw)
    beam_model = build_beam(beam_data_structured)
    if add_combos is not None:
        if add_combos.lower() == "nbcc_2020":
            load_combos = nbcc_2020_combos()
        elif add_combos.lower() == "nbcc_2015":
            load_combos = nbcc_2015_combos()
        elif add_combos.lower() == "nbcc_2010":
            load_combos = nbcc_2010_combos()
        else:
            raise ValueError(f"Add combos must be one of 'nbcc_2020', 'nbcc_2015', or 'nbcc_2010', not {add_combos}")
        for combo_name, combo_factors in load_combos.items():
            beam_model.add_load_combo(combo_name, combo_factors)
    return beam_model


def extract_arrays_all_combos(solved_beam_model: FEModel3D, result_type: str, direction: Optional[str] = None, n_points: int=200) -> dict:
    """
    Returns a dictionary keyed by load combo name that contains the resulting arrays of the 'solved_beam_model',
    for the given 'result_type' and 'direction' with 'n_points' as the number of values in the array.

    'solved_beam_model': A PyNite.FEModel3D object that contains one member and has been successfully analyzed
    'result_type': str, one of {'shear', 'moment', 'deflection', 'axial', 'torque'}
    'direction': str that corresponds to the 'result_type':
        'shear': {'Fx', 'Fy', 'Fz'}
        'moment': {'Mx', 'My', 'Mz'}
        'deflection': {'dx', 'dy', 'dz'}
    'n_points': the number of values in the resulting arrays

    The keys in the resulting dictionary represent the names of all of the load combos in the model. The
    values are (n_points, 2)-shaped arrays that contain an x-array (of beam locations) and a y-array (of results).
    """
    all_combos = {}
    beam_name = list(solved_beam_model.Members.keys())[0]
    member = solved_beam_model.Members[beam_name]
    if result_type == "shear":
        method = member.shear_array
    elif result_type == "moment":
        method = member.moment_array
    elif result_type == "deflection":
        method = member.deflection_array
    elif result_type == "axial":
        method = member.axial_array
    elif result_type == "torque":
        method = member.torque_array
    
    for combo_name in solved_beam_model.LoadCombos:
        if result_type in ('shear', 'moment', 'deflection'):
            all_combos[combo_name] = method(direction, n_points, combo_name)
        else:
            all_combos[combo_name] = method(n_points, combo_name)
    return all_combos


def calc_shear_modulus(E: float, nu: float) -> float:
    """
    Returns the shear modulus based on Poisson's ratio and
    the elastic modulus

    'E': Elastic modulus
    'nu': Poisson's ratio
    """
    G = E / (2 * (1 + nu))
    return G


def euler_buckling_load(
    l: float,
    E: float, 
    I: float,
    k: float
) -> float:
    """
    Returns the critical Euler buckling load for an axially
    loaded member.

    'l': length
    'E': Elastic modulus
    'I': Moment of inertia
    'k': Effective length factor

    This function assumes consistent units are used for each parameter.
    """
    P_cr = math.pi**2 * E * I / (k * l)**2
    return P_cr


def beam_reactions_ss_cant(
    w: float,
    b: float,
    a: float,
) -> float:
    """
    Returns the reactions "R1" and "R2" for a simply-supported beam with 
    a cantilever (over-hang) on one end. "R1" is the cantilever reaction and
    "R2" is the backspan reaction.

    'w': the magnitude of the distributed load (+ve is in the gravity direction,
        -ve is in the uplift direction)
    'b': the length of the back-span
    'a': the length of the cantilever

    This function assumes consitent units are used for each parameter.
    """
    r2 = w / (2*b) * (b**2 - a**2)
    r1 = w / (2*b) * (b + a)**2
    return r1, r2


def fe_model_ss_cant(
    w: float,
    b: float,
    a: float,
    E: float = 1.,
    I: float = 1.,
    A: float = 1.,
    J: float = 1.,
    nu: float = 1.,
    rho: float = 1.,
) -> FEModel3D:
    """
    Returns a PyNite FEModel3D for a simply-supported beam with a cantilever
    (over-hang) on one end.

    'w': the magnitude of the distributed load (+ve is in the gravity direction,
        -ve is in the uplift direction)
    'b': the length of the back-span
    'a': the length of the cantilever
    'E': elastic modulus of beam material
    'I': moment of inertia (in same direction was 'w') of beam member
    'A': area of beam member
    'J': polar moment of inertia of beam member
    'nu': poisson's ratio of beam material
    'rho': Specific gravity of beam material
    """
    G = calc_shear_modulus(nu, E)
    beam_model = FEModel3D()
    beam_model.add_material("beam_material", E, G, nu, rho)
    beam_model.add_node("N0", 0, 0, 0)
    beam_model.add_node("N1", b, 0, 0)
    beam_model.add_node("N2", b + a, 0, 0)
    
    beam_model.add_member("M0", "N0", "N2", "beam_material", Iy=1, Iz=I, J=J, A=A)

    beam_model.def_support(
        "N0",
        support_DX=1,
        support_DY=1,
        support_DZ=1,
        support_RX=1,
        support_RY=1,
        support_RZ=0
    )
    beam_model.def_support(
        "N1",
        support_DX=0,
        support_DY=1,
        support_DZ=0,
        support_RX=0,
        support_RY=0,
        support_RZ=0
    )

    beam_model.add_member_dist_load(
        "M0", "Fy", w, w, 0, case="Case 1"
    )
    return beam_model


def get_nodes_location_fixity(beam_model: FEModel3D) -> dict[str, tuple[str, float]]:
    """
    Returns the node location and fixity for all nodes in 'beam_model'
    """
    acc = {}
    for node_name, node_data in beam_model.Nodes.items():
        x_coord = node_data.X
        all_supports = [
            node_data.support_DX,
            node_data.support_DY,
            node_data.support_DZ,
            node_data.support_RX,
            node_data.support_RY,
            node_data.support_RZ
        ]
        if sum(all_supports) == 6:
            acc.update({node_name, (x_coord, "F")})
        elif sum(all_supports) >= 2:
            acc.update({node_name, (x_coord, "P")})
        elif sum(all_supports) == 1:
            acc.update({node_name, (x_coord, "R")})
        elif sum(all_supports) == 0:
            acc.update({node_name, (x_coord, None)})
    return acc


def get_all_node_reactions(beam_model: FEModel3D) -> dict[str, dict[str, dict[str, float]]]:
    """
    Returns a nested dictionary of node reactions for all load combos keyed
    by node name, direction, and then load combo name.
    """
    acc = {}
    for node_name, node_data in beam_model.Nodes.items():
        acc[node_name] = {}
        for attr in ('RxnFX', 'RxnFY', 'RxnFZ', 'RxnMX', 'RxnMY', 'RxnMZ'):
            acc[node_name][attr.replace("Rxn","")] = getattr(node_data, attr)
    return acc


def get_max_reaction(all_node_reactions: dict, node_name: str, direction: str) -> tuple[float, str]:
    """
    Returns the maximum reaction and combo name for the reactions in 'all_node_reactions'
    for the given 'node_name' and 'direction'.
    """
    reactions = all_node_reactions[node_name][direction]
    max_reaction = max(reactions.values())
    combos_by_reactions = {v: k for k, v in reactions.items()}
    return max_reaction, combos_by_reactions[max_reaction]


def get_min_reaction(all_node_reactions: dict, node_name: str, direction: str) -> tuple[float, str]:
    """
    Returns the minimum reaction and combo name for the reactions in 'all_node_reactions'
    for the given 'node_name' and 'direction'.
    """
    reactions = all_node_reactions[node_name][direction]
    min_reaction = min(reactions.values())
    combos_by_reactions = {v: k for k, v in reactions.items()}
    return min_reaction, combos_by_reactions[min_reaction]


def is_member_direction_empty(beam_model: FEModel3D, direction: str) -> bool:
    """
    Returns True if there are no loads or reactions present on 'beam_model'
    for the given direction.
    """
    all_node_reactions = get_all_node_reactions(beam_model)
    
    reaction_indicators = []
    for node_name, reaction_dict in all_node_reactions.items():
        reactions = reaction_dict[direction.upper()]
        sum_reactions = sum(reactions.values())
        if not math.isclose(sum_reactions, 0):
            reaction_indicators.append(1)
        else:
            reaction_indicators.append(0)
    return not sum(reaction_indicators)

def get_nodes_location_fixity(beam_model: FEModel3D) -> dict[str, tuple[str, float]]:
    """
    Returns the node location and fixity for all nodes in 'beam_model'
    """
    acc = {}
    for node_name, node_data in beam_model.Nodes.items():
        x_coord = node_data.X
        all_supports = [
            node_data.support_DX,
            node_data.support_DY,
            node_data.support_DZ,
            node_data.support_RX,
            node_data.support_RY,
            node_data.support_RZ
        ]
        if sum(all_supports) == 6:
            acc.update({x_coord, "F"})
        elif sum(all_supports) >= 2:
            acc.update({x_coord, "P"})
        elif sum(all_supports) == 1:
            acc.update({x_coord, "R"})
        elif sum(all_supports) == 0:
            acc.update({x_coord, None})
        
    return acc