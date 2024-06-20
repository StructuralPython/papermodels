import csv
import math
from typing import Optional

from PyNite import FEModel3D
from papermodels.fileio.utils import str_to_float, str_to_int


import csv
from typing import Optional
import math
from PyNite import FEModel3D
from papermodels.fileio.utils import (
    str_to_int,
    str_to_float,
    read_csv_file,
    convert_unit_string,
    parse_unit_system,
)

DEFAULT_ATTRS = {
    "L": 1.0,
    "E": 1.0,
    "Iz": 1.0,
    "Iy": 1.0,
    "A": 1.0,
    "J": 1.0,
    "nu": 1.0,
    "rho": 1.0,
}


def get_structured_beam_data(raw_data: list[list[str]]) -> dict:
    """
    Returns a dictionary that has string keys describing the attributes of a beam for analysis.
    """
    unit_designation = None
    if "units" in raw_data[-1][0].lower():
        unit_designation = raw_data.pop()[0].split(":")[-1].replace(" ", "")
    beam_attributes = parse_beam_attributes(raw_data[1], unit_designation)
    numeric_beam_data = convert_to_numeric(raw_data[2:])
    if unit_designation is not None:
        numeric_beam_data = convert_to_units(numeric_beam_data, unit_designation)
    beam_name = raw_data[0][0]
    supports = numeric_beam_data[0]
    loads = numeric_beam_data[1:]
    structured_data = {}
    structured_data["Name"] = beam_name
    structured_data.update(beam_attributes)
    structured_data["Supports"] = parse_supports(supports, unit_designation)
    structured_data["Loads"] = parse_loads(loads)
    return structured_data


def parse_beam_attributes(
    beam_attributes: list[str], unit_designation: Optional[str] = None
) -> dict[str, float]:
    """
    Returns a dictionary of beam attributes according to which attributes are present
    in 'beam_attributes' in the order according to the beam file format BEAM_FORMAT.md,
    Workbook_04 edition. The order of attributes are as follows: Length,E,Iz,[Iy,A,J,nu,rho]
    """
    beam_attrs = map_beam_attributes(parse_kwargs(beam_attributes, unit_designation))
    return DEFAULT_ATTRS | beam_attrs


def parse_kwargs(
    attr_items: list[str], unit_designation: Optional[str] = None
) -> dict[str, str | float]:
    """
    Returns a dictionary of any of the elements in 'attr_items' that
    contain an "=" character as key, value pairs.
    """
    acc = {}
    for attr_item in attr_items:
        if "=" in attr_item:
            attr, item = attr_item.split("=", 1)
            attr = attr.rstrip(" ")
            item = item.strip(" ")
            as_float = str_to_float(item)
            if unit_designation is not None:
                unit_system = parse_unit_system(unit_designation)
                as_float = convert_unit_string(as_float, unit_system)
            acc.update({attr: as_float})
    return acc


def map_beam_attributes(parsed_attrs: dict[str]) -> dict[str]:
    """
    Returns a list of str data that represents the beam attributes in the
    original attributes line with any kwargs remove and placed in the expected
    positional order.
    """
    attr_map = {
        "l": "L",
        "e": "E",
        "Ix": "Iz",
        "ix": "Iz",
        "iz": "Iz",
        "iy": "Iy",
        "a": "A",
        "j": "J",
    }
    acc = {}
    for attr, value in parsed_attrs.items():
        new_attr = attr_map.get(attr, attr)
        acc.update({new_attr: value})
    return acc


def parse_supports(
    supports: list[str], unit_designation: Optional[str] = None
) -> dict[float, str]:
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
        as_float = str_to_float(sup_loc)
        if unit_designation is not None:
            unit_system = parse_unit_system(unit_designation)
            as_float = convert_unit_string(as_float, unit_system)
        parsed.update({str_to_float(as_float): sup_type})
    return parsed


def convert_to_units(
    raw_data: list[list[str]], unit_designation: str
) -> list[list[float]]:
    """
    Returns a nested list of floats representing all of the numeric string data in
    'raw_data' being converted into a floats.

    If the data cannot be converted, a ValueError will be raised.
    """
    unit_system = parse_unit_system(unit_designation)
    outer_acc = []
    for line in raw_data:
        inner_acc = []
        for element in line:
            inner_acc.append(convert_unit_string(element, unit_system))
        outer_acc.append(inner_acc)
    return outer_acc


def parse_loads(loads_data: list[str | float]) -> list[dict]:
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
                    "Case": load_case,
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
                    "Case": load_case,
                }
            )

        elif load_type == "TRIBW":
            magnitude = load[1]
            width = load[2]
            start_loc = load[3]
            end_loc = load[4]
            udl_start, udl_end = calculate_line_load_from_trib_width(magnitude, width)
            parsed.append(
                {
                    "Type": "Dist",
                    "Direction": load_dir.title(),
                    "Start Magnitude": udl_start,
                    "End Magnitude": udl_end,
                    "Start Location": start_loc,
                    "End Location": end_loc,
                    "Case": load_case,
                }
            )

        elif load_type == "TRIBVW":
            magnitude = load[1]
            width_start = load[2]
            width_end = load[3]
            start_loc = load[4]
            end_loc = load[5]
            line_start, line_end = calculate_line_load_from_trib_width(
                magnitude, width_start, width_end
            )
            parsed.append(
                {
                    "Type": "Dist",
                    "Direction": load_dir.title(),
                    "Start Magnitude": line_start,
                    "End Magnitude": line_end,
                    "Start Location": start_loc,
                    "End Location": end_loc,
                    "Case": load_case,
                }
            )
    return parsed


def calculate_line_load_from_trib_width(
    area_load: float, trib_width_start: float, trib_width_end: Optional[float] = None
) -> tuple[float, float]:
    """
    Returns the equivalent line load for an 'area_load'for the given
    trib_width. Units of quantities provided must be consistent for correct results.

    'area_load' - a quantity in force/area
    'trib_width_start' - a quantity of length in the same units as area. If only 'trib_width_start'
        is provided, then the line load is assumed to be uniform.
    'trib_width_end' - a quantity of length in the same units as area.
    """
    if trib_width_end is not None:
        return area_load * trib_width_start, area_load * trib_width_end
    return (area_load * trib_width_start,) * 2


def remove_comments(raw_data: list[list[str]]) -> list[list[str]]:
    """
    Returns 'raw_data' but with any comments removed. A comment is represented
    the same as in Python with a "#" followed by arbitrary characters.
    """
    outer_acc = []
    for line in raw_data:
        inner_acc = []
        for elem in line:
            if "#" in elem:
                not_comment, comment = elem.replace(" ", "").split("#")
                inner_acc.append(not_comment)
            else:
                inner_acc.append(elem)
        outer_acc.append(inner_acc)
    return outer_acc


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
            node_data.support_RZ,
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


def get_all_node_reactions(
    beam_model: FEModel3D,
) -> dict[str, dict[str, dict[str, float]]]:
    """
    Returns a nested dictionary of node reactions for all load combos keyed
    by node name, direction, and then load combo name.
    """
    acc = {}
    for node_name, node_data in beam_model.Nodes.items():
        acc[node_name] = {}
        for attr in ("RxnFX", "RxnFY", "RxnFZ", "RxnMX", "RxnMY", "RxnMZ"):
            acc[node_name][attr.replace("Rxn", "")] = getattr(node_data, attr)
    return acc


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
