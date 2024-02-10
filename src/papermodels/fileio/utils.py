import csv
from enum import Enum
from fractions import Fraction
from typing import Optional
import forallpeople as si
import parse

si.environment("structural")
import forallpeople.tuplevector as vec


class UnitSystem(Enum):
    kPa = 0
    MPa = 1
    GPa = 2
    psf = 3
    ksf = 4
    psi = 5
    ksi = 6


def str_to_int(s: str) -> int | str:
    """
    Returns the string 's' converted into an integer if possible.
    Returns the original string, 's', otherwise.
    """
    try:
        return int(s)
    except ValueError:
        return s


def str_to_float(s: str) -> float | str:
    """
    Returns the string 's' converted into a float if possible.
    Returns the original string, 's', otherwise.
    """
    try:
        return float(s)
    except ValueError:
        return s


def read_csv_file(filename: str) -> str:
    """
    Returns data contained in the file, 'filename' as a list of lists.
    It is assumed that the data in the file is "csv-ish", meaning with
    comma-separated values.
    """
    with open(filename, "r") as file:
        csv_data = list(csv.reader(file))
    return csv_data


def parse_unit_system(unit_designation: str) -> UnitSystem:
    """
    Returns a UnitSystem enum based on the provided 'unit_designation'.

    'unit_designation': one of {'kPa', 'MPa', 'GPa', 'psf', 'ksf', 'psi', 'ksi'}

    The unit designation is considered enough to describe units of both force and
    distance in a single three-character string. Their meanings are as follows:

    kPa = kN / m2
    MPa = N / mm2
    GPa = kN / mm2
    psf = lbs / ft2
    psi = lbs / inch2
    ksf = kips / ft2
    ksi = kips / inch2
    """
    if unit_designation.lower() == "kpa":
        return UnitSystem.kPa
    elif unit_designation.lower() == "mpa":
        return UnitSystem.MPa
    elif unit_designation.lower() == "gpa":
        return UnitSystem.GPa
    elif unit_designation.lower() == "psf":
        return UnitSystem.psf
    elif unit_designation.lower() == "psi":
        return UnitSystem.psi
    elif unit_designation.lower() == "ksf":
        return UnitSystem.ksf
    elif unit_designation.lower() == "ksi":
        return UnitSystem.ksi


def parse_unit_string(unit_string: str | float) -> float | si.Physical:
    """
    Returns a Physical to represent the physical quantity described by 'unit_string'

    'unit_string': a string formatted as "{number} {unit}" where unit must match a variable
    name within the forallpeople.environment dictionary. If the unit is not found, a
    KeyError is raised.
    """
    if isinstance(unit_string, (float, int)):
        return unit_string
    elif isinstance(unit_string, str) and ":" in unit_string:
        return unit_string
    try:
        magnitude_value, unit_value = unit_string.split(" ")
    except ValueError:
        return unit_string

    unit_value, exponent = parse_exponent(unit_value)
    magnitude = str_to_float(magnitude_value)
    unit = getattr(si, unit_value)
    return magnitude * unit**exponent


def parse_exponent(unit_value: str) -> tuple[str, int]:
    """
    Returns a tuple that includes the exponent described on the 'unit_value'.

    e.g. unit_value = "mm4" -> ("mm", 4)
         unit_value = "kN" -> ("kN", 1)
    """
    if unit_value[-1].isnumeric():
        if "-" in unit_value:
            return (unit_value[:-2], -str_to_int(unit_value[-1]))
        else:
            return (unit_value[:-1], str_to_int(unit_value[-1]))
    else:
        return (unit_value, 1)


def parse_arch_notation(dimension_string: str) -> str:
    """
    Returns a string representing the nominal value of 'dimension_string' in inches.
    e.g.
        1'6 -> "18.0 inch"
        0'4 5/8 -> "4.625 inch"
    """
    if isinstance(dimension_string, (float, int)):
        return dimension_string
    if "'" in dimension_string:
        ft_value, inch_value = dimension_string.split("'")
    elif '"' in dimension_string:
        ft_value = "0"
        inch_value = dimension_string.replace('"', "")
    else:
        return dimension_string
    if "/" in inch_value:
        whole_inch, fractional_inch = inch_value.split(" ")
        whole_inch_magnitude = str_to_float(whole_inch)
        fractional_inch_magnitude = float(Fraction(fractional_inch))
        inch_magnitude = whole_inch_magnitude + fractional_inch_magnitude
        ft_magnitude = str_to_float(ft_value)
    else:
        ft_magnitude = str_to_float(ft_value)
        inch_magnitude = str_to_float(inch_value)
    return f"{ft_magnitude * 12 + inch_magnitude} inch"


def convert_unit_string(unit_string: str, unit_system: UnitSystem) -> float:
    """
    Returns a float of the 'unit_string' after it has been normalized into the
    units of 'unit_system'.
    """
    unit_string = parse_arch_notation(unit_string)
    quantity = parse_unit_string(unit_string)
    quantity_type = get_quantity_type(quantity)
    if quantity_type == "length":
        if unit_system in (UnitSystem.GPa, UnitSystem.MPa):
            scaled_quantity = quantity.si().prefix("m")
        elif unit_system == UnitSystem.kPa:
            scaled_quantity = quantity.si().prefix("unity")
        elif unit_system in (UnitSystem.psf, UnitSystem.ksf):
            scaled_quantity = quantity.to("ft")
        else:
            scaled_quantity = quantity.to("inch")
    elif quantity_type == "force":
        if unit_system in (UnitSystem.GPa, UnitSystem.kPa):
            scaled_quantity = quantity.si().prefix("k")
        elif unit_system == UnitSystem.MPa:
            scaled_quantity = quantity.si().prefix("unity")
        elif unit_system in (UnitSystem.psf, UnitSystem.psi):
            scaled_quantity = quantity.to("lb")
        else:
            scaled_quantity = quantity.to("kip")
    elif quantity_type == "pressure":
        if unit_system == UnitSystem.GPa:
            scaled_quantity = quantity.si().prefix("G")
        elif unit_system == UnitSystem.MPa:
            scaled_quantity = quantity.si().prefix("M")
        elif unit_system == UnitSystem.kPa:
            scaled_quantity = quantity.si().prefix("k")
        elif unit_system == UnitSystem.psi:
            scaled_quantity = quantity.to("psi")
        elif unit_system == UnitSystem.psf:
            scaled_quantity = quantity.to("psf")
        elif unit_system == UnitSystem.ksf:
            scaled_quantity = quantity.to("ksf")
        elif unit_system == UnitSystem.ksi:
            scaled_quantity = quantity.to("ksi")
    elif quantity_type == "line":
        if unit_system == UnitSystem.GPa:
            scaled_quantity = quantity.si().prefix("M")
        elif unit_system == UnitSystem.MPa:
            scaled_quantity = quantity.si().prefix("k")
        elif unit_system == UnitSystem.kPa:
            scaled_quantity = quantity.si().prefix("k")
        elif unit_system == UnitSystem.psi:
            scaled_quantity = quantity.to("lb_in")
        elif unit_system == UnitSystem.psf:
            scaled_quantity = quantity.to("lb_ft")
        elif unit_system == UnitSystem.ksf:
            scaled_quantity = quantity.to("kip_ft")
        elif unit_system == UnitSystem.ksi:
            scaled_quantity = quantity.to("kip_in")
    elif quantity_type == "moment":
        if unit_system == UnitSystem.GPa:
            scaled_quantity = quantity.si().prefix("G")
        elif unit_system == UnitSystem.MPa:
            scaled_quantity = quantity.si().prefix("M")
        elif unit_system == UnitSystem.kPa:
            scaled_quantity = quantity.si().prefix("k")
        elif unit_system == UnitSystem.psi:
            scaled_quantity = quantity.to("lbin")
        elif unit_system == UnitSystem.psf:
            scaled_quantity = quantity.to("lbft")
        elif unit_system == UnitSystem.ksf:
            scaled_quantity = quantity.to("kipft")
        elif unit_system == UnitSystem.ksi:
            scaled_quantity = quantity.to("kipin")
    else:
        return quantity
    return float(scaled_quantity)


def get_quantity_type(p: si.Physical) -> Optional[str]:
    """
    Returns one of
    {'length', 'force', 'pressure', 'line', 'torque', None} in correspondence
    with the .dimensions attribute of 'p'.

    Dimensions of 'p' that are powers of the above quantity types (such as
    area or volume) are all considered whatever the "base" type is.

    e.g. 'length' applies to lengths, areas, volumes, moments of area, etc.
    """
    if not isinstance(p, si.Physical):
        return None
    if p.dimensions in [
        (0, 1, 0, 0, 0, 0, 0),
        (0, 2, 0, 0, 0, 0, 0),
        (0, 3, 0, 0, 0, 0, 0),
        (0, 4, 0, 0, 0, 0, 0),
    ]:
        return "length"
    elif p.dimensions == (1, 1, -2, 0, 0, 0, 0):
        return "force"
    elif p.dimensions == (1, -1, -2, 0, 0, 0, 0):
        return "pressure"
    elif p.dimensions == (1, 0, -2, 0, 0, 0, 0):
        return "line"
    elif p.dimensions == (1, 2, -2, 0, 0, 0, 0):
        return "moment"
    else:
        return None
