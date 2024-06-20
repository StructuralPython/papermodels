from collections import UserDict
from dataclasses import dataclass
from typing import Optional, Any
from shapely import Polygon

from papermodels.paper.annotations import Annotation
from papermodels.units.parsing import parse_unit_system, convert_unit_string, UnitSystem


@dataclass
class LoadElement:
    """
    Represents a loading area extracted from a markup
    """

    load_polygon: Polygon
    occupancy: Optional[str] = None
    project_occupancies: Optional[dict] = None
    load_components: Optional[dict] = None
    unit_system: Optional[str | UnitSystem] = None

    def __post_init__(self):
        if self.unit_system is not None and isinstance(self.unit_system, str):
            self.unit_system = parse_unit_system(self.unit_system)

        if self.occupancy is not None and self.project_occupancies is not None:
            self.load_components = self.project_occupancies[self.occupancy]

        if self.load_components is not None:
            self.load_components = self.parse_load_components()

    @classmethod
    def from_parsed_annotation(
        cls,
        annotation_props: dict,
        project_occupancies: Optional[dict] = None,
        unit_system: Optional[str] = None,
    ):
        """
        Returns a LoadElement based on the data from 'annot' and 'annotation_props'
        which is a dictionary of the values resulting from the function
        papermodels.paper.annotations.parse_annotations.

        'annot': the Annotation
        'annotation_props': dict with keys "type" and "geometry", where "type"
            corresponds to the legend data associated with the annotation that
            would have been extracted during the parse_annotations function.

            For loads, the value for "type" should be a str formatted like the following:
            "Load
            [occupancy=<occupancy_name>]
            [[load_component=load_magnitude]
             [load_component=load_magnitude], etc. ]

             e.g.
             Load
             occupancy=Residential

             e.g.
             Load
             D=34 psf
             L=4.8
             S=1.2 kPa

             The load component must be on the left-hand side of an "="
             and the load magnitude can either be a str or float representing
             a magnitude. If not convertible to a float, the magnitude will need
             to be parsed.
        'project_occupancies': a dict describine occupancy names and their load components
        'unit_system': one of {'psf', 'ksf', 'psi', 'ksi', 'kPa', 'MPa', 'GPa'}
        """
        annotation_type = annotation_props["type"]
        parsed_load_text = parse_load_text(annotation_type)
        occupancy = parsed_load_text.pop("occupancy", None)
        load_components = parsed_load_text.copy()
        return cls(
            load_polygon=annotation_props["geometry"],
            occupancy=occupancy,
            project_occupancies=project_occupancies,
            load_components=load_components,
            unit_system=unit_system,
        )

    def parse_load_components(self) -> dict:
        """
        Returns a copy of self.load_components where the values of the keys have
        been converted into a float with appropriate scaling as per their units
        and the overall unit system.
        """
        return {
            component: convert_unit_string(magnitude, self.unit_system)
            for component, magnitude in self.load_components.items()
        }


@dataclass
class LoadArray(UserDict):
    """
    A type alias to represent LoadElements extracted from the paper model.

    The actual dictionary is also stored under the 'data' attribute, as per the
    conventional UserDict behaviour.

    The dictionary stored within should be in the following structure.
        {
            "page_number_1": [
                load_element_1,
                load_element_2,
                load_element_3
            ],
            "page_number_2": [
                load_element_3,
                load_element_4
            ]
        }
    """

    @classmethod
    def from_annotations(cls, parsed_load_annotations: dict[Annotation]):
        """
        Returns a LoadArray-style dict populated with LoadElements.

        parsed_load_annotations is a dict with Annotation objects as keys
        and dict of values. The values dict has keys, "type" and "geometry".
        The value of the "type" key will match the annotation's corresponding
        legend entry so the use of the user-defined legend will persist.
        """
        for annot, annot_props in parsed_load_annotations.items():
            annot_type = annot_props["type"]
            annot_geom = annot_props["geometry"]


def parse_load_text(text: str) -> dict:
    """
    Returns a dict representing the parsed load data present in 'text'
    """
    text_by_line = text.split("\n")
    parsed_load_data = {}
    for line in text_by_line:
        if "occupancy" in line.lower() and "=" in line.lower():
            lhs, rhs = line.split("=")
            rhs = rhs.strip(" ")
            parsed_load_data["occupancy"] = rhs
        elif "=" in line.lower():
            lhs, rhs = line.split("=")
            lhs = lhs.rstrip(" ").strip(" ")
            rhs = rhs.strip(" ").rstrip(" ")
            parsed_load_data[lhs] = rhs
    return parsed_load_data
