from collections import UserDict
from dataclasses import dataclass
from typing import Optional, Any
from shapely import Polygon

            
@dataclass
class LoadElement:
    """
    Represents a loading area extracted from a markup
    """
    load_polygon: Polygon
    text_description: str
    occupancy: Optional[str] = None
    project_occupancies: Optional[str] = None
    load_components: Optional[dict] = None
    boundary: Optional[tuple[float]] = None # might be useful; remove if not

    def __post_init__(self):
        self.occupancy = self.parse_occupancy(self.text_description)
        if self.occupancy is not None and self.project_occupancies is not None:
            self.load_components = self.load_components_from_occupancies(self.occupancy, self.project_occupancies)
        else:
            self.load_components = self.parse_load_components(self.text_description)
        self.boundary = self.geometry.boundary


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
    pass

    def __post_init__(self):
        pass # TODO: This is where I left off


        


