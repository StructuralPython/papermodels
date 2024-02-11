from dataclasses import dataclass
from typing import Optional, Any
from shapely import Polygon

@dataclass
class LoadParameters:
    occupancy: Optional[str] = None
    project_occupancies: Optional[dict] = None
    load_components: Optional[dict] = None

    def __post_init__(self):
        if self.occupancy is not None and self.project_occupancies is not None:
            self.load_components = self.get(self.occupancy, None)


@dataclass
class TribLoad:
    """
    A class to represent a load that is generated from the intersection of an Element
    trib region and a load polygon.
    """
    load_polygon: Polygon
    load_parameters: dict
    trib_region: Polygon
    element_id: str

