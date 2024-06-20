from __future__ import annotations
from dataclasses import dataclass
from enum import IntEnum
import pathlib
from typing import Optional, Any
from shapely import Polygon
from papermodels.datatypes.element import Element
from papermodels.datatypes.analysis_models import PyNiteFEModel
from papermodels.geometry import geom_ops as go
from papermodels.fileio.model_files import get_structured_beam_data
from papermodels.fileio.utils import read_csv_file


class ReactionType(IntEnum):
    MOMENT = 0
    POINT = 1
    UNIFORM_LINE = 2
    VARIABLE_LINE = 3
    UNIFORM_AREA = 4
    VARIABLE_AREA = 5


class AnalysisModelType(IntEnum):
    PYNITEFEA: 0
    SAP2000: 0


@dataclass
class ElementModel:
    element: Optional[Element] = None
    reaction_type: Optional[ReactionType] = None
    model: Optional[Any] = None
    structured_element_data: Optional[dict] = None
    section: dict = None
    trib_region: Optional[Polygon] = None

    @classmethod
    def from_element(cls, element: Element):
        return cls(element=element)

    @classmethod
    def from_file(cls, filepath: pathlib.Path | str):
        raise NotImplemented

    def determine_loads(self, loads: list[dict]) -> None:
        """
        Returns None. Performs a geometric intersection with each load present
        in 'loads' to determine if any of the loads intersect with the
        self.trib_area of the element.
        """
        if self.trib_area is None:
            return
        for load in loads:
            if load["geometry"].intersects(self.trib_region):
                applied_load = go.get_applied_load(load["geometry"], self.trib_region)
                if self.structured_element_data is None:
                    self.structured_element_data = {}
                self.structured_element_data["Loads"].append(applied_load)

    def get_reactions(self) -> tuple[float]:
        raise NotImplemented

    def create_model(self) -> None:
        raise NotImplemented

    def analyze_model(self) -> None:
        raise NotImplemented


@dataclass
class BeamModel(ElementModel):

    def from_file(cls, filepath: pathlib.Path):
        raw_data = read_csv_file(filepath)
        structured_beam_data = get_structured_beam_data(raw_data)
        analysis_model = PyNiteFEModel(structured_beam_data)
        analysis_model.create_model()
        beam_model = cls(
            structured_element_data=structured_beam_data,
            reaction_type=ReactionType.POINT,
        )
        return beam_model


@dataclass
class ColumnModel(ElementModel):
    pass


@dataclass
class WallModel(ElementModel):
    pass


class JoistArrayModel(ElementModel):
    """
    Models a spread of joists over a region where the distance
    between the supports may vary linearly.
    """

    initial_offset: float | int = (0.0,)
    joist_at_start: bool = (True,)
    joist_at_end: bool = (False,)
    cantilever_tolerance: float = (1e-2,)

    def __post_init__(self):
        joist_supports = [inter[2] for inter in self.element.intersections]
        joist_prototype = self.element.geometry
        self.id = self.element.tag
        self.spacing = 400  # Need to include this in the legend and thus, the Element
        self.initial_offset = float(self.initial_offset)
        self._joist_prototype = joist_prototype
        self._cantilever_tolerance = self.cantilever_tolerance
        self._extents = go.get_joist_extents(joist_prototype, joist_supports)
        self._supports = go.determine_support_order(joist_prototype, joist_supports)
        self._cantilevers = go.get_cantilever_segments(joist_prototype, self._supports)
        self.vector_parallel = go.get_direction_vector(joist_prototype)
        self.vector_normal = go.rotate_90(self.vector_parallel, ccw=False)
        self.joist_at_start = float(self.joist_at_start)
        self.joist_at_end = float(self.joist_at_end)
        self.joist_locations = go.get_joist_locations(
            self.get_extent_edge("start"),
            self.get_extent_edge("end"),
            self.spacing,
            self.initial_offset,
            self.joist_at_start,
        )
        self.joists = [
            self.generate_joist(idx) for idx, _ in enumerate(self.joist_locations)
        ]
        self.joist_trib_widths = [
            self.get_joist_trib_widths(idx)
            for idx, _ in enumerate(self.joist_locations)
        ]
        self.joist_trib_areas = [
            self.generate_trib_area(idx) for idx, _ in enumerate(self.joist_locations)
        ]

    # def __repr__(self):
    #     return class_representation(self)

    @classmethod
    def from_element(
        cls,
        element: Optional[Element],
        initial_offset: float | int = 0.0,
        joist_at_start: bool = True,
        joist_at_end: bool = False,
        cantilever_tolerance: float = 1e-2,
    ) -> JoistArrayModel:
        return cls(
            element, initial_offset, joist_at_start, joist_at_end, cantilever_tolerance
        )

    def generate_joist(self, index: int):
        """
        Returns i, j coordinates of the joist in the JoistArray at the position
        of 'index'. Raises IndexError if 'index' is not within the joist array
        extents given the spacing.

        'index': joists are numbered from 0 (first joist, at joist extent) and
            go to n, the last joist in the array.
        """
        start_centroid = self.get_extent_edge("start").centroid
        try:
            joist_distance = self.joist_locations[index]
        except IndexError as e:
            raise IndexError(
                f"Joist index {index} is beyond the extent of the joist array for {self.id}. "
                f"Last index is {len(self.joist_locations) - 1} @ {self.joist_locations[-1]}"
            ) from None

        if index != 0 and index != len(self.joist_locations) - 1:
            new_centroid = go.project_node(
                start_centroid, -self.vector_normal, joist_distance
            )

            system_bounds = go.get_system_bounds(
                self._joist_prototype, list(self._supports.values())
            )
            projection_distance = go.get_magnitude(system_bounds)
            ray_aj = go.project_node(
                new_centroid, -self.vector_parallel, projection_distance
            )
            ray_a = go.create_linestring([new_centroid, ray_aj])
            ray_bj = go.project_node(
                new_centroid, self.vector_parallel, projection_distance
            )
            ray_b = go.create_linestring([new_centroid, ray_bj])
            support_a_loc = ray_a.intersection(self._supports["A"])
            support_b_loc = ray_b.intersection(self._supports["B"])

            end_a = support_a_loc
            end_b = support_b_loc

        # These clauses req'd to deal with floating point error possible
        # on the end joists (occurs after performing project_node)
        elif index == 0:
            end_a = support_a_loc = self._extents["A"][0]
            end_b = support_b_loc = self._extents["B"][0]
        elif index == len(self.joist_locations) - 1:
            end_a = support_a_loc = self._extents["A"][1]
            end_b = support_b_loc = self._extents["B"][1]

        if self._cantilevers["A"]:
            end_a = go.project_node(
                support_a_loc, -self.vector_parallel, self._cantilevers["A"]
            )
        if self._cantilevers["B"]:
            end_b = go.project_node(
                support_b_loc, self.vector_parallel, self._cantilevers["B"]
            )
        return go.create_linestring([end_a, end_b])

    def get_extent_edge(self, edge: str = "start"):
        """
        Gets the "joist" that would exist at the edge of the array

        'edge': one of {'start', 'end'}
        """
        if edge == "start":
            node_i = self._extents["A"][0]
            node_j = self._extents["B"][0]
        elif edge == "end":
            node_i = self._extents["A"][1]
            node_j = self._extents["B"][1]
        return go.create_linestring([node_i, node_j])

    def get_joist_trib_widths(self, index) -> tuple[float, float]:
        """
        Returns the trib widths of the the joist at 'index'. The trib
        widths are a tuple representing the left and right width,
        respectively.
        """
        if index < 0:
            # Convert -ve index lookup to a +ve index lookup
            index = len(self.joist_locations) + index
        if index == 0:  # The first joist
            trib_widths = (0.0, self.joist_locations[1] / 2.0)
        elif index == len(self.joist_locations) - 1:  # The last joist
            spacing_left = self.joist_locations[-1] - self.joist_locations[-2]
            trib_widths = (spacing_left / 2.0, 0.0)
        else:
            spacing_left = self.joist_locations[index] - self.joist_locations[index - 1]
            spacing_right = (
                self.joist_locations[index + 1] - self.joist_locations[index]
            )
            trib_widths = (spacing_left / 2.0, spacing_right / 2.0)
        return trib_widths

    def generate_trib_area(self, index: int) -> Polygon:
        """
        Returns a tuple of Polygon representing the tributary area of the 'joist' based on the
        given 'trib_widths'
        """
        joist = self.joists[index]
        trib_widths = self.joist_trib_widths[index]
        i_node, j_node = joist.boundary.geoms  # Point, Point
        trib_left, trib_right = trib_widths  # float, float

        # Left
        if trib_left != 0.0:
            i_left = go.project_node(i_node, self.vector_normal, trib_left)
            j_left = go.project_node(j_node, self.vector_normal, trib_left)
            trib_area_left = go.get_convex_hull(
                go.create_multipolygon([i_left, j_left, j_node, i_node])
            )
        else:
            trib_area_left = go.create_polygon([])

        # Right
        if trib_right != 0.0:
            i_right = go.project_node(i_node, -self.vector_normal, trib_right)
            j_right = go.project_node(j_node, -self.vector_normal, trib_right)
            trib_area_right = go.create_convex_hull(
                go.create_multipoint([i_right, j_right, j_node, i_node])
            )
        else:
            trib_area_right = go.create_polygon([])
        return trib_area_left | trib_area_right
