from dataclasses import dataclass
from decimal import Decimal
from typing import Any



@dataclass(frozen=True)
class Annotation:
    page: int
    object_type: str
    text: str
    vertices: list[Decimal]
    line_color: tuple[Decimal]
    fill_color: tuple[Decimal]
    line_type: Any
    line_weight: Decimal
    line_opacity: Decimal
    fill_opacity: Decimal
    matrix: list[Decimal]
    local_id: int

    def __eq__(self, other):
        return all([
            self.page==other.page,
            self.object_type==other.object_type,
            self.text==other.text,
            self.vertices==other.vertices,
            self.matrix==other.matrix,
            self.line_color==other.line_color,
            self.fill_color==other.fill_color,
            self.line_type==other.line_type,
            self.line_weight==other.line_weight,
            self.line_opacity==other.line_opacity,
            self.fill_opacity==other.fill_opacity,
            self.matrix==other.matrix,
        ])



A0 = Annotation(
    object_type="Rectangle",
    vertices=[
        35.149490138888886,
        9.475074888888887,
        35.149490138888886,
        9.827993777777777,
        35.50240902777777,
        9.827993777777777,
        35.50240902777777,
        9.475074888888887,
    ],
    page="0",
    text=None,
    line_color=(0.0, 0.0, 0.0),
    line_weight=3.0,
    fill_color=(1, 1, 1),
    line_type=None,
    line_opacity=None,
    fill_opacity=None,
    matrix=[1.0, 0.0, 0.0, 1.0, -1992.727, -537.1696],
    local_id=3
)

A1 = Annotation(
    object_type="PolyLine",
    vertices=[
        35.149490138888886,
        9.475074888888887,
        35.149490138888886,
        9.827993777777777,
        35.50240902777777,
        9.827993777777777,
        35.50240902777777,
        9.475074888888887,
    ],
    page="1",
    text=None,
    line_color=(0.0, 0.0, 0.0),
    line_weight=3.0,
    fill_color=(1, 1, 1),
    line_type=None,
    line_opacity=None,
    fill_opacity=None,
    matrix=[1.0, 0.0, 0.0, 1.0, -1992.727, -537.1696],
    local_id=0
)
