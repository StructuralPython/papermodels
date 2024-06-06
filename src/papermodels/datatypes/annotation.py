from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy_utils import ColorType
from sqlalchemy_utils import ScalarListType

from dataclasses import dataclass
from decimal import Decimal
from typing import Any
from papermodels.datatypes.utils import class_representation


@dataclass(frozen=True)
class Annotation:
    page: int
    object_type: str
    text: str
    vertices: list[Decimal]
    matrix: list[Decimal]
    line_color: tuple[Decimal]
    fill_color: tuple[Decimal]
    line_type: Any
    line_weight: Decimal
    line_opacity: Decimal
    fill_opacity: Decimal
    matrix: list[Decimal]


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
)
