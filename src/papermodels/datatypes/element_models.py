from dataclasses import dataclass
from enum import Enum
from typing import Optional
from papermodels.datatypes.element import Element

class ReactionType(Enum):
    MOMENT = 0
    POINT = 1
    UNIFORM_LINE = 2
    VARIABLE_LINE = 3
    UNIFORM_AREA = 4
    VARIABLE_AREA = 5

@dataclass
class ElementModel:
    element: Optional[Element] = None
    reaction_type: Optional[ReactionType] = None

    @classmethod
    def from_element(cls, element: Element) -> None:
        return cls(element)

    def reactions(self) -> tuple[float]:
        raise NotImplemented
    

    