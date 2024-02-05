from __future__ import annotations

from dataclasses import dataclass
from papermodels.datatypes.element import Element

@dataclass
class WallModel:
    element: Element

    @classmethod
    def from_element(cls, element: Element) -> WallModel:
        return cls(element)