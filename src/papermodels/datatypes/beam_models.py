from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
from papermodels.datatypes.element import Element

@dataclass
class BeamModel:
    element: Optional[Element] = None

    @classmethod
    def from_element(cls, element: Element) -> BeamModel:
        return cls(element)