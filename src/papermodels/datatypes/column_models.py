from __future__ import annotations
from dataclasses import dataclass
from papermodels.datatypes.element import Element

@dataclass
class ColumnModel:
    element: Element

    @classmethod
    def from_element(cls, element: Element) -> ColumnModel:
        return cls(element)