from dataclasses import dataclass
from typing import Any, Optional

@dataclass
class Joist:
    """
    Models a joist with a uniform load of 
    'w' on all spans of the joist that exist.

                 w
    ||||||||||||||||||||||||||||
    ----------------------------
        ^                  ^
        R1                 R2
    < a ><      span      >< b >
    """
    span: float | Any
    a: float | Any = 0.
    b: float | Any = 0.

    def get_r1(self):
        w = 1
        r1_a = w * self.a / (2 * self.span) * (2 * self.span + self.a)
        r1_span = w * self.span / 2
        r1_b = -w * self.a**2 / (2 * self.span)
        total_r1 = sum([r1_a, r1_span, r1_b])
        total_load = self.get_total_load()
        return total_r1 / total_load
    
    def get_r2(self):
        w = 1
        r2_b = w * self.b / (2 * self.span) * (2 * self.span + self.b)
        r2_span = w * self.span / 2
        r2_a = -w * self.b**2 / (2 * self.span)
        total_r2 = sum([r2_a, r2_span, r2_b])
        total_load = self.get_total_load()
        return total_r2 / total_load

    def get_total_load(self):
        w = 1
        total_load_a = w * self.a
        total_load_span = w * self.span
        total_load_b = w * self.b
        total_load = sum([total_load_a, total_load_span, total_load_b])
        return total_load

