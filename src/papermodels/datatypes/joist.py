from dataclasses import dataclass
from typing import Any, Optional
import pycba as cba


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
    
    def __post_init__(self):
        L=[self.a, self.span, self.b]
        EI=[1e3, 1e3, 1e3]
        R=[0., 0., -1., 0., -1., 0., 0., 0.,]

        if self.a == 0:
            L.pop(0)
            EI.pop(0)
            R.pop(0)
            R.pop(0)
        if self.b == 0:
            L.pop()
            EI.pop()
            R.pop()
            R.pop()

        self._pycba_model = cba.BeamAnalysis(
            L,
            EI,
            R,
        )
        for idx, _ in enumerate(L):
            self._pycba_model.add_udl(idx + 1, 1) # 1-based idx


    def get_r1(self):
        self._pycba_model.analyze()
        total_r1 = self._pycba_model._beam_results.R[0]
        total_load = self.get_total_load()
        return round(total_r1 / total_load, 9)
    
    def get_r2(self):
        self._pycba_model.analyze()
        total_r2 = self._pycba_model._beam_results.R[1]
        total_load = self.get_total_load()
        return round(total_r2 / total_load, 9)

    def get_total_load(self):
        w = 1
        total_load_a = w * self.a
        total_load_span = w * self.span
        total_load_b = w * self.b
        total_load = sum([total_load_a, total_load_span, total_load_b])
        return total_load

