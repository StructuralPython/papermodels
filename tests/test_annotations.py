from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt

def test__annotation_to_wkt():
    wkt_0 = _annotation_to_wkt(A0)
    wkt_1 = _annotation_to_wkt(A1)
    # assert wkt_0 == 

