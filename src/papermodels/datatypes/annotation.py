from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy_utils import ColorType
from sqlalchemy_utils import ScalarListType

from papermodels.datatypes.utils import class_representation

Base = declarative_base()


class Annotation(Base):
    """
    A data type to represent a PDF annotation exported through an FDF file.
    """

    __tablename__ = "annotations"

    annot_id = Column(Integer, primary_key=True)
    object_type = Column(
        String
    )  # One of "Line", "Circle", "Rectangle", "Square", "PolyLine", "Polygon"
    page = Column(String)  # "0" or integer as a string
    text = Column(String)  # Optional[str]
    vertices = Column(ScalarListType)  # [x1, y1, x2, y2, x3, y3, ..., xn, yn]; auto-closed depending on object_type
    matrix = Column(ScalarListType)  # [1., 0., 0., 1., 0., 0.]
    line_color = Column(ColorType)  # (r, g, b)
    line_weight = Column(Float)  # 0 <= x
    line_type = Column(ScalarListType)
    line_opacity = Column(Float)  # 0 <= x <= 1
    fill_color = Column(ColorType)  # (r, g, b)
    fill_opacity = Column(Float)  # 0 <= x <= 1

    def __repr__(self):
        return class_representation(self)


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
