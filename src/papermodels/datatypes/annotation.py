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
    page = Column(String)
    text = Column(String)
    vertices = Column(ScalarListType)  # [x1, y1, x2, y2, x3, y3, ..., xn, yn]
    matrix = Column(ScalarListType)  # [0., 0., 0., 0., 0., 0.]
    line_color = Column(ColorType)  # (r, g, b)
    line_weight = Column(Float)
    line_type = Column(ScalarListType)
    line_opacity = Column(Float)
    fill_color = Column(ColorType)
    fill_opacity = Column(Float)

    def __repr__(self):
        return class_representation(self)


A0 = Annotation(
    object_type="Polygon",
    vertices=[
        10.128380527777777,
        22.09484263888889,
        22.983913194444444,
        22.09484263888889,
        22.983913194444444,
        11.688489097222222,
        22.983913194444444,
        9.36691498611111,
        10.128369944444445,
        9.366916749999998,
        10.126076888888889,
        11.636233888888889,
        11.542573152777777,
        13.496487555555555,
        11.542573152777777,
        16.50542623611111,
        10.126076888888889,
        17.912785555555555,
    ],
    page="0",
    text=None,
    line_color=(1.0, 0.0, 0.0),
    line_weight=1.0,
    fill_color=(0.5019608, 1.0, 0.0),
    line_type=None,
    line_opacity=None,
    fill_opacity=None,
    matrix=[1.0, 0.0, 0.0, 1.0, -568.5768, -525.5376],
)
