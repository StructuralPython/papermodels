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
    __tablename__ = 'annotations'
    
    annot_id = Column(Integer, primary_key=True)
    object_type = Column(String)
    page = Column(String)
    text = Column(String)
    vertices = Column(ScalarListType)
    matrix = Column(ScalarListType)
    line_color = Column(ColorType)
    line_weight = Column(Float)
    line_type = Column(ScalarListType)
    line_opacity = Column(Float)
    fill_color = Column(ColorType)
    fill_opacity = Column(Float)
    
    def __repr__(self):
        return class_representation(self)
