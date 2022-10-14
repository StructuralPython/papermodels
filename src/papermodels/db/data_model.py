from typing import Any
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, ForeignKey
from sqlalchemy_utils import ColorType
from sqlalchemy_utils import ScalarListType


Base = declarative_base()


annot_geom = Table(
    "annot_geom",
    Base.metadata,
    Column("annot_id", Integer, ForeignKey("annotations.annot_id")),
    Column("geom_id", Integer, ForeignKey("geometry.geom_id")),
)

annot_label = Table(
    "annot_label",
    Base.metadata,
    Column("annot_id", Integer, ForeignKey("annotations.annot_id")),
    Column("label_id", Integer, ForeignKey("labels.label_id")),
)


class Annotation(Base):
    """
    A data type to represent a PDF annotation exported through an FDF file.
    """

    __tablename__ = "annotations"

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


# Examples
# A0 = Annotation(
#     annot_id=0,
#     object_type="rect",
#     page=
# )


class Geometry(Base):

    __tablename__ = "geometry"

    geom_id = Column(Integer, primary_key=True)
    wkt = Column(String)
    label = Column(String)

    def __repr__(self):
        return class_representation(self)


class Labels(Base):

    __tablename__ = "labels"

    label_id = Column(Integer, primary_key=True)
    label = Column(String)
    assigned_to = Column(String)

    def __repr__(self):
        return class_representation(self)


class LegendEntry(Base):

    __tablename__ = "legend"

    legend_id = Column(Integer, primary_key=True)
    key = Column(String)
    value = Column(String)

    def __repr__(self):
        return class_representation(self)


def class_representation(object: Any) -> str:
    """
    Returns a generic repr string for a given 'object' where 'object' is a
    "data class"-like object.
    """
    class_name = str(type(object)).split(".")[-1].replace("'", "").replace(">", "")
    attrs = []
    for k, v in object.__dict__.items():
        if k.startswith("_"):
            continue
        if isinstance(v, str):
            v = f"'{v}'"
        attrs.append(f"{k}={v}")
    return f"{class_name}({', '.join(attrs)})"


def db_connect():
    """
    Returns a database connection for use. Currently set to use
    sqlite in memory for testing
    """
    engine = create_engine("sqlite+pysqlite:///:memory:", echo=True, future=True)
    return engine
