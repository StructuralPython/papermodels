from __future__ import annotations
from decimal import Decimal
from shapely import LineString, Polygon
from typing import Union, Any, Optional
from papermodels.paper.annotations import Annotation
from papermodels.paper.pdf import load_pdf_annotations

def shapely_to_annotation(
    geom: Union[LineString, Polygon],
    page: int = 0,
    text: str = "",
    line_color: tuple[Decimal] = (Decimal("0.0"), Decimal("0.0"), Decimal("0.0")),
    fill_color: Optional[tuple[Decimal]] = None,
    line_type: Optional[Any] = None,
    line_weight: Decimal = Decimal("1.0"),
    line_opacity: Decimal = Decimal("1.0"),
    fill_opacity: Decimal = Decimal("1.0"),
    matrix: list[Decimal] = [Decimal("1.0"), Decimal("0.0"), Decimal("0.0"), Decimal("1.0"), Decimal("0.0"), Decimal("0.0")],
    local_id: Optional[int] = None
) -> Annotation:
    """
    Returns a valid Annotation object based on the provided 'geom' and the selected
    options.
    """
    if geom.geom_type == "Polygon":
        annot_type = "Polygon"
        shape = geom.exterior
    elif geom.geom_type == "LineString":
        annot_type = "Line"
        shape = geom
    
    vertices = []
    for coord in shape.coords:
        x, y = coord
        vertices.append(Decimal(x))
        vertices.append(Decimal(y))
    
    return Annotation(
        page=page,
        object_type=annot_type,
        text=text,
        vertices=vertices,
        line_color=line_color,
        fill_color=fill_color,
        line_type=line_type,
        line_weight=line_weight,
        line_opacity=line_opacity,
        fill_opacity=fill_opacity,
        matrix=matrix,
        local_id=local_id
    )