from __future__ import annotations
from copy import deepcopy
from typing import Optional, Any, Union
import pikepdf as pike
import re
from decimal import Decimal
from papermodels.datatypes.annotation import Annotation
from papermodels.paper.annotations import scale_annotations
from papermodels.geometry import geom_ops
import pathlib
import parse
import numpy as np


def load_pdf_annotations(pdf_path: pathlib.Path | str, show_skipped: bool = False) -> list[Annotation]:
    """
    Returns a lists of pdf annotations keyed by page index.

    'pdf_path': Path-like object representing the path to the PDF file to open.
    """
    pdf_path = pathlib.Path(pdf_path)
    with pike.open(pdf_path) as pdf_obj:
        annots_in_pdf = []
        skipped_annots = []
        for page_num, page_data in enumerate(pdf_obj.pages):
            rotate = page_data.get("/Rotate", None)
            for annot_idx, annot in enumerate(page_data.obj.Annots):
                pm_annot = pike_annotation_to_pm_annotation(annot, annot_idx, page_num, rotate)
                if pm_annot is not None:
                    annots_in_pdf.append(pm_annot)
                else:
                    skipped_annots.append(annot)
        if show_skipped:
            print(skipped_annots)
        return annots_in_pdf


def update_pdf_annotations(
        pdf_path: pathlib.Path, 
        parsed_annotations: dict[Annotation, dict],
        export_path: pathlib.Path,
        append_text: bool = False,
        ) -> None:
    """
    Returns None. Creates a copy of the file at 'pdf_path' with all of annotations representing 
    structural elements (collector and transfer elements, not loads or legend entries) having
    their text field updated with the assigned tag of the structural element.
    """
    pdf_path = pathlib.Path(pdf_path).resolve()
    with pike.open(pdf_path,) as pdf_obj:
        for page_idx, page_data in enumerate(pdf_obj.pages):
            for annot_idx, annot in enumerate(page_data.obj.Annots):
                for parsed_annotation, annot_attrs in parsed_annotations.items():
                    annotations_equal = compare_annotations(parsed_annotation, annot, page_idx)
                    if annotations_equal:
                        text_to_add = f"tag: {annot_attrs['tag']}"
                        replace_text = not append_text
                        updated_pike_annot = update_annotation_text_field(annot, text_to_add, replace_text)
                        page_data.obj.Annots[annot_idx] = updated_pike_annot
        pdf_obj.save(export_path)


def update_annotation_text_field(pike_annot: pike.Annotation, text: str, replace: bool = False):
    """
    Returns a copy of the 'pike_annot' with the text field updated with 'text'. If 'replace'
    is True then the text field is completely replaced with 'text'. If not, then 'text'
    is appended on a new line to the existing 'text'.
    """
    orig_text = pike_annot.get("/Contents")
    if "tag:" in str(orig_text) and not replace: # Don't append a tag if there is one there
        return pike_annot
    elif replace or orig_text is None:
       pike_annot['/Contents'] = pike.String(text)
    else:
       pike_annot['/Contents'] = pike.String(f"{str(orig_text)}\r\n{text}")
    return pike_annot


def compare_annotations(pm_annot: Annotation, pike_annot: pike.Annotation, page_num: int) -> bool:
    """
    Returns True if the `pm_annot` is nominally equal to the 'pike_annot' by converting
    the 'pike_annot' to the 'pm_annot' and seeing if they are equal. Will return True
    for duplicate annotations.
    """
    converted = pike_annotation_to_pm_annotation(pike_annot, pm_annot.local_id, page_num)
    if converted is None:
        return False
    return pm_annot == converted


def pike_annotation_to_pm_annotation(annot, annot_idx: int, page_idx: int, rotate: Optional[int] = None) -> Optional[Annotation]:
    """
    Returns either an Annotation object or None. None is returned if:
        - The annotation has a "parent" key
        - Has the word "ignore" in the text field
        - If that annotation is a type that cannot be read yet
        - The annotation does not have a Subj key

    Note: Numerical data within pike annotation objects are of type Decimal
    """
    if "/Parent" in annot:
        return None
    # Content stream for Bluebeam PDFs is linked as an XObject
    # The below gets the linked XObject content stream
    content_stream = annot.get("/AP", {}).get("/N", None)
    if content_stream is not None:
        utf_stream = pike.Stream.read_bytes(content_stream).decode("utf-8")
        stream_dict = parse_content_stream(utf_stream)
    else:
        stream_dict = {}
    if "/Subj" in annot:
        annot_type = str(annot["/Subj"])
    else:
        return None
    if annot_type.lower() in ("polygon", "polyline", "circle", "ellipse"):
        vertices = tuple(annot.get("/Vertices", tuple()))
    elif annot_type.lower() in (
        "rectangle",
        "square",
        "rectangle sketch to scale",
    ):
        bbox = tuple(annot.get("/Rect", tuple()))
        buffers = tuple(annot.get("/RD", (0, 0, 0, 0)))
        x1, y1, x2, y2 = bbox
        b1, b2, b3, b4 = buffers
        vertices = (
            x1,
            y1,
            x1,
            y2 - b4 - b2,
            x2 - b3 - b1,
            y2 - b4 - b2,
            x2 - b3 - b1,
            y1,
        )
    elif annot_type == "Line":
        vertices = tuple(annot.get("/L", tuple()))
    else:
        print(f"Cannot read (yet): {annot_type}")
        return None
    if rotate == 90:
        vertex_array = geom_ops.vertices_to_array(vertices)
        rotated_vertices = geom_ops.rotate_90_coords(vertex_array, ccw=True)
        vertices = geom_ops.flatten_vertex_array(rotated_vertices)
    elif rotate == 270:
        vertex_array = geom_ops.vertices_to_array(vertices)
        rotated_vertices = geom_ops.rotate_90_coords(vertex_array, ccw=False)
        vertices = geom_ops.flatten_vertex_array(rotated_vertices)
    text = str(annot.get("/Contents", ""))
    rc = annot.get("/RC")
    text = parse_html_text_content(str(rc)) if rc is not None else text
    if "ignore" in text:
        return None
    line_color = tuple(annot.get("/C", stream_dict.get("RG", (0, 0, 0))))
    fill_color = tuple(annot.get("/IC", stream_dict.get("rg", (1, 1, 1))))
    # fill_opacity and global_opacity are different things and should be separated
    # as of 2025-03-31 I am treating global opacity as being the same as fill
    # opacity for a quick fix for editors that do not allow control over only
    # the fill opacity.
    fill_opacity = annot.get("/FillOpacity", Decimal("1.0"))
    line_opacity = Decimal("1.0")
    global_opacity = annot.get("/CA", Decimal("1.0"))
    if global_opacity != Decimal("1.0") and fill_opacity == Decimal("1.0"):
        fill_opacity = global_opacity
    if global_opacity != Decimal("1.0") and line_opacity == Decimal("1.0"):
        line_opacity = global_opacity
    line_width = stream_dict.get("w", (1,))[0]
    line_type = None
    matrix = (1, 0, 0, 1, 0, 0)

    annotation = Annotation(
        object_type=annot_type,
        page=page_idx,
        vertices=tuple(vertices),
        text=text,
        line_color=line_color,
        fill_color=fill_color,
        fill_opacity=fill_opacity,
        line_weight=line_width,
        line_type=line_type,
        line_opacity=line_opacity,
        matrix=matrix,
        local_id=annot_idx
    )
    return annotation

def parse_content_stream(stream: str) -> dict[str, list]:
    """
    Parses data in the content stream into an ordered dictionary.
    Keys are operators. Values are operands.

    Currently only parses non-path data reliably. Path data
    is not parsed separately so, in the output dict, repetitive
    operators are overwritten in the dict (e.g. multiple 'l'
    operators will only have one entry).
    """
    commands = {}
    operand_with_operator = re.compile(r"[0-9.\s]+[a-zA-Z]+")
    operators = re.compile(f"[a-zA-Z]+")
    operands = re.compile(r"[\d.]+")

    matches = operand_with_operator.findall(stream)
    for match in matches:
        operator = operators.findall(match)[0]
        operand = operands.findall(match)
        numerical_operands = []
        for element in operand:
            if "." in element:
                elem = Decimal(element)
            else:
                elem = int(element)
            numerical_operands.append(elem)
        commands.update({operator: numerical_operands})
    return commands


def parse_html_text_content(xml_code: str) -> str:
    """
    Returns 'Python-formatted' text for html text in the annotation popup
    """
    if xml_code is None:
        return ""
    else:
        format_code = "<body{xml_code}>{html_code}</body>"
        results = parse.search(format_code, xml_code)
        html_code = results.named["html_code"]
        python_text = (
            html_code.replace("</p><p>", "\n").replace("<p>", "").replace("</p>", "")
        )
        return python_text