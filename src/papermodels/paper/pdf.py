from __future__ import annotations
from typing import Optional, Any, Union
import pikepdf as pike
import re
from decimal import Decimal
from papermodels.datatypes.annotation import Annotation
import pathlib
import parse


def load_pdf_annotations(pdf_path: pathlib.Path | str) -> list[Annotation]:
    """
    Returns a lists of pdf annotations keyed by page index.

    'pdf_path': Path-like object representing the path to the PDF file to open.
    """
    pdf_path = pathlib.Path(pdf_path)
    pdf_obj = pike.open(pdf_path)
    annots_in_pdf = []

    # Numeric data extracted from pikepdf is type Decimal
    for page_num, page_data in enumerate(pdf_obj.pages):
        for annot in page_data.obj.Annots:
            if "/Parent" in annot:
                continue
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
            if annot_type.lower() in ("polygon", "polyline", "circle", "ellipse"):
                vertices = list(annot.get("/Vertices", []))
            elif annot_type.lower() in (
                "rectangle",
                "square",
                "rectangle sketch to scale",
            ):
                bbox = list(annot.get("/Rect", []))
                buffers = list(annot.get("/RD", [0, 0, 0, 0]))
                x1, y1, x2, y2 = bbox
                b1, b2, b3, b4 = buffers
                vertices = [
                    x1,
                    y1,
                    x1,
                    y2 - b4 - b2,
                    x2 - b3 - b1,
                    y2 - b4 - b2,
                    x2 - b3 - b1,
                    y1,
                ]
            elif annot_type == "Line":
                vertices = list(annot.get("/L", []))
            else:
                print(f"Cannot read (yet): {annot_type}")

            text = str(annot.get("/Contents", ""))
            rc = annot.get("/RC")
            text = parse_html_text_content(str(rc)) if rc is not None else ""
            line_color = tuple(annot.get("/C", stream_dict.get("RG", (0, 0, 0))))
            fill_color = tuple(annot.get("/IC", stream_dict.get("rg", (1, 1, 1))))
            fill_opacity = annot.get("/FillOpacity", Decimal("1.0"))
            line_width = stream_dict.get("w", [1])[0]
            line_type = None
            line_opacity = 1
            matrix = (1, 0, 0, 1, 0, 0)
            annotation = Annotation(
                object_type=annot_type,
                page=page_num,
                vertices=tuple(vertices),
                text=text,
                line_color=line_color,
                fill_color=fill_color,
                fill_opacity=fill_opacity,
                line_weight=line_width,
                line_type=line_type,
                line_opacity=line_opacity,
                matrix=matrix,
            )
            annots_in_pdf.append(annotation)
    return annots_in_pdf


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
