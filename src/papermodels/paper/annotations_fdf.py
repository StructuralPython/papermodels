from typing import Optional, Any, Union
import pathlib

from colour import Color
import numpy as np
import parse

from papermodels.db.data_model import Annotation


def read_fdf_file(file_path: pathlib.Path) -> list[str]:
    """
    Reads the FDF file and returns a list of strings.
    """
    with open(file_path, 'rb') as file:
        acc = []
        for line in file:
            try:
                acc.append(line.decode('utf-8'))
            except:
                pass
    return acc


def get_annotations_from_fdf(fdf_str: str) -> list[Annotation]:
    """
    Separates FDF data by objects
    """
    annotations = []
    in_stream_data = False
    stream_data = None
    stream_properties = {}
    annotation_properties = {}
    annot_type, vertices = None, None
    annotation = None
    
    # This upcoming for/if/elif block assumes the following:
    #   1. The FDF file contains geometric annotations, which we want, combined with
    #      other annotation data (e.g. bounding boxes) that are related to the geometric
    #      annotations but which we do not want.
    #   2. The general format of an FDF geometric annotation is as follows:
    #     a) The geometric annotation in 'obj<<' format as an object
    #     b) Some other information in the 'obj<<' format as other objects (bounding boxes, etc.)
    #     c) The stream data that contains many of the geometry's properties
    #   3. We are going to get the vertices information (geometry) from the 
    #     'obj<<' formatted data and most of its properties from the stream
    #     which means we need to retain the geometric data as we iterate before
    #     we finally get to that object's applicable stream data.
    #  A visual inspection of an FDF file with geometric markup should be able to inform
    #  the general approach taken.
        
    for line in fdf_str:
        if "endstream" in line and stream_data:
            stream_properties = extract_stream_properties(stream_data)
            stream_data = None
        elif in_stream_data == True:
            stream_data = line
            in_stream_data = False
            continue
        elif not parse.search("{} 0 obj<<", line) and "stream" not in line:
            continue
        elif "stream" in line:
            in_stream_data = True
            continue       
        type_and_vertices = extract_type_and_vertices(line)
        object_properties = extract_object_properties(line)
        if annot_type and vertices and stream_properties:
            annotation_properties.update(stream_properties)
            annotation = Annotation(
                object_type=annot_type,
                vertices=vertices,
                page=page,
                text=text,
                line_color=Color(rgb=stream_properties["line_color"]),
                line_weight=stream_properties["line_weight"],
                fill_color=Color(rgb=stream_properties["fill_color"]),
                line_type=stream_properties["line_type"],
                line_opacity=object_properties["line_opacity"],
                fill_opacity=object_properties["fill_opacity"],
                matrix=matrix,
            )
            stream_properties = {}
            type_and_vertices = None
            annotations.append(annotation)
        elif type_and_vertices:
            annot_type, vertices = type_and_vertices
            text = object_properties.pop("label")
            page = object_properties.pop("page")
            annotation_properties.update(object_properties)
        elif object_properties['matrix'] is not None:
            matrix = object_properties.pop('matrix')
    return annotations


def extract_type_and_vertices(line: str) -> Optional[tuple[str, str]]:
    """
    Returns a tuple of two strings representing the annotation type and
    a string of vertices in the annotation type extracted from 'line',
    a str of FDF data.
    If 'line' does not include an annotation with vertices, None is returned.
    """
    possible_annotation = parse.search("obj<</Subj({})", line)
    if possible_annotation:
        annot_type = possible_annotation[0]
        vertices = None
        if annot_type == "Line":
            vertices = parse.search("/L[{}]", line)[0]
            vertices = [float(vertex) for vertex in vertices.split(' ')]
            return (annot_type, vertices)
        elif annot_type in ("Circle", "PolyLine", "Polygon"):
            vertices = parse.search("/Vertices[{}]", line)[0]
            vertices = [float(vertex) for vertex in vertices.split(' ')]
            return (annot_type, vertices)
        elif annot_type in ("Rectangle", "Square"):
            bbox = parse.search("/Rect[{}]", line)[0]
            x1, y1, x2, y2 = bbox.split()
            vertices = " ".join([x1,y1, x1, y2, x2, y2, x2, y1])
            vertices = [float(vertex) for vertex in vertices.split(' ')]
            return (annot_type, vertices)

        
def extract_object_properties(line: str) -> Optional[dict]:
    object_properties = {}
    object_properties.update(extract_object_opacity(line))
    object_properties.update(extract_label(line))
    object_properties.update(extract_page(line))
    object_properties.update(extract_matrix(line))
    return object_properties
        
        
def extract_object_opacity(line: str) -> Optional[dict]:
    fill_opacity = parse.search("/FillOpacity {:g}/", line)
    line_opacity = parse.search("/LineOpacity {:g}/", line)
    if fill_opacity: fill_opacity = fill_opacity[0]
    if line_opacity: line_opacity = line_opacity[0]
    return {"fill_opacity": fill_opacity, "line_opacity": line_opacity}


def extract_label(line: str) -> Optional[dict]:
    label = parse.search("/Contents({})/", line)
    if label: label = label[0]
    return {"label": label}


def extract_page(line: str) -> Optional[dict]:
    page = parse.search("/Page {}", line)
    if page: 
        page = page[0]
    return {"page": page}


def extract_matrix(line: str) -> Optional[dict]:
    matrix = parse.search("/Matrix[{}]/", line)
    if matrix: 
        matrix = matrix[0]
        matrix = [float(num) for num in matrix.split()]
    return {"matrix": matrix}


def extract_stream_properties(stream_line: str) -> dict:
    """
    Returns a dict of properties which are available from the stream data:
    'line_color', 'fill_color', 'line_weight'.
    """
    line_color = parse_line_color(stream_line)
    fill_color = parse_fill_color(stream_line)
    line_weight = parse_line_weight(stream_line)
    line_type = parse_line_type(stream_line)
    return {
        "line_color": line_color, 
        "fill_color": fill_color, 
        "line_weight": line_weight,
        "line_type": line_type
    }


def parse_line_color(stream_line: str) -> tuple[int]:
    """
    Returns a tuple representing the parsed line color specification contained
    within 'stream_line', a line of text representing the FDF data stream.
    Returns None if no line color data is found in 'stream_line'
    
    The returned tuple is in the format of (R, G, B) where each of R, G, B are a float
    from 0.0 to 1.0
    """
    line_color_result = parse.search("{:g} {:g} {:g} RG", stream_line, case_sensitive=True)
    if line_color_result:
        return tuple(line_color_result)
    else:
        return (0, 0, 0)
    
    
def parse_fill_color(stream_line: str) -> tuple[int]:
    """
    Returns a tuple representing the parsed line color specification contained
    within 'stream_line', a line of text representing the FDF data stream.
    Returns None if no fill color data is found in 'stream_line'
    
    The returned tuple is in the format of (R, G, B) where each of R, G, B are a float
    from 0.0 to 1.0
    """
    fill_color_result = parse.search("{:g} {:g} {:g} rg", stream_line, case_sensitive=True)
    if fill_color_result:
        return tuple(fill_color_result)
    else:
        return (1, 1, 1)

    
def parse_line_weight(stream_line: str) -> float:
    """
    Returns a float representing the parsed line weight specification contained
    within 'stream_line', a line of text representing the FDF data stream.
    Returns None if no line weight is found in 'stream_line'
    
    The returned value represents a line weight in points (1 point = 1/72 of an inch)
    """
    line_weight_result = parse.search(" {:g} w", stream_line, case_sensitive=True)
    if line_weight_result:
        return line_weight_result[0]
    

def parse_line_type(stream_line: str) -> tuple[float, tuple]:
    """
    Returns a tuple representing the parsed line type specification contained
    within 'stream_line', a line of text representing the FDF data stream.
    Returns None if no line type data is found in 'stream_line'
    """
    line_type_data = parse.search(" [{}] {:g} d", stream_line)
    if line_type_data:
        acc = []
        for line in line_type_data[0].split(" "):
            acc.append(float(line))
        line_type = (line_type_data[1], tuple(acc))
        return line_type
    
    
def scale_object(annot: Annotation, scale: float) -> str:
    """
    Scale the annotation. Each vertex in 'annot' will be multiplied
    by 'scale'
    """
    return


# def group_vertices(vertices: str, close = False) -> list[str]:
#     """
#     Returns a list of (x, y) tuples from a string of vertices in the format of:
#     'x1 y1 x2 y2 x3 y3 ... xn yn'
#     """
#     acc = []
#     coordinates = []
#     for idx, ordinate in enumerate(vertices.split(" ")):
#         if idx % 2:
#             coordinates.append(ordinate)
#             acc.append(" ".join(coordinates))
#             coordinates = []
#         else:
#             coordinates.append(ordinate)
#     if close:
#         acc.append(acc[0])
#     return ", ".join(acc)
