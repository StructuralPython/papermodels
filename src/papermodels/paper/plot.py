from typing import Optional, Any, Union
import pathlib

from colour import Color
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import parse

from papermodels.db.data_model import Annotation

def plot_annotations(annots: list[Annotation], size: Optional[float], dpi: Optional[float]) -> None:
    """
    Plots annotations with matplotlib
    """
    fig, ax = plt.subplots()
    for idx, annot in enumerate(annots):
        if annot.object_type == "Polygon" or annot.object_type == "Rectangle":
            xy = xy_vertices(annot.vertices)
            ax.add_patch(
                Polygon(
                    xy=xy.T,
                    closed=True,
                    linestyle=annot.line_type,
                    linewidth=annot.line_weight,
                    ec=annot.line_color,
                    fc=annot.fill_color,
                    alpha=annot.fill_opacity,
                    zorder=idx,
                )
            )
        elif annot.object_type == "Line" or annot.object_type == "PolyLine":
            xy = xy_vertices(annot.vertices)
            ax.plot(
                xy[0],
                xy[1],
                linestyle=annot.line_type,
                linewidth=annot.line_weight,
                color = annot.line_color,
                alpha = annot.line_opacity,
                zorder=idx
            )

    plt.axis("scaled")
    if size:
        fig.set_size_inches(size, size)
    if dpi:
        fig.set_dpi(dpi)
    plt.show()
    

def xy_vertices(vertices: str, close = False) -> list[list[float]]:
    """
    Returns a list of lists of floats to emulate a 2d numpy array of x, y values
    """
    x = []
    y = []
    for idx, ordinate in enumerate(vertices):
        if idx % 2:
            y.append(float(ordinate))
        else:
            x.append(float(ordinate))
    return np.asarray([x, y])