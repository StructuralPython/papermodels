from __future__ import annotations
from typing import Optional, Any, Union
import pathlib

from colour import Color
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
import numpy as np
import parse
from ..datatypes.annotation import Annotation
from shapely.ops import polylabel

def plot_annotations(
    annots: list[Annotation] | dict[Annotation, dict],
    figsize: int | float | tuple[int | float, int | float] = (17, 11),
    dpi: float = 100,
    plot_tags: bool = False,
) -> Figure:
    """
    Plots that annotations, 'annots' in matplotlib. Size and dpi can be adjusted
    to make the plot bigger/smaller. Size is in inches and dpi stands for
    "dots per inch". For a biggish plot, values of size=12, dpi=200 gives
    good results.
    """
    if isinstance(figsize, (int, float)):
        figsize = (figsize, figsize)

    fig = Figure(figsize=figsize, dpi=dpi, tight_layout=True)
    ax = fig.gca()
    annotation_dict = isinstance(annots, dict)
    has_tags = False
    min_extent = np.array([float("-inf"), float("-inf")])
    max_extent = np.array([float("inf"), float("inf")])
    for idx, annot in enumerate(annots):
        if annotation_dict:
            has_tags = "tag" in annots[annot]

        xy = xy_vertices(annot.vertices, dpi)
        if sum(max_extent) == float("inf"):
            min_extent = np.maximum(min_extent, np.max(xy, axis=1))
            max_extent = np.minimum(max_extent, np.min(xy, axis=1))
        else:
            min_extent = np.minimum(min_extent, np.min(xy, axis=1))
            max_extent = np.maximum(max_extent, np.max(xy, axis=1))

        if annot.object_type.lower() in (
            "polygon",
            "square",
            "rectangle",
            "rectangle sketch to scale",
        ):
            xy = xy_vertices(annot.vertices, dpi)
            if sum(max_extent) == float("inf"):
                min_extent = np.maximum(max_extent, np.max(xy, axis=1))
                max_extent = np.minimum(min_extent, np.min(xy, axis=1))
            else:
                min_extent = np.minimum(min_extent, np.min(xy, axis=1))
                max_extent = np.maximum(max_extent, np.max(xy, axis=1))
            ax.add_patch(
                Polygon(
                    xy=xy.T,
                    closed=True,
                    linestyle=annot.line_type,
                    linewidth=float(annot.line_weight),
                    ec=tuple(float(elem) for elem in annot.line_color),
                    fc=tuple(float(elem) for elem in annot.fill_color),
                    alpha=float(annot.fill_opacity),
                    zorder=idx,
                )
            )
        elif annot.object_type.lower() in ("line", "polyline"):
            xy = xy_vertices(annot.vertices, dpi)
            ax.plot(
                xy[0],
                xy[1],
                linestyle=annot.line_type,
                linewidth=float(annot.line_weight),
                color=tuple(float(elem) for elem in annot.line_color),
                alpha=float(annot.line_opacity),
                zorder=idx,
            )
        if annotation_dict and has_tags and plot_tags:
            tag = annots[annot]["tag"]
            geom = annots[annot]["geometry"]
            rep_point = np.array(geom.representative_point().coords[0])
            centroid_point = np.array(geom.centroid.coords[0])
            plot_point = (rep_point + centroid_point) / 2
            ax.annotate(
                tag,
                plot_point * dpi / 72,
                zorder=100 * len(annots),
            )

    ax.set_aspect("equal")
    plot_margin_metric = np.linalg.norm(
        max_extent - min_extent
    )  # Distance between bot-left and top-right
    ax.set_xlim(
        min_extent[0] - plot_margin_metric * 0.05,
        max_extent[0] + plot_margin_metric * 0.05,
    )
    ax.set_ylim(
        min_extent[1] - plot_margin_metric * 0.05,
        max_extent[1] + plot_margin_metric * 0.05,
    )
    return fig


def xy_vertices(vertices: str, dpi: float, close=False) -> list[list[float]]:
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
    scaled_vertices = np.asarray([x, y]) * dpi / 72
    return scaled_vertices
