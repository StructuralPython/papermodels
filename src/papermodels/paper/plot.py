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
import textalloc as ta
from ..datatypes.element import Element


def plot_elements(
    elements: list[Element],
    figsize: int | float | tuple[int | float, int | float] = (17, 11),
    dpi: float = 100,
    plot_trib_areas: bool = False,
    plot_extent_polygons: bool = False,
    plot_tags: bool = False,
) -> Figure:
    """
    Plots the elements in matplotlib. Size and dpi can be adjusted
    to make the plot bigger/smaller. Size is in inches and dpi stands for
    "dots per inch". For a biggish plot, values of size=12, dpi=200 gives
    good results.

    Colors are based on the kind of element that is being plotted
    """
    if isinstance(figsize, (int, float)):
        figsize = (figsize, figsize)

    fig = Figure(figsize=figsize, dpi=dpi, tight_layout=True)
    ax = fig.gca()
    has_tags = False
    min_extent = np.array([float("-inf"), float("-inf")])
    max_extent = np.array([float("inf"), float("inf")])
    text_annotations = []
    initial_positions_x = []
    initial_positions_y = []
    tags = []
    lines_x = []
    lines_y = []
    linear_polygons = []
    point_polygons = []
    extent_polygons = []
    trib_areas = []
    collectors = []
    transfers = []
    plot_colors = {
        "poly_edge": "k",
        "poly_face": (0.5, 0.5, 0.5),
        "linear_poly": "green",
        "collector_lines": (0.7, 0.7, 0.7),
        "transfer_lines": (0.3, 0.3, 0.3),
        "trib_edge": "red",
        "trib_face": (0.8, 0.1, 0.0),
        "extent_edge": "yellow",
        "extent_face": (0.8, 0.7, 0.0),
        "prototype_line": "yellow",
    }
    plot_objs = []
    for element in elements:
        plot_attrs = get_element_plotting_attributes(element, is_subelement=True)
        plot_objs.append(plot_attrs)
        if element.subelements:
            for subelem in element.subelements:
                plot_attrs = get_element_plotting_attributes(
                    subelem, is_subelement=True
                )
                plot_objs.append(plot_attrs)

    for po in plot_objs:
        xy = po["xy"]
        if sum(max_extent) == float("inf"):
            min_extent = np.maximum(min_extent, np.max(xy, axis=1))
            max_extent = np.minimum(max_extent, np.min(xy, axis=1))
        else:
            min_extent = np.minimum(min_extent, np.min(xy, axis=1))
            max_extent = np.maximum(max_extent, np.max(xy, axis=1))

        # For tagging
        initial_positions_x.append(po["anchor_point"][0])
        initial_positions_y.append(po["anchor_point"][1])
        tags.append(po["tag"])

        if po["is_poly"]:
            ax.add_patch(
                Polygon(
                    xy=xy.T,
                    closed=True,
                    # linestyle=annot.line_type,
                    # linewidth=float(annot.line_weight),
                    ec=plot_colors["poly_edge"],
                    fc=plot_colors["poly_face"],
                    alpha=0.6,
                    # zorder=idx,
                )
            )
        elif po["is_line"]:
            lines_x.append(xy[0])
            lines_y.append(xy[1])
            if po["rank"] == 0:
                line_color = plot_colors["collector_lines"]
                z_order = None
            else:
                line_color = plot_colors["transfer_lines"]
                z_order = None
            if po["extent_polygon_xy"] is not None:
                line_color = plot_colors["prototype_line"]
                z_order = 1e6
            ax.plot(
                xy[0],
                xy[1],
                # linestyle=annot.line_type,
                # linewidth=float(annot.line_weight),
                color=line_color,
                alpha=1.0,
                zorder=z_order,
            )

        if po["trib_area_xy"] is not None and plot_trib_areas:
            ax.add_patch(
                Polygon(
                    xy=po["trib_area_xy"].T,
                    closed=True,
                    # linestyle=annot.line_type,
                    # linewidth=float(annot.line_weight),
                    ec=plot_colors["trib_edge"],
                    fc=plot_colors["trib_face"],
                    alpha=0.3,
                    zorder=1000,
                )
            )

        elif po["extent_polygon_xy"] is not None and plot_extent_polygons:
            ax.add_patch(
                Polygon(
                    xy=po["extent_polygon_xy"].T,
                    closed=True,
                    # linestyle=annot.line_type,
                    # linewidth=float(annot.line_weight),
                    ec=plot_colors["extent_edge"],
                    fc=plot_colors["extent_face"],
                    alpha=0.3,
                    zorder=1000,
                )
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
    ta.allocate(
        ax=ax,
        x=initial_positions_x,
        y=initial_positions_y,
        text_list=tags,
        # x_lines=lines_x,
        # y_lines=lines_y,
        textsize=8,
        textcolor="k",
        linecolor="k",
        avoid_label_lines_overlap=True,
        avoid_crossing_label_lines=True,
    )
    return fig


def get_element_plotting_attributes(
    element: Element, is_subelement: bool = False
) -> dict[str, Any]:
    """
    Returns a dict of plotting attributes for the element.

    e.g.
    {
        "is_poly": True,
        "is_line": False,
        "is_linear": False,
        "is_subelement": False,
        "tag": str,
        "rank": int,
        "xy": np.array([[...], [...]]),
        "anchor_point": np.array([..., ...]),
        "extent_polygon_xy": np.array([[...], [...]]) | None,
        "trib_area_xy": np.array([[...], [...]]) | None,

    }
    """
    is_poly = False
    is_line = False
    is_linear = False
    if element.geometry.geom_type == "Polygon":
        coords = element.geometry.exterior.coords
        is_poly = True
    elif element.geometry.geom_type == "LineString":
        coords = element.geometry.coords
        is_line = True

    trib_xy = None
    extent_xy = None
    if element.reaction_type == "linear":
        is_linear = True
    if element.trib_area:
        trib_coords = element.trib_area.exterior.coords
        trib_xy = np.array(list(zip(*trib_coords)))
    if element.extent_polygon:
        extent_coords = element.extent_polygon.exterior.coords
        extent_xy = np.array(list(zip(*extent_coords)))

    xy = np.array(list(zip(*coords)))
    tag = element.tag
    rank = element.rank
    geom = element.geometry
    rep_point = np.array(geom.representative_point().coords[0])
    centroid_point = np.array(geom.centroid.coords[0])
    anchor_point = (rep_point + centroid_point) / 2
    plotting_attributes = {
        "is_poly": is_poly,
        "is_line": is_line,
        "is_linear": is_linear,
        "is_subelement": is_subelement,
        "tag": tag,
        "rank": rank,
        "xy": xy,
        "anchor_point": anchor_point,
        "trib_area_xy": trib_xy,
        "extent_polygon_xy": extent_xy,
    }
    return plotting_attributes


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
    text_annotations = []
    initial_positions_x = []
    initial_positions_y = []
    lines_x = []
    lines_y = []
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
            geom = annots[annot]["geometry"]
            minx, miny, maxx, maxy = geom.bounds
            lines_x.append([minx * dpi / 72, maxx * dpi / 72])
            lines_y.append([miny * dpi / 72, maxy * dpi / 72])

            xy = xy_vertices(annot.vertices, dpi)
            if sum(max_extent) == float("inf"):
                min_extent = np.maximum(max_extent, np.max(xy, axis=1))
                max_extent = np.minimum(min_extent, np.min(xy, axis=1))
            else:
                min_extent = np.minimum(min_extent, np.min(xy, axis=1))
                max_extent = np.maximum(max_extent, np.max(xy, axis=1))
            face_color = (
                tuple(float(elem) for elem in annot.fill_color)
                if annot.fill_color is not None
                else None
            )
            ax.add_patch(
                Polygon(
                    xy=xy.T,
                    closed=True,
                    # linestyle=annot.line_type,
                    linewidth=float(annot.line_weight),
                    ec=tuple(float(elem) for elem in annot.line_color),
                    fc=face_color,
                    alpha=0.2 * float(annot.fill_opacity),
                    zorder=idx,
                )
            )
        elif annot.object_type.lower() in ("line", "polyline"):
            xy = xy_vertices(annot.vertices, dpi)
            lines_x.append(xy[0])
            lines_y.append(xy[1])
            ax.plot(
                xy[0],
                xy[1],
                # linestyle=annot.line_type,
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
            initial_positions_x.append(plot_point[0] * dpi / 72)
            initial_positions_y.append(plot_point[1] * dpi / 72)
            text_annotations.append(tag)
            # ax.annotate(
            #     tag,
            #     plot_point * dpi / 72,
            #     zorder=100 * len(annots),
            # )

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
    ta.allocate(
        ax=ax,
        x=initial_positions_x,
        y=initial_positions_y,
        text_list=text_annotations,
        x_lines=lines_x,
        y_lines=lines_y,
        textsize=8,
        textcolor="k",
        linecolor="k",
        avoid_label_lines_overlap=True,
        avoid_crossing_label_lines=True,
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
