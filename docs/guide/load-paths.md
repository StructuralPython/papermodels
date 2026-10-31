# How load paths are built

## Elements

Every recognised markup becomes an {class}`~papermodels.datatypes.element.Element`:
one 2D geometry plus its connections.

- **Intersections** connect elements on the *same* plane, e.g. a joist crossing
  the beam that supports it.
- **Correspondents** connect polygons on *adjacent* planes, e.g. a column on
  level 2 sitting on the column below it.

Load flows from lower rank to higher rank across intersections, and downward
across correspondents.

## The geometry graph

{class}`~papermodels.datatypes.geometry_graph.GeometryGraph` is a
{class}`networkx.DiGraph` whose nodes are element tags. An edge `A → B` means
*B supports A*. From this graph papermodels identifies:

collector elements
: Rank-0 elements with nothing above them. They collect area load directly.

transfer elements
: Elements that receive load from other elements.

orphaned elements
: Lines with fewer than two supports, or polygons with none. These are reported
  as warnings because they usually point to a markup that does not quite touch
  its support.

Plot the graph with {meth}`~papermodels.datatypes.geometry_graph.GeometryGraph.plot_connectivity`.

## Transfer rules for polygons

Walls and columns that transfer load out (rank 0) follow these rules:

1. A **point** reaction (a column) transfers to exactly one element on the
   plane: the frame element with the largest overlap, or, failing that, the
   polygon with the largest overlapping area.
2. A **linear** reaction (a wall) transfers only to polygons of the same family
   (tags starting with the same letter, so `WLL0.0` → `W0.1` but not
   `WLL0.0` → `CT0.2`) and to beams.
3. A polygon corresponds with at most one polygon of the same family on the
   plane below.

## Collector behaviour

Collectors start as single prototype elements. Calling
{meth}`~papermodels.datatypes.geometry_graph.GeometryGraph.assign_collector_behaviour`
with a collector model, typically
{class}`~papermodels.datatypes.joist_models.JoistArrayModel`, expands each
prototype into its joists and trib areas:

```python
from papermodels.datatypes.joist_models import JoistArrayModel
from papermodels.datatypes.element import create_element_filter

graph.assign_collector_behaviour(JoistArrayModel, spacing=0.4)
```

Pass a `filter_function` (see
{func}`~papermodels.datatypes.element.create_element_filter`) to apply
different models or spacings to different collectors.

Joists that would land where a support is missing are moved to the nearest valid
station, or dropped with a warning. The trib areas tile the array region with no
gaps or overlaps.

## Loaded elements

{meth}`~papermodels.datatypes.geometry_graph.GeometryGraph.create_loaded_elements`
walks the graph from top to bottom and returns a
{class}`~papermodels.datatypes.element.LoadedElement` for every element. Each
one combines:

- its supports (reaction locations),
- the point and distributed loads delivered by the elements it supports, and
- the area loads that fall on its trib area.
