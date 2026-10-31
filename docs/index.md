# papermodels

**A load path designer for structural engineers.**

FEA is the primary analysis tool in structural engineering, yet many common
structures, such as a wood-framed house, are not worth a finite element model.
For those, engineers track loads by hand: simplify the structure, then follow
each load from where it lands down to the foundation. Past a handful of
members, that becomes a data-management problem.

papermodels automates hand load-tracking. You mark up a PDF (or DXF) drawing
with ordinary annotations, and papermodels:

1. reads the markups and classifies them against a legend you draw on the page,
2. builds a directed **load path graph** of how every member supports every other,
   across floors,
3. generates joist arrays and tributary areas,
4. distributes area loads onto each member, and
5. exports each member as a structural model (JSON) for your design tools.

```python
from decimal import Decimal
from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.joist_models import JoistArrayModel

graph = GeometryGraph.from_pdf_file("sketch_to_scale.pdf", scale=Decimal(1 / 72 * 4))
graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
loaded_elements = graph.create_loaded_elements()

loaded_elements["FB2.0"].model()
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Getting started
:link: getting-started
:link-type: doc
Install papermodels and run your first marked-up drawing.
:::

:::{grid-item-card} Marking up drawings
:link: guide/markup
:link-type: doc
Legends, ranks, reaction types, tags, trib areas and loads.
:::

:::{grid-item-card} Examples
:link: examples/index
:link-type: doc
Worked notebooks, each with its marked-up PDF.
:::

:::{grid-item-card} API reference
:link: api/index
:link-type: doc
Every module, class and function.
:::
::::

```{toctree}
:hidden:
:caption: Guide

getting-started
guide/markup
guide/load-paths
guide/exporting
```

```{toctree}
:hidden:
:caption: Examples

examples/index
```

```{toctree}
:hidden:
:caption: Reference

api/index
design/index
```
