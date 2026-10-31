# Getting started

## Installation

papermodels is on PyPI:

```bash
pip install papermodels
```

To plot the load path graph with Graphviz layouts, install the optional extra
(requires Graphviz on your system):

```bash
pip install "papermodels[graphviz]"
```

## The workflow

A papermodels session has four steps.

### 1. Mark up a drawing

Open a plan drawing in any PDF editor that creates standard annotations
(Bluebeam, Adobe Acrobat, PDF-XChange, and so on) and draw:

- a **legend**: one sample annotation of each kind of element, with text describing it,
- the **structural elements**: lines for joists and beams, rectangles for walls and columns,
- **loading areas** (polygons with an occupancy).

Each page is one plane (floor level) of the structure, ordered from the lowest
level up. See {doc}`guide/markup` for the full rules.

### 2. Build the graph

```python
from decimal import Decimal
from papermodels.datatypes.geometry_graph import GeometryGraph

graph = GeometryGraph.from_pdf_file(
    "my_house.pdf",
    scale=Decimal(1 / 72 * 4),  # convert PDF points to model units
)
```

PDF coordinates are in points (1/72 inch). The `scale` turns them into model
units. For standard metric drawing scales, use the constants in
{mod}`papermodels.scales`, e.g. `papermodels.scales.SCALE_1_TO_100` (millimetres)
or `SCALE_1_TO_100_M` (metres).

Inspect what was recognised:

```python
graph.plot_annotations(0)   # markups on page 0, with their tags
graph.plot_connectivity()   # the load path graph
```

### 3. Generate joists and distribute loads

Joists are drawn once, as a **prototype** line, and expanded into a full joist
array by a collector model:

```python
from papermodels.datatypes.joist_models import JoistArrayModel

graph.assign_collector_behaviour(JoistArrayModel, spacing=1.0)
loaded_elements = graph.create_loaded_elements()
```

`create_loaded_elements()` returns a dict of
{class}`~papermodels.datatypes.element.LoadedElement` keyed by tag. Each one
knows its supports, the loads delivered to it from above, and the area loads on
its trib area.

### 4. Export the member models

```python
loaded_elements["FB2.0"].model()   # a dict describing the member
```

To write every member to disk, see {doc}`guide/exporting`.

## Next steps

- {doc}`guide/markup`: how to draw markups papermodels understands.
- {doc}`guide/load-paths`: how the load path graph is built.
- {doc}`examples/index`: worked examples with PDFs.
