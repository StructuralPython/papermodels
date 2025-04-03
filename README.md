# papermodels

papermodels is a lot of things but, ultimately, is intended to represent a class of structural engineering software that does not currently exist: **a load path designer**.

FEA has taken hold as the primary analysis tool in structural engineering yet there are many common structures that need to be designed where an adequate FE model would be considered ineffectual, inappropriate for the use case, or generally "not worth it". A typical example is a wood-framed house.

When an FE model is not used, what do engineers do? We commonly track loads by hand. When performed by hand, we make simplifications and design the structure to conform with those simplifications. However, load tracking a structure with any more than a small amount of members quickly becomes a challenge in data management.

papermodels is a software tool to "automate" the process of tracking loads by hand.

## Road map

### Pre-processing

- [x] Import PDF annotations as structural markup
- [x] Generate an inter-page geometry graph of an annotated PDF document
- [x] Plot annotations
- [x] Filter annotations by properties
- [x] Scale annotations
- [x] Categorize annotations according to a legend (e.g. green line == "beam")
- [x] Auto-tag annotations as structural members
- [x] Convert the geometry graph into a load path graph

### Load distribution model
- [x] Convert arbitrary load polygon into exact distributed load
- [x] One-way arbitrary load distribution as simply-supported joists
- [ ] One-way arbitrary load distribution as continuous joists 

### Load path graphs
- [x] Create load path graph from scratch
- [x] Create load path graph from geometry graph
- [x] Create individual members from scratch
- [ ] Create individual members from text file
- [x] Save individual members to text file
- [ ] Save load path graph to text or binary file
- [ ] Reload load path graph from text or binary file

### Arbitrary load factoring
- [ ] Data model to handle creation of arbitrary load factoring approaches

## Demo Usage

**Please see the Examples directory**

The below will parse annotations from the PDF (provided that there are legend entries setup correctly) and will create the `LoadedElement`s.

```python
from decimal import Decimal
from papermodels.paper.pdf import load_pdf_annotations
from papermodels.datatypes.geometry_graph import GeometryGraph
from papermodels.datatypes.joist_models import JoistArrayModel

graph = GeometryGraph.from_pdf_file("sketch_to_scale.pdf", scale=Decimal(1 / 72 * 4))
graph.generate_subelements(JoistArrayModel.create_subelements, spacing=1.3333)
les = graph.create_loaded_elements()

graph.plot()
```

The `LoadedElement`s can then be dumped to disk for further processing:

```python
import pathlib
output_dir = pathlib.Path("model_files")

file_format = "toml" # also try "json"

for loaded_element in les:
    if file_format == "toml":
        output_toml = output_dir / "toml"
        output_toml.mkdir(parents=True, exist_ok=True)
        with open(output_dir / "toml" / f"{loaded_element.tag}.toml", "wb") as file:
            loaded_element.dump_toml(file)
    elif file_format == "json":
        output_json = output_dir / "json"
        output_json.mkdir(parents=True, exist_ok=True)
        with open(output_dir / "json" / f"{loaded_element.tag}.json", "w") as file:
            loadeded_element.dump_json(file)
```