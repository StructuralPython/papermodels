# Exporting member models

Each {class}`~papermodels.datatypes.element.LoadedElement` describes one member
as a plain dictionary:

```python
member = loaded_elements["FB2.0"]
member.model()             # dict, numbers rounded to 3 decimal places
member.model(precision=4)
```

```{versionchanged} 0.8.0
`model` is a method, not a property. Call `.model()`.
```

## Writing files

{meth}`~papermodels.datatypes.element.LoadedElement.dump_json` writes a member's
model to an open text file:

```python
import pathlib

output_dir = pathlib.Path("model_files")
output_dir.mkdir(parents=True, exist_ok=True)

for tag, member in loaded_elements.items():
    with open(output_dir / f"{tag}.json", "w") as file:
        member.dump_json(file, precision=4)
```

```{note}
Earlier versions also had a `dump_toml` method. It has been removed. To write
TOML, pass `member.model()` to a
TOML writer such as [tomli-w](https://pypi.org/project/tomli-w/).
```

## Load transfer diagrams with beamz9000

{meth}`~papermodels.datatypes.element.LoadedElement.to_beamz` returns a
horizontal member as a beam in the
[beamz9000](https://beamz9000.structuralpython.com) scheme format, for a
free-body diagram of its symbolic load transfer: which members and loading
areas load it, and where. The loads have no magnitudes (each is a placeholder
of 1); each is labelled with its source instead:

```python
import beamz9000 as bz   # uv add beamz9000

beam = loaded_elements["FB3.1"].to_beamz(node_labels="auto")
bz.draw(bz.parse_scheme({"FB3.1": beam}), "FB3.1.svg")
```

- A point load from a member above (a joist of an array, `J3.0-5`, or a
  column, `CB3.0`) is a point load labelled with its tag.
- A distributed load from a member above (joists `J3.0`, a wall `WB2.0`) is a
  distributed load labelled with its tag.
- A loading area is a hatched pressure load labelled with its occupancy, one
  per stretch of the member it covers.
- The supports are a pin, then rollers. An end that bears on a support
  starts the beam at that support's centre; an end that runs past the
  outermost support is a cantilever. Loads are trimmed to the beam's ends.
- With `length_unit="ft"` (or `"m"`), the spans carry that unit; add a
  forallpeople environment to the scheme:
  `bz.parse_scheme({"_environment": "structural", "FB3.1": beam})`.
- Any other beamz9000 beam keys pass through, e.g. `node_labels="auto"` or
  `display={"legend": "top"}`.

## Tagged PDFs

To keep a record of which tag belongs to which markup, write the tags back into
the PDF:

```python
graph.export_tagged_pdf()                       # writes <name>-tagged.pdf
graph.export_tagged_pdf("markup-v2-tagged.pdf")
```
