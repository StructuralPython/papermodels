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

## Beam diagrams with beamz9000

{meth}`~papermodels.datatypes.element.LoadedElement.to_beamz` returns a
horizontal member as a beam in the
[beamz9000](https://beamz9000.structuralpython.com) scheme format, so it can be
drawn as a free-body diagram and analysed. Give the pressures of the
occupancies its loading areas use:

```python
import beamz9000 as bz   # uv add "beamz9000[pycba]"

occupancies = {"office": {"D": 0.05, "L": 0.05}, "corridor": {"D": 0.05, "L": 0.1}}
beam = loaded_elements["FB2.0"].to_beamz(occupancies, results=["forces", "reactions"])
bz.draw(bz.parse_scheme({"FB2.0": beam}), "FB2.0.svg")
```

- The supports are a pin, then rollers. An end that bears on a support
  starts the beam at that support's centre; an end that runs past the
  outermost support is a cantilever.
- Each loading area becomes a `pressure` load named by its occupancy, with
  trib widths that carry its share of the area, trimmed to the beam's ends.
- Loads transferred from members above aren't included: their magnitudes are
  only known once the load path is solved, e.g. with loadtracer.
- Any other beamz9000 beam keys pass through, e.g. `e=`, `i=`,
  `node_labels="auto"` or `load_combos=[...]`.
- With `length_unit="ft"` (or `"m"`), lengths and trib widths carry that unit.
  Write the pressures with units too (`"50 psf"`) and add a forallpeople
  environment to the scheme:

  ```python
  beam = member.to_beamz({"office": {"D": "50 psf", "L": "50 psf"}}, length_unit="ft")
  scheme = bz.parse_scheme({"_environment": "structural", "FB2.0": beam})
  ```

## Tagged PDFs

To keep a record of which tag belongs to which markup, write the tags back into
the PDF:

```python
graph.export_tagged_pdf()                       # writes <name>-tagged.pdf
graph.export_tagged_pdf("markup-v2-tagged.pdf")
```
