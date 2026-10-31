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

## Tagged PDFs

To keep a record of which tag belongs to which markup, write the tags back into
the PDF:

```python
graph.export_tagged_pdf()                       # writes <name>-tagged.pdf
graph.export_tagged_pdf("markup-v2-tagged.pdf")
```
