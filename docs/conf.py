"""Sphinx configuration for the papermodels documentation."""

import json
import pathlib
import shutil

import papermodels

DOCS_DIR = pathlib.Path(__file__).parent
REPO_ROOT = DOCS_DIR.parent

# -- Project information -----------------------------------------------------

project = "papermodels"
author = "Connor Ferster"
copyright = "2024-%Y, Connor Ferster"
release = papermodels.__version__
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_design",
    "myst_nb",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "deflist", "fieldlist"]
myst_heading_anchors = 3

# Notebooks are rendered with their saved outputs; they are not re-executed
# (they need the example PDFs next to them and produce files on disk).
nb_execution_mode = "off"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_mock_imports = ["pygraphviz"]
autoclass_content = "class"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}

# -- HTML output -------------------------------------------------------------

html_theme = "shibuya"
html_title = "papermodels"
html_baseurl = "https://papermodels.structuralpython.com/"
html_copy_source = False

html_theme_options = {
    "github_url": "https://github.com/structuralpython/papermodels",
    "accent_color": "teal",
    "nav_links": [
        {"title": "Guide", "url": "getting-started"},
        {"title": "Examples", "url": "examples/index"},
        {"title": "API", "url": "api/index"},
    ],
}

html_context = {
    "source_type": "github",
    "source_user": "structuralpython",
    "source_repo": "papermodels",
    "source_version": "main",
    "source_docs_path": "/docs/",
}

# -- Example notebooks -------------------------------------------------------
# The notebooks live in Examples/ next to their PDFs. Sphinx can only read
# files inside the docs source directory, so copy them into
# docs/examples/notebooks/ (git-ignored) and give each a page title.

EXAMPLES = [
    ("Example 1/Example_1.ipynb", "example-1", "Example 1: From sketch to load path"),
    (
        "Example 2/Example_2.ipynb",
        "example-2",
        "Example 2: Manual tags and tagged PDFs",
    ),
    ("Example 3/Example.ipynb", "example-3", "Example 3: Vertical load transfers"),
    ("Example 4/Example_3.ipynb", "example-4", "Example 4: Joist prototypes"),
    (
        "Example 4/Annotation Legend Entries.ipynb",
        "legend-entries",
        "Columns vs walls: reaction types",
    ),
    (
        "Example 5/Example_5.ipynb",
        "example-5",
        "Example 5: Multiple plans on one sheet",
    ),
    ("Example 6/Example_6.ipynb", "example-6", "Example 6: User-drawn trib areas"),
]


def _copy_example_notebooks(app):
    out_dir = DOCS_DIR / "examples" / "notebooks"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    for src, slug, title in EXAMPLES:
        src_path = REPO_ROOT / "Examples" / src
        nb = json.loads(src_path.read_text(encoding="utf-8"))
        source_link = (
            "https://github.com/structuralpython/papermodels/tree/main/Examples/"
            + src_path.parent.name.replace(" ", "%20")
        )
        nb["cells"].insert(
            0,
            {
                "cell_type": "markdown",
                "id": "papermodels-docs-title",
                "metadata": {},
                "source": [
                    f"# {title}\n",
                    "\n",
                    f"The notebook and its PDF are in [`Examples/{src_path.parent.name}`]"
                    f"({source_link}) in the repository.",
                ],
            },
        )
        (out_dir / f"{slug}.ipynb").write_text(json.dumps(nb), encoding="utf-8")


# -- Docstring conversion ----------------------------------------------------
# papermodels docstrings describe parameters as
#
#     'name': Description that may wrap onto
#         indented continuation lines.
#
# which reStructuredText reads as broken definition lists. Convert those lines
# into ``:param name:`` fields, and separate any other indented blocks with
# blank lines, so the source keeps its style and the site renders cleanly.

import inspect
import re

_QUOTED_PARAM = re.compile(
    r"^(?P<names>'[^']+'(?:\s*(?:,|and)\s*'[^']+')*)\s*(?::|\s-\s)\s*(?P<text>.*)$"
)
_BARE_PARAM = re.compile(r"^(?P<name>[A-Za-z_][\w/]*)\s*:\s*(?P<text>.*)$")


def _param_names(obj) -> set[str]:
    try:
        return set(inspect.signature(obj).parameters)
    except (TypeError, ValueError):
        return set()


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip())


def _escape(name: str) -> str:
    return name.replace("*", r"\*")


def _literalize(line: str) -> str:
    # |a -> b| is prose notation in these docstrings, not an RST substitution.
    line = re.sub(r"\|([^|]+)\|", r"``\1``", line)
    stripped = line.strip()
    if stripped.startswith("def "):
        return line[: _indent(line)] + f"``{stripped}``"
    return line


def _convert_docstring(lines: list[str], known: set[str]) -> list[str]:
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        names = None
        if _indent(line) == 0:
            m = _QUOTED_PARAM.match(line)
            if m:
                names = re.findall(r"'([^']+)'", m.group("names"))
                text = m.group("text")
            else:
                m = _BARE_PARAM.match(line)
                if m and m.group("name").split("/")[0] in known:
                    names = [m.group("name")]
                    text = m.group("text")
        if names is None:
            # Separate an indentation change from the previous line with a blank line.
            if (
                out
                and out[-1].strip()
                and line.strip()
                and _indent(line) != _indent(out[-1])
            ):
                out.append("")
            out.append(_literalize(line))
            i += 1
            continue

        # Collect the body: following lines that are blank or indented.
        body = []
        i += 1
        while i < n and (not lines[i].strip() or _indent(lines[i]) > 0):
            body.append(lines[i])
            i += 1
        while body and not body[-1].strip():
            body.pop()
        cont = min((_indent(l) for l in body if l.strip()), default=0)
        field_body = [text] if text else []
        prev_indent = cont
        for l in body:
            if not l.strip():
                field_body.append("")
                continue
            rel = _indent(l) - cont
            if field_body and field_body[-1].strip() and rel != prev_indent - cont:
                field_body.append("")
            prev_indent = _indent(l)
            field_body.append(" " * rel + _literalize(l.strip()))
        if out and out[-1].strip() and not out[-1].startswith(":param"):
            out.append("")
        for name in names:
            out.append(
                f":param {_escape(name)}: " + (field_body[0] if field_body else "")
            )
            out.extend(("   " + fb) if fb else "" for fb in field_body[1:])
        if (
            i < n
            and lines[i].strip()
            and not (_QUOTED_PARAM.match(lines[i]) or _BARE_PARAM.match(lines[i]))
        ):
            out.append("")
    return out


def _process_docstring(app, what, name, obj, options, lines):
    known = _param_names(obj) if callable(obj) else set()
    lines[:] = _convert_docstring(list(lines), known)


def _setup_docstrings(app):
    app.connect("autodoc-process-docstring", _process_docstring)


def setup(app):
    app.connect("builder-inited", _copy_example_notebooks)
    _setup_docstrings(app)
