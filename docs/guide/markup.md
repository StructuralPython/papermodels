# Marking up drawings

papermodels reads standard PDF annotations (or DXF entities, where the layer
name plays the role of annotation text). Nothing about the drawing itself is
read, only the markups you place on it.

## Pages are planes

Each page of the PDF is one plane of the structure, ordered from the **lowest
level first**. Elements on one page transfer load to elements on the page
below wherever their geometries overlap.

### Several plans on one sheet

If one sheet holds several plans, draw a **page** annotation (a rectangle with
the text `page`) around each plan, and an **origin** annotation (a small
rectangle with the text `origin`) at the same reference point on each plan,
such as the same grid intersection. Plans are ordered left to right, lowest
level on the left. Each virtual page is then treated exactly like a separate PDF
page, including how its elements are tagged.

Pages that your PDF software has rotated are un-rotated automatically.

```{note}
The word `origin` cannot appear in any other legend entry (in a Type, Rank or
Occupancy value).
```

## The legend

papermodels does not guess what a markup means. It matches each markup to a
**legend entry**: an annotation with the same appearance (line color, fill
color, line weight, dash pattern, opacity and shape) whose text starts with the
legend identifier (`Legend` by default). Legend entries need only appear on one
page of the document.

A legend entry's text is a list of `Key: Value` lines:

```text
Legend
Type: Floor Beam
Rank: 1
```

Keys are case-insensitive and spaces in them become underscores, so
`Reaction Type` becomes `reaction_type`. Any extra key you add is carried
through to the element.

### Structural elements

Structural elements need a **Type** and a **Rank**.

`Type`
: A descriptive name. Its initials become the tag prefix, so `Floor Beam` →
  `FB` and `Joist` → `J`.

`Rank`
: An integer giving the element's position in the load path. Lower ranks
  are supported by higher ranks, so joists might be rank 0, beams rank 1,
  and walls and columns a high rank such as 50. Elements only transfer load
  to elements of a higher rank.

`Reaction Type` (optional)
: `Point` (default) or `Linear`. This matters for polygons, which can
  represent both columns and walls. A column delivers a point load and a
  wall delivers a line load, and papermodels cannot tell which from an aspect
  ratio, so you state it.

Draw beams and joists as **lines** (single segments) and walls and columns as
**rectangles**.

```text
Legend          Legend
Type: Column    Type: Wall
Rank: 50        Rank: 50
                Reaction Type: Linear
```

### Loading areas

A loading area is a polygon whose legend entry has an **Occupancy** field:

```text
Legend
Occupancy: Residential
```

Loading areas can be any shape. Each one is intersected with the trib areas of
the elements below it and converted into exact distributed loads with the
[load-distribution](https://pypi.org/project/load-distribution/) package.

### Trib areas

To give an element a tributary area by hand, rather than having one generated,
draw a polygon whose legend entry has `Type: Trib Area` and place it over the
element.

```{warning}
If an element with a hand-drawn trib area has nothing transferring load to it,
it becomes a collector. Do not assign a collector behaviour to it, or the
generated trib area will overwrite yours. See
{doc}`/examples/notebooks/example-6`.
```

### Joist prototypes, extents and containers

A joist array is drawn as a single **prototype** joist. By default the array
spans the region where all of the prototype's supports exist. To control the
extent, draw one of:

- an **extent** line (`Type` containing `extent`), which sets how far the array
  runs, with the prototype's cantilevers held constant, or
- a **joist container** polygon (`Type` containing `container`), the outline the
  joists fill. Each joist is the container cut along its station line, so
  cantilevers and backspans can vary along the array.

See {class}`~papermodels.datatypes.joist_models.JoistArrayModel` and the
{doc}`/design/index` notes for how stations are placed.

## Tags

Every element gets a unique tag of the form `{INITIALS}{PAGE}.{INDEX}`, for
example `FB2.1` for the second floor beam on page 2. Joists generated from a
prototype `J2.0` are tagged `J2.0-0`, `J2.0-1`, and so on.

To fix a tag yourself, put a `tag:` line in the annotation's text:

```text
tag: FB2.5
```

Manual tags are kept, and automatic tags are numbered around them.
{meth}`GeometryGraph.export_tagged_pdf()
<papermodels.datatypes.geometry_graph.GeometryGraph.export_tagged_pdf>` writes
every assigned tag back into a copy of the PDF. Use that copy as your new markup
source and tags will stay stable as you edit. See
{doc}`/examples/notebooks/example-2`.
