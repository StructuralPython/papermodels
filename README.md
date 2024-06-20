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
- [ ] Convert the geometry graph into a load path graph

### Load distribution model
- [x] Convert arbitrary load polygon into exact distributed load
- [x] One-way arbitrary load distribution as simply-supported joists
- [ ] One-way arbitrary load distribution as continuous joists 

### Load path graphs
- [ ] Create load path graph from scratch
- [ ] Create load path graph from geometry graph
- [ ] Create individual members from scratch
- [ ] Create individual members from text file
- [ ] Save individual members to text file
- [ ] Save load path graph to text or binary file
- [ ] Reload load path graph from text or binary file

### Arbitrary load factoring
- [ ] Data model to handle creation of arbitrary load factoring approaches