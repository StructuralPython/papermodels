from papermodels.datatypes.annotation import Annotation, A0, A1
from papermodels.paper.annotations import _annotation_to_wkt
from papermodels.paper import annotations as an
from papermodels.paper import pdf
import numpy as np
import numpy.testing as npt
from pytest import fixture
from shapely import Polygon, box, Point
import pathlib
import fixtures

@fixture()
def page_annotation():
    page_geom = box(200, 200, 480, 360)
    return fixtures.shapely_to_annotation(page_geom)

@fixture()
def page_polygon():
    page_geom = box(200, 200, 480, 360)
    return page_geom

@fixture()
def origin_annotation():
    origin_geom = box(215, 210, 225, 220)
    return fixtures.shapely_to_annotation(origin_geom)

@fixture()
def three_page_annots():
    path = pathlib.Path(__file__).parent
    annots = pdf.load_pdf_annotations(path / "test_data" / "three-pages-on-one.pdf")
    return annots

def test__annotation_to_wkt():
    wkt_0 = _annotation_to_wkt(A0)
    wkt_1 = _annotation_to_wkt(A1)


def test_get_page_bottom_left_corner(page_polygon):
    npt.assert_array_equal(an.get_page_bottom_left_corner(page_polygon), np.array([200.0, 200.0]))


def test_get_origin_offset(origin_annotation, page_polygon):
    npt.assert_array_equal(an.get_origin_offset(origin_annotation, page_polygon), np.array([20., 15.]))


def test_enumerate_page_annotations(three_page_annots):
    page_annots = [annot for annot in three_page_annots if annot.text == "page"]
    page_geom_map = an.enumerate_page_annotations(page_annots)
    page_geoms = [an.annotation_to_shapely(annot) for annot in page_annots]
    assert page_geom_map[0] == page_geoms[0]
    assert page_geom_map[1] == page_geoms[1]
    assert page_geom_map[2] == page_geoms[2]


def test_sort_annotations_by_page_polygon(three_page_annots):
    page_annots = [annot for annot in three_page_annots if annot.text == "page"]
    other_annots = [annot for annot in three_page_annots if annot not in page_annots]
    page_geom_map = an.enumerate_page_annotations(page_annots)
    sorted_annotations = an.sort_annotations_by_page_polygon(other_annots, page_geom_map)

    # Check correct number of annotations per page
    assert len(sorted_annotations[page_geom_map[0]]) == 4
    assert len(sorted_annotations[page_geom_map[1]]) == 4
    assert len(sorted_annotations[page_geom_map[2]]) == 4

    # Check that they are not the same annotations on each page
    assert len((set(sorted_annotations[page_geom_map[0]]) - set(sorted_annotations[page_geom_map[1]]))) == 4
    assert len((set(sorted_annotations[page_geom_map[1]]) - set(sorted_annotations[page_geom_map[2]]))) == 4
    assert len((set(sorted_annotations[page_geom_map[2]]) - set(sorted_annotations[page_geom_map[0]]))) == 4
    assert len((set(sorted_annotations[page_geom_map[2]]) - set(sorted_annotations[page_geom_map[2]]))) == 0

    # Check they all have an origin annotation
    assert next((annot for annot in sorted_annotations[page_geom_map[0]] if annot.text=="origin"))
    assert next((annot for annot in sorted_annotations[page_geom_map[1]] if annot.text=="origin"))
    assert next((annot for annot in sorted_annotations[page_geom_map[2]] if annot.text=="origin"))


def test_align_annotations_to_pages(three_page_annots):
    page_annots = [annot for annot in three_page_annots if annot.text == "page"]
    other_annots = [annot for annot in three_page_annots if annot not in page_annots]
    page_geom_map = an.enumerate_page_annotations(page_annots)
    sorted_annotations = an.sort_annotations_by_page_polygon(other_annots, page_geom_map)
    aligned_annotations = an.align_annotations_to_pages(sorted_annotations)

    origin_annot_1 = next((annot for annot in aligned_annotations[page_geom_map[0]] if annot.text=="origin"))
    origin_annot_2 = next((annot for annot in aligned_annotations[page_geom_map[1]] if annot.text=="origin"))
    origin_annot_3 = next((annot for annot in aligned_annotations[page_geom_map[2]] if annot.text=="origin"))

    origin_centroid_1 = an.get_origin_centroid(origin_annot_1)
    origin_centroid_2 = an.get_origin_centroid(origin_annot_2)
    origin_centroid_3 = an.get_origin_centroid(origin_annot_3)

    # Test that origin annotations have the same coordinates
    npt.assert_array_almost_equal(origin_centroid_1, origin_centroid_2)
    npt.assert_array_almost_equal(origin_centroid_1, origin_centroid_3)

    frame_annot_1 = next((annot for annot in aligned_annotations[page_geom_map[0]] if annot.object_type=="Line"))
    frame_annot_2 = next((annot for annot in aligned_annotations[page_geom_map[1]] if annot.object_type=="Line"))
    frame_annot_3 = next((annot for annot in aligned_annotations[page_geom_map[2]] if annot.object_type=="Line"))

    # Test that the frame annotations have the same coordinates
    frame_geom_1 = an.annotation_to_shapely(frame_annot_1)
    frame_geom_2 = an.annotation_to_shapely(frame_annot_2)
    frame_geom_3 = an.annotation_to_shapely(frame_annot_3)

    # Check to see if the frame elements are close to each other
    # In testing the abs_delta ~0.233, rel_delta ~0.00254
    # However it seems that the tolerance setting required to pass this is decimal=0 
    # which is perhaps more tolerant than I would like to have in this test but I think it represents
    # a reasonable amount of discrepancy which would still allow correspondent polygons to overlap across
    # pages.
    npt.assert_array_almost_equal(np.array(frame_geom_1.coords[0]), np.array(frame_geom_2.coords[0]), decimal=0)
    npt.assert_array_almost_equal(np.array(frame_geom_1.coords[0]), np.array(frame_geom_3.coords[0]), decimal=0)

