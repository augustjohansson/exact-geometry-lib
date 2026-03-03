# Tests for IntersectionConstruction - converted from C++ tests + new cases


import math
import pytest
from simpex import Point, IntersectionConstruction as IC


def point_norm(p):
    return math.sqrt(p.x()**2 + p.y()**2 + p.z()**2)


def point_dist(p, q):
    d = Point(p.x() - q.x(), p.y() - q.y(), p.z() - q.z())
    return point_norm(d)


# ---------------------------------------------------------------------------
# Point-point intersections
# ---------------------------------------------------------------------------

class TestPointPoint1D:
    def test_same_point(self):
        res = IC.intersection_point_point_1d(1.0, 1.0)
        assert len(res) == 1
        assert abs(res[0] - 1.0) < 1e-15

    def test_different_points(self):
        res = IC.intersection_point_point_1d(0.0, 1.0)
        assert len(res) == 0


class TestPointPoint2D:
    def test_same_point(self):
        p = Point(1.0, 2.0, 0.0)
        res = IC.intersection_point_point_2d(p, p)
        assert len(res) == 1
        assert point_dist(res[0], p) < 1e-15

    def test_different_points(self):
        p = Point(0.0, 0.0, 0.0)
        q = Point(1.0, 0.0, 0.0)
        res = IC.intersection_point_point_2d(p, q)
        assert len(res) == 0


class TestPointPoint3D:
    def test_same_point(self):
        p = Point(1.0, 2.0, 3.0)
        res = IC.intersection_point_point_3d(p, p)
        assert len(res) == 1
        assert point_dist(res[0], p) < 1e-15

    def test_different_points(self):
        p = Point(0.0, 0.0, 0.0)
        q = Point(1.0, 1.0, 1.0)
        res = IC.intersection_point_point_3d(p, q)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Segment-point intersections
# ---------------------------------------------------------------------------

class TestSegmentPoint1D:
    def test_point_inside(self):
        res = IC.intersection_segment_point_1d(0.0, 1.0, 0.5)
        assert len(res) == 1
        assert abs(res[0] - 0.5) < 1e-15

    def test_point_at_endpoint(self):
        res = IC.intersection_segment_point_1d(0.0, 1.0, 0.0)
        assert len(res) == 1

    def test_point_outside(self):
        res = IC.intersection_segment_point_1d(0.0, 1.0, 2.0)
        assert len(res) == 0


class TestSegmentPoint2D:
    def test_point_on_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.5, 0.0, 0.0)
        res = IC.intersection_segment_point_2d(p0, p1, q0)
        assert len(res) == 1
        assert point_dist(res[0], q0) < 1e-15

    def test_point_off_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.5, 1.0, 0.0)
        res = IC.intersection_segment_point_2d(p0, p1, q0)
        assert len(res) == 0


class TestSegmentPoint3D:
    def test_point_on_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 1.0, 1.0)
        q0 = Point(0.5, 0.5, 0.5)
        res = IC.intersection_segment_point_3d(p0, p1, q0)
        assert len(res) == 1
        assert point_dist(res[0], q0) < 1e-15

    def test_point_off_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.5, 0.0, 1.0)
        res = IC.intersection_segment_point_3d(p0, p1, q0)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Triangle-point intersections
# ---------------------------------------------------------------------------

class TestTrianglePoint2D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)

    def test_point_inside(self):
        q0 = Point(0.25, 0.25, 0.0)
        res = IC.intersection_triangle_point_2d(self.p0, self.p1, self.p2, q0)
        assert len(res) == 1
        assert point_dist(res[0], q0) < 1e-15

    def test_point_outside(self):
        q0 = Point(1.0, 1.0, 0.0)
        res = IC.intersection_triangle_point_2d(self.p0, self.p1, self.p2, q0)
        assert len(res) == 0


class TestTrianglePoint3D:
    def test_point_on_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.25, 0.25, 0.0)
        res = IC.intersection_triangle_point_3d(p0, p1, p2, q0)
        assert len(res) == 1
        assert point_dist(res[0], q0) < 1e-15

    def test_point_off_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.25, 0.25, 1.0)
        res = IC.intersection_triangle_point_3d(p0, p1, p2, q0)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Tetrahedron-point intersections
# ---------------------------------------------------------------------------

class TestTetrahedronPoint3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)
        self.p3 = Point(0.0, 0.0, 1.0)

    def test_point_inside(self):
        q0 = Point(0.1, 0.1, 0.1)
        res = IC.intersection_tetrahedron_point_3d(
            self.p0, self.p1, self.p2, self.p3, q0)
        assert len(res) == 1
        assert point_dist(res[0], q0) < 1e-15

    def test_point_outside(self):
        q0 = Point(1.0, 1.0, 1.0)
        res = IC.intersection_tetrahedron_point_3d(
            self.p0, self.p1, self.p2, self.p3, q0)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Segment-segment intersections
# ---------------------------------------------------------------------------

class TestSegmentSegment1D:
    def test_overlapping(self):
        res = IC.intersection_segment_segment_1d(0.0, 2.0, 1.0, 3.0)
        assert len(res) == 2

    def test_touching_endpoint(self):
        res = IC.intersection_segment_segment_1d(0.0, 1.0, 1.0, 2.0)
        assert len(res) >= 1

    def test_non_overlapping(self):
        res = IC.intersection_segment_segment_1d(0.0, 1.0, 2.5, 3.5)
        assert len(res) == 0

    def test_contained(self):
        res = IC.intersection_segment_segment_1d(0.0, 3.0, 1.0, 2.0)
        assert len(res) == 2


class TestSegmentSegment2D:
    def test_crossing(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 1.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(0.0, 1.0, 0.0)
        res = IC.intersection_segment_segment_2d(p0, p1, q0, q1)
        assert len(res) == 1
        expected = Point(0.5, 0.5, 0.0)
        assert point_dist(res[0], expected) < 1e-14

    def test_non_crossing(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.0, 1.0, 0.0)
        q1 = Point(1.0, 1.0, 0.0)
        res = IC.intersection_segment_segment_2d(p0, p1, q0, q1)
        assert len(res) == 0

    def test_t_junction(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(2.0, 0.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(1.0, 1.0, 0.0)
        res = IC.intersection_segment_segment_2d(p0, p1, q0, q1)
        assert len(res) == 1
        expected = Point(1.0, 0.0, 0.0)
        assert point_dist(res[0], expected) < 1e-14


# ---------------------------------------------------------------------------
# Triangle-segment intersections
# ---------------------------------------------------------------------------

class TestTriangleSegment2D:
    def test_segment_crossing_edge(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.5, -0.5, 0.0)
        q1 = Point(0.5, 0.5, 0.0)
        res = IC.intersection_triangle_segment_2d(p0, p1, p2, q0, q1)
        assert len(res) >= 1

    def test_segment_inside(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.1, 0.1, 0.0)
        q1 = Point(0.2, 0.1, 0.0)
        res = IC.intersection_triangle_segment_2d(p0, p1, p2, q0, q1)
        assert len(res) == 2

    def test_segment_outside(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(2.0, 2.0, 0.0)
        q1 = Point(3.0, 2.0, 0.0)
        res = IC.intersection_triangle_segment_2d(p0, p1, p2, q0, q1)
        assert len(res) == 0


class TestTriangleSegment3D:
    def test_failing_case_from_cpp(self):
        """Converted from C++ test: failing case triangle_segment_3d"""
        p0 = Point(0.333333333333333315, 0, 0.333333333333333315)
        p1 = Point(0.66666666666666663, 0, 0.333333333333333315)
        p2 = Point(0.66666666666666663, 0.333333333333333315, 0.333333333333333315)
        q0 = Point(0.5, 0.25, 0.25)
        q1 = Point(0.75, 0.25, 0.5)
        ref_res = Point(0.583333333333333259, 0.25, 0.333333333333333315)
        res = IC.intersection_triangle_segment_3d(p0, p1, p2, q0, q1)
        assert len(res) == 1
        assert point_dist(res[0], ref_res) < 1e-15

    def test_segment_piercing_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.25, 0.25, -0.5)
        q1 = Point(0.25, 0.25, 0.5)
        res = IC.intersection_triangle_segment_3d(p0, p1, p2, q0, q1)
        assert len(res) == 1
        expected = Point(0.25, 0.25, 0.0)
        assert point_dist(res[0], expected) < 1e-14

    def test_segment_outside_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(2.0, 2.0, -1.0)
        q1 = Point(2.0, 2.0, 1.0)
        res = IC.intersection_triangle_segment_3d(p0, p1, p2, q0, q1)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Tetrahedron-segment intersections
# ---------------------------------------------------------------------------

class TestTetrahedronSegment3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)
        self.p3 = Point(0.0, 0.0, 1.0)

    def test_segment_inside(self):
        q0 = Point(0.1, 0.1, 0.1)
        q1 = Point(0.2, 0.1, 0.1)
        res = IC.intersection_tetrahedron_segment_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1)
        assert len(res) == 2

    def test_segment_outside(self):
        q0 = Point(2.0, 2.0, 2.0)
        q1 = Point(3.0, 3.0, 3.0)
        res = IC.intersection_tetrahedron_segment_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Triangle-triangle intersections
# ---------------------------------------------------------------------------

class TestTriangleTriangle2D:
    def test_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(2.0, 0.0, 0.0)
        p2 = Point(0.0, 2.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(1.0, 2.0, 0.0)
        res = IC.intersection_triangle_triangle_2d(p0, p1, p2, q0, q1, q2)
        assert len(res) >= 1

    def test_non_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(3.0, 0.0, 0.0)
        q1 = Point(4.0, 0.0, 0.0)
        q2 = Point(3.0, 1.0, 0.0)
        res = IC.intersection_triangle_triangle_2d(p0, p1, p2, q0, q1, q2)
        assert len(res) == 0


class TestTriangleTriangle3D:
    def test_coplanar_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(2.0, 0.0, 0.0)
        p2 = Point(0.0, 2.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(1.0, 2.0, 0.0)
        res = IC.intersection_triangle_triangle_3d(p0, p1, p2, q0, q1, q2)
        assert len(res) >= 1

    def test_non_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(5.0, 0.0, 0.0)
        q1 = Point(6.0, 0.0, 0.0)
        q2 = Point(5.0, 1.0, 0.0)
        res = IC.intersection_triangle_triangle_3d(p0, p1, p2, q0, q1, q2)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Tetrahedron-triangle intersections
# ---------------------------------------------------------------------------

class TestTetrahedronTriangle3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)
        self.p3 = Point(0.0, 0.0, 1.0)

    def test_triangle_inside(self):
        q0 = Point(0.1, 0.1, 0.1)
        q1 = Point(0.2, 0.1, 0.1)
        q2 = Point(0.1, 0.2, 0.1)
        res = IC.intersection_tetrahedron_triangle_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1, q2)
        assert len(res) == 3

    def test_triangle_outside(self):
        q0 = Point(2.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(2.0, 1.0, 0.0)
        res = IC.intersection_tetrahedron_triangle_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1, q2)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# Tetrahedron-tetrahedron intersections
# ---------------------------------------------------------------------------

class TestTetrahedronTetrahedron3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)
        self.p3 = Point(0.0, 0.0, 1.0)

    def test_identical(self):
        res = IC.intersection_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            self.p0, self.p1, self.p2, self.p3)
        assert len(res) == 4

    def test_non_overlapping(self):
        q0 = Point(2.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(2.0, 1.0, 0.0)
        q3 = Point(2.0, 0.0, 1.0)
        res = IC.intersection_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            q0, q1, q2, q3)
        assert len(res) == 0

    def test_partially_overlapping(self):
        q0 = Point(0.5, 0.0, 0.0)
        q1 = Point(1.5, 0.0, 0.0)
        q2 = Point(0.5, 1.0, 0.0)
        q3 = Point(0.5, 0.0, 1.0)
        res = IC.intersection_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            q0, q1, q2, q3)
        assert len(res) >= 1
