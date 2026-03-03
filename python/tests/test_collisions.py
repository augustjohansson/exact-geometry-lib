# Tests for CollisionPredicates - covering all collision pairs from point to tet


import pytest
from geometry import Point, CollisionPredicates as CP


# ---------------------------------------------------------------------------
# Segment-point collisions
# ---------------------------------------------------------------------------

class TestSegmentPoint1D:
    def test_point_inside_segment(self):
        assert CP.collides_segment_point_1d(0.0, 1.0, 0.5)

    def test_point_at_left_endpoint(self):
        assert CP.collides_segment_point_1d(0.0, 1.0, 0.0)

    def test_point_at_right_endpoint(self):
        assert CP.collides_segment_point_1d(0.0, 1.0, 1.0)

    def test_point_outside_left(self):
        assert not CP.collides_segment_point_1d(0.0, 1.0, -0.5)

    def test_point_outside_right(self):
        assert not CP.collides_segment_point_1d(0.0, 1.0, 1.5)


class TestSegmentPoint2D:
    def test_point_on_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        pt = Point(0.5, 0.0, 0.0)
        assert CP.collides_segment_point_2d(p0, p1, pt)

    def test_point_at_endpoint(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        assert CP.collides_segment_point_2d(p0, p1, p0)
        assert CP.collides_segment_point_2d(p0, p1, p1)

    def test_point_off_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        pt = Point(0.5, 0.1, 0.0)
        assert not CP.collides_segment_point_2d(p0, p1, pt)

    def test_point_collinear_outside(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        pt = Point(2.0, 0.0, 0.0)
        assert not CP.collides_segment_point_2d(p0, p1, pt)

    def test_diagonal_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 1.0, 0.0)
        pt = Point(0.5, 0.5, 0.0)
        assert CP.collides_segment_point_2d(p0, p1, pt)


class TestSegmentPoint3D:
    def test_point_on_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 1.0, 1.0)
        pt = Point(0.5, 0.5, 0.5)
        assert CP.collides_segment_point_3d(p0, p1, pt)

    def test_point_off_segment(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        pt = Point(0.5, 0.0, 0.1)
        assert not CP.collides_segment_point_3d(p0, p1, pt)

    def test_point_at_endpoint(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        assert CP.collides_segment_point_3d(p0, p1, p0)
        assert CP.collides_segment_point_3d(p0, p1, p1)


# ---------------------------------------------------------------------------
# Segment-segment collisions
# ---------------------------------------------------------------------------

class TestSegmentSegment1D:
    def test_overlapping(self):
        assert CP.collides_segment_segment_1d(0.0, 1.0, 0.5, 1.5)

    def test_touching_at_endpoint(self):
        assert CP.collides_segment_segment_1d(0.0, 1.0, 1.0, 2.0)

    def test_non_overlapping(self):
        # Note: the 1D check uses a tolerance of min(len_p, len_q).
        # To ensure non-overlap: need gap > max length of shorter segment.
        assert not CP.collides_segment_segment_1d(0.0, 1.0, 2.5, 3.5)

    def test_one_inside_other(self):
        assert CP.collides_segment_segment_1d(0.0, 2.0, 0.5, 1.5)

    def test_identical(self):
        assert CP.collides_segment_segment_1d(0.0, 1.0, 0.0, 1.0)


class TestSegmentSegment2D:
    def test_crossing(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 1.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(0.0, 1.0, 0.0)
        assert CP.collides_segment_segment_2d(p0, p1, q0, q1)

    def test_parallel_non_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.0, 1.0, 0.0)
        q1 = Point(1.0, 1.0, 0.0)
        assert not CP.collides_segment_segment_2d(p0, p1, q0, q1)

    def test_collinear_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(2.0, 0.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        assert CP.collides_segment_segment_2d(p0, p1, q0, q1)

    def test_non_crossing(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.5, 0.5, 0.0)
        q1 = Point(0.5, 2.0, 0.0)
        assert not CP.collides_segment_segment_2d(p0, p1, q0, q1)

    def test_t_junction(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.5, 0.0, 0.0)
        q1 = Point(0.5, 1.0, 0.0)
        assert CP.collides_segment_segment_2d(p0, p1, q0, q1)


class TestSegmentSegment3D:
    def test_crossing(self):
        p0 = Point(0.0, 0.5, 0.0)
        p1 = Point(1.0, 0.5, 0.0)
        q0 = Point(0.5, 0.0, 0.0)
        q1 = Point(0.5, 1.0, 0.0)
        assert CP.collides_segment_segment_3d(p0, p1, q0, q1)

    def test_skew(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        q0 = Point(0.0, 0.0, 1.0)
        q1 = Point(0.0, 1.0, 1.0)
        assert not CP.collides_segment_segment_3d(p0, p1, q0, q1)


# ---------------------------------------------------------------------------
# Triangle-point collisions
# ---------------------------------------------------------------------------

class TestTrianglePoint2D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)

    def test_point_inside(self):
        pt = Point(0.25, 0.25, 0.0)
        assert CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, pt)

    def test_point_on_vertex(self):
        assert CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, self.p0)
        assert CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, self.p1)
        assert CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, self.p2)

    def test_point_on_edge(self):
        pt = Point(0.5, 0.0, 0.0)
        assert CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, pt)

    def test_point_outside(self):
        pt = Point(1.0, 1.0, 0.0)
        assert not CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, pt)

    def test_point_far_outside(self):
        pt = Point(5.0, 5.0, 0.0)
        assert not CP.collides_triangle_point_2d(self.p0, self.p1, self.p2, pt)


class TestTrianglePoint3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)

    def test_point_on_triangle(self):
        pt = Point(0.25, 0.25, 0.0)
        assert CP.collides_triangle_point_3d(self.p0, self.p1, self.p2, pt)

    def test_point_above_triangle(self):
        pt = Point(0.25, 0.25, 1.0)
        assert not CP.collides_triangle_point_3d(self.p0, self.p1, self.p2, pt)

    def test_point_on_vertex(self):
        assert CP.collides_triangle_point_3d(self.p0, self.p1, self.p2, self.p0)

    def test_point_off_plane_but_inside_projection(self):
        pt = Point(0.25, 0.25, 0.5)
        assert not CP.collides_triangle_point_3d(self.p0, self.p1, self.p2, pt)


# ---------------------------------------------------------------------------
# Triangle-segment collisions
# ---------------------------------------------------------------------------

class TestTriangleSegment2D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)

    def test_segment_crossing_triangle(self):
        q0 = Point(0.5, -0.5, 0.0)
        q1 = Point(0.5, 0.5, 0.0)
        assert CP.collides_triangle_segment_2d(self.p0, self.p1, self.p2, q0, q1)

    def test_segment_inside_triangle(self):
        q0 = Point(0.1, 0.1, 0.0)
        q1 = Point(0.2, 0.2, 0.0)
        assert CP.collides_triangle_segment_2d(self.p0, self.p1, self.p2, q0, q1)

    def test_segment_outside_triangle(self):
        q0 = Point(2.0, 2.0, 0.0)
        q1 = Point(3.0, 3.0, 0.0)
        assert not CP.collides_triangle_segment_2d(self.p0, self.p1, self.p2, q0, q1)

    def test_segment_on_edge(self):
        q0 = Point(0.0, 0.0, 0.0)
        q1 = Point(0.5, 0.0, 0.0)
        assert CP.collides_triangle_segment_2d(self.p0, self.p1, self.p2, q0, q1)


class TestTriangleSegment3D:
    def test_segment_piercing_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.25, 0.25, -0.5)
        q1 = Point(0.25, 0.25, 0.5)
        assert CP.collides_triangle_segment_3d(p0, p1, p2, q0, q1)

    def test_segment_missing_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(2.0, 2.0, -0.5)
        q1 = Point(2.0, 2.0, 0.5)
        assert not CP.collides_triangle_segment_3d(p0, p1, p2, q0, q1)

    def test_segment_in_plane_of_triangle(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.1, 0.1, 0.0)
        q1 = Point(0.2, 0.1, 0.0)
        assert CP.collides_triangle_segment_3d(p0, p1, p2, q0, q1)


# ---------------------------------------------------------------------------
# Triangle-triangle collisions
# ---------------------------------------------------------------------------

class TestTriangleTriangle2D:
    def test_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.5, 0.0, 0.0)
        q1 = Point(1.5, 0.0, 0.0)
        q2 = Point(0.5, 1.0, 0.0)
        assert CP.collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2)

    def test_non_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(2.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(2.0, 1.0, 0.0)
        assert not CP.collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2)

    def test_sharing_edge(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(1.0, 0.0, 0.0)
        q1 = Point(0.0, 0.0, 0.0)
        q2 = Point(0.0, -1.0, 0.0)
        assert CP.collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2)

    def test_one_inside_other(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(4.0, 0.0, 0.0)
        p2 = Point(0.0, 4.0, 0.0)
        q0 = Point(0.5, 0.5, 0.0)
        q1 = Point(1.5, 0.5, 0.0)
        q2 = Point(0.5, 1.5, 0.0)
        assert CP.collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2)


class TestTriangleTriangle3D:
    def test_coplanar_overlapping(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(0.5, 0.0, 0.0)
        q1 = Point(1.5, 0.0, 0.0)
        q2 = Point(0.5, 1.0, 0.0)
        assert CP.collides_triangle_triangle_3d(p0, p1, p2, q0, q1, q2)

    def test_intersecting_in_3d(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.5, 1.0, 0.0)
        q0 = Point(0.5, 0.5, -0.5)
        q1 = Point(0.5, 0.5, 0.5)
        q2 = Point(0.5, -0.5, 0.0)
        assert CP.collides_triangle_triangle_3d(p0, p1, p2, q0, q1, q2)

    def test_non_colliding(self):
        p0 = Point(0.0, 0.0, 0.0)
        p1 = Point(1.0, 0.0, 0.0)
        p2 = Point(0.0, 1.0, 0.0)
        q0 = Point(5.0, 0.0, 0.0)
        q1 = Point(6.0, 0.0, 0.0)
        q2 = Point(5.0, 1.0, 0.0)
        assert not CP.collides_triangle_triangle_3d(p0, p1, p2, q0, q1, q2)


# ---------------------------------------------------------------------------
# Tetrahedron-point collision
# ---------------------------------------------------------------------------

class TestTetrahedronPoint3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)
        self.p3 = Point(0.0, 0.0, 1.0)

    def test_point_inside(self):
        pt = Point(0.1, 0.1, 0.1)
        assert CP.collides_tetrahedron_point_3d(self.p0, self.p1, self.p2, self.p3, pt)

    def test_point_at_vertex(self):
        assert CP.collides_tetrahedron_point_3d(
            self.p0, self.p1, self.p2, self.p3, self.p0)
        assert CP.collides_tetrahedron_point_3d(
            self.p0, self.p1, self.p2, self.p3, self.p3)

    def test_point_on_face(self):
        pt = Point(0.25, 0.25, 0.0)
        assert CP.collides_tetrahedron_point_3d(self.p0, self.p1, self.p2, self.p3, pt)

    def test_point_outside(self):
        pt = Point(1.0, 1.0, 1.0)
        assert not CP.collides_tetrahedron_point_3d(
            self.p0, self.p1, self.p2, self.p3, pt)

    def test_point_far_outside(self):
        pt = Point(5.0, 0.0, 0.0)
        assert not CP.collides_tetrahedron_point_3d(
            self.p0, self.p1, self.p2, self.p3, pt)


# ---------------------------------------------------------------------------
# Tetrahedron-segment collision
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
        assert CP.collides_tetrahedron_segment_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1)

    def test_segment_piercing(self):
        q0 = Point(0.1, 0.1, -0.5)
        q1 = Point(0.1, 0.1, 0.5)
        assert CP.collides_tetrahedron_segment_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1)

    def test_segment_outside(self):
        q0 = Point(2.0, 2.0, 2.0)
        q1 = Point(3.0, 3.0, 3.0)
        assert not CP.collides_tetrahedron_segment_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1)

    def test_segment_on_face(self):
        q0 = Point(0.1, 0.0, 0.0)
        q1 = Point(0.5, 0.0, 0.0)
        assert CP.collides_tetrahedron_segment_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1)


# ---------------------------------------------------------------------------
# Tetrahedron-triangle collision
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
        assert CP.collides_tetrahedron_triangle_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1, q2)

    def test_triangle_piercing(self):
        q0 = Point(0.1, 0.1, -0.5)
        q1 = Point(0.2, 0.1, -0.5)
        q2 = Point(0.1, 0.1, 0.5)
        assert CP.collides_tetrahedron_triangle_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1, q2)

    def test_triangle_outside(self):
        q0 = Point(2.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(2.0, 1.0, 0.0)
        assert not CP.collides_tetrahedron_triangle_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1, q2)

    def test_triangle_on_face(self):
        q0 = Point(0.1, 0.1, 0.0)
        q1 = Point(0.3, 0.1, 0.0)
        q2 = Point(0.1, 0.3, 0.0)
        assert CP.collides_tetrahedron_triangle_3d(
            self.p0, self.p1, self.p2, self.p3, q0, q1, q2)


# ---------------------------------------------------------------------------
# Tetrahedron-tetrahedron collision
# ---------------------------------------------------------------------------

class TestTetrahedronTetrahedron3D:
    def setup_method(self):
        self.p0 = Point(0.0, 0.0, 0.0)
        self.p1 = Point(1.0, 0.0, 0.0)
        self.p2 = Point(0.0, 1.0, 0.0)
        self.p3 = Point(0.0, 0.0, 1.0)

    def test_overlapping(self):
        q0 = Point(0.5, 0.0, 0.0)
        q1 = Point(1.5, 0.0, 0.0)
        q2 = Point(0.5, 1.0, 0.0)
        q3 = Point(0.5, 0.0, 1.0)
        assert CP.collides_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            q0, q1, q2, q3)

    def test_non_overlapping(self):
        q0 = Point(2.0, 0.0, 0.0)
        q1 = Point(3.0, 0.0, 0.0)
        q2 = Point(2.0, 1.0, 0.0)
        q3 = Point(2.0, 0.0, 1.0)
        assert not CP.collides_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            q0, q1, q2, q3)

    def test_sharing_face(self):
        # Second tet shares face p0-p1-p2 with first tet
        q0 = Point(0.0, 0.0, 0.0)
        q1 = Point(1.0, 0.0, 0.0)
        q2 = Point(0.0, 1.0, 0.0)
        q3 = Point(0.0, 0.0, -1.0)
        assert CP.collides_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            q0, q1, q2, q3)

    def test_identical(self):
        assert CP.collides_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            self.p0, self.p1, self.p2, self.p3)

    def test_one_inside_other(self):
        q0 = Point(0.1, 0.1, 0.1)
        q1 = Point(0.2, 0.1, 0.1)
        q2 = Point(0.1, 0.2, 0.1)
        q3 = Point(0.1, 0.1, 0.2)
        assert CP.collides_tetrahedron_tetrahedron_3d(
            self.p0, self.p1, self.p2, self.p3,
            q0, q1, q2, q3)
