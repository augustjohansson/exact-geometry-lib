# Tests for ConvexTriangulation - converted from C++ tests and extended


import math
import pytest
from simpex import Point, ConvexTriangulation, CollisionPredicates


def orient3d_val(a, b, c, d):
    """Compute the orientation of 4 points (sign of determinant)."""
    ax, ay, az = a.x() - d.x(), a.y() - d.y(), a.z() - d.z()
    bx, by, bz = b.x() - d.x(), b.y() - d.y(), b.z() - d.z()
    cx, cy, cz = c.x() - d.x(), c.y() - d.y(), c.z() - d.z()
    return (ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx))


def triangulation_volume_3d(triangulation):
    """Compute total volume of a 3D triangulation."""
    vol = 0.0
    for tet in triangulation:
        x0, x1, x2, x3 = tet[0], tet[1], tet[2], tet[3]
        v = (x0.x() * (x1.y()*x2.z() + x3.y()*x1.z() + x2.y()*x3.z()
                       - x2.y()*x1.z() - x1.y()*x3.z() - x3.y()*x2.z())
             - x1.x() * (x0.y()*x2.z() + x3.y()*x0.z() + x2.y()*x3.z()
                         - x2.y()*x0.z() - x0.y()*x3.z() - x3.y()*x2.z())
             + x2.x() * (x0.y()*x1.z() + x3.y()*x0.z() + x1.y()*x3.z()
                         - x1.y()*x0.z() - x0.y()*x3.z() - x3.y()*x1.z())
             - x3.x() * (x0.y()*x1.z() + x1.y()*x2.z() + x2.y()*x0.z()
                         - x1.y()*x0.z() - x2.y()*x1.z() - x0.y()*x2.z()))
        vol += abs(v)
    return vol / 6.0


def pure_triangular(triangulation, dim):
    """Check all simplices have dim+1 vertices."""
    return all(len(s) == dim + 1 for s in triangulation)


def triangulation_selfintersects_3d(triangulation):
    """
    Check if a 3D triangulation self-intersects.
    Mirrors the local helper from the C++ test (respects shared vertices).
    """
    from simpex import CollisionPredicates
    dim = 3
    for i in range(len(triangulation)):
        t1 = triangulation[i]
        for j in range(i + 1, len(triangulation)):
            t2 = triangulation[j]
            # Count shared vertices
            t1_shared = set()
            t2_shared = set()
            shared = 0
            for v1 in range(dim + 1):
                for v2 in range(dim + 1):
                    if t1[v1] == t2[v2]:
                        shared += 1
                        t1_shared.add(v1)
                        t2_shared.add(v2)
            if shared == 0:
                if CollisionPredicates.collides_tetrahedron_tetrahedron_3d(
                        t1[0], t1[1], t1[2], t1[3],
                        t2[0], t2[1], t2[2], t2[3]):
                    return True
            else:
                for a in range(dim + 1):
                    if a not in t1_shared:
                        if CollisionPredicates.collides_tetrahedron_point_3d(
                                t2[0], t2[1], t2[2], t2[3], t1[a]):
                            return True
                    if a not in t2_shared:
                        if CollisionPredicates.collides_tetrahedron_point_3d(
                                t1[0], t1[1], t1[2], t1[3], t2[a]):
                            return True
    return False

def has_degenerate_3d(triangulation):
    """Check for degenerate (zero-volume) tetrahedra."""
    for tet in triangulation:
        if abs(orient3d_val(tet[0], tet[1], tet[2], tet[3])) == 0:
            return True
    return False

class TestConvexTriangulationGrahamScan3D:
    def test_trivial_single_tet(self):
        """Converted from C++: test trivial case (single tetrahedron)."""
        pts = [Point(0, 0, 0), Point(0, 0, 1),
               Point(0, 1, 0), Point(1, 0, 0)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert len(tri) == 1
        assert abs(triangulation_volume_3d(tri) - 1.0/6.0) < 1e-14

    def test_trivial_cube(self):
        """Converted from C++: test trivial case 2 (unit cube)."""
        pts = [Point(0, 0, 0), Point(0, 0, 1), Point(0, 1, 0), Point(0, 1, 1),
               Point(1, 0, 0), Point(1, 0, 1), Point(1, 1, 0), Point(1, 1, 1)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert pure_triangular(tri, 3)
        assert not has_degenerate_3d(tri)
        assert not triangulation_selfintersects_3d(tri)
        assert abs(triangulation_volume_3d(tri) - 1.0) < 1e-14

    def test_coplanar_points(self):
        """Converted from C++: test coplanar points."""
        pts = [Point(0,   0,   0), Point(0,   0,   1), Point(0,   1,   0),
               Point(0,   1,   1), Point(1,   0,   0), Point(1,   0,   1),
               Point(1,   1,   0), Point(1,   1,   1), Point(0.1, 0.1, 0)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert pure_triangular(tri, 3)
        assert not has_degenerate_3d(tri)
        assert not triangulation_selfintersects_3d(tri)
        assert abs(triangulation_volume_3d(tri) - 1.0) < 1e-14

    def test_coplanar_colinear_points(self):
        """Converted from C++: test coplanar colinear points."""
        pts = [Point(0, 0,   0), Point(0, 0,   1), Point(0, 1,   0),
               Point(0, 1,   1), Point(1, 0,   0), Point(1, 0,   1),
               Point(1, 1,   0), Point(1, 1,   1), Point(0, 0.1, 0)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert pure_triangular(tri, 3)
        assert not has_degenerate_3d(tri)
        assert not triangulation_selfintersects_3d(tri)
        assert abs(triangulation_volume_3d(tri) - 1.0) < 1e-14

    def test_failing_case_1(self):
        """Converted from C++: test failing case."""
        pts = [Point(0.7, 0.6, 0.1),
               Point(0.7, 0.6, 0.5),
               Point(0.1, 0.1, 0.1),
               Point(0.8333333333333333, 0.8333333333333333, 0),
               Point(0.1, 0.15, 0.1),
               Point(0.1, 0.45, 0.1),
               Point(0.16, 0.15, 0.1),
               Point(0.61, 0.525, 0.1),
               Point(0.46, 0.6, 0.100000000000000006)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert pure_triangular(tri, 3)
        assert not has_degenerate_3d(tri)
        assert not triangulation_selfintersects_3d(tri)

    def test_failing_case_2(self):
        """Converted from C++: test failing case 2."""
        pts = [Point(0.7, 0.6, 0.5),
               Point(0.7, 0.1, 0.1),
               Point(0.8, 0, 0),
               Point(0.1, 0.1, 0.1),
               Point(0.16, 0.1, 0.1),
               Point(0.592, 0.1, 0.1),
               Point(0.16, 0.1, 0.14),
               Point(0.52, 0.1, 0.38),
               Point(0.7, 0.1, 0.38)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert pure_triangular(tri, 3)
        assert not has_degenerate_3d(tri)
        assert not triangulation_selfintersects_3d(tri)

    def test_failing_case_3(self):
        """Converted from C++: test failing case 3."""
        pts = [Point(0.495926, 0.512037, 0.144444),
               Point(0.376482, 0.519121, 0.284321),
               Point(0.386541, 0.599783, 0.0609262),
               Point(0.388086, 0.60059, 0.0607155),
               Point(0.7, 0.6, 0.5),
               Point(0.504965, 0.504965, 0.0447775),
               Point(0.833333, 0.833333, 0)]
        tri = ConvexTriangulation.triangulate_graham_scan_3d(pts)
        assert pure_triangular(tri, 3)
        assert not has_degenerate_3d(tri)
        assert not triangulation_selfintersects_3d(tri)


# ---------------------------------------------------------------------------
# 2D triangulation tests (Graham scan 2D)
# ---------------------------------------------------------------------------

class TestConvexTriangulationGrahamScan2D:
    def test_single_triangle(self):
        pts = [Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0)]
        tri = ConvexTriangulation.triangulate_graham_scan_2d(pts)
        assert len(tri) == 1
        assert len(tri[0]) == 3

    def test_square(self):
        pts = [Point(0, 0, 0), Point(1, 0, 0),
               Point(1, 1, 0), Point(0, 1, 0)]
        tri = ConvexTriangulation.triangulate_graham_scan_2d(pts)
        assert len(tri) == 2
        assert all(len(s) == 3 for s in tri)


# ---------------------------------------------------------------------------
# Generic triangulate() wrapper
# ---------------------------------------------------------------------------

class TestConvexTriangulationWrapper:
    def test_triangulate_1d(self):
        pts = [Point(0, 0, 0), Point(1, 0, 0)]
        tri = ConvexTriangulation.triangulate(pts, 2, 1)
        assert len(tri) == 1
        assert len(tri[0]) == 2

    def test_triangulate_2d(self):
        pts = [Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0)]
        tri = ConvexTriangulation.triangulate(pts, 2, 2)
        assert len(tri) >= 1
        assert all(len(s) == 3 for s in tri)

    def test_triangulate_3d(self):
        pts = [Point(0, 0, 0), Point(1, 0, 0),
               Point(0, 1, 0), Point(0, 0, 1)]
        tri = ConvexTriangulation.triangulate(pts, 3, 3)
        assert len(tri) == 1
        assert len(tri[0]) == 4
