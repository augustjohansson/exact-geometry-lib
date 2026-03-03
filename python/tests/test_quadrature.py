# Tests for SimplexQuadrature


import math
import pytest
from geometry import Point, SimplexQuadrature


# ---------------------------------------------------------------------------
# Helper: compute area/volume analytically
# ---------------------------------------------------------------------------

def triangle_area(p0, p1, p2):
    ax = p1.x() - p0.x(); ay = p1.y() - p0.y()
    bx = p2.x() - p0.x(); by = p2.y() - p0.y()
    return 0.5 * abs(ax*by - ay*bx)


def tet_volume(p0, p1, p2, p3):
    ax = p1.x()-p0.x(); ay = p1.y()-p0.y(); az = p1.z()-p0.z()
    bx = p2.x()-p0.x(); by = p2.y()-p0.y(); bz = p2.z()-p0.z()
    cx = p3.x()-p0.x(); cy = p3.y()-p0.y(); cz = p3.z()-p0.z()
    return abs(ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx)) / 6.0


# ---------------------------------------------------------------------------
# Interval quadrature (1D)
# ---------------------------------------------------------------------------

class TestSimplexQuadratureInterval:
    def test_unit_interval_weight_sum(self):
        """Weights on unit interval [0,1] must sum to its length."""
        sq = SimplexQuadrature(1, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_interval(coords, 1)
        assert abs(sum(wts) - 1.0) < 1e-14

    def test_scaled_interval_weight_sum(self):
        sq = SimplexQuadrature(1, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(3.0, 0.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_interval(coords, 1)
        assert abs(sum(wts) - 3.0) < 1e-14

    def test_exact_polynomial_integration_order2(self):
        """Integrate x over [0,1]: exact answer is 0.5."""
        sq = SimplexQuadrature(1, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_interval(coords, 1)
        n = len(wts)
        result = sum(pts[i]*wts[i] for i in range(n))
        assert abs(result - 0.5) < 1e-14

    def test_quadrature_rule_wrapper_interval(self):
        sq = SimplexQuadrature(1, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule(coords, 1)
        assert abs(sum(wts) - 1.0) < 1e-14

    def test_higher_order_interval(self):
        sq = SimplexQuadrature(1, 5)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_interval(coords, 1)
        assert abs(sum(wts) - 1.0) < 1e-14


# ---------------------------------------------------------------------------
# Triangle quadrature (2D)
# ---------------------------------------------------------------------------

class TestSimplexQuadratureTriangle:
    def test_unit_triangle_weight_sum(self):
        """Weights on unit triangle must sum to its area (0.5)."""
        sq = SimplexQuadrature(2, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
        assert abs(sum(wts) - 0.5) < 1e-14

    def test_scaled_triangle_weight_sum(self):
        sq = SimplexQuadrature(2, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(0.0, 2.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
        expected_area = triangle_area(coords[0], coords[1], coords[2])
        assert abs(sum(wts) - expected_area) < 1e-14

    def test_exact_constant_integration(self):
        """Integrate f=1 over unit triangle: result = area = 0.5."""
        sq = SimplexQuadrature(2, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
        assert abs(sum(wts) - 0.5) < 1e-14

    def test_quadrature_rule_wrapper_triangle(self):
        sq = SimplexQuadrature(2, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule(coords, 2)
        assert abs(sum(wts) - 0.5) < 1e-14

    def test_higher_order_triangle(self):
        sq = SimplexQuadrature(2, 5)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
        assert abs(sum(wts) - 0.5) < 1e-14

    def test_points_layout(self):
        """Check that pts has gdim coordinates per point."""
        sq = SimplexQuadrature(2, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
        assert len(pts) == len(wts) * 2  # gdim=2 coords per point


# ---------------------------------------------------------------------------
# Tetrahedron quadrature (3D)
# ---------------------------------------------------------------------------

class TestSimplexQuadratureTetrahedron:
    def test_unit_tet_weight_sum(self):
        """Weights on unit tet must sum to its volume (1/6)."""
        sq = SimplexQuadrature(3, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                  Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)]
        pts, wts = sq.compute_quadrature_rule_tetrahedron(coords, 3)
        assert abs(sum(wts) - 1.0/6.0) < 1e-14

    def test_scaled_tet_weight_sum(self):
        sq = SimplexQuadrature(3, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0),
                  Point(0.0, 2.0, 0.0), Point(0.0, 0.0, 2.0)]
        pts, wts = sq.compute_quadrature_rule_tetrahedron(coords, 3)
        expected_vol = tet_volume(coords[0], coords[1], coords[2], coords[3])
        assert abs(sum(wts) - expected_vol) < 1e-14

    def test_quadrature_rule_wrapper_tet(self):
        sq = SimplexQuadrature(3, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                  Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)]
        pts, wts = sq.compute_quadrature_rule(coords, 3)
        assert abs(sum(wts) - 1.0/6.0) < 1e-14

    def test_higher_order_tet(self):
        sq = SimplexQuadrature(3, 4)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                  Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)]
        pts, wts = sq.compute_quadrature_rule_tetrahedron(coords, 3)
        assert abs(sum(wts) - 1.0/6.0) < 1e-13

    def test_points_layout(self):
        """Check that pts has gdim=3 coordinates per point."""
        sq = SimplexQuadrature(3, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                  Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)]
        pts, wts = sq.compute_quadrature_rule_tetrahedron(coords, 3)
        assert len(pts) == len(wts) * 3  # gdim=3 coords per point


# ---------------------------------------------------------------------------
# Compress
# ---------------------------------------------------------------------------

class TestSimplexQuadratureCompress:
    def test_compress_returns_indices(self):
        sq = SimplexQuadrature(2, 2)
        coords = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)]
        pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
        indices = SimplexQuadrature.compress((pts, wts), 2, 2)
        # Result is either empty (no compression needed) or a subset of indices
        assert isinstance(indices, list)
