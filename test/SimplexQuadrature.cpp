// Unit tests for SimplexQuadrature, ported from python/tests/test_quadrature.py

#include <cmath>
#include <numeric>
#include "../geometry/SimplexQuadrature.h"
#include "../geometry/Point.h"
#include "catch/catch.hpp"

using namespace simpex;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static double sum_weights(const std::pair<std::vector<double>, std::vector<double>>& qr)
{
  const auto& wts = qr.second;
  return std::accumulate(wts.begin(), wts.end(), 0.0);
}

// ---------------------------------------------------------------------------
// Interval quadrature (1D)
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: interval")
{
  SECTION("unit interval weight sum")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 1.0) < 1e-14);
  }

  SECTION("scaled interval weight sum")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(3.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 3.0) < 1e-14);
  }

  SECTION("exact polynomial integration: integrate x over [0,1] = 0.5")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    const auto& pts = qr.first;
    const auto& wts = qr.second;
    double result = 0.0;
    for (std::size_t i = 0; i < wts.size(); ++i)
      result += pts[i] * wts[i];
    CHECK(std::abs(result - 0.5) < 1e-14);
  }

  SECTION("wrapper compute_quadrature_rule for interval")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 1.0) < 1e-14);
  }

  SECTION("higher order interval")
  {
    SimplexQuadrature sq(1, 5);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 1.0) < 1e-14);
  }
}

// ---------------------------------------------------------------------------
// Triangle quadrature (2D)
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: triangle")
{
  SECTION("unit triangle weight sum")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-14);
  }

  SECTION("scaled triangle weight sum")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(0.0, 2.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    // area = 0.5 * 2 * 2 = 2.0
    CHECK(std::abs(sum_weights(qr) - 2.0) < 1e-14);
  }

  SECTION("wrapper compute_quadrature_rule for triangle")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-14);
  }

  SECTION("higher order triangle")
  {
    SimplexQuadrature sq(2, 5);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-14);
  }

  SECTION("points layout: gdim=2 coords per point")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(qr.first.size() == qr.second.size() * 2);
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron quadrature (3D)
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: tetrahedron")
{
  SECTION("unit tet weight sum")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-14);
  }

  SECTION("scaled tet weight sum")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0),
      Point(0.0, 2.0, 0.0), Point(0.0, 0.0, 2.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    // volume = (2^3)/6 = 8/6 = 4/3
    CHECK(std::abs(sum_weights(qr) - 4.0/3.0) < 1e-14);
  }

  SECTION("wrapper compute_quadrature_rule for tet")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-14);
  }

  SECTION("higher order tet")
  {
    SimplexQuadrature sq(3, 4);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-13);
  }

  SECTION("points layout: gdim=3 coords per point")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(qr.first.size() == qr.second.size() * 3);
  }
}
