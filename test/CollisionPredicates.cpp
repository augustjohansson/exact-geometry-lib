// Unit tests for CollisionPredicates, ported from python/tests/test_collisions.py

#include "../geometry/CollisionPredicates.h"
#include "../geometry/Point.h"
#include "catch/catch.hpp"

using namespace simpex;

// ---------------------------------------------------------------------------
// Segment-point collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: segment-point 1D")
{
  SECTION("point inside segment")
  {
    CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 0.5));
  }
  SECTION("point at left endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 0.0));
  }
  SECTION("point at right endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 1.0));
  }
  SECTION("point outside left")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, -0.5));
  }
  SECTION("point outside right")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 1.5));
  }
}

TEST_CASE("CollisionPredicates: segment-point 2D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);

  SECTION("point on segment")
  {
    CHECK(CollisionPredicates::collides_segment_point_2d(p0, p1, Point(0.5, 0.0, 0.0)));
  }
  SECTION("point at endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_2d(p0, p1, p0));
    CHECK(CollisionPredicates::collides_segment_point_2d(p0, p1, p1));
  }
  SECTION("point off segment")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_2d(p0, p1, Point(0.5, 0.1, 0.0)));
  }
  SECTION("point collinear outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_2d(p0, p1, Point(2.0, 0.0, 0.0)));
  }
  SECTION("diagonal segment midpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 0.0), Point(0.5, 0.5, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: segment-point 3D")
{
  SECTION("point on segment")
  {
    CHECK(CollisionPredicates::collides_segment_point_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), Point(0.5, 0.5, 0.5)));
  }
  SECTION("point off segment")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.5, 0.0, 0.1)));
  }
  SECTION("point at endpoint")
  {
    const Point p0(0.0, 0.0, 0.0);
    const Point p1(1.0, 0.0, 0.0);
    CHECK(CollisionPredicates::collides_segment_point_3d(p0, p1, p0));
    CHECK(CollisionPredicates::collides_segment_point_3d(p0, p1, p1));
  }
}

// ---------------------------------------------------------------------------
// Segment-segment collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: segment-segment 1D")
{
  SECTION("overlapping")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 0.5, 1.5));
  }
  SECTION("touching at endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 1.0, 2.0));
  }
  SECTION("non-overlapping")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 2.5, 3.5));
  }
  SECTION("one inside other")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 2.0, 0.5, 1.5));
  }
  SECTION("identical")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 0.0, 1.0));
  }
}

TEST_CASE("CollisionPredicates: segment-segment 2D")
{
  SECTION("crossing")
  {
    CHECK(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 0.0),
      Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)));
  }
  SECTION("parallel non-overlapping")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(1.0, 1.0, 0.0)));
  }
  SECTION("collinear overlapping")
  {
    CHECK(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0),
      Point(1.0, 0.0, 0.0), Point(3.0, 0.0, 0.0)));
  }
  SECTION("T-junction")
  {
    CHECK(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.5, 0.0, 0.0), Point(0.5, 1.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: segment-segment 3D")
{
  SECTION("crossing in same plane")
  {
    CHECK(CollisionPredicates::collides_segment_segment_3d(
      Point(0.0, 0.5, 0.0), Point(1.0, 0.5, 0.0),
      Point(0.5, 0.0, 0.0), Point(0.5, 1.0, 0.0)));
  }
  SECTION("skew segments")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_segment_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 0.0, 1.0), Point(0.0, 1.0, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Triangle-point collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: triangle-point 2D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);

  SECTION("point inside")
  {
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(0.25, 0.25, 0.0)));
  }
  SECTION("point on vertex")
  {
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, p0));
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, p1));
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, p2));
  }
  SECTION("point on edge")
  {
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(0.5, 0.0, 0.0)));
  }
  SECTION("point outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(1.0, 1.0, 0.0)));
  }
  SECTION("point far outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(5.0, 5.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: triangle-point 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);

  SECTION("point on triangle plane")
  {
    CHECK(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, Point(0.25, 0.25, 0.0)));
  }
  SECTION("point above triangle")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, Point(0.25, 0.25, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Triangle-triangle collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: triangle-triangle 2D")
{
  SECTION("overlapping triangles")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(0.25, 0.25, 0.0), Point(1.25, 0.25, 0.0), Point(0.25, 1.25, 0.0)));
  }
  SECTION("non-overlapping triangles")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_triangle_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(2.0, 0.0, 0.0), Point(3.0, 0.0, 0.0), Point(2.0, 1.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: triangle-triangle 3D")
{
  SECTION("intersecting triangles")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(0.25, 0.25, -0.5), Point(0.25, 0.25, 0.5), Point(0.5, 0.5, 0.0)));
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron-point collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: tetrahedron-point 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);
  const Point p3(0.0, 0.0, 1.0);

  SECTION("point inside tetrahedron")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(
      p0, p1, p2, p3, Point(0.1, 0.1, 0.1)));
  }
  SECTION("point at vertex")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(p0, p1, p2, p3, p0));
  }
  SECTION("point outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_point_3d(
      p0, p1, p2, p3, Point(1.0, 1.0, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron-tetrahedron collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: tetrahedron-tetrahedron 3D")
{
  SECTION("overlapping tetrahedra")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0),
      Point(0.1, 0.1, 0.1), Point(1.1, 0.1, 0.1),
      Point(0.1, 1.1, 0.1), Point(0.1, 0.1, 1.1)));
  }
  SECTION("non-overlapping tetrahedra")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0),
      Point(3.0, 0.0, 0.0), Point(4.0, 0.0, 0.0),
      Point(3.0, 1.0, 0.0), Point(3.0, 0.0, 1.0)));
  }
}
