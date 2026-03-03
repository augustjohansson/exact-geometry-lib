// CGALComparison.cpp
//
// Catch2 tests that verify CGAL (EPICK) and Shewchuk predicates give
// identical results on deterministic random inputs.
//
// Newly added tests to notify the user of:
//   * "CGAL vs Shewchuk: orient2d sign agreement"        (10 000 random 2D points)
//   * "CGAL vs Shewchuk: orient3d sign agreement"        (10 000 random 3D points)
//   * "CGAL vs Shewchuk: collides_triangle_point_3d"     (10 000 random pairs)
//   * "CGAL vs Shewchuk: collides_tetrahedron_point_3d"  (5 000 random pairs)
//   * "CGAL vs Shewchuk: collides_triangle_triangle_3d"  (random pairs + canonical cases)
//   * "CGAL vs Shewchuk: collides_tetrahedron_tetrahedron_3d" (random pairs + canonical cases)
//
// All random generators use seed 42 for reproducibility.

#ifdef GEOMETRY_WITH_CGAL

#include <array>
#include <cmath>
#include <random>

#include "../geometry/CGALPredicates.h"
#include "../geometry/CollisionPredicates.h"
#include "../geometry/predicates.h"
#include "../geometry/Point.h"
#include "catch/catch.hpp"

using namespace simpex;

// ---------------------------------------------------------------------------
// Deterministic RNG — seed 42
// ---------------------------------------------------------------------------
static std::mt19937 rng_cmp(42);
static std::uniform_real_distribution<double> dist_cmp(0.0, 1.0);

static double r01() { return dist_cmp(rng_cmp); }

static Point rpt3() { return Point(r01(), r01(), r01()); }
static Point rpt2() { return Point(r01(), r01(), 0.0); }

static std::array<Point, 3> rtri3d()
{
  while (true)
  {
    Point p0 = rpt3(), p1 = rpt3(), p2 = rpt3();
    if ((p1 - p0).cross(p2 - p0).norm() > 1e-8)
      return {p0, p1, p2};
  }
}

static std::array<Point, 4> rtet()
{
  while (true)
  {
    Point p0 = rpt3(), p1 = rpt3(), p2 = rpt3(), p3 = rpt3();
    if (std::abs(orient3d(p0, p1, p2, p3)) > 1e-8)
      return {p0, p1, p2, p3};
  }
}

static int sgn(double v) { return (v > 0) ? 1 : (v < 0) ? -1 : 0; }

// ---------------------------------------------------------------------------
// orient2d
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: orient2d sign agreement", "[cgal]")
{
  const int N = 10000;
  int n_mismatch = 0;
  for (int i = 0; i < N; ++i)
  {
    Point a = rpt2(), b = rpt2(), c = rpt2();
    if (sgn(orient2d(a, b, c)) != sgn(cgal_orient2d(a, b, c)))
      ++n_mismatch;
  }
  CHECK(n_mismatch == 0);
}

// ---------------------------------------------------------------------------
// orient3d
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: orient3d sign agreement", "[cgal]")
{
  const int N = 10000;
  int n_mismatch = 0;
  for (int i = 0; i < N; ++i)
  {
    Point a = rpt3(), b = rpt3(), c = rpt3(), d = rpt3();
    if (sgn(orient3d(a, b, c, d)) != sgn(cgal_orient3d(a, b, c, d)))
      ++n_mismatch;
  }
  CHECK(n_mismatch == 0);
}

// ---------------------------------------------------------------------------
// collides_triangle_point_3d  (NEW — comparison test added to track agreement)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: collides_triangle_point_3d agreement", "[cgal]")
{
  SECTION("random triangle-point pairs")
  {
    const int N = 10000;
    int n_mismatch = 0;
    for (int i = 0; i < N; ++i)
    {
      auto t = rtri3d();
      Point q = rpt3();
      bool s = CollisionPredicates::collides_triangle_point_3d(t[0], t[1], t[2], q);
      bool c = cgal_collides_triangle_point_3d(t[0], t[1], t[2], q);
      if (s != c) ++n_mismatch;
    }
    CHECK(n_mismatch == 0);
  }

  SECTION("point on vertex")
  {
    Point p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0);
    CHECK(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, p0)
          == cgal_collides_triangle_point_3d(p0, p1, p2, p0));
    CHECK(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, p1)
          == cgal_collides_triangle_point_3d(p0, p1, p2, p1));
  }

  SECTION("point clearly outside triangle plane")
  {
    Point p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0);
    Point q(5, 5, 5);
    CHECK(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, q)
          == cgal_collides_triangle_point_3d(p0, p1, p2, q));
  }
}

// ---------------------------------------------------------------------------
// collides_tetrahedron_point_3d  (NEW — comparison test added to track agreement)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: collides_tetrahedron_point_3d agreement", "[cgal]")
{
  SECTION("random tet-point pairs")
  {
    const int N = 5000;
    int n_mismatch = 0;
    for (int i = 0; i < N; ++i)
    {
      auto t = rtet();
      Point q = rpt3();
      bool s = CollisionPredicates::collides_tetrahedron_point_3d(t[0], t[1], t[2], t[3], q);
      bool c = cgal_collides_tetrahedron_point_3d(t[0], t[1], t[2], t[3], q);
      if (s != c) ++n_mismatch;
    }
    CHECK(n_mismatch == 0);
  }

  SECTION("point at tet vertex")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0), p3(0,0,1);
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,p0)
          == cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,p0));
  }

  SECTION("point inside unit tet")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0), p3(0,0,1);
    Point q(0.1, 0.1, 0.1);
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,q)
          == cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,q));
  }

  SECTION("point clearly outside unit tet")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0), p3(0,0,1);
    Point q(5,5,5);
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,q)
          == cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,q));
  }
}

// ---------------------------------------------------------------------------
// collides_triangle_triangle_3d
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: collides_triangle_triangle_3d agreement", "[cgal]")
{
  SECTION("random triangle pairs")
  {
    const int N = 300;   // 44 850 pairs
    std::vector<std::array<Point, 3>> tris(N);
    for (auto& t : tris) t = rtri3d();

    int n_mismatch = 0;
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j)
      {
        bool s = CollisionPredicates::collides_triangle_triangle_3d(
          tris[i][0], tris[i][1], tris[i][2],
          tris[j][0], tris[j][1], tris[j][2]);
        bool c = cgal_collides_triangle_triangle_3d(
          tris[i][0], tris[i][1], tris[i][2],
          tris[j][0], tris[j][1], tris[j][2]);
        if (s != c) ++n_mismatch;
      }
    CHECK(n_mismatch == 0);
  }

  SECTION("known intersecting: triangles crossing in 3D")
  {
    // t1 lies in z=0 plane; t2 passes through it
    Point p0(0,0,0), p1(2,0,0), p2(0,2,0);
    Point q0(1,0,-1), q1(1,0,1), q2(0,1,0);
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2)
          == cgal_collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2));
  }

  SECTION("known intersecting: shared vertex")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0);
    Point q0(0,0,0), q1(0,0,1), q2(1,1,1);
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2)
          == cgal_collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2));
  }

  SECTION("known non-intersecting: well separated")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0);
    Point q0(5,5,5), q1(6,5,5), q2(5,6,5);
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2)
          == cgal_collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2));
  }
}

// ---------------------------------------------------------------------------
// collides_tetrahedron_tetrahedron_3d
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: collides_tetrahedron_tetrahedron_3d agreement", "[cgal]")
{
  SECTION("random tetrahedron pairs")
  {
    const int N = 150;   // 11 175 pairs
    std::vector<std::array<Point, 4>> tets(N);
    for (auto& t : tets) t = rtet();

    int n_mismatch = 0;
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j)
      {
        bool s = CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
          tets[i][0], tets[i][1], tets[i][2], tets[i][3],
          tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
        bool c = cgal_collides_tetrahedron_tetrahedron_3d(
          tets[i][0], tets[i][1], tets[i][2], tets[i][3],
          tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
        if (s != c) ++n_mismatch;
      }
    CHECK(n_mismatch == 0);
  }

  SECTION("known: unit tet vs itself (full overlap)")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0), p3(0,0,1);
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
            p0,p1,p2,p3, p0,p1,p2,p3)
          == cgal_collides_tetrahedron_tetrahedron_3d(
               p0,p1,p2,p3, p0,p1,p2,p3));
  }

  SECTION("known: well separated tetrahedra")
  {
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0), p3(0,0,1);
    Point q0(5,0,0), q1(6,0,0), q2(5,1,0), q3(5,0,1);
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
            p0,p1,p2,p3, q0,q1,q2,q3)
          == cgal_collides_tetrahedron_tetrahedron_3d(
               p0,p1,p2,p3, q0,q1,q2,q3));
  }

  SECTION("known: partially overlapping tetrahedra")
  {
    // Unit tet and another shifted by (0.5, 0, 0) should overlap
    Point p0(0,0,0), p1(1,0,0), p2(0,1,0), p3(0,0,1);
    Point q0(0.5,0,0), q1(1.5,0,0), q2(0.5,1,0), q3(0.5,0,1);
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
            p0,p1,p2,p3, q0,q1,q2,q3)
          == cgal_collides_tetrahedron_tetrahedron_3d(
               p0,p1,p2,p3, q0,q1,q2,q3));
  }
}

#endif // GEOMETRY_WITH_CGAL
