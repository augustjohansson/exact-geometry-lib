// Performance benchmark: Shewchuk's exact predicates vs CGAL
//
// Compares timing of:
//   1. orient2d  - 2D orientation predicate
//   2. orient3d  - 3D orientation predicate
//   3. Triangle-triangle collision (3D)
//   4. Tetrahedron-tetrahedron collision (3D)
//
// Usage: ./geometry_bench

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

// Shewchuk predicates and geometry library
#include "../geometry/predicates.h"
#include "../geometry/Point.h"
#include "../geometry/CollisionPredicates.h"

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2   CGALPoint2;
typedef K::Point_3   CGALPoint3;
typedef K::Triangle_3 CGALTriangle3;
typedef K::Tetrahedron_3 CGALTetrahedron3;

using namespace dolfin;

// ---------------------------------------------------------------------------
// Timing helper
// ---------------------------------------------------------------------------
using Clock = std::chrono::steady_clock;
using Ms    = std::chrono::duration<double, std::milli>;

static double elapsed_ms(const Clock::time_point& start,
                          const Clock::time_point& end)
{
  return std::chrono::duration_cast<Ms>(end - start).count();
}

// ---------------------------------------------------------------------------
// Degeneracy threshold for random simplex generation
// ---------------------------------------------------------------------------
static constexpr double DEGENERACY_THRESHOLD = 1e-8;

static std::mt19937 rng(42);
static std::uniform_real_distribution<double> dist(0.0, 1.0);

static double rand01() { return dist(rng); }

// Generate a non-degenerate random 3D tetrahedron (vertices in [0,1]^3)
// Retry if degenerate (nearly zero volume).
static std::vector<Point> random_tet()
{
  while (true)
  {
    Point p0(rand01(), rand01(), rand01());
    Point p1(rand01(), rand01(), rand01());
    Point p2(rand01(), rand01(), rand01());
    Point p3(rand01(), rand01(), rand01());
    if (std::abs(orient3d(p0, p1, p2, p3)) > DEGENERACY_THRESHOLD)
      return {p0, p1, p2, p3};
  }
}

// Generate a non-degenerate random 3D triangle (vertices in [0,1]^3)
static std::vector<Point> random_tri3d()
{
  while (true)
  {
    Point p0(rand01(), rand01(), rand01());
    Point p1(rand01(), rand01(), rand01());
    Point p2(rand01(), rand01(), rand01());
    // Check that triangle is not degenerate (cross-product norm > eps)
    Point d1 = p1 - p0;
    Point d2 = p2 - p0;
    Point cross = d1.cross(d2);
    if (cross.norm() > DEGENERACY_THRESHOLD)
      return {p0, p1, p2};
  }
}

// Generate a random 2D point (x, y in [0,1])
static Point random_point2d() { return Point(rand01(), rand01(), 0.0); }

// Generate a random 3D point (x, y, z in [0,1])
static Point random_point3d() { return Point(rand01(), rand01(), rand01()); }

// ---------------------------------------------------------------------------
// Print benchmark result
// ---------------------------------------------------------------------------
static void print_result(const std::string& label,
                          long long n_ops,
                          double shewchuk_ms,
                          double cgal_ms)
{
  double speedup = cgal_ms / shewchuk_ms;
  std::printf("  %-52s  %8.1f ms   %8.1f ms   %5.2fx\n",
              label.c_str(), shewchuk_ms, cgal_ms, speedup);
}

// ---------------------------------------------------------------------------
// Benchmark: orient2d
// ---------------------------------------------------------------------------
static void bench_orient2d(int n_calls)
{
  // Pre-generate random 2D points
  std::vector<Point> pts(3 * n_calls);
  for (auto& p : pts) p = random_point2d();

  // --- Shewchuk ---
  volatile double sink = 0.0; // prevent dead-code elimination
  auto t0 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
  {
    const Point& a = pts[3*i];
    const Point& b = pts[3*i + 1];
    const Point& c = pts[3*i + 2];
    sink += orient2d(a, b, c);
  }
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);

  // --- CGAL ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
  {
    const Point& a = pts[3*i];
    const Point& b = pts[3*i + 1];
    const Point& c = pts[3*i + 2];
    CGALPoint2 ca(a.x(), a.y());
    CGALPoint2 cb(b.x(), b.y());
    CGALPoint2 cc(c.x(), c.y());
    cgal_sink += (int)CGAL::orientation(ca, cb, cc);
  }
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);

  (void)sink; (void)cgal_sink;

  char label[80];
  std::snprintf(label, sizeof(label), "orient2d (%d calls)", n_calls);
  print_result(label, n_calls, shewchuk_ms, cgal_ms);
}

// ---------------------------------------------------------------------------
// Benchmark: orient3d
// ---------------------------------------------------------------------------
static void bench_orient3d(int n_calls)
{
  std::vector<Point> pts(4 * n_calls);
  for (auto& p : pts) p = random_point3d();

  // --- Shewchuk ---
  volatile double sink = 0.0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
  {
    const Point& a = pts[4*i];
    const Point& b = pts[4*i + 1];
    const Point& c = pts[4*i + 2];
    const Point& d = pts[4*i + 3];
    sink += orient3d(a, b, c, d);
  }
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);

  // --- CGAL ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
  {
    const Point& a = pts[4*i];
    const Point& b = pts[4*i + 1];
    const Point& c = pts[4*i + 2];
    const Point& d = pts[4*i + 3];
    CGALPoint3 ca(a.x(), a.y(), a.z());
    CGALPoint3 cb(b.x(), b.y(), b.z());
    CGALPoint3 cc(c.x(), c.y(), c.z());
    CGALPoint3 cd(d.x(), d.y(), d.z());
    cgal_sink += (int)CGAL::orientation(ca, cb, cc, cd);
  }
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);

  (void)sink; (void)cgal_sink;

  char label[80];
  std::snprintf(label, sizeof(label), "orient3d (%d calls)", n_calls);
  print_result(label, n_calls, shewchuk_ms, cgal_ms);
}

// ---------------------------------------------------------------------------
// Benchmark: triangle-triangle collision (3D)
// ---------------------------------------------------------------------------
static void bench_tri_tri_3d(int n_tris)
{
  // Generate triangles
  std::vector<std::vector<Point>> tris(n_tris);
  for (auto& t : tris) t = random_tri3d();

  long long n_pairs = (long long)n_tris * (n_tris - 1) / 2;

  // --- Shewchuk ---
  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = i + 1; j < n_tris; ++j)
      sink += (int)CollisionPredicates::collides_triangle_triangle_3d(
        tris[i][0], tris[i][1], tris[i][2],
        tris[j][0], tris[j][1], tris[j][2]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);

  // Build CGAL triangles
  std::vector<CGALTriangle3> cgal_tris(n_tris);
  for (int i = 0; i < n_tris; ++i)
  {
    const auto& t = tris[i];
    cgal_tris[i] = CGALTriangle3(
      CGALPoint3(t[0].x(), t[0].y(), t[0].z()),
      CGALPoint3(t[1].x(), t[1].y(), t[1].z()),
      CGALPoint3(t[2].x(), t[2].y(), t[2].z()));
  }

  // --- CGAL ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = i + 1; j < n_tris; ++j)
      cgal_sink += (int)CGAL::do_intersect(cgal_tris[i], cgal_tris[j]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);

  (void)sink; (void)cgal_sink;

  char label[80];
  std::snprintf(label, sizeof(label),
                "triangle-triangle 3D (%d tris, %lld pairs)",
                n_tris, n_pairs);
  print_result(label, n_pairs, shewchuk_ms, cgal_ms);
}

// ---------------------------------------------------------------------------
// Benchmark: tetrahedron-tetrahedron collision (3D)
// ---------------------------------------------------------------------------
static void bench_tet_tet_3d(int n_tets)
{
  // Generate tetrahedra
  std::vector<std::vector<Point>> tets(n_tets);
  for (auto& t : tets) t = random_tet();

  long long n_pairs = (long long)n_tets * (n_tets - 1) / 2;

  // --- Shewchuk ---
  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = i + 1; j < n_tets; ++j)
      sink += (int)CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);

  // Build CGAL tetrahedra
  std::vector<CGALTetrahedron3> cgal_tets(n_tets);
  for (int i = 0; i < n_tets; ++i)
  {
    const auto& t = tets[i];
    cgal_tets[i] = CGALTetrahedron3(
      CGALPoint3(t[0].x(), t[0].y(), t[0].z()),
      CGALPoint3(t[1].x(), t[1].y(), t[1].z()),
      CGALPoint3(t[2].x(), t[2].y(), t[2].z()),
      CGALPoint3(t[3].x(), t[3].y(), t[3].z()));
  }

  // --- CGAL ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = i + 1; j < n_tets; ++j)
      cgal_sink += (int)CGAL::do_intersect(cgal_tets[i], cgal_tets[j]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);

  (void)sink; (void)cgal_sink;

  char label[80];
  std::snprintf(label, sizeof(label),
                "tetrahedron-tetrahedron 3D (%d tets, %lld pairs)",
                n_tets, n_pairs);
  print_result(label, n_pairs, shewchuk_ms, cgal_ms);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main()
{
  std::printf("\n");
  std::printf("=== Performance Benchmark: Shewchuk vs CGAL (EPICK) ===\n\n");
  std::printf("  %-52s  %12s   %12s   %7s\n",
              "Test", "Shewchuk", "CGAL", "Speedup");
  std::printf("  %s\n", std::string(90, '-').c_str());

  bench_orient2d(10000000);
  bench_orient3d(10000000);
  bench_tri_tri_3d(1000);
  bench_tet_tet_3d(500);

  std::printf("\n");
  std::printf("Notes:\n");
  std::printf("  Speedup > 1 means Shewchuk is faster than CGAL.\n");
  std::printf("  CGAL uses Exact_predicates_inexact_constructions_kernel (EPICK).\n");
  std::printf("  Shewchuk orient2d/orient3d are adaptive-precision predicates.\n");
  std::printf("\n");

  return 0;
}
