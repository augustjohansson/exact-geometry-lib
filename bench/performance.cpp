// Performance benchmark: Shewchuk's exact predicates vs CGAL
//
// For each test:
//   1. Times Shewchuk implementation
//   2. Times CGAL implementation (using CGAL::do_intersect / CGAL::orientation)
//   3. Runs a correctness check comparing Shewchuk vs CGAL on the same data
//
// Random points use a fixed seed (42) for reproducibility.
//
// Exit code: 0 if all correctness checks pass, 1 if any mismatches are detected.
//
// Usage: ./simpex_bench

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

// Shewchuk predicates and geometry library
#include "../geometry/predicates.h"
#include "../geometry/Point.h"
#include "../geometry/CollisionPredicates.h"

// CGAL predicates (requires -DGEOMETRY_WITH_CGAL)
#include "../geometry/CGALPredicates.h"

using namespace simpex;

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
// Deterministic RNG (seed 42) — same points every run
// ---------------------------------------------------------------------------
static constexpr double DEGENERACY_THRESHOLD = 1e-8;

static std::mt19937 rng(42);
static std::uniform_real_distribution<double> dist(0.0, 1.0);

static double rand01() { return dist(rng); }

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

static std::vector<Point> random_tri3d()
{
  while (true)
  {
    Point p0(rand01(), rand01(), rand01());
    Point p1(rand01(), rand01(), rand01());
    Point p2(rand01(), rand01(), rand01());
    Point d1 = p1 - p0;
    Point d2 = p2 - p0;
    if (d1.cross(d2).norm() > DEGENERACY_THRESHOLD)
      return {p0, p1, p2};
  }
}

static Point random_point2d() { return Point(rand01(), rand01(), 0.0); }
static Point random_point3d() { return Point(rand01(), rand01(), rand01()); }

static std::vector<Point> random_seg3d()
{
  return { random_point3d(), random_point3d() };
}

// ---------------------------------------------------------------------------
// Result struct
// ---------------------------------------------------------------------------
struct BenchResult
{
  std::string label;
  long long   n_ops;
  double      shewchuk_ms;
  double      cgal_ms;
  long long   n_mismatch;
};

static void print_timing_row(const BenchResult& r)
{
  double speedup = r.cgal_ms / r.shewchuk_ms;
  std::printf("  %-52s  %8.1f ms   %8.1f ms   %5.2fx\n",
              r.label.c_str(), r.shewchuk_ms, r.cgal_ms, speedup);
}

// ---------------------------------------------------------------------------
// Helper: compare orientation sign
// ---------------------------------------------------------------------------
static int sgn(double v) { return (v > 0) ? 1 : (v < 0) ? -1 : 0; }

// ---------------------------------------------------------------------------
// Benchmark: orient2d
// ---------------------------------------------------------------------------
static BenchResult bench_orient2d(int n_calls)
{
  // Pre-generate random 2D points (deterministic order)
  std::vector<Point> pts(3 * n_calls);
  for (auto& p : pts) p = random_point2d();

  // --- Shewchuk timing ---
  volatile double sink = 0.0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
    sink += orient2d(pts[3*i], pts[3*i+1], pts[3*i+2]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  // --- CGAL timing ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
    cgal_sink += (int)cgal_orient2d(pts[3*i], pts[3*i+1], pts[3*i+2]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  // --- Correctness check: compare signs ---
  long long n_mismatch = 0;
  for (int i = 0; i < n_calls; ++i)
  {
    double s = orient2d(pts[3*i], pts[3*i+1], pts[3*i+2]);
    double c = cgal_orient2d(pts[3*i], pts[3*i+1], pts[3*i+2]);
    if (sgn(s) != sgn(c)) ++n_mismatch;
  }

  char label[80];
  std::snprintf(label, sizeof(label), "orient2d (%d calls)", n_calls);
  return {label, n_calls, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Benchmark: orient3d
// ---------------------------------------------------------------------------
static BenchResult bench_orient3d(int n_calls)
{
  std::vector<Point> pts(4 * n_calls);
  for (auto& p : pts) p = random_point3d();

  // --- Shewchuk timing ---
  volatile double sink = 0.0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
    sink += orient3d(pts[4*i], pts[4*i+1], pts[4*i+2], pts[4*i+3]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  // --- CGAL timing ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_calls; ++i)
    cgal_sink += (int)cgal_orient3d(pts[4*i], pts[4*i+1], pts[4*i+2], pts[4*i+3]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  // --- Correctness check: compare signs ---
  long long n_mismatch = 0;
  for (int i = 0; i < n_calls; ++i)
  {
    double s = orient3d(pts[4*i], pts[4*i+1], pts[4*i+2], pts[4*i+3]);
    double c = cgal_orient3d(pts[4*i], pts[4*i+1], pts[4*i+2], pts[4*i+3]);
    if (sgn(s) != sgn(c)) ++n_mismatch;
  }

  char label[80];
  std::snprintf(label, sizeof(label), "orient3d (%d calls)", n_calls);
  return {label, n_calls, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Benchmark: triangle-triangle collision (3D)
// ---------------------------------------------------------------------------
static BenchResult bench_tri_tri_3d(int n_tris)
{
  std::vector<std::vector<Point>> tris(n_tris);
  for (auto& t : tris) t = random_tri3d();

  long long n_pairs = (long long)n_tris * (n_tris - 1) / 2;

  // --- Shewchuk timing ---
  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = i + 1; j < n_tris; ++j)
      sink += (int)CollisionPredicates::collides_triangle_triangle_3d(
        tris[i][0], tris[i][1], tris[i][2],
        tris[j][0], tris[j][1], tris[j][2]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  // --- CGAL timing ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = i + 1; j < n_tris; ++j)
      cgal_sink += (int)cgal_collides_triangle_triangle_3d(
        tris[i][0], tris[i][1], tris[i][2],
        tris[j][0], tris[j][1], tris[j][2]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  // --- Correctness check ---
  long long n_mismatch = 0;
  for (int i = 0; i < n_tris; ++i)
    for (int j = i + 1; j < n_tris; ++j)
    {
      bool s = CollisionPredicates::collides_triangle_triangle_3d(
        tris[i][0], tris[i][1], tris[i][2],
        tris[j][0], tris[j][1], tris[j][2]);
      bool c = cgal_collides_triangle_triangle_3d(
        tris[i][0], tris[i][1], tris[i][2],
        tris[j][0], tris[j][1], tris[j][2]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "triangle-triangle 3D (%d tris, %lld pairs)",
                n_tris, n_pairs);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Benchmark: tetrahedron-tetrahedron collision (3D)
// ---------------------------------------------------------------------------
static BenchResult bench_tet_tet_3d(int n_tets)
{
  std::vector<std::vector<Point>> tets(n_tets);
  for (auto& t : tets) t = random_tet();

  long long n_pairs = (long long)n_tets * (n_tets - 1) / 2;

  // --- Shewchuk timing ---
  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = i + 1; j < n_tets; ++j)
      sink += (int)CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  // --- CGAL timing ---
  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = i + 1; j < n_tets; ++j)
      cgal_sink += (int)cgal_collides_tetrahedron_tetrahedron_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  // --- Correctness check ---
  long long n_mismatch = 0;
  for (int i = 0; i < n_tets; ++i)
    for (int j = i + 1; j < n_tets; ++j)
    {
      bool s = CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
      bool c = cgal_collides_tetrahedron_tetrahedron_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tets[j][0], tets[j][1], tets[j][2], tets[j][3]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "tetrahedron-tetrahedron 3D (%d tets, %lld pairs)",
                n_tets, n_pairs);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Benchmark: triangle-segment collision (3D)
// ---------------------------------------------------------------------------
static BenchResult bench_tri_seg_3d(int n_tris, int n_segs)
{
  std::vector<std::vector<Point>> tris(n_tris);
  for (auto& t : tris) t = random_tri3d();
  std::vector<std::vector<Point>> segs(n_segs);
  for (auto& s : segs) s = random_seg3d();

  long long n_pairs = (long long)n_tris * n_segs;

  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = 0; j < n_segs; ++j)
      sink += (int)CollisionPredicates::collides_triangle_segment_3d(
        tris[i][0], tris[i][1], tris[i][2],
        segs[j][0], segs[j][1]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = 0; j < n_segs; ++j)
      cgal_sink += (int)cgal_collides_triangle_segment_3d(
        tris[i][0], tris[i][1], tris[i][2],
        segs[j][0], segs[j][1]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  long long n_mismatch = 0;
  for (int i = 0; i < n_tris; ++i)
    for (int j = 0; j < n_segs; ++j)
    {
      bool s = CollisionPredicates::collides_triangle_segment_3d(
        tris[i][0], tris[i][1], tris[i][2],
        segs[j][0], segs[j][1]);
      bool c = cgal_collides_triangle_segment_3d(
        tris[i][0], tris[i][1], tris[i][2],
        segs[j][0], segs[j][1]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "triangle-segment 3D (%d tris, %d segs)", n_tris, n_segs);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Benchmark: tetrahedron-segment collision (3D)
// ---------------------------------------------------------------------------
static BenchResult bench_tet_seg_3d(int n_tets, int n_segs)
{
  std::vector<std::vector<Point>> tets(n_tets);
  for (auto& t : tets) t = random_tet();
  std::vector<std::vector<Point>> segs(n_segs);
  for (auto& s : segs) s = random_seg3d();

  long long n_pairs = (long long)n_tets * n_segs;

  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_segs; ++j)
      sink += (int)CollisionPredicates::collides_tetrahedron_segment_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        segs[j][0], segs[j][1]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_segs; ++j)
      cgal_sink += (int)cgal_collides_tetrahedron_segment_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        segs[j][0], segs[j][1]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  long long n_mismatch = 0;
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_segs; ++j)
    {
      bool s = CollisionPredicates::collides_tetrahedron_segment_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        segs[j][0], segs[j][1]);
      bool c = cgal_collides_tetrahedron_segment_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        segs[j][0], segs[j][1]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "tetrahedron-segment 3D (%d tets, %d segs)", n_tets, n_segs);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Benchmark: tetrahedron-triangle collision (3D)
// ---------------------------------------------------------------------------
static BenchResult bench_tet_tri_3d(int n_tets, int n_tris)
{
  std::vector<std::vector<Point>> tets(n_tets);
  for (auto& t : tets) t = random_tet();
  std::vector<std::vector<Point>> tris(n_tris);
  for (auto& t : tris) t = random_tri3d();

  long long n_pairs = (long long)n_tets * n_tris;

  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_tris; ++j)
      sink += (int)CollisionPredicates::collides_tetrahedron_triangle_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tris[j][0], tris[j][1], tris[j][2]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_tris; ++j)
      cgal_sink += (int)cgal_collides_tetrahedron_triangle_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tris[j][0], tris[j][1], tris[j][2]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  long long n_mismatch = 0;
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_tris; ++j)
    {
      bool s = CollisionPredicates::collides_tetrahedron_triangle_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tris[j][0], tris[j][1], tris[j][2]);
      bool c = cgal_collides_tetrahedron_triangle_3d(
        tets[i][0], tets[i][1], tets[i][2], tets[i][3],
        tris[j][0], tris[j][1], tris[j][2]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "tetrahedron-triangle 3D (%d tets, %d tris)", n_tets, n_tris);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}

// ---------------------------------------------------------------------------
// Near-miss helpers: shift all q-points by 'delta' along z so entities are
// close but (usually) not intersecting.
// ---------------------------------------------------------------------------
static BenchResult bench_tri_tri_3d_near(int n_tris, double delta)
{
  std::vector<std::vector<Point>> trisA(n_tris);
  for (auto& t : trisA) t = random_tri3d();

  // Build near-miss triangles: translate each in z by delta
  std::vector<std::vector<Point>> trisB(n_tris);
  for (int i = 0; i < n_tris; ++i)
  {
    const Point off(0.0, 0.0, delta);
    trisB[i] = { trisA[i][0] + off, trisA[i][1] + off, trisA[i][2] + off };
  }

  long long n_pairs = (long long)n_tris * n_tris;

  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = 0; j < n_tris; ++j)
      sink += (int)CollisionPredicates::collides_triangle_triangle_3d(
        trisA[i][0], trisA[i][1], trisA[i][2],
        trisB[j][0], trisB[j][1], trisB[j][2]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tris; ++i)
    for (int j = 0; j < n_tris; ++j)
      cgal_sink += (int)cgal_collides_triangle_triangle_3d(
        trisA[i][0], trisA[i][1], trisA[i][2],
        trisB[j][0], trisB[j][1], trisB[j][2]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  long long n_mismatch = 0;
  for (int i = 0; i < n_tris; ++i)
    for (int j = 0; j < n_tris; ++j)
    {
      bool s = CollisionPredicates::collides_triangle_triangle_3d(
        trisA[i][0], trisA[i][1], trisA[i][2],
        trisB[j][0], trisB[j][1], trisB[j][2]);
      bool c = cgal_collides_triangle_triangle_3d(
        trisA[i][0], trisA[i][1], trisA[i][2],
        trisB[j][0], trisB[j][1], trisB[j][2]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "tri-tri 3D near-miss delta=%.0e (%d tris)", delta, n_tris);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}

static BenchResult bench_tet_tet_3d_near(int n_tets, double delta)
{
  std::vector<std::vector<Point>> tetsA(n_tets);
  for (auto& t : tetsA) t = random_tet();

  std::vector<std::vector<Point>> tetsB(n_tets);
  for (int i = 0; i < n_tets; ++i)
  {
    const Point off(0.0, 0.0, delta);
    tetsB[i] = { tetsA[i][0] + off, tetsA[i][1] + off,
                 tetsA[i][2] + off, tetsA[i][3] + off };
  }

  long long n_pairs = (long long)n_tets * n_tets;

  volatile int sink = 0;
  auto t0 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_tets; ++j)
      sink += (int)CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
        tetsA[i][0], tetsA[i][1], tetsA[i][2], tetsA[i][3],
        tetsB[j][0], tetsB[j][1], tetsB[j][2], tetsB[j][3]);
  auto t1 = Clock::now();
  double shewchuk_ms = elapsed_ms(t0, t1);
  (void)sink;

  volatile int cgal_sink = 0;
  auto t2 = Clock::now();
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_tets; ++j)
      cgal_sink += (int)cgal_collides_tetrahedron_tetrahedron_3d(
        tetsA[i][0], tetsA[i][1], tetsA[i][2], tetsA[i][3],
        tetsB[j][0], tetsB[j][1], tetsB[j][2], tetsB[j][3]);
  auto t3 = Clock::now();
  double cgal_ms = elapsed_ms(t2, t3);
  (void)cgal_sink;

  long long n_mismatch = 0;
  for (int i = 0; i < n_tets; ++i)
    for (int j = 0; j < n_tets; ++j)
    {
      bool s = CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
        tetsA[i][0], tetsA[i][1], tetsA[i][2], tetsA[i][3],
        tetsB[j][0], tetsB[j][1], tetsB[j][2], tetsB[j][3]);
      bool c = cgal_collides_tetrahedron_tetrahedron_3d(
        tetsA[i][0], tetsA[i][1], tetsA[i][2], tetsA[i][3],
        tetsB[j][0], tetsB[j][1], tetsB[j][2], tetsB[j][3]);
      if (s != c) ++n_mismatch;
    }

  char label[80];
  std::snprintf(label, sizeof(label),
                "tet-tet 3D near-miss delta=%.0e (%d tets)", delta, n_tets);
  return {label, n_pairs, shewchuk_ms, cgal_ms, n_mismatch};
}


int main()
{
  std::printf("\n");
  std::printf("=== Performance Benchmark: Shewchuk vs CGAL (EPICK) ===\n");
  std::printf("    Random seed: 42 (deterministic)\n\n");
  std::printf("  %-52s  %12s   %12s   %7s\n",
              "Test", "Shewchuk", "CGAL", "Speedup");
  std::printf("  %s\n", std::string(90, '-').c_str());

  std::vector<BenchResult> results;
  results.push_back(bench_orient2d(10000000));
  print_timing_row(results.back());
  results.push_back(bench_orient3d(10000000));
  print_timing_row(results.back());
  results.push_back(bench_tri_tri_3d(1000));
  print_timing_row(results.back());
  results.push_back(bench_tri_seg_3d(500, 500));
  print_timing_row(results.back());
  results.push_back(bench_tet_tet_3d(500));
  print_timing_row(results.back());
  results.push_back(bench_tet_seg_3d(300, 300));
  print_timing_row(results.back());
  results.push_back(bench_tet_tri_3d(300, 300));
  print_timing_row(results.back());

  std::printf("\n=== Near-miss benchmarks (close entities, 1e-20 separation) ===\n\n");
  std::printf("  %-52s  %12s   %12s   %7s\n",
              "Test", "Shewchuk", "CGAL", "Speedup");
  std::printf("  %s\n", std::string(90, '-').c_str());

  results.push_back(bench_tri_tri_3d_near(500, 1e-20));
  print_timing_row(results.back());
  results.push_back(bench_tet_tet_3d_near(300, 1e-20));
  print_timing_row(results.back());

  // ---------------------------------------------------------------------------
  // Correctness summary
  // ---------------------------------------------------------------------------
  std::printf("\n=== Correctness Summary: Shewchuk vs CGAL ===\n\n");
  std::printf("  %-52s  %6s   %s\n", "Test", "Status", "Mismatches");
  std::printf("  %s\n", std::string(72, '-').c_str());

  bool all_pass = true;
  for (const auto& r : results)
  {
    const char* status = (r.n_mismatch == 0) ? "PASS" : "FAIL";
    if (r.n_mismatch > 0) all_pass = false;
    std::printf("  %-52s  %6s   %lld / %lld\n",
                r.label.c_str(), status, r.n_mismatch, r.n_ops);
  }

  std::printf("\n  Overall: %s\n\n", all_pass ? "ALL PASS" : "SOME FAILURES DETECTED");

  std::printf("Notes:\n");
  std::printf("  Speedup > 1 means Shewchuk is faster than CGAL.\n");
  std::printf("  CGAL uses Exact_predicates_inexact_constructions_kernel (EPICK).\n");
  std::printf("  Shewchuk orient2d/orient3d are adaptive-precision predicates.\n");
  std::printf("\n");

  return all_pass ? 0 : 1;
}
