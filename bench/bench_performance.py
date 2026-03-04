"""
Performance benchmark: Shewchuk's exact predicates vs CGAL (Python)

For each test:
  1. Times Shewchuk implementation
  2. Times CGAL implementation (using cgal_orient2d / cgal_collides_*)
  3. Runs a correctness check comparing Shewchuk vs CGAL on the same data

Random points use a fixed seed (42) for reproducibility.

Usage:
  python bench/bench_performance.py
  # or from the build/python directory:
  SIMPEX_BUILD_DIR=build/python python bench/bench_performance.py
"""

import sys
import os
import math
import random
import time

# Locate the built simpex module
_script_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root  = os.path.join(_script_dir, '..')
_build_dir  = os.environ.get(
    'SIMPEX_BUILD_DIR',
    os.path.join(_repo_root, 'build', 'python')
)
if _build_dir not in sys.path:
    sys.path.insert(0, _build_dir)

import simpex
from simpex import Point, CollisionPredicates

CGAL_AVAILABLE = getattr(simpex, 'CGAL_AVAILABLE', False)

# ---------------------------------------------------------------------------
# Deterministic RNG — seed 42
# ---------------------------------------------------------------------------
DEGENERACY_THRESHOLD = 1e-8
rng = random.Random(42)


def rand01():
    return rng.random()


def random_point2d():
    return Point(rand01(), rand01(), 0.0)


def random_point3d():
    return Point(rand01(), rand01(), rand01())


def _cross3d(d1, d2):
    return Point(d1.y()*d2.z() - d1.z()*d2.y(),
                 d1.z()*d2.x() - d1.x()*d2.z(),
                 d1.x()*d2.y() - d1.y()*d2.x())


def random_tri3d():
    while True:
        p0, p1, p2 = random_point3d(), random_point3d(), random_point3d()
        d1 = Point(p1.x()-p0.x(), p1.y()-p0.y(), p1.z()-p0.z())
        d2 = Point(p2.x()-p0.x(), p2.y()-p0.y(), p2.z()-p0.z())
        if _cross3d(d1, d2).norm() > DEGENERACY_THRESHOLD:
            return (p0, p1, p2)


def random_tet():
    while True:
        p0, p1, p2, p3 = (random_point3d(), random_point3d(),
                          random_point3d(), random_point3d())
        if abs(simpex.orient3d(p0, p1, p2, p3)) > DEGENERACY_THRESHOLD:
            return (p0, p1, p2, p3)


# ---------------------------------------------------------------------------
# Timing helper
# ---------------------------------------------------------------------------
def elapsed_ms(start, end):
    return (end - start) * 1000.0


def _sgn(v):
    return 1 if v > 0 else (-1 if v < 0 else 0)


# ---------------------------------------------------------------------------
# Result struct
# ---------------------------------------------------------------------------
class BenchResult:
    def __init__(self, label, n_ops, shewchuk_ms, cgal_ms, n_mismatch):
        self.label        = label
        self.n_ops        = n_ops
        self.shewchuk_ms  = shewchuk_ms
        self.cgal_ms      = cgal_ms
        self.n_mismatch   = n_mismatch


def _print_timing_row(r):
    speedup = r.cgal_ms / r.shewchuk_ms if r.shewchuk_ms > 0 else float('inf')
    print(f"  {r.label:<52}  {r.shewchuk_ms:>10.1f} ms   {r.cgal_ms:>10.1f} ms   {speedup:>5.2f}x")


# ---------------------------------------------------------------------------
# Benchmark: orient2d
# ---------------------------------------------------------------------------
def bench_orient2d(n_calls):
    pts = [(random_point2d(), random_point2d(), random_point2d())
           for _ in range(n_calls)]

    t0 = time.perf_counter()
    for a, b, c in pts:
        simpex.orient2d(a, b, c)
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        t2 = time.perf_counter()
        for a, b, c in pts:
            simpex.cgal_orient2d(a, b, c)
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)

        n_mismatch = sum(
            1 for a, b, c in pts
            if _sgn(simpex.orient2d(a, b, c)) != _sgn(simpex.cgal_orient2d(a, b, c))
        )
    else:
        cgal_ms, n_mismatch = float('nan'), 0

    return BenchResult(f"orient2d ({n_calls} calls)",
                       n_calls, shewchuk_ms, cgal_ms, n_mismatch)


# ---------------------------------------------------------------------------
# Benchmark: orient3d
# ---------------------------------------------------------------------------
def bench_orient3d(n_calls):
    pts = [(random_point3d(), random_point3d(),
            random_point3d(), random_point3d())
           for _ in range(n_calls)]

    t0 = time.perf_counter()
    for a, b, c, d in pts:
        simpex.orient3d(a, b, c, d)
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        t2 = time.perf_counter()
        for a, b, c, d in pts:
            simpex.cgal_orient3d(a, b, c, d)
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)

        n_mismatch = sum(
            1 for a, b, c, d in pts
            if _sgn(simpex.orient3d(a, b, c, d)) != _sgn(simpex.cgal_orient3d(a, b, c, d))
        )
    else:
        cgal_ms, n_mismatch = float('nan'), 0

    return BenchResult(f"orient3d ({n_calls} calls)",
                       n_calls, shewchuk_ms, cgal_ms, n_mismatch)


# ---------------------------------------------------------------------------
# Benchmark: triangle-triangle collision (3D)
# ---------------------------------------------------------------------------
def bench_tri_tri_3d(n_tris):
    tris = [random_tri3d() for _ in range(n_tris)]
    n_pairs = n_tris * (n_tris - 1) // 2

    t0 = time.perf_counter()
    for i in range(n_tris):
        for j in range(i + 1, n_tris):
            CollisionPredicates.collides_triangle_triangle_3d(
                tris[i][0], tris[i][1], tris[i][2],
                tris[j][0], tris[j][1], tris[j][2])
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        t2 = time.perf_counter()
        for i in range(n_tris):
            for j in range(i + 1, n_tris):
                simpex.cgal_collides_triangle_triangle_3d(
                    tris[i][0], tris[i][1], tris[i][2],
                    tris[j][0], tris[j][1], tris[j][2])
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)

        n_mismatch = sum(
            1 for i in range(n_tris) for j in range(i + 1, n_tris)
            if CollisionPredicates.collides_triangle_triangle_3d(
                   tris[i][0], tris[i][1], tris[i][2],
                   tris[j][0], tris[j][1], tris[j][2])
               != simpex.cgal_collides_triangle_triangle_3d(
                   tris[i][0], tris[i][1], tris[i][2],
                   tris[j][0], tris[j][1], tris[j][2])
        )
    else:
        cgal_ms, n_mismatch = float('nan'), 0

    return BenchResult(f"triangle-triangle 3D ({n_tris} tris, {n_pairs} pairs)",
                       n_pairs, shewchuk_ms, cgal_ms, n_mismatch)


# ---------------------------------------------------------------------------
# Benchmark: tetrahedron-tetrahedron collision (3D)
# ---------------------------------------------------------------------------
def bench_tet_tet_3d(n_tets):
    tets = [random_tet() for _ in range(n_tets)]
    n_pairs = n_tets * (n_tets - 1) // 2

    t0 = time.perf_counter()
    for i in range(n_tets):
        for j in range(i + 1, n_tets):
            CollisionPredicates.collides_tetrahedron_tetrahedron_3d(
                tets[i][0], tets[i][1], tets[i][2], tets[i][3],
                tets[j][0], tets[j][1], tets[j][2], tets[j][3])
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        t2 = time.perf_counter()
        for i in range(n_tets):
            for j in range(i + 1, n_tets):
                simpex.cgal_collides_tetrahedron_tetrahedron_3d(
                    tets[i][0], tets[i][1], tets[i][2], tets[i][3],
                    tets[j][0], tets[j][1], tets[j][2], tets[j][3])
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)

        n_mismatch = sum(
            1 for i in range(n_tets) for j in range(i + 1, n_tets)
            if CollisionPredicates.collides_tetrahedron_tetrahedron_3d(
                   tets[i][0], tets[i][1], tets[i][2], tets[i][3],
                   tets[j][0], tets[j][1], tets[j][2], tets[j][3])
               != simpex.cgal_collides_tetrahedron_tetrahedron_3d(
                   tets[i][0], tets[i][1], tets[i][2], tets[i][3],
                   tets[j][0], tets[j][1], tets[j][2], tets[j][3])
        )
    else:
        cgal_ms, n_mismatch = float('nan'), 0

    return BenchResult(f"tetrahedron-tetrahedron 3D ({n_tets} tets, {n_pairs} pairs)",
                       n_pairs, shewchuk_ms, cgal_ms, n_mismatch)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    print("=== Performance Benchmark: Shewchuk vs CGAL (EPICK) [Python] ===")
    print("    Random seed: 42 (deterministic)")
    print()
    print(f"  {'Test':<52}  {'Shewchuk':>12}   {'CGAL':>12}   {'Speedup':>7}")
    print(f"  {'-'*90}")

    results = []
    results.append(bench_orient2d(1_000_000))
    _print_timing_row(results[-1])
    results.append(bench_orient3d(1_000_000))
    _print_timing_row(results[-1])
    results.append(bench_tri_tri_3d(400))
    _print_timing_row(results[-1])
    results.append(bench_tet_tet_3d(300))
    _print_timing_row(results[-1])

    # Correctness summary
    print()
    print("=== Correctness Summary: Shewchuk vs CGAL ===")
    print()
    print(f"  {'Test':<52}  {'Status':>6}   Mismatches")
    print(f"  {'-'*72}")

    all_pass = True
    for r in results:
        status = "PASS" if r.n_mismatch == 0 else "FAIL"
        if r.n_mismatch > 0:
            all_pass = False
        print(f"  {r.label:<52}  {status:>6}   {r.n_mismatch} / {r.n_ops}")

    print()
    print(f"  Overall: {'ALL PASS' if all_pass else 'SOME FAILURES DETECTED'}")
    print()
    print("Notes:")
    print("  Speedup > 1 means Shewchuk is faster than CGAL.")
    if not CGAL_AVAILABLE:
        print("  CGAL not available: rebuild the Python module with CGAL.")
        print("  See README.md for instructions on compiling with CGAL.")
    else:
        print("  CGAL uses Exact_predicates_inexact_constructions_kernel (EPICK).")
    print("  Python function-call overhead dominates for single-predicate tests.")
    print("  See bench/performance.cpp for a lower-overhead C++ version.")
    print()
