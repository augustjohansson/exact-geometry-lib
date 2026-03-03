"""
Performance benchmark: Shewchuk's exact predicates vs CGAL (Python)

Mirrors the C++ benchmark in bench/performance.cpp.

Compares timing of:
  1. orient2d  - 2D orientation predicate
  2. orient3d  - 3D orientation predicate
  3. Triangle-triangle collision (3D)
  4. Tetrahedron-tetrahedron collision (3D)

Usage:
  python bench_performance.py
  # or from the repo root:
  python bench/bench_performance.py
"""

import sys
import os
import math
import random
import time

# Locate the built geometry module
_script_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root  = os.path.join(_script_dir, '..')
_build_dir  = os.environ.get(
    'GEOMETRY_BUILD_DIR',
    os.path.join(_repo_root, 'build', 'python')
)
if _build_dir not in sys.path:
    sys.path.insert(0, _build_dir)

import geometry
from geometry import Point, CollisionPredicates

# ---------------------------------------------------------------------------
# Degeneracy threshold for random simplex generation
# ---------------------------------------------------------------------------
DEGENERACY_THRESHOLD = 1e-8

# ---------------------------------------------------------------------------
# Check CGAL availability
# ---------------------------------------------------------------------------
CGAL_AVAILABLE = getattr(geometry, 'CGAL_AVAILABLE', False)


# ---------------------------------------------------------------------------
# Random data generators (seeded for reproducibility)
# ---------------------------------------------------------------------------
rng = random.Random(42)

def rand01():
    return rng.random()

def random_point2d():
    return Point(rand01(), rand01(), 0.0)

def random_point3d():
    return Point(rand01(), rand01(), rand01())

def cross3d(d1, d2):
    return Point(d1.y()*d2.z() - d1.z()*d2.y(),
                 d1.z()*d2.x() - d1.x()*d2.z(),
                 d1.x()*d2.y() - d1.y()*d2.x())

def random_tri3d():
    """Generate a non-degenerate random 3D triangle."""
    while True:
        p0 = random_point3d()
        p1 = random_point3d()
        p2 = random_point3d()
        d1 = Point(p1.x()-p0.x(), p1.y()-p0.y(), p1.z()-p0.z())
        d2 = Point(p2.x()-p0.x(), p2.y()-p0.y(), p2.z()-p0.z())
        cp = cross3d(d1, d2)
        if cp.norm() > DEGENERACY_THRESHOLD:
            return (p0, p1, p2)

def random_tet():
    """Generate a non-degenerate random tetrahedron."""
    while True:
        p0 = random_point3d()
        p1 = random_point3d()
        p2 = random_point3d()
        p3 = random_point3d()
        if abs(geometry.orient3d(p0, p1, p2, p3)) > DEGENERACY_THRESHOLD:
            return (p0, p1, p2, p3)


# ---------------------------------------------------------------------------
# Timing helper
# ---------------------------------------------------------------------------
def elapsed_ms(start, end):
    return (end - start) * 1000.0


# ---------------------------------------------------------------------------
# Print benchmark result
# ---------------------------------------------------------------------------
HEADER_PRINTED = False

def print_header():
    global HEADER_PRINTED
    if not HEADER_PRINTED:
        print()
        print("=== Performance Benchmark: Shewchuk vs CGAL (EPICK) [Python] ===")
        print()
        print(f"  {'Test':<52}  {'Shewchuk':>12}   {'CGAL':>12}   {'Speedup':>7}")
        print(f"  {'-'*90}")
        HEADER_PRINTED = True

def print_result(label, shewchuk_ms, cgal_ms):
    speedup = cgal_ms / shewchuk_ms if shewchuk_ms > 0 else float('inf')
    print(f"  {label:<52}  {shewchuk_ms:>10.1f} ms   {cgal_ms:>10.1f} ms   {speedup:>5.2f}x")


# ---------------------------------------------------------------------------
# Benchmark: orient2d
# ---------------------------------------------------------------------------
def bench_orient2d(n_calls):
    # Pre-generate random 2D points
    pts = [(random_point2d(), random_point2d(), random_point2d())
           for _ in range(n_calls)]

    # --- Shewchuk ---
    t0 = time.perf_counter()
    for a, b, c in pts:
        geometry.orient2d(a, b, c)
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        # --- CGAL ---
        t2 = time.perf_counter()
        for a, b, c in pts:
            geometry.cgal_orient2d(a, b, c)
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)
    else:
        cgal_ms = float('nan')

    label = f"orient2d ({n_calls} calls)"
    print_result(label, shewchuk_ms, cgal_ms)


# ---------------------------------------------------------------------------
# Benchmark: orient3d
# ---------------------------------------------------------------------------
def bench_orient3d(n_calls):
    pts = [(random_point3d(), random_point3d(),
            random_point3d(), random_point3d())
           for _ in range(n_calls)]

    # --- Shewchuk ---
    t0 = time.perf_counter()
    for a, b, c, d in pts:
        geometry.orient3d(a, b, c, d)
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        # --- CGAL ---
        t2 = time.perf_counter()
        for a, b, c, d in pts:
            geometry.cgal_orient3d(a, b, c, d)
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)
    else:
        cgal_ms = float('nan')

    label = f"orient3d ({n_calls} calls)"
    print_result(label, shewchuk_ms, cgal_ms)


# ---------------------------------------------------------------------------
# Benchmark: triangle-triangle collision (3D)
# ---------------------------------------------------------------------------
def bench_tri_tri_3d(n_tris):
    tris = [random_tri3d() for _ in range(n_tris)]
    n_pairs = n_tris * (n_tris - 1) // 2

    # --- Shewchuk ---
    t0 = time.perf_counter()
    for i in range(n_tris):
        for j in range(i + 1, n_tris):
            CollisionPredicates.collides_triangle_triangle_3d(
                tris[i][0], tris[i][1], tris[i][2],
                tris[j][0], tris[j][1], tris[j][2])
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        # --- CGAL ---
        t2 = time.perf_counter()
        for i in range(n_tris):
            for j in range(i + 1, n_tris):
                geometry.cgal_collides_triangle_triangle_3d(
                    tris[i][0], tris[i][1], tris[i][2],
                    tris[j][0], tris[j][1], tris[j][2])
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)
    else:
        cgal_ms = float('nan')

    label = f"triangle-triangle 3D ({n_tris} tris, {n_pairs} pairs)"
    print_result(label, shewchuk_ms, cgal_ms)


# ---------------------------------------------------------------------------
# Benchmark: tetrahedron-tetrahedron collision (3D)
# ---------------------------------------------------------------------------
def bench_tet_tet_3d(n_tets):
    tets = [random_tet() for _ in range(n_tets)]
    n_pairs = n_tets * (n_tets - 1) // 2

    # --- Shewchuk ---
    t0 = time.perf_counter()
    for i in range(n_tets):
        for j in range(i + 1, n_tets):
            CollisionPredicates.collides_tetrahedron_tetrahedron_3d(
                tets[i][0], tets[i][1], tets[i][2], tets[i][3],
                tets[j][0], tets[j][1], tets[j][2], tets[j][3])
    t1 = time.perf_counter()
    shewchuk_ms = elapsed_ms(t0, t1)

    if CGAL_AVAILABLE:
        # --- CGAL ---
        t2 = time.perf_counter()
        for i in range(n_tets):
            for j in range(i + 1, n_tets):
                geometry.cgal_collides_tetrahedron_tetrahedron_3d(
                    tets[i][0], tets[i][1], tets[i][2], tets[i][3],
                    tets[j][0], tets[j][1], tets[j][2], tets[j][3])
        t3 = time.perf_counter()
        cgal_ms = elapsed_ms(t2, t3)
    else:
        cgal_ms = float('nan')

    label = f"tetrahedron-tetrahedron 3D ({n_tets} tets, {n_pairs} pairs)"
    print_result(label, shewchuk_ms, cgal_ms)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print_header()

    # Use fewer iterations than C++ for orient2d/orient3d due to Python
    # function-call overhead, but use similar simplex counts for collision tests.
    bench_orient2d(1_000_000)
    bench_orient3d(1_000_000)
    bench_tri_tri_3d(400)
    bench_tet_tet_3d(300)

    print()
    print("Notes:")
    print("  Speedup > 1 means Shewchuk is faster than CGAL.")
    if not CGAL_AVAILABLE:
        print("  CGAL not available: rebuild the Python module with CGAL.")
    else:
        print("  CGAL uses Exact_predicates_inexact_constructions_kernel (EPICK).")
    print("  Python function-call overhead dominates for single-predicate tests.")
    print("  See bench/performance.cpp for a lower-overhead C++ version.")
    print()
