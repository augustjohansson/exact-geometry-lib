# simpex — exact geometry library

A standalone C++ geometry library providing:

- **Collision detection** using `orient2d` and `orient3d` exact arithmetic predicates
- **Intersection construction** for simplices (points, segments, triangles, tetrahedra)
- **Convex triangulation** via Graham scan
- **Simplex quadrature** generation (1D, 2D, 3D)
- **Python bindings** via pybind11
- **Performance benchmarks** comparing Shewchuk's predicates vs CGAL (EPICK)

The predicates are based on Jonathan Richard Shewchuk's robust geometric predicates.

## Building

Requires:
- CMake >= 3.10
- C++11 compiler (C++14 for bench/Python bindings)
- Eigen3 (header-only; install via `apt-get install libeigen3-dev`)
- pybind11 (optional, for Python bindings; install via `pip install pybind11`)
- CGAL (optional; install via `apt-get install libcgal-dev`)

### Basic build (no Python, no CGAL)

```bash
mkdir build && cd build
cmake .. -DBUILD_PYTHON=OFF -DBUILD_BENCHMARKS=OFF
make
ctest
```

### Build with Python bindings

```bash
mkdir build && cd build
cmake .. -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
make
ctest
```

### Compiling with CGAL

CGAL enables:
- `simpex_bench` — C++ performance benchmark (Shewchuk vs CGAL timing + correctness)
- `simpex_cgal_tests` — C++ Catch2 tests verifying CGAL and Shewchuk agree on every predicate
- CGAL Python bindings (`simpex.cgal_orient2d`, `simpex.cgal_collides_*`, etc.)

Install CGAL first:

```bash
# Debian/Ubuntu
sudo apt-get install libcgal-dev

# macOS (Homebrew)
brew install cgal
```

Then build:

```bash
mkdir build && cd build
cmake .. -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
make
ctest          # runs simpex_tests + simpex_cgal_tests + python_tests
```

CMake automatically detects CGAL. You can check whether it was found:

```bash
cmake .. 2>&1 | grep -i cgal
```

If CGAL is found you will see the `simpex_bench` and `simpex_cgal_tests` targets built.
If CGAL is not found, those targets are silently skipped.

## Running Tests

```bash
cd build
ctest -V
```

To run only the Python tests directly:

```bash
SIMPEX_BUILD_DIR=build/python python3 -m pytest python/tests -v
```

To run CGAL vs Shewchuk comparison tests (requires CGAL):

```bash
./build/test/simpex_cgal_tests
```

## Python Bindings

After building, the Python module `simpex` is available in the `build/python/` directory.

```python
import sys
sys.path.insert(0, "build/python")
from simpex import Point, CollisionPredicates, IntersectionConstruction
from simpex import ConvexTriangulation, SimplexQuadrature

# Create a point
p = Point(1.0, 2.0, 3.0)

# Collision detection
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(1.0, 0.0, 0.0)
pt = Point(0.5, 0.0, 0.0)
print(CollisionPredicates.collides_segment_point_2d(p0, p1, pt))  # True

# CGAL predicates (when built with CGAL)
import simpex
if simpex.CGAL_AVAILABLE:
    print(simpex.cgal_orient2d(p0, p1, pt))
    print(simpex.cgal_collides_triangle_triangle_3d(p0, p1, pt, p0, p1, pt))
    print(simpex.cgal_collides_tetrahedron_point_3d(p0, p1, pt, p, pt))
```

## Library Contents

| File | Description |
|------|-------------|
| `Point.h` | 3D point class |
| `predicates.h/cpp` | orient2d, orient3d exact predicates (Shewchuk) |
| `CGALPredicates.h` | CGAL EPICK predicates (orient2d/3d, collides_*, requires CGAL) |
| `CollisionPredicates.h/cpp` | Pairwise collision detection |
| `IntersectionConstruction.h/cpp` | Intersection computation |
| `ConvexTriangulation.h/cpp` | Convex triangulation |
| `SimplexQuadrature.h/cpp` | Quadrature rules for simplices |
| `GeometryPredicates.h/cpp` | Degeneracy checks |
| `GeometryTools.h` | Geometry utilities |
| `GeometryDebugging.h/cpp` | Debugging utilities |
| `CGALExactArithmetic.h` | Optional CGAL verification (disabled by default) |
| `python/bindings.cpp` | pybind11 Python bindings |
| `python/tests/` | Python test suite (138 tests) |
| `test/` | C++ Catch2 test suite |
| `test/CGALComparison.cpp` | C++ tests comparing CGAL vs Shewchuk (requires CGAL) |
| `bench/performance.cpp` | C++ performance benchmark (Shewchuk vs CGAL) |
| `bench/bench_performance.py` | Python performance benchmark |

## Performance Benchmarks

The `bench/` directory contains benchmarks comparing Shewchuk's exact predicates against
CGAL's `Exact_predicates_inexact_constructions_kernel` (EPICK).  Each benchmark:
1. Times both implementations on the same deterministic random data (seed 42)
2. Checks correctness — verifies Shewchuk and CGAL give identical results
3. Prints a per-test PASS/FAIL summary at the end

**C++ benchmark** (built automatically when CGAL is available):

```bash
./build/bench/simpex_bench
```

**Python benchmark**:

```bash
SIMPEX_BUILD_DIR=build/python python3 bench/bench_performance.py
```

## CGAL vs Shewchuk Agreement

The `simpex_cgal_tests` Catch2 executable (built when CGAL is found) verifies that the
Shewchuk-based predicates and CGAL EPICK predicates give identical results on:

| Test | Description |
|------|-------------|
| `orient2d` | 10 000 random 2D point triples |
| `orient3d` | 10 000 random 3D point quadruples |
| `collides_triangle_point_3d` | 10 000 random triangle-point pairs |
| `collides_tetrahedron_point_3d` | 5 000 random tet-point pairs |
| `collides_triangle_triangle_3d` | ~45 000 random triangle pairs + canonical cases |
| `collides_tetrahedron_tetrahedron_3d` | ~11 000 random tet pairs + canonical cases |

All random generators are seeded with 42 for full reproducibility.

## Optional CGAL Verification

To enable CGAL-based verification of collision predicates (for debugging),
compile with `-DGEOMETRY_ENABLE_CGAL_DEBUGGING=ON`.

## License

LGPL v3 (inherited from DOLFIN/FEniCS project)
