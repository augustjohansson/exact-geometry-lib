# exact-geometry-lib

A standalone C++ geometry library providing:

- **Collision detection** using `orient2d` and `orient3d` exact arithmetic predicates
- **Intersection construction** for simplices (points, segments, triangles, tetrahedra)
- **Convex triangulation** via Graham scan
- **Simplex quadrature** generation (1D, 2D, 3D)
- **Python bindings** via pybind11

The predicates are based on Jonathan Richard Shewchuk's robust geometric predicates.

## Building

Requires:
- CMake >= 3.5
- C++11 compiler
- Eigen3 (header-only library; install system-wide via e.g. `apt-get install libeigen3-dev`)
- pybind11 (optional, for Python bindings; install via `pip install pybind11`)

```bash
mkdir build && cd build
cmake .. -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
make
ctest
```

To build without Python bindings:

```bash
cmake .. -DBUILD_PYTHON=OFF
make
ctest
```

## Running Tests

```bash
cd build
ctest -V
```

To run only the Python tests directly:

```bash
cd build/python
python3 -m pytest ../../python/tests -v
```

## Python Bindings

After building, the Python module `geometry` is available in the `build/python/` directory.

```python
import sys
sys.path.insert(0, "build/python")
from geometry import Point, CollisionPredicates, IntersectionConstruction
from geometry import ConvexTriangulation, SimplexQuadrature

# Create a point
p = Point(1.0, 2.0, 3.0)

# Collision detection
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(1.0, 0.0, 0.0)
pt = Point(0.5, 0.0, 0.0)
print(CollisionPredicates.collides_segment_point_2d(p0, p1, pt))  # True

# Intersection construction
res = IntersectionConstruction.intersection_triangle_segment_3d(...)

# Convex triangulation
tri = ConvexTriangulation.triangulate_graham_scan_3d(points)

# Simplex quadrature
sq = SimplexQuadrature(2, 3)  # triangle, order 3
pts, wts = sq.compute_quadrature_rule_triangle(coords, 2)
```

## Library Contents

| File | Description |
|------|-------------|
| `Point.h/cpp` | 3D point class |
| `predicates.h/cpp` | orient2d, orient3d exact predicates |
| `CollisionPredicates.h/cpp` | Pairwise collision detection |
| `IntersectionConstruction.h/cpp` | Intersection computation |
| `ConvexTriangulation.h/cpp` | Convex triangulation |
| `SimplexQuadrature.h/cpp` | Quadrature rules for simplices |
| `GeometryPredicates.h/cpp` | Degeneracy checks |
| `GeometryTools.h` | Geometry utilities |
| `GeometryDebugging.h/cpp` | Debugging utilities |
| `CGALExactArithmetic.h` | Optional CGAL verification (disabled by default) |
| `python/bindings.cpp` | pybind11 Python bindings |
| `python/tests/` | Python test suite |

## Optional CGAL Verification

To enable CGAL-based verification of collision predicates (for debugging),
compile with `-DGEOMETRY_ENABLE_CGAL_DEBUGGING`.

## License

LGPL v3 (inherited from DOLFIN/FEniCS project)

