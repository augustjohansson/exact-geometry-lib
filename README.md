# exact-geometry-lib

A standalone C++ geometry library providing:

- **Collision detection** using `orient2d` and `orient3d` exact arithmetic predicates
- **Intersection construction** for simplices (points, segments, triangles, tetrahedra)
- **Convex triangulation** via Graham scan
- **Simplex quadrature** generation (1D, 2D, 3D)

The predicates are based on Jonathan Richard Shewchuk's robust geometric predicates.

## Building

Requires:
- CMake >= 3.5
- C++11 compiler
- Eigen3 (header-only library; install system-wide via e.g. `apt-get install libeigen3-dev`)

```bash
mkdir build && cd build
cmake ..
make
ctest
```

## Running Tests

```bash
cd build
ctest -V
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

## Optional CGAL Verification

To enable CGAL-based verification of collision predicates (for debugging),
compile with `-DGEOMETRY_ENABLE_CGAL_DEBUGGING`.

## License

LGPL v3 (inherited from DOLFIN/FEniCS project)
