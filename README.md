# simpex — exact geometry library

C++ library for collision detection and intersection of simplices using Shewchuk's exact arithmetic predicates.

## Building

```bash
# Debian/Ubuntu dependencies
sudo apt-get install cmake libeigen3-dev
pip install pybind11   # optional, for Python bindings

mkdir build && cd build
cmake ..
make
ctest
```

To build with Python bindings:

```bash
cmake .. -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
```

## Compiling with CGAL

CGAL enables the performance benchmark (`simpex_bench`), CGAL comparison tests (`simpex_cgal_tests`), and Python CGAL bindings.

```bash
sudo apt-get install libcgal-dev   # Debian/Ubuntu
brew install cgal                  # macOS
```

Then build normally — CMake detects CGAL automatically.

To enable CGAL-based internal verification of collision predicates:

```bash
cmake .. -DSIMPEX_ENABLE_CGAL_DEBUGGING=ON
```

## CMake Options

| Option                        | Default | Description                                                           |
|-------------------------------|---------|-----------------------------------------------------------------------|
| `SIMPEX_NATIVE_ARCH`          | ON      | Enable `-march=native` (not portable across machines).                |
| `SIMPEX_WARNINGS`             | ON      | Enable strict compiler warnings.                                      |
| `SIMPEX_ENABLE_CGAL_DEBUGGING`| ON      | Enable CGAL-based verification of geometry predicates (requires CGAL).|
| `BUILD_TESTING`               | ON      | Build unit tests.                                                     |
| `BUILD_PYTHON`                | ON      | Build Python bindings (requires pybind11).                            |
| `BUILD_BENCHMARKS`            | ON      | Build C++ performance benchmark (requires CGAL).                      |

## Tolerances

The library uses the following hard-coded tolerance values.  Users integrating
simpex into applications with extreme coordinate magnitudes should be aware of
these thresholds.

| Location | Value | Purpose |
|---|---|---|
| `geometry/ConvexTriangulation.cpp` (three places) | `3.0e-16` | Deduplication threshold: two points are considered distinct if they differ by more than this value in any coordinate dimension (`unique_points`). |
| `geometry/ConvexTriangulation.cpp` | `1e-14` | Orientation test threshold: a face normal is accepted only if `|orient3d| > 1e-14`. |
| `geometry/ConvexTriangulation.cpp` | `1e-14` | Volume verification: the reconstructed triangulation volume must match the reference volume to within this tolerance. |
| `geometry/GeometryPredicates.cpp` | `3.0e-16` | Collinearity test: two unit vectors are considered non-collinear if `-(|u·v| - 1) > 3.0e-16`. |
| `geometry/IntersectionConstruction.cpp` | `3.0e-16` | Degenerate segment test: a segment endpoint is considered to lie on the splitting plane if `|orient| < 3.0e-16`. |
| `geometry/Point.h` | `3.0e-16` | Assertion: the axis vector passed to `Point::rotate` must be a unit vector to within `3.0e-16`. |
| `geometry/CGALExactArithmetic.h` | `1e-15` | Coincident-point guard: two simplex vertices are treated as identical when converting to CGAL exact arithmetic if their distance is less than `1e-15`. |

All collision and orientation predicates (orient2d, orient3d, in-circle,
in-sphere) use Shewchuk's adaptive exact arithmetic, which carries no fixed
tolerance: results are computed to full floating-point precision with an
adaptively chosen working precision for near-degenerate cases.

## Running

```bash
ctest -V                          # all tests
./build/bench/simpex_bench        # C++ performance benchmark (requires CGAL)
SIMPEX_BUILD_DIR=build/python python3 bench/bench_performance.py
```

## Python Bindings

```python
import sys; sys.path.insert(0, "build/python")
from simpex import Point, CollisionPredicates, IntersectionConstruction

p0, p1, pt = Point(0,0,0), Point(1,0,0), Point(0.5,0,0)
print(CollisionPredicates.collides_segment_point_2d(p0, p1, pt))  # True

import simpex
if simpex.CGAL_AVAILABLE:
    print(simpex.cgal_collides_triangle_triangle_3d(p0, p1, pt, p0, p1, pt))
```

## License

The code for this library is largely directly copied from DOLFIN (https://bitbucket.org/fenics-project/dolfin/).

DOLFIN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

DOLFIN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with DOLFIN. If not, see <https://www.gnu.org/licenses/>.
