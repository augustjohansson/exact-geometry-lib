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

To enable CGAL-based internal debugging of collision predicates:

```bash
cmake .. -DGEOMETRY_ENABLE_CGAL_DEBUGGING=ON
```

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
