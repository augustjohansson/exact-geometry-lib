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

## Files

| File | Description |
|------|-------------|
| `geometry/CGALPredicates.h` | CGAL EPICK collision predicates (requires CGAL) |
| `geometry/CGALExactArithmetic.h` | CGAL exact-arithmetic verification (requires `-DGEOMETRY_ENABLE_CGAL_DEBUGGING=ON`) |
| `geometry/CollisionPredicates.h/cpp` | Collision detection (Shewchuk) |
| `geometry/predicates.h/cpp` | orient2d, orient3d |
| `test/CGALComparison.cpp` | Catch2 tests: CGAL vs Shewchuk agreement |
| `bench/performance.cpp` | C++ timing + correctness benchmark |

## License

LGPL v3
