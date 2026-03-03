cmake_minimum_required(VERSION 3.10)
project(exact-geometry-lib VERSION 1.0.0 LANGUAGES CXX)

# Options
option(BUILD_TESTING "Build tests" ON)
option(BUILD_PYTHON "Build Python bindings" ON)
option(BUILD_BENCHMARKS "Build performance benchmarks (requires CGAL)" ON)

# Optimization flags for GCC and Clang
option(GEOMETRY_NATIVE_ARCH "Enable -march=native for CPU-specific optimizations (not portable across machines)" ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  add_compile_options(-O3)
  if(GEOMETRY_NATIVE_ARCH)
    add_compile_options(-march=native)
  endif()
endif()

# Add geometry library
add_subdirectory(geometry)

# Add Python bindings
if(BUILD_PYTHON)
  add_subdirectory(python)
endif()

# Add benchmarks
if(BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()

# Add tests
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)

  # Add Python tests if Python bindings are built
  if(BUILD_PYTHON)
    find_program(PYTEST_EXECUTABLE pytest)
    if(PYTEST_EXECUTABLE)
      add_test(
        NAME python_tests
        COMMAND ${PYTEST_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/tests -v
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
      )
    endif()
  endif()
endif()
