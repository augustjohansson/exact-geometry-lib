cmake_minimum_required(VERSION 3.10)
project(simpex VERSION 1.0.0 LANGUAGES CXX)

# Options
option(BUILD_TESTING "Build tests" ON)
option(BUILD_PYTHON "Build Python bindings" ON)
option(BUILD_BENCHMARKS "Build performance benchmarks (requires CGAL)" ON)
option(SIMPEX_NATIVE_ARCH "Enable -march=native for CPU-specific optimizations (not portable across machines)" ON)
option(SIMPEX_WARNINGS "Enable strict compiler warnings" ON)

# Optimization and warning flags for GCC and Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  add_compile_options(-O3)
  if(SIMPEX_NATIVE_ARCH)
    add_compile_options(-march=native)
  endif()
  if(SIMPEX_WARNINGS)
    add_compile_options(
      -Wall
      -Wextra
      -Wshadow
      -Wconversion
      -Wsign-conversion
      -Wdouble-promotion
      -Wunused
      -pedantic
    )
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
    find_program(PYTHON3_EXECUTABLE NAMES python3 python)
    if(PYTHON3_EXECUTABLE)
      add_test(
        NAME python_tests
        COMMAND ${PYTHON3_EXECUTABLE} -m pytest ${CMAKE_SOURCE_DIR}/python/tests -v
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
      )
      set_tests_properties(python_tests PROPERTIES
        ENVIRONMENT "SIMPEX_BUILD_DIR=${CMAKE_BINARY_DIR}/python")
    endif()
  endif()
endif()
