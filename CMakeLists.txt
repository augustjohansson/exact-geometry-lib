cmake_minimum_required(VERSION 3.5)
project(exact-geometry-lib VERSION 1.0.0 LANGUAGES CXX)

# Options
option(BUILD_TESTING "Build tests" ON)
option(BUILD_PYTHON "Build Python bindings" ON)

# Add geometry library
add_subdirectory(geometry)

# Add Python bindings
if(BUILD_PYTHON)
  add_subdirectory(python)
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
