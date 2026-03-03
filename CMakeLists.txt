cmake_minimum_required(VERSION 3.5)
project(exact-geometry-lib VERSION 1.0.0 LANGUAGES CXX)

# Options
option(BUILD_TESTING "Build tests" ON)

# Add geometry library
add_subdirectory(geometry)

# Add tests
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()
