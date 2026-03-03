cmake_minimum_required(VERSION 3.5)

find_package(CGAL QUIET)

if(NOT CGAL_FOUND)
  message(STATUS "CGAL not found: skipping bench/geometry_bench")
  return()
endif()

add_executable(geometry_bench performance.cpp)
target_link_libraries(geometry_bench PRIVATE geometry_lib CGAL::CGAL)
target_compile_features(geometry_bench PRIVATE cxx_std_14)
