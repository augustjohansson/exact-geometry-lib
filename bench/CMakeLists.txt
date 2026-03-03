cmake_minimum_required(VERSION 3.10)

find_package(CGAL QUIET)

if(NOT CGAL_FOUND)
  message(STATUS "CGAL not found: skipping bench/simpex_bench")
  return()
endif()

add_executable(simpex_bench performance.cpp)
target_link_libraries(simpex_bench PRIVATE simpex_lib CGAL::CGAL)
target_compile_features(simpex_bench PRIVATE cxx_std_14)
