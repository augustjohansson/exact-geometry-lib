cmake_minimum_required(VERSION 3.5)

find_package(pybind11 REQUIRED)

pybind11_add_module(geometry bindings.cpp)
target_link_libraries(geometry PRIVATE geometry_lib)
target_compile_features(geometry PRIVATE cxx_std_14)

# Optionally link CGAL to expose CGAL predicates in Python
find_package(CGAL QUIET)
if(CGAL_FOUND)
  target_compile_definitions(geometry PRIVATE GEOMETRY_WITH_CGAL)
  target_link_libraries(geometry PRIVATE CGAL::CGAL)
endif()
