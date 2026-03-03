cmake_minimum_required(VERSION 3.10)

find_package(pybind11 REQUIRED)

pybind11_add_module(simpex bindings.cpp)
target_link_libraries(simpex PRIVATE simpex_lib)
target_compile_features(simpex PRIVATE cxx_std_14)

# Optionally link CGAL to expose CGAL predicates in Python
find_package(CGAL QUIET)
if(CGAL_FOUND)
  target_compile_definitions(simpex PRIVATE GEOMETRY_WITH_CGAL)
  target_link_libraries(simpex PRIVATE CGAL::CGAL)
endif()
