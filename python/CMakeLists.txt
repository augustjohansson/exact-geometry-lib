cmake_minimum_required(VERSION 3.5)

find_package(pybind11 REQUIRED)

pybind11_add_module(geometry bindings.cpp)
target_link_libraries(geometry PRIVATE geometry_lib)
target_compile_features(geometry PRIVATE cxx_std_11)
