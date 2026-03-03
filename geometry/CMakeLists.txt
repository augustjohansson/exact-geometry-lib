cmake_minimum_required(VERSION 3.10)
option(GEOMETRY_ENABLE_CGAL_DEBUGGING "Enable CGAL-based verification of geometry predicates" OFF)

set(GEOMETRY_SOURCES
  predicates.cpp
  CollisionPredicates.cpp
  ConvexTriangulation.cpp
  GeometryDebugging.cpp
  GeometryPredicates.cpp
  IntersectionConstruction.cpp
  SimplexQuadrature.cpp
)

set(GEOMETRY_HEADERS
  CGALExactArithmetic.h
  CollisionPredicates.h
  ConvexTriangulation.h
  GeometryDebugging.h
  GeometryPredicates.h
  GeometryTools.h
  IntersectionConstruction.h
  Point.h
  SimplexQuadrature.h
  predicates.h
)

add_library(geometry_lib STATIC ${GEOMETRY_SOURCES})
target_include_directories(geometry_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(geometry_lib PUBLIC cxx_std_11)
set_target_properties(geometry_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(GEOMETRY_ENABLE_CGAL_DEBUGGING)
  target_compile_definitions(geometry_lib PUBLIC DOLFIN_ENABLE_GEOMETRY_DEBUGGING)
  find_package(CGAL REQUIRED)
  target_link_libraries(geometry_lib PUBLIC CGAL::CGAL)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(geometry_lib PUBLIC Eigen3::Eigen)

install(TARGETS geometry_lib ARCHIVE DESTINATION lib)
install(FILES ${GEOMETRY_HEADERS} DESTINATION include/geometry)
