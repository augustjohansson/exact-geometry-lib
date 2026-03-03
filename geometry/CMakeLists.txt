cmake_minimum_required(VERSION 3.5)

# Optional CGAL verification (for debugging)
option(GEOMETRY_ENABLE_CGAL_DEBUGGING "Enable CGAL-based verification of geometry predicates" OFF)

set(GEOMETRY_SOURCES
  Point.cpp
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

add_library(geometry STATIC ${GEOMETRY_SOURCES})
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(geometry PUBLIC cxx_std_11)

if(GEOMETRY_ENABLE_CGAL_DEBUGGING)
  target_compile_definitions(geometry PUBLIC DOLFIN_ENABLE_GEOMETRY_DEBUGGING)
  find_package(CGAL REQUIRED)
  target_link_libraries(geometry PUBLIC CGAL::CGAL)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(geometry PUBLIC Eigen3::Eigen)

install(TARGETS geometry ARCHIVE DESTINATION lib)
install(FILES ${GEOMETRY_HEADERS} DESTINATION include/geometry)
