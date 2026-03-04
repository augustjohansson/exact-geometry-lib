cmake_minimum_required(VERSION 3.10)
option(SIMPEX_ENABLE_CGAL_DEBUGGING "Enable CGAL-based verification of geometry predicates" ON)
set(GEOMETRY_SOURCES
  predicates.cpp
  CollisionPredicates.cpp
  ConvexTriangulation.cpp
  GeometryPredicates.cpp
  IntersectionConstruction.cpp
  SimplexQuadrature.cpp
  SimplexQuadratureCompression.cpp
  SimplexQuadratureTables.cpp
)

set(GEOMETRY_HEADERS
  CGALExactArithmetic.h
  CollisionPredicates.h
  ConvexTriangulation.h
  GeometryPredicates.h
  GeometryTools.h
  IntersectionConstruction.h
  Point.h
  SimplexQuadrature.h
  SimplexQuadratureCompression.h
  predicates.h
)

add_library(simpex_lib STATIC ${GEOMETRY_SOURCES})
target_include_directories(simpex_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(simpex_lib PUBLIC cxx_std_11)
set_target_properties(simpex_lib PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME simpex)

if(SIMPEX_ENABLE_CGAL_DEBUGGING)
  target_compile_definitions(simpex_lib PUBLIC SIMPEX_ENABLE_GEOMETRY_DEBUGGING)
  find_package(CGAL REQUIRED)
  target_link_libraries(simpex_lib PUBLIC CGAL::CGAL)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(simpex_lib PUBLIC Eigen3::Eigen)

install(TARGETS simpex_lib ARCHIVE DESTINATION lib)
install(FILES ${GEOMETRY_HEADERS} DESTINATION include/simpex)
