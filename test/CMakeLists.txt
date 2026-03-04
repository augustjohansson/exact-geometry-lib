cmake_minimum_required(VERSION 3.10)

set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_CollisionPredicates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_ConvexTriangulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_IntersectionConstruction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_SimplexQuadrature.cpp
)

add_executable(simpex_tests ${TEST_SOURCES})
target_include_directories(simpex_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simpex_tests PRIVATE simpex_lib)
target_compile_features(simpex_tests PRIVATE cxx_std_11)

enable_testing()
add_test(NAME simpex_tests COMMAND simpex_tests)

# ---------------------------------------------------------------------------
# CGAL comparison tests (optional — only built when CGAL is found)
# Verifies that CGAL (EPICK) and Shewchuk predicates agree on random inputs.
# ---------------------------------------------------------------------------
find_package(CGAL QUIET)
if(CGAL_FOUND)
  add_executable(simpex_cgal_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_CGALComparison.cpp
  )
  target_include_directories(simpex_cgal_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(simpex_cgal_tests PRIVATE simpex_lib CGAL::CGAL)
  target_compile_features(simpex_cgal_tests PRIVATE cxx_std_14)
  target_compile_definitions(simpex_cgal_tests PRIVATE GEOMETRY_WITH_CGAL)
  add_test(NAME simpex_cgal_tests COMMAND simpex_cgal_tests)
else()
  message(STATUS "CGAL not found: skipping CGAL comparison tests (simpex_cgal_tests)")
endif()
