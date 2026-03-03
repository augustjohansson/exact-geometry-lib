cmake_minimum_required(VERSION 3.10)

set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CollisionPredicates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ConvexTriangulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/IntersectionConstruction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/SimplexQuadrature.cpp
)

add_executable(simpex_tests ${TEST_SOURCES})
target_include_directories(simpex_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simpex_tests PRIVATE simpex_lib)
target_compile_features(simpex_tests PRIVATE cxx_std_11)

enable_testing()
add_test(NAME simpex_tests COMMAND simpex_tests)
