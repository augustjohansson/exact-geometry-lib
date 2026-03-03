cmake_minimum_required(VERSION 3.5)

# Test executable
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ConvexTriangulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/IntersectionConstruction.cpp
)

add_executable(geometry_tests ${TEST_SOURCES})
target_include_directories(geometry_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(geometry_tests PRIVATE geometry_lib)
target_compile_features(geometry_tests PRIVATE cxx_std_11)

enable_testing()
add_test(NAME geometry_tests COMMAND geometry_tests)
