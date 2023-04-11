cmake_minimum_required(VERSION 3.10)
project(delaunay-triangulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(include)

# Main target
add_executable(delaunay-triangulation
    src/main.cpp
    src/delaunay.cpp
)

# Tests target
# add_executable(test_delaunay
#     # test/test_delaunay.cpp
#     src/delaunay.cpp
# )

enable_testing()

# Test case (replace "/path/to" with the path to the test_delaunay executable)
# add_test(NAME test_delaunay COMMAND /path/to/test_delaunay)