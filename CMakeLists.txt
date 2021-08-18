cmake_minimum_required(VERSION 3.5)

project(NISTPP VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${CMAKE_BINARY_DIR}/conan_paths.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/nistpp/include)

add_subdirectory(nistpp)

enable_testing()
add_subdirectory(unit_tests)

