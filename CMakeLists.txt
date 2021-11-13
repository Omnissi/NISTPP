cmake_minimum_required(VERSION 3.5)

project(NISTPP VERSION 1.0.0)

include(${CMAKE_BINARY_DIR}/conan_paths.cmake)

include(cmake/cxx_flags.cmake)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Sprout/)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/nistpp/include)

add_subdirectory(nistpp)

enable_testing()
add_subdirectory(unit_tests)

