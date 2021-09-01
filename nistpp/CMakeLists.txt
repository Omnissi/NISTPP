cmake_minimum_required(VERSION 3.5)

project(nistpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(boost COMPONENTS math REQUIRED)

file(GLOB_RECURSE PROJ_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/*.h*)
file(GLOB_RECURSE PROJ_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/*.c*)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(${PROJECT_NAME} STATIC ${PROJ_HEADER} ${PROJ_SOURCE})

target_link_libraries(${PROJECT_NAME} PUBLIC Boost::boost)
