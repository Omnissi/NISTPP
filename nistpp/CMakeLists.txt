cmake_minimum_required(VERSION 3.5)

project(nistpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (WIN32)
    find_package(boost COMPONENTS math REQUIRED)
ELSE()
    find_package(Boost COMPONENTS math REQUIRED)
ENDIF()
find_package(kissfft REQUIRED)

file(GLOB_RECURSE PROJ_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/*.h*)
file(GLOB_RECURSE PROJ_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/*.c*)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(${PROJECT_NAME} STATIC ${PROJ_HEADER} ${PROJ_SOURCE})

target_link_libraries(${PROJECT_NAME} PUBLIC Boost::boost kissfft::kissfft-double-openmp)

set(HEADERS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(${PROJECT_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${HEADERS_PATH}>
        $<INSTALL_INTERFACE:include>)

#install(DIRECTORY ${CMAKE_SOURCE_DIR}/nistpp/include/
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

file(GLOB_RECURSE PUBLIC_INCLUDE_FILES "include/*.h*")

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${PUBLIC_INCLUDE_FILES}")

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include/nistpp
        )


install(EXPORT ${PROJECT_NAME} NAMESPACE ${PROJECT_NAME}:: DESTINATION cmake)