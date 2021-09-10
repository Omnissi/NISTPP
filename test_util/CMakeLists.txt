cmake_minimum_required(VERSION 3.5)

project(test_util LANGUAGES CXX)

add_definitions(-DGTEST_LANGUAGE_CXX14)

find_package(Boost COMPONENTS filesystem REQUIRED)
find_package(Threads REQUIRED)

add_executable(${PROJECT_NAME} main.cpp ${PROJ_TESTS})

target_compile_definitions(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_defs})

target_link_libraries(${PROJECT_NAME} PRIVATE
        nistpp
        Boost::filesystem
        Threads::Threads)