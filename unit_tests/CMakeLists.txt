cmake_minimum_required(VERSION 3.5)

project(nistpp_tests LANGUAGES CXX)

add_definitions(-DGTEST_LANGUAGE_CXX14)

find_package(Threads REQUIRED)
find_package(Boost COMPONENTS filesystem REQUIRED)
find_package(GTest REQUIRED)

file(GLOB_RECURSE PROJ_TESTS
   tst_frequency.h
   tst_block_frequency.h)

add_executable(${PROJECT_NAME} main.cpp ${PROJ_TESTS})

target_compile_definitions(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_defs})

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PRIVATE
        nistpp
        Threads::Threads
        GTest::GTest)

