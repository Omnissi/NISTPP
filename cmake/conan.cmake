include(conan_cmake_wrapper)

if(NOT EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    conan_cmake_run(CONANFILE conanfile.py
        BUILD missing
    )
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(KEEP_RPATHS)
conan_define_targets()

conan_set_std()
conan_set_libcxx()
if(WIN32)
    conan_set_vs_runtime()
endif()

conan_check_compiler()

string(REPLACE "." ";" VERSION_LIST ${CONAN_PACKAGE_VERSION})
list(GET VERSION_LIST 0 MAJOR)
list(GET VERSION_LIST 1 MINOR)
list(GET VERSION_LIST 2 PATCH)
list(GET VERSION_LIST 3 TWEAK)

set(${PROJECT_NAME}_VERSION_MAJOR ${MAJOR})
set(${PROJECT_NAME}_VERSION_MINOR ${MINOR})
set(${PROJECT_NAME}_VERSION_PATCH ${PATCH})
set(${PROJECT_NAME}_VERSION_TWEAK ${TWEAK})
set(${PROJECT_NAME}_VERSION ${CONAN_PACKAGE_VERSION})

include_directories(${CONAN_INCLUDE_DIRS})
link_directories(${CONAN_LIBS})
