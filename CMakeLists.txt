cmake_minimum_required(VERSION 3.14)
project(gslib_cpp VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(GSLIB_BUILD_TESTS "Build tests" ON)
option(GSLIB_BUILD_PYTHON "Build Python bindings" OFF)

# Core static library
add_library(gslib_core STATIC
    src/setrot.cpp
    src/sqdist.cpp
    src/cova3.cpp
    src/gamv.cpp
    src/kt3d.cpp
    src/declus.cpp
    src/sort.cpp
)
target_include_directories(gslib_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/extern
)

# OpenMP (optional)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(gslib_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# Python bindings
if(GSLIB_BUILD_PYTHON)
    find_package(pybind11 QUIET)
    if(NOT pybind11_FOUND)
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt)
            add_subdirectory(extern/pybind11)
        else()
            message(FATAL_ERROR "pybind11 not found. Install it or add as submodule in extern/pybind11")
        endif()
    endif()
    pybind11_add_module(gslib_cpp python/bindings.cpp)
    target_link_libraries(gslib_cpp PRIVATE gslib_core)
endif()

# Tests
if(GSLIB_BUILD_TESTS)
    enable_testing()

    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_executable(test_gslib
        tests/test_sqdist.cpp
        tests/test_cova3.cpp
        tests/test_gamv.cpp
        tests/test_kt3d.cpp
        tests/test_declus.cpp
        tests/test_numerical.cpp
    )
    target_link_libraries(test_gslib PRIVATE gslib_core GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(test_gslib)
endif()
