# Copyright 2016 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.23)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 23)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

project(
    runko
    VERSION 5.0.0
    DESCRIPTION "Modern C++/python CPU/GPU plasma toolbox."
    HOMEPAGE_URL "https://github.com/hel-astro-lab"
    LANGUAGES CXX
)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

add_library(runko_includes INTERFACE)
target_include_directories(runko_includes INTERFACE ${PROJECT_SOURCE_DIR})

add_library(runko_warnings INTERFACE)
target_compile_options(runko_warnings
    INTERFACE
    # FIXME: make sure that these are actually used
    # (they might not be with hipcc nor with apple clang)
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic -Wno-missing-braces>
)

add_subdirectory(external/corgi)
add_subdirectory(external/tyvi)
add_subdirectory(bindings)
add_subdirectory(docs)

if (PROJECT_IS_TOP_LEVEL AND NOT DEFINED SKBUILD_PROJECT_NAME)
    message(STATUS "Enabling tests...")
    message(STATUS "Runko tests are run with python: ${PYTHON_EXECUTABLE}")

    include(CTest)
    add_subdirectory(tests)
endif ()
