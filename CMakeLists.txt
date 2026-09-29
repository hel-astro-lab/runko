# Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.23)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 26)
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

# Pika defines option(BUILD_SHARED_LIBS ...) which also makes tyvi a shared library.
# However, for the python module we want them to be included into runko_cpp_bindings.so
# otherwise it would not know where to find libtyvi.so or libpika.so.
#
# This has to be set with option(BUILD_SHARED_LIBS ...) as cmake docs state:
#
# > Note that if bringing external dependencies directly into the build,
# > such as with FetchContent or a direct call to add_subdirectory(),
# > and one of those dependencies has such a call to option(BUILD_SHARED_LIBS ...),
# > the top level project must also call option(BUILD_SHARED_LIBS ...)
# > before bringing in its dependencies. Failure to do so can lead to
# > different behavior between the first and subsequent CMake runs.
option(BUILD_SHARED_LIBS "Don't touch! Otherwise breaks python install." OFF)

if(${tyvi_BACKEND} STREQUAL "hip")
    message(STATUS "runko: enabling hip")
    enable_language(HIP)
endif()

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

find_package(Python COMPONENTS Interpreter Development REQUIRED)

if(NOT DEFINED ENV{pybind11_DIR})
    # Make pybind11 python install automatically findable.
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

find_package(pybind11 CONFIG REQUIRED)

add_library(runko_includes INTERFACE)
target_include_directories(runko_includes INTERFACE ${PROJECT_SOURCE_DIR}/src)

add_library(runko_warnings INTERFACE)
target_compile_options(runko_warnings
    INTERFACE
    -Wconversion
    -Wall
    -Wextra
    -Wpedantic
    -Werror
)

add_subdirectory(external/tyvi)
add_subdirectory(external/entt)

# Required for shared libraries.
set_target_properties(tyvi PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

# Wrap dependencies so that they can be marked SYSTEM libraries (suppresses warnings from them).
add_library(runko_dependencies INTERFACE)
target_link_libraries(runko_dependencies
    INTERFACE
    tyvi
    EnTT::EnTT
)
target_include_directories(runko_dependencies
    SYSTEM INTERFACE
    $<TARGET_PROPERTY:tyvi,INTERFACE_INCLUDE_DIRECTORIES>
)


add_subdirectory(src)
add_subdirectory(docs)

# If building with scikit-build-core not it editable mode disable tests.
if(DEFINED SKBUILD_PROJECT_NAME AND NOT ${SKBUILD_STATE} MATCHES "editable")
    set(RUNKO_ENABLE_TESTS OFF)
else()
    set(RUNKO_ENABLE_TESTS ON)
endif()

if(PROJECT_IS_TOP_LEVEL AND ${RUNKO_ENABLE_TESTS})
    message(STATUS "Enabling runko tests...")
    message(STATUS "Runko tests are run with python: ${Python_EXECUTABLE}")

    include(CTest)
    add_subdirectory(tests)
else()
    message(STATUS "Skipping runko tests.")
endif()
