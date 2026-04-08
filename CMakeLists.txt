cmake_minimum_required (VERSION 3.23)

set (CMAKE_CXX_STANDARD 23)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project (
    runko
    VERSION 5.0.0
    DESCRIPTION "Modern C++/python CPU/GPU plasma toolbox."
    HOMEPAGE_URL "https://github.com/hel-astro-lab"
    LANGUAGES CXX
)

# Add the cmake folder to get FindPkg files
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

##################################################
# Dependencies

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# HDF5: try standard detection first, fall back to manual search via HDF5_ROOT.
# hipcc's --offload-arch flag breaks FindHDF5's try_link tests on Cray.
find_package(HDF5 QUIET)
if(NOT HDF5_FOUND OR NOT HDF5_LIBRARIES)
    find_path(HDF5_INCLUDE_DIR hdf5.h HINTS ENV HDF5_ROOT PATH_SUFFIXES include REQUIRED)
    find_library(HDF5_LIBRARY hdf5 HINTS ENV HDF5_ROOT PATH_SUFFIXES lib REQUIRED)
    find_library(HDF5_HL_LIBRARY hdf5_hl HINTS ENV HDF5_ROOT PATH_SUFFIXES lib)
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    set(HDF5_LIBRARIES ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY})
    set(HDF5_FOUND TRUE)
endif()
if(NOT TARGET hdf5::hdf5)
    add_library(hdf5::hdf5 INTERFACE IMPORTED)
    target_include_directories(hdf5::hdf5 INTERFACE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(hdf5::hdf5 INTERFACE ${HDF5_LIBRARIES})
endif()

##################################################

add_subdirectory(external/corgi)
add_subdirectory(external/tyvi)

##################################################
# Targets

# Use phony target for warnings and include directories.
add_library(runko_includes INTERFACE)
target_include_directories (runko_includes INTERFACE ${PROJECT_SOURCE_DIR})

add_library(runko_warnings INTERFACE)
target_compile_options(runko_warnings
    INTERFACE
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic -Wno-missing-braces>
)


# Code Coverage Configuration
add_library(runko_coverage_config INTERFACE)

option(RUNKO_CODE_COVERAGE "Enable coverage reporting" OFF)

if(RUNKO_CODE_COVERAGE AND NOT (CMAKE_BUILD_TYPE MATCHES "Debug"))
    message(WARNING "Code coverage enabled in non-debug build!")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(runko_coverage_config INTERFACE --coverage)
    target_link_options(runko_coverage_config INTERFACE --coverage -fprofile-arcs)
endif()

add_subdirectory(bindings)
add_subdirectory(docs)

#--------------------------------------------------
# unit tests

message(STATUS "Runko tests ran with python: ${Python_EXECUTABLE}")

include(CTest)
add_subdirectory(tests/multirank_tests)

add_custom_target(check-runko ALL
                  ${Python_EXECUTABLE} -m unittest discover -s ../tests/ -v
                  DEPENDS pycorgi runko_cpp_bindings
                  VERBATIM
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests
)


