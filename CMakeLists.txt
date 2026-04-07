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
find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(external/corgi)

# Cray GTL library for GPU-aware MPI on LUMI (set via preset)
set(GTL_LIB_PATH "" CACHE PATH "Path to Cray GTL library (set via preset on LUMI)")
if(TARGET mpi4cpp AND GTL_LIB_PATH AND EXISTS "${GTL_LIB_PATH}")
    target_link_directories(mpi4cpp PUBLIC ${GTL_LIB_PATH})
    target_link_libraries(mpi4cpp PUBLIC mpi_gtl_hsa)
endif()

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


