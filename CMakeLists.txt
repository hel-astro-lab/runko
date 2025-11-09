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
# optimization flags
#
# FIXME: move these to cmake presets

#--------------------------------------------------
# debug mode (default for code development)

#set(CMAKE_CXX_FLAGS_DEBUG "-g2 -O2 -DDEBUG -march=native -ftree-vectorize -fopt-info-all-all=all.all -Wall -Wextra -pedantic -Wno-error=unknown-pragmas -Wno-error=unused-parameter -save-temps=obj")
#set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O2 -DDEBUG -march=native -ftree-vectorize -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-fopenmp -g -Og -DDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -DDEBUG -march=native -ftree-vectorize")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -DDEBUG -march=znver3 -mtune=znver3 -ftree-vectorize")

#--------------------------------------------------
# release mode flags

set(CMAKE_CXX_FLAGS_RELEASE "-fopenmp -O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -flto -ffp=4 -march=znver3 -mtune=znver3 -fopenmp -fsave-loopmark") // cray compiler flags
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp") # common gcc
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -flto -ffp=4 -march=znver3 -mtune=znver3 -fopenmp -fsave-loopmark") # aggressive cray compiler flags
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=znver3 -mtune=znver3 -ftree-vectorize") # cray compiler flags
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -flto=auto -fomit-frame-pointer -fopt-info-all-all=all.all") #rome w/ output

#--------------------------------------------------
# other compiliation modes

# AddressSanitizer
#set(CMAKE_CXX_FLAGS "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g3 -O2 -march=native -ftree-vectorize")
set(CMAKE_CXX_FLAGS_ASAN "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g3 -O2 -march=native -ftree-vectorize")

# LeakSanitizer
set(CMAKE_CXX_FLAGS_LSAN "-fsanitize=leak -fno-omit-frame-pointer -g3 -O2")

# MemorySanitizer
set(CMAKE_CXX_FLAGS_MSAN "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g3 -O2")

#--------------------------------------------------

##################################################
# Build type

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Flags (debug): ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    message(STATUS "Flags (release): ${CMAKE_CXX_FLAGS_RELEASE}")
endif()


##################################################
# Dependencies

find_package (Python COMPONENTS Interpreter Development REQUIRED)

# -fopenmp added manually.
# This is required to fix missing symbol errors on LUMI.
# Symbols in question: __kmpc_for_static_*
# Potentially relevant: https://stackoverflow.com/questions/60567645/undefined-symbol-in-clang-compiled-python-extension-using-openmp
# find_package (OpenMP REQUIRED)

find_package (HDF5 REQUIRED)
# This is not directly needed in runko as we need it only for tests,
# but mpi4cpp will anyways need this, so we require it here.
find_package (MPI REQUIRED)

# Corgi has git submodules pybind11 and MPI/mpi4cpp.
# Mpi4cpp should find MPI, which added to targets that link with corgi target.
#
# If corgi happens to be replaced, then pybind11 and MPI # have to be added explicitly here.
add_subdirectory(external/corgi)

add_subdirectory(external/tyvi)

##################################################
# Targets

# Use phony target for warnings and include directories.
add_library(runko_includes INTERFACE)
target_include_directories (runko_includes INTERFACE ${PROJECT_SOURCE_DIR})

add_library(runko_warnings INTERFACE)

# FIXME: actually add the warnings to runko_warnings and fix them.
# target_compile_options(runko_warnings
#     INTERFACE
#     $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
#     $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic -Wno-missing-braces>
# )


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


