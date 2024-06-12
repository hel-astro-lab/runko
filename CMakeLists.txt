cmake_minimum_required (VERSION 3.5)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

project (runko LANGUAGES CXX C)

# https://cliutils.gitlab.io/modern-cmake/chapters/packages/CUDA.html
include(CheckLanguage)
check_language(CUDA)

# Add the cmake folder to get FindPkg files
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
 
##################################################
# optimization flags

# debug mode (default for code development)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O2 -DDEBUG -march=native -ftree-vectorize -fopt-info-all-all=all.all -Wall -Wextra -pedantic -Wno-error=unknown-pragmas -Wno-error=unused-parameter -save-temps=obj")
#set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O2 -DDEBUG -march=native -ftree-vectorize -Wall -Wextra -pedantic")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -DDEBUG -march=native -ftree-vectorize")

#--------------------------------------------------
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -ftree-vectorize -fopt-info-all-all=all.all") general
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -mtune=native -ftree-vectorize")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -flto=auto -fomit-frame-pointer -fopt-info-all-all=all.all")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=znver1 -mtune=znver1 -flto=auto -fomit-frame-pointer -fopt-info-all-all=all.all")

# epyc rome
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer -fopt-info-all-all=all.all")

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


#--------------------------------------------------
# librarires

# set location of python libraries
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# set numerical precision variable
set (FP_PRECISION DP)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



##################################################
# dependencies

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) # include main repository folder so that header links work


include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})


find_package(Python COMPONENTS Interpreter Development REQUIRED)

# fmt preinstallation fails with intel compilers
#FIND_PACKAGE (FMT REQUIRED)
#INCLUDE_DIRECTORIES ( "$ENV{FMT_INCLUDE_DIR}" )

# hdf5 for i/o
#FIND_PACKAGE (HDF5 COMPONENTS CXX REQUIRED)
#FIND_PACKAGE (HDF5 COMPONENTS CXX)
FIND_PACKAGE (HDF5)

# hpc stuff 
#FIND_PACKAGE( FFTW3 ) # this is not needed
#find_package (MPI) # assumed to be provided by the compiler
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/external/corgi/mpi4cpp/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/external/corgi/mpi4cpp/tools/optional-lite/include )

##################################################
# openMP
#see:
# - https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
OPTION (USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)

# Find OpenMP
if(APPLE AND USE_OpenMP)
    set (OpenMP_C_FLAGS "-fopenmp-simd")
    set (OpenMP_CXX_FLAGS "-fopenmp-simd")
    set (OpenMP_EXE_LINKER_FLAGS "-fopenmp-simd")
endif()

if(USE_OpenMP)
  find_package(OpenMP)
endif(USE_OpenMP)

if (USE_OpenMP AND OPENMP_FOUND)
    #include_directories("${OPENMP_INCLUDES}")
    #link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
endif()



##################################################
# targets build

#-------------------------------------------------- 
# Warning flags based on the compiler
#set (WARNING_FLAGS
#    $<$<CXX_COMPILER_ID:GNU>:
#        -Wall -Wextra -Wpedantic
#    >
#    $<$<CXX_COMPILER_ID:Clang>:
#        -Wall -Wextra -Wpedantic -Wno-missing-braces -g
#    >
#)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -h gnu")

# add code coverage testing for GGC
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
#endif()

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage -fprofile-arcs)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage -fprofile-arcs)
  endif()
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")



#-------------------------------------------------- 
add_subdirectory(external/corgi)
add_subdirectory(bindings)
add_subdirectory(docs)

#-------------------------------------------------- 
# unit tests

if(DEFINED Python_EXECUTABLE AND NOT DEFINED PYTHON_EXECUTABLE)
  set(PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
endif()

message(STATUS "Runko tests ran with python: ${PYTHON_EXECUTABLE}")

add_custom_target(check-runko ALL
                  ${PYTHON_EXECUTABLE} -m unittest discover -s ../tests/ -v
                  DEPENDS pycorgi pyrunko 
                  VERBATIM
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
                  )


