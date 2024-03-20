cmake_minimum_required (VERSION 3.12)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# https://cliutils.gitlab.io/modern-cmake/chapters/packages/CUDA.html
project (runko LANGUAGES CXX C)

include(CheckLanguage)
check_language(CUDA)

# Add the cmake folder to get FindPkg files
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
 
##################################################
# optimization flags

#--------------------------------------------------
# debug mode (default)
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -DDEBUG -march=native -ftree-vectorize -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS "-g -O2 -DDEBUG -march=native -ftree-vectorize")

#--------------------------------------------------
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -march=native -ftree-vectorize -fopt-info-all-all=all.all") general
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -march=native -ftree-vectorize -fopt-info-all-all=all.all")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -flto=auto -fomit-frame-pointer -fopt-info-all-all=all.all")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -march=znver1 -mtune=znver1 -flto=auto -fomit-frame-pointer -fopt-info-all-all=all.all")

# epyc rome
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer -fopt-info-all-all=all.all")

#--------------------------------------------------
# other compiliation modes

# AddressSanitizer
set(CMAKE_CXX_FLAGS_ASAN "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1")

# LeakSanitizer
set(CMAKE_CXX_FLAGS_LSAN "-fsanitize=leak -fno-omit-frame-pointer -g -O1")

# MemorySanitizer
set(CMAKE_CXX_FLAGS_MSAN "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2")

#--------------------------------------------------

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Flags (debug): ${CMAKE_CXX_FLAGS}")
message(STATUS "Flags RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")



#--------------------------------------------------
# librarires

# set location of python libraries
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# set numerical precision variable
set (FP_PRECISION DP)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



##################################################
# dependencies

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# fmt preinstallation fails with intel compilers
#FIND_PACKAGE (FMT REQUIRED)
#INCLUDE_DIRECTORIES ( "$ENV{FMT_INCLUDE_DIR}" )

# hdf5 for i/o
#FIND_PACKAGE (HDF5 COMPONENTS CXX REQUIRED)
#FIND_PACKAGE (HDF5 COMPONENTS CXX)
FIND_PACKAGE (HDF5)

#FIND_PACKAGE( FFTW3 ) # this does not exists

# hpc stuff 
#find_package (MPI)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/corgi/mpi4cpp/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/corgi/mpi4cpp/tools/optional-lite/include )

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
add_subdirectory(corgi)
add_subdirectory(bindings)
add_subdirectory(docs)

#-------------------------------------------------- 
# unit tests
add_custom_target(check-runko ALL
                  ${PYTHON_EXECUTABLE} -m unittest discover -s ../tests/ -v
                  DEPENDS pycorgi pyrunko 
                  VERBATIM
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
                  )


