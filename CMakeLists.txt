cmake_minimum_required (VERSION 3.12)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

project (runko LANGUAGES CUDA CXX C)
# https://cliutils.gitlab.io/modern-cmake/chapters/packages/CUDA.html

include(CheckLanguage)
check_language(CUDA)

# Add the cmake folder to get FindPkg files
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
 
# optimization flags
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -fprofile-generate -fprofile-use")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -ftree-vectorize -fopt-info-vec")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=haswell -mtune=haswell -g") #latest optimization
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -mtune=native") #general optimization


# AddressSanitizer
set(CMAKE_CXX_FLAGS_ASAN "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1")

# LeakSanitizer
set(CMAKE_CXX_FLAGS_LSAN "-fsanitize=leak -fno-omit-frame-pointer -g -O1")

# MemorySanitizer
set(CMAKE_CXX_FLAGS_MSAN "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2")



# set location of python libraries
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# set numerical precision variable
set (FP_PRECISION DP)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



##################################################
# dependencies

#FIND_PACKAGE (EIGEN3 REQUIRED NO_MODULE)
#if(APPLE)
#    INCLUDE_DIRECTORIES (SYSTEM "$ENV{EIGEN3_INCLUDE_DIR}" )
#endif()
include_directories(tools/eigen/)


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
    set (OpenMP_C_FLAGS "-fopenmp")
    set (OpenMP_CXX_FLAGS "-fopenmp")
    set (OpenMP_EXE_LINKER_FLAGS "-fopenmp")
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


#add_subdirectory(vlasov)
#add_subdirectory(em-fields)
#add_subdirectory(pic)
#add_subdirectory(io)
#add_subdirectory(tests)

add_subdirectory(bindings)
add_subdirectory(docs)




#-------------------------------------------------- 
# Documentation Target
#
#option(BUILD_DOC "Build Doxygen documentation" ON)
#
## check if Doxygen is installed
#find_package(Doxygen)
#if (DOXYGEN_FOUND)
#    # set input and output files for autogeneration
#    set(DOXYGEN_IN ${MAINFOLDER}/docs/Doxyfile.in)
#    set(DOXYGEN_OUT ${MAINFOLDER}/docs/Doxyfile)
#
#    # generate the Doxyfile
#    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#
#    # generate mainpage.md
#    configure_file(${MAINFOLDER}/docs/mainpage.md.in 
#                   ${MAINFOLDER}/docs/mainpage.md
#                  )
#    
#    add_custom_target( doc_doxygen ALL
#                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#                       COMMENT "Generating API documentation with Doxygen"
#                       VERBATIM )
#else (DOXYGEN_FOUND)
#    message(
#        "Doxygen not found, skipping generation of documentation")
#endif (DOXYGEN_FOUND)



#-------------------------------------------------- 
# unit tests
#
add_custom_target(check-runko ALL
                  ${PYTHON_EXECUTABLE} -m unittest discover -s ../tests/ -v
                  DEPENDS pycorgi pyrunko 
                  VERBATIM
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
                  )


