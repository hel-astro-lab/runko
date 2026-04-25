# Cray MPICH GPU Transport Layer linkage for LUMI.
#
# Wired in by the `lumi-gpu-release` preset via CMAKE_PROJECT_runko_INCLUDE,
# which CMake invokes immediately after `project(runko)`. Keeps the top-level
# CMakeLists.txt and bindings/CMakeLists.txt vendor-neutral.
#
# Why not CMAKE_*_LINKER_FLAGS in the preset:
#   Cray's CC wrapper drops `-Wl,-rpath` and `-l` items from those preset
#   variables, so the GTL never reaches the link line. Calling
#   `link_libraries(<full .so path>)` puts the library on the link line as a
#   positional argument that the wrapper preserves verbatim.

if(NOT DEFINED ENV{CRAY_MPICH_ROOTDIR})
    message(FATAL_ERROR
        "lumi-gtl.cmake: CRAY_MPICH_ROOTDIR is unset. Load the cray-mpich "
        "module before running `cmake --preset lumi-gpu-release`.")
endif()

find_library(MPI_GTL_HSA_LIB
    NAMES mpi_gtl_hsa
    PATHS "$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib"
    NO_DEFAULT_PATH
    REQUIRED
)

# Top-level link_libraries propagates to every target in subdirs added after
# this point, which is all of them. Static libs (corgi, tyvi) ignore it; the
# pybind11 shared module (runko_cpp_bindings) picks it up.
link_libraries("${MPI_GTL_HSA_LIB}")
