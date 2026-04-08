# Minimal CPU-only rocThrust config.
# Provides the roc::rocthrust INTERFACE target with thrust headers from the
# rocm-libraries source tree. No HIP, rocprim, or GPU compiler needed.

# Locate the thrust headers relative to this config file:
#   cmake/rocthrust-cpu/rocthrust-config.cmake  (this file)
#   external/tyvi/rocm-libraries/projects/rocthrust/thrust/  (headers)
get_filename_component(_ROCTHRUST_CPU_CONFIG_DIR "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)
get_filename_component(_ROCTHRUST_CPU_ROOT "${_ROCTHRUST_CPU_CONFIG_DIR}/../../external/tyvi/rocm-libraries/projects/rocthrust" ABSOLUTE)

if(NOT EXISTS "${_ROCTHRUST_CPU_ROOT}/thrust/device_vector.h")
    message(FATAL_ERROR
        "rocThrust headers not found at ${_ROCTHRUST_CPU_ROOT}/thrust/\n"
        "Run: cd external/tyvi && git clone --no-checkout --depth=1 --filter=tree:0 "
        "https://github.com/ROCm/rocm-libraries.git && cd rocm-libraries && "
        "git sparse-checkout init --cone && git sparse-checkout set projects/rocthrust && "
        "git checkout develop")
endif()

if(NOT TARGET roc::rocthrust)
    add_library(roc::rocthrust INTERFACE IMPORTED)
    set_target_properties(roc::rocthrust PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_ROCTHRUST_CPU_ROOT}"
    )
endif()

set(rocthrust_FOUND TRUE)
set(rocthrust_INCLUDE_DIR "${_ROCTHRUST_CPU_ROOT}")
set(rocthrust_INCLUDE_DIRS "${_ROCTHRUST_CPU_ROOT}")

unset(_ROCTHRUST_CPU_CONFIG_DIR)
unset(_ROCTHRUST_CPU_ROOT)
