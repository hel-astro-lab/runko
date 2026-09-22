# Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

set (TOOLS_FILES
     ./runko/bindings/pytools.c++
     ./runko/bindings/pyactions.c++
     ./runko/tools/config_parser.c++
     ./runko/tools/system.c++
     ./runko/simulation_context.c++
     ./runko/communication_common.c++
     ./runko/comm/local.c++
     ./runko/comm/external.c++
     ./runko/runtime.c++
     ./runko/coords.c++
     ./runko/actions/env.c++
     ./runko/actions/sim_env.c++
     ./runko/actions/emf.c++
     ./runko/actions/pic.c++
     )

set (EMF_FILES
     ./runko/bindings/pyemf.c++
     ./runko/emf/tile.c++
     ./runko/emf/virtual_tile.c++
     ./runko/emf/yee_lattice.c++
     ./runko/emf/yee_lattice_fdtd2.c++
     ./runko/emf/yee_lattice_stencil.c++
     ./runko/emf/yee_lattice_current_filter_binomial2.c++
     )

set (PIC_FILES
     ./runko/bindings/pypic.c++
     ./runko/pic/tile.c++
     ./runko/pic/tile_communication.c++
     ./runko/pic/virtual_tile.c++
     ./runko/pic/particle.c++
     ./runko/pic/particle_current_zigzag_1st.c++
     ./runko/pic/reflector_wall.c++
     )

set (IO_FILES
     ./runko/io/snapshots/mpiio_fields.c++
     ./runko/io/snapshots/mpiio_particles.c++
     ./runko/io/snapshots/mpiio_spectra.c++
    )

set (RUNKO_CPP_BINDINGS_FILES
     "./runko/bindings/runko_cpp_bindings.c++"
     ${TOOLS_FILES}
     ${EMF_FILES}
     ${PIC_FILES}
     ${IO_FILES}
    )


add_library(runko_cpp_bindings MODULE)
target_sources(runko_cpp_bindings PRIVATE ${RUNKO_CPP_BINDINGS_FILES})
target_link_libraries(runko_cpp_bindings PRIVATE pybind11::module)

pybind11_extension(runko_cpp_bindings)
if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS.
    pybind11_strip(runko_cpp_bindings)
endif()

set_target_properties(runko_cpp_bindings PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                         HIP_VISIBILITY_PRESET "hidden"
                                         CUDA_VISIBILITY_PRESET "hidden")

target_link_libraries(runko_cpp_bindings
    PRIVATE
    runko_includes
    runko_warnings
    runko_dependencies
)


if(${tyvi_BACKEND} STREQUAL "hip")
     foreach(file ${RUNKO_CPP_BINDINGS_FILES})
         set_source_files_properties("${file}" PROPERTIES LANGUAGE HIP)
     endforeach()
     set_property(TARGET runko_cpp_bindings PROPERTY HIP_STANDARD 26)
endif()

# Install runko_cpp_bindings.
install(TARGETS runko_cpp_bindings DESTINATION .)
# Vendor-in corgi if we are building the python package:
install(TARGETS pycorgi DESTINATION .)


