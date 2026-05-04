# Copyright 2016 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: GPL-3.0-or-later

set (TOOLS_FILES
     pytools.c++
     ../tools/config_parser.c++
     ../tools/system.c++
     )

set (EMF_FILES
     pyemf.c++
     ../core/emf/tile.c++
     ../core/emf/virtual_tile.c++
     ../core/emf/yee_lattice.c++
     ../core/emf/yee_lattice_FDTD2.c++
     ../core/emf/yee_lattice_stencil.c++
     ../core/emf/yee_lattice_current_filter_binomial2.c++
     )

set (PIC_FILES
     pypic.c++
     ../core/pic/tile.c++
     ../core/pic/tile_communication.c++
     ../core/pic/virtual_tile.c++
     ../core/pic/particle.c++
     ../core/pic/particle_current_zigzag_1st.c++
     ../core/pic/reflector_wall.c++
     )

set (IO_FILES
     ../io/snapshots/mpiio_fields.c++
     ../io/snapshots/mpiio_particles.c++
     ../io/snapshots/mpiio_spectra.c++
    )


pybind11_add_module(runko_cpp_bindings
    runko_cpp_bindings.c++
    ${TOOLS_FILES}
    ${EMF_FILES}
    ${PIC_FILES}
    ${IO_FILES}
)

set_target_properties(tyvi PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

target_link_libraries(runko_cpp_bindings
    PRIVATE
    runko_includes
    runko_warnings
    corgi
    tyvi
)

# Install runko_cpp_bindings.
install(TARGETS runko_cpp_bindings DESTINATION .)
# Vendor-in corgi if we are building the python package:
install(TARGETS pycorgi DESTINATION .)


