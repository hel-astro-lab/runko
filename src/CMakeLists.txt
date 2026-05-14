# Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

set (TOOLS_FILES
     ./runko/bindings/pytools.c++
     ./runko/tools/config_parser.c++
     ./runko/tools/system.c++
     )

set (EMF_FILES
     ./runko/bindings/pyemf.c++
     ./runko/emf/tile.c++
     ./runko/emf/virtual_tile.c++
     ./runko/emf/yee_lattice.c++
     ./runko/emf/yee_lattice_FDTD2.c++
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


pybind11_add_module(runko_cpp_bindings
    ./runko/bindings/runko_cpp_bindings.c++
    ${TOOLS_FILES}
    ${EMF_FILES}
    ${PIC_FILES}
    ${IO_FILES}
)

target_link_libraries(runko_cpp_bindings
    PRIVATE
    runko_includes
    runko_warnings
    runko_dependencies
)

# Install runko_cpp_bindings.
install(TARGETS runko_cpp_bindings DESTINATION .)
# Vendor-in corgi if we are building the python package:
install(TARGETS pycorgi DESTINATION .)


