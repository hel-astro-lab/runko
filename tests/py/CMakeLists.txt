# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

function(add_python_test name)
    add_test(NAME py-${name} COMMAND
        ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test_${name}.py"
    )
    set_property(TEST py-${name}
        PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src:${PROJECT_SOURCE_DIR}
    )
endfunction(add_python_test)

add_python_test(auto_outdir)
add_python_test(configuration)
add_python_test(emf_antenna)
add_python_test(emf_antenna_time_evolution)
add_python_test(emf_current_filter_binomial2)
add_python_test(emf_edge_bc)
add_python_test(emf_FDTD2)
add_python_test(emf)
add_python_test(emf_stencil)
add_python_test(hilbert)
add_python_test(module_meta)
add_python_test(mpiio_fields)
add_python_test(mpiio_particles)
add_python_test(mpiio_spectra)
add_python_test(pic_current_depositer_zigzag_1st_atomic)
add_python_test(pic_current_depositer_zigzag_1st)
add_python_test(pic_moving_injector)
add_python_test(pic_particle_pusher)
add_python_test(pic_particle_sorting)
add_python_test(pic)
add_python_test(pic_reflector_wall)
add_python_test(ram_usage)
add_python_test(simulation)
add_python_test(tile_grid)
add_python_test(virtual_tile_specializations)
