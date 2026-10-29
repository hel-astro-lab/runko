# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

function(add_python_mpi_test name no_mpi_proc)
    set(test_name py-multirank-${name})
    if(ARGC EQUAL 3)
        list(GET ARGV 2 extra_postfix)
        set(test_name ${test_name}${extra_postfix})
    endif()
    add_test(NAME ${test_name} COMMAND
        ${MPIEXEC} --map-by :OVERSUBSCRIBE
        ${MPIEXEC_NUMPROC_FLAG} ${no_mpi_proc}
        ${Python_EXECUTABLE} ${py_test_wrapper}
        "${CMAKE_CURRENT_SOURCE_DIR}/test_${name}.py"
    )
    set_property(TEST ${test_name}
        PROPERTY ENVIRONMENT
        "PYTHONPATH=${CMAKE_BINARY_DIR}/src:${PROJECT_SOURCE_DIR}"
        "PIKA_COMMANDLINE_OPTIONS=--pika:threads=2"
    )
endfunction(add_python_mpi_test)

add_python_mpi_test(testing_setup 2)
add_python_mpi_test(emf_simulation 3 "_uneven")
add_python_mpi_test(emf_simulation 4 "_even")
add_python_mpi_test(runko_simulation 2)
add_python_mpi_test(pic_simulation 4)
add_python_mpi_test(mpiio_fields 4)
add_python_mpi_test(mpiio_particles 4)
add_python_mpi_test(mpiio_spectra 4)
