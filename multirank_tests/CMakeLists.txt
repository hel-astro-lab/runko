function(add_python_mpi_test name no_mpi_proc)
    add_test(NAME ${name} COMMAND
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${no_mpi_proc}
        ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test_${name}.py"
    )
endfunction(add_python_mpi_test)

add_python_mpi_test(testing_setup 2)
add_python_mpi_test(emf_simulation 4)
add_python_mpi_test(runko_simulation 2)
add_python_mpi_test(pic_simulation 4)
