# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

function(add_cpp_mpi_test name no_mpi_proc)
    add_executable("${name}.test" "test_${name}.c++")
    if(${tyvi_BACKEND} STREQUAL "hip")
        set_source_files_properties(
          "test_${name}.c++"
          PROPERTIES
          LANGUAGE HIP
          HIP_STANDARD 26
        )
    endif()

    target_include_directories(
        "${name}.test"
        SYSTEM PRIVATE
        # Share the same boost.ut with cpp tests.
        "${CMAKE_CURRENT_SOURCE_DIR}/../cpp"
    )

    target_link_libraries(
        "${name}.test"
        PRIVATE
        runko_includes
        runko_warnings
        runko_dependencies
    )

    add_test(NAME "cpp-multirank-${name}" COMMAND
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${no_mpi_proc} "${name}.test"
    )
endfunction(add_cpp_mpi_test)

add_cpp_mpi_test(testing_setup 2)

