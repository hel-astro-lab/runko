# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

set(TESTS_TO_RUN
    unit_testing
    parse_atom_args
    qed_random
    qed2
    qed1
    pic_particle
)

foreach(test ${TESTS_TO_RUN})
    add_executable("${test}.test" "test_${test}.c++")

    target_include_directories(
        "${test}.test"
        SYSTEM PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )

    target_link_libraries(
        "${test}.test"
        PRIVATE
        runko_includes
        runko_warnings
        runko_dependencies
    )

    add_test(NAME "cpp-${test}" COMMAND "${test}.test")
endforeach()

# The container is compiled into the python module only, so its test builds it itself
# (its headers include pybind11).
target_sources(pic_particle.test
               PRIVATE "${PROJECT_SOURCE_DIR}/src/runko/pic/particle.c++")
target_link_libraries(pic_particle.test PRIVATE pybind11::pybind11)

# Provide a simple smoke test to make sure that the CLI works and can display a --help message.
add_test(NAME cpp-unit_test_cli COMMAND unit_testing.test --help)
set_property(TEST cpp-unit_test_cli PROPERTY PASS_REGULAR_EXPRESSION "help")

# Setup HIP for tests.
if(${tyvi_BACKEND} STREQUAL "hip")
    foreach(test ${TESTS_TO_RUN})
         set_source_files_properties("test_${test}.c++" PROPERTIES LANGUAGE HIP HIP_STANDARD 23)
    endforeach()
    set_source_files_properties("${PROJECT_SOURCE_DIR}/src/runko/pic/particle.c++"
                                PROPERTIES LANGUAGE HIP HIP_STANDARD 23)
endif()
