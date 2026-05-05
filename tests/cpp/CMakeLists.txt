# Copyright 2026 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: GPL-3.0-or-later

set(TESTS_TO_RUN
    unit_testing
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

# Provide a simple smoke test to make sure that the CLI works and can display a --help message.
add_test(NAME cpp-unit_test_cli COMMAND unit_testing.test --help)
set_property(TEST cpp-unit_test_cli PROPERTY PASS_REGULAR_EXPRESSION "help")
