# Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(Sphinx QUIET)

if (Sphinx_FOUND)

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)

    message(STATUS "Enabling target 'docs' using sphinx: ${SPHINX_EXECUTABLE}")

    add_custom_target(docs
                      COMMAND ${SPHINX_EXECUTABLE} -b html
                      ${SPHINX_SOURCE} ${SPHINX_BUILD}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")
else()
    message(STATUS "Target 'docs' disabled (requires sphinx)")
endif ()

