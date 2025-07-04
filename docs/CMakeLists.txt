find_package(Sphinx)

if (Sphinx_FOUND)

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)

    add_custom_target(docs
                      COMMAND ${SPHINX_EXECUTABLE} -b html
                      ${SPHINX_SOURCE} ${SPHINX_BUILD}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")
else()
    message(STATUS "Target 'docs' disabled (requires sphinx)")
endif ()

