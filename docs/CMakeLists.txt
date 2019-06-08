find_package(Doxygen)
if (NOT DOXYGEN)
  message(STATUS "Target 'doxygen' disabled (requires doxygen)")
  return ()
endif ()

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR})


set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
 
#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
 
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
#                   DEPENDS pycorgi pyplasmabox
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs")
 
add_custom_target(doxygen DEPENDS ${DOXYGEN_INDEX_FILE})


find_package(Sphinx)
#if (NOT Sphinx)
#  message(STATUS "Target 'docs' disabled (requires sphinx)")
#  return ()
#endif ()
 
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)
 
add_custom_target(docs 
                  COMMAND ${SPHINX_EXECUTABLE} -b html
                  -Dbreathe_projects.PlasmaBox=${DOXYGEN_OUTPUT_DIR}/xml
                  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	              DEPENDS 
		          ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
		          ${DOXYGEN_INDEX_FILE}
                  COMMENT "Generating documentation with Sphinx")

