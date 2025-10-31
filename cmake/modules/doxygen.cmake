# make doc_doxygen optional if someone does not have / like doxygen

if(NOT DEFINED BUILD_DOC)
    option(BUILD_DOC "Build documentation using Doxygen" ON)
endif()

find_package(Doxygen REQUIRED)

if(BUILD_DOC)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs)

    add_custom_command(
            OUTPUT ${DOXYGEN_OUT}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
    )

    add_custom_target(doc
            DEPENDS ${DOXYGEN_OUT}/html/index.html
            COMMENT "Documentation generated in ${DOXYGEN_OUT}"
    )

    set_target_properties(doc PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD YES)

    add_custom_target(doxygen DEPENDS doc)

    message(STATUS "Doxygen target 'doc' (or 'doxygen') added. Run 'make doc' to generate docs in ${DOXYGEN_OUT}.")
else()
    message(STATUS "Doxygen target disabled (BUILD_DOC=OFF).")
endif()