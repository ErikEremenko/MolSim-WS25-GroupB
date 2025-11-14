option(BUILD_DOC "Build documentation using Doxygen" ON)

find_package(Doxygen QUIET COMPONENTS dot) # doxygen not required

if(BUILD_DOC AND Doxygen_FOUND)
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")
    set(DOXYGEN_FILE_PATTERNS "*.h;*.hpp;*.c;*.cc;*.cpp;*.md")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doxys_documentation")

    doxygen_add_docs(doc_doxygen
            "${PROJECT_SOURCE_DIR}/src"
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/README.md"
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            COMMENT "Generate API documentation with Doxygen"
    )

    set_target_properties(doc_doxygen PROPERTIES
            EXCLUDE_FROM_ALL TRUE
            EXCLUDE_FROM_DEFAULT_BUILD TRUE
    )

    message(STATUS "Doxygen target 'doc_doxygen' added. Run 'make doc_doxygen' to generate docs in ${DOXYGEN_OUTPUT_DIRECTORY}.")
else()
    message(STATUS "Doxygen target disabled (BUILD_DOC=OFF).")
endif()