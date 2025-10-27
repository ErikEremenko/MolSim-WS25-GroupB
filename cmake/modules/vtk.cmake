option(ENABLE_VTK_OUTPUT "Enable VTK output" OFF)

if(ENABLE_VTK_OUTPUT)
    message(STATUS "VTK output enabled")
    # VTK Library, only add required components for less unnecessary includes
    find_package(VTK REQUIRED COMPONENTS
            CommonCore
            CommonDataModel
            IOXML
    )

    if(VTK_FOUND)
        message (STATUS "Found VTK Version: ${VTK_VERSION}")
    else ()
        message(FATAL_ERROR "VTK not found")
    endif ()

    if(VTK_VERSION VERSION_GREATER_EQUAL 8.9)
         include_directories(${VTK_INCLUDE_DIRS})
    else()
        include(${VTK_USE_FILE})
    endif ()

    target_link_libraries(MolSim
            PRIVATE
            ${VTK_LIBRARIES}
    )
    target_compile_definitions(MolSim PRIVATE ENABLE_VTK_OUTPUT)
endif()
