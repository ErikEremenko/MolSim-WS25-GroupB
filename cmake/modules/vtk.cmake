# cmake/modules/vtk.cmake
option(ENABLE_VTK_OUTPUT "Enable VTK output" OFF)

function(molsim_enable_vtk tgt)
    if(ENABLE_VTK_OUTPUT)
        message(STATUS "VTK output enabled")
        find_package(VTK REQUIRED COMPONENTS
                CommonCore
                CommonDataModel
                IOXML
        )
        # Link imported module targets
        target_link_libraries(${tgt} PUBLIC
                VTK::CommonCore
                VTK::CommonDataModel
                VTK::IOXML
        )
        # Propagate the feature macro to consumers (MolSim, tests)
        target_compile_definitions(${tgt} PUBLIC ENABLE_VTK_OUTPUT)
    endif()
endfunction()
