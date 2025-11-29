find_package(OpenMP REQUIRED COMPONENTS CXX)

function(molsim_enable_OpenMP TGT)
    message(STATUS "OpenMP enabled for target ${TGT}")
    target_link_libraries(${TGT} PUBLIC OpenMP::OpenMP_CXX)
endfunction()