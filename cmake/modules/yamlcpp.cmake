include(FetchContent)

FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG master
)
FetchContent_MakeAvailable(yaml-cpp)

function(molsim_enable_yamlcpp TGT)
    message(STATUS "yaml-cpp enabled for target ${TGT}")
    target_link_libraries(${TGT} PUBLIC yaml-cpp::yaml-cpp)
endfunction()