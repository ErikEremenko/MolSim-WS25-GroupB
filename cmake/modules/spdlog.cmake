include(FetchContent)

FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(spdlog)

function(molsim_enable_spdlog TGT)
    message(STATUS "spdlog enabled for target ${TGT}")
    target_link_libraries(${TGT} PUBLIC spdlog::spdlog)
endfunction()