enable_testing()
include(FetchContent)
include(GoogleTest)

function(molsim_enable_testing)
    message(STATUS "Google Test enabled")
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG 52eb8108c5bdec04579160ae17225d66034bd723
    )
    FetchContent_MakeAvailable(googletest)

    file(GLOB_RECURSE TEST_FILES "${CMAKE_SOURCE_DIR}/test/*.cpp")

    add_executable(MolSimTests
            ${TEST_FILES}
    )

    target_link_libraries(MolSimTests PUBLIC
            molsim_core
            GTest::gtest_main
    )

    gtest_discover_tests(MolSimTests)
    message(STATUS "Test executable 'MolSimTests' created.")
endfunction()
