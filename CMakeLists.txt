cmake_minimum_required(VERSION 3.14)
project(MolSim_GroupB VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# --- Project Modules ---
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(vtk)            # Defines molsim_enable_vtk()
include(clang-format)
option(BUILD_DOC "Build documentation using Doxygen" ON)
include(doxygen)

# --- Core Library ---
# Contains all shared application logic
file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/Simulation.cpp") # Exclude main()

add_library(molsim_core ${SRC_FILES})
target_include_directories(molsim_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_options(molsim_core PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

# Apply VTK to the core library if enabled
molsim_enable_vtk(molsim_core)

# --- Main Application Executable ---
add_executable(MolSim src/Simulation.cpp) # The file containing main()
target_link_libraries(MolSim PRIVATE molsim_core)

# --- Testing ---
enable_testing()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571
)
FetchContent_MakeAvailable(googletest)

# Test executable
add_executable(MolSimTests
        test/ParticleContainerTest.cpp
)

target_link_libraries(MolSimTests PUBLIC
        molsim_core      # This link now provides the include directory
        GTest::gtest_main
)

# Discover tests for CTest
include(GoogleTest)
gtest_discover_tests(MolSimTests)
