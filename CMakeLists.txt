cmake_minimum_required(VERSION 3.14)
project(MolSim_GroupB VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# Project Modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(vtk)
include(clang-format)
include(doxygen)
include(spdlog)
include(gtest)
include(openmp)

# Core Library
# Contains all shared application logic
file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp") # Exclude main()

add_library(molsim_core ${SRC_FILES})
target_include_directories(molsim_core
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_options(molsim_core PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

# Apply dependencies to the core library
molsim_enable_vtk(molsim_core)
molsim_enable_spdlog(molsim_core)
molsim_enable_OpenMP(molsim_core)

# Main Application Executable
add_executable(MolSim src/main.cpp)
target_link_libraries(MolSim PRIVATE molsim_core)

molsim_enable_testing()
