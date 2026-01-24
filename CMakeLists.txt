cmake_minimum_required(VERSION 3.14)
project(MolSim_GroupB VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# Enable Link-Time Optimization (LTO) for Release builds (allows cross-translation-unit inlining and optimization)
include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)
option(ENABLE_LTO "Enable Link-Time Optimization for Release builds" ON)

if(ENABLE_LTO AND LTO_SUPPORTED)
    message(STATUS "LTO: Link-Time Optimization enabled for Release builds")
else()
    if(ENABLE_LTO)
        message(WARNING "LTO: Not supported by compiler - ${LTO_ERROR}")
    endif()
endif()

# Project Modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(vtk)
include(clang-format)
include(doxygen)
include(spdlog)
include(gtest)
include(openmp)
include(yamlcpp)

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

# Enable LTO for core library in Release mode
if(ENABLE_LTO AND LTO_SUPPORTED)
    set_property(TARGET molsim_core PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

# -march=native enables the best SIMD instructions for the current CPU (AVX2 locally, AVX-512 on cluster)
option(ENABLE_NATIVE_ARCH "Enable -march=native for SIMD optimization" ON)

target_compile_options(molsim_core PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)
# Add SIMD flags for Release builds
# -ffast-math for aggressive floating-point optimization
# -funroll-loops for loop unrolling
if(ENABLE_NATIVE_ARCH)
    target_compile_options(molsim_core PRIVATE
        $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Release>>:-march=native -mtune=native -ffast-math -funroll-loops>
        $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CONFIG:Release>>:-march=native -mtune=native -ffast-math -funroll-loops>
        $<$<AND:$<CXX_COMPILER_ID:Intel>,$<CONFIG:Release>>:-xHost -fp-model fast=2 -unroll>
    )
    message(STATUS "SIMD: Using native architecture optimization (-march=native). Fast math and loop unrolling enabled.")
endif()

# Apply dependencies to the core library
molsim_enable_vtk(molsim_core)
molsim_enable_spdlog(molsim_core)
molsim_enable_yamlcpp(molsim_core)
molsim_enable_OpenMP(molsim_core)

# Main Application Executable
add_executable(MolSim src/main.cpp
)
target_link_libraries(MolSim PRIVATE molsim_core)

# Enable LTO for the main executable in Release mode
if(ENABLE_LTO AND LTO_SUPPORTED)
    set_property(TARGET MolSim PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

# Define custom macro for project source directory
target_compile_definitions(MolSim PRIVATE PROJ_SRC_DIR="${CMAKE_SOURCE_DIR}")

molsim_enable_testing()
