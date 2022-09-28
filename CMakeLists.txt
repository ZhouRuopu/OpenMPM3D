cmake_minimum_required(VERSION 3.12)
project(MPM3D CXX)

# Disallow in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "${PROJECT_NAME} requires an out-of-source Build. Please create a seperate binary directory and run CMake there.")
endif()

# Allow more readable "if...then...else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

#################### project version ####################
set(${PROJECT_NAME}_MAJOR_VERSION 1)
set(${PROJECT_NAME}_MINOR_VERSION 0)
set(${PROJECT_NAME}_PATCH_LEVEL 0)
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION})
set(${PROJECT_NAME}_VERSION_FULL ${${PROJECT_NAME}_VERSION}.${${PROJECT_NAME}_PATCH_LEVEL})

#################### set default cmake build type to Release ####################
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE PATH "MPM3D build type" FORCE)
endif()

#################### build MPM3D with selected solvers ####################
option(MPM_SOLVER "Build MPM solver in MPM3D." ON)
option(SGMPM_SOLVER "Build SGMPM solver in MPM3D." OFF)
option(HFEMPM_SOLVER "Build Hybrid FE-MPM solver in MPM3D. " OFF)
option(FVM_SOLVER "Build FVM solver in MPM3D." OFF)

if(MPM_SOLVER)
    add_definitions(-DMPM_SOLVER)
endif()

if(SGMPM_SOLVER)
    add_definitions(-DSGMPM_SOLVER)
endif()

if(HFEMPM_SOLVER)
    add_definitions(-DHFEMPM_SOLVER)
endif()

if(FVM_SOLVER)
    add_definitions(-DFVM_SOLVER)
endif()

#################### build MPM3D with double precision ####################
if(NOT MPM3D_USE_DOUBLE)
    option(MPM3D_USE_DOUBLE "Build double precision version." ON)
endif()

if(MPM3D_USE_DOUBLE)
    add_definitions(-D_MPM_DOUBLE)
endif()

#################### build MPM3D with massive particles ####################
if(NOT MPM3D_MASSIVE_PARTICLE)
    option(MPM3D_MASSIVE_PARTICLE "Build massive particles version" OFF)
endif()

if(MPM3D_MASSIVE_PARTICLE)
    add_definitions(-D_MPM_MASSIVE_PARTICLE)
endif()

#################### CMake warning level option ####################
option(MPM3D_NO_WARNING "Set warning level." ON)
if(MPM3D_NO_WARNING)
    if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake|VCExpress|MSBuild)")
        set(CMAKE_CXX_WARNING_LEVEL 0)
        if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
            string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        endif()
    endif()
endif()

#################### VTK support ####################
option(MPM3D_USE_VTKDATA "Build VTK unstructured grid data format support. This requires a precompiled VTK." ON)

if(MPM3D_USE_VTKDATA)
    find_package(VTK)
    if(VTK_USE_FILE)
        include(${VTK_USE_FILE})
        include_directories(${VTK_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "Please set VTK_DIR.")
    endif()
else()
    add_definitions(-D_MPM_NOVTKDATA)
endif()

add_subdirectory(src)