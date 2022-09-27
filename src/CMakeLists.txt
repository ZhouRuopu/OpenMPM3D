project(MPM3D_SRC)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

#################### set source files ####################

set(SRC_LIST)
#################### set include files ####################

set(INC_LIST)
#################### compile procedure ####################
set(MPM3D_BIN "MPM3D")

add_executable(${MPM3D_BIN} ${SRC_LIST} ${INC_LIST})

if(MPM3D_USE_VTKDATA)
    target_link_libraries(${MPM3D_BIN} vtkCommon vtkFiltering vtkIO)
endif()

if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake|VCExpress|MSBuild)")
    set_target_properties(${MPM3D_BIN} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT.lib")
endif()

#################### set include directories ####################
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

#################### INSTALL ####################
set(INSTALL_LIB_DIR lib)
set(INSTALL_BIN_DIR bin)

install(TARGETS ${MPM3D_BIN}
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR})
include(CPack)