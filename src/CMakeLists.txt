project(MPM3D_SRC)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

#################### set source files ####################
#------------------- main -----------------------------------------------#
aux_source_directory(main                       SRCS_MPM3D)

source_group(Sources\ Files\\MAIN               FILES ${SRCS_MPM3D})
#------------------- Material -------------------------------------------#
aux_source_directory(material                   SRCS_MATERIAL)
aux_source_directory(material/strength          SRCS_STRENGTH)
aux_source_directory(material/eos               SRCS_EOS)
aux_source_directory(material/failure           SRCS_FAILURE)

source_group(Sources\ Files\\MATERIAL           FILES ${SRCS_MATERIAL})
source_group(Sources\ Files\\MATERIAL\\STRENGTH FILES ${SRCS_STRENGTH})
source_group(Sources\ Files\\MATERIAL\\EOS      FILES ${SRCS_EOS})
source_group(Sources\ Files\\MATERIAL\\FAILURE  FILES ${SRCS_FAILURE})
#------------------- body -----------------------------------------------#

#------------------- grid -----------------------------------------------#

#------------------- solver ---------------------------------------------#

#------------------- utility --------------------------------------------#

set(SRC_LIST
    ${SRCS_MPM3D}
    ${SRCS_MATERIAL}
    ${SRCS_STRENGTH}
    ${SRCS_EOS}
    ${SRCS_FAILURE})
#################### set include files ####################
#------------------- main ----------------------------------------------#
file(GLOB INCS_MPM3D                            main/*.h*)

source_group(Header\ Files\\MAIN                FILES ${INCS_MPM3D})
#------------------- Material ------------------------------------------#
file(GLOB INCS_MATERIAL                         material/*.h*)
file(GLOB INCS_STRENGTH                         material/strength/*.h*)
file(GLOB INCS_EOS                              material/eos/*.h*)
file(GLOB INCS_FAILURE                          material/failure/*.h*)

source_group(Header\ Files\\MATERIAL            FILES ${INCS_MATERIAL})
source_group(Header\ Files\\MATERIAL\\STRENGTH  FILES ${INCS_STRENGTH})
source_group(Header\ Files\\MATERIAL\\EOS       FILES ${INCS_EOS})
source_group(Header\ Files\\MATERIAL\\FAILURE   FILES ${INCS_FAILURE})
#------------------- body ----------------------------------------------#

#------------------- grid ----------------------------------------------#

#------------------- solver --------------------------------------------#

#------------------- utility -------------------------------------------#

set(INC_LIST
    ${INCS_MPM3D}
    ${INCS_MATERIAL}
    ${INCS_STRENGTH}
    ${INCS_EOS}
    ${INCS_FAILURE})
#################### compile procedure ####################
set(MPM3D_BIN "MPM3D")

add_executable(${MPM3D_BIN} ${SRC_LIST} ${INC_LIST})

if(MPM3D_USE_VTKDATA)
    target_link_libraries(${MPM3D_BIN} ${VTK_LIBRARIES})
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