CMAKE_MINIMUM_REQUIRED (VERSION 3.0)

include_directories(${XERCESC_INCLUDE_DIR})
link_directories(${XERCESC_LIB_DIR})

include_directories("Common")
include_directories("Common/Headers")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_definitions(-DMPM_CODE)
add_definitions(-DCONST_ARRAYS)
add_definitions(-DPOROELASTICITY)
add_definitions(-DRESTART_OPTION)
add_definitions(-DADJUST_COPIED_PK)
add_definitions(-D_OPENMP)

#THESE ARE FOR XML READ ONLY
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/NairnMPM/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/NairnFEA/src)
ADD_SUBDIRECTORY(NairnMPM)
ADD_SUBDIRECTORY(Common)

ADD_EXECUTABLE (Nairn_MPM "Common/System/main.cpp")
TARGET_LINK_LIBRARIES (Nairn_MPM 
                                 NairnMPMClass  
 
                                 xerces-c
                                 )
