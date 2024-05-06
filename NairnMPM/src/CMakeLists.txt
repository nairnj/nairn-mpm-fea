
include_directories("src")
include_directories("..") #NARINMPM
include_directories("../../Common")
include_directories("../..")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


FILE(GLOB_RECURSE BC_SRCS "Boundary_Conditions/*.h" "Boundary_Conditions/*.cpp")
FILE(GLOB_RECURSE BC_SRCS "Cracks/*.h" "Cracks/*.cpp")
FILE(GLOB_RECURSE CUSTOM_TASKS "Custom_Tasks/*.h" "Custom_Tasks/*.cpp")
FILE(GLOB_RECURSE Elements_SRCS "Elements/*.h" "Elements/*.cpp")
FILE(GLOB_RECURSE EXCEPT_SRCS "Exceptions/*.h" "Exceptions/*.cpp")
FILE(GLOB_RECURSE GLOBAL_QTS "Global_Quantities/*.h" "Global_Quantities/*.cpp")

FILE(GLOB_RECURSE NODES_SRCS "Nodes/*.h" "Nodes/*.cpp")
FILE(GLOB_RECURSE MATERIALS_SRCS "Materials/*.h"  "Materials/*.cpp")

FILE(GLOB_RECURSE MPMCLASSES_SRCS "MPM_Classes/*.h" "MPM_Classes/*.cpp")
FILE(GLOB_RECURSE NAIRNMPMCLASS_SRCS "NairnMPM_Class/*.h" "NairnMPM_Class/*.cpp")
FILE(GLOB_RECURSE PATCHES_SRCS "Patches/*.h" "Patches/*.cpp")
FILE(GLOB_RECURSE READ_MPM_SRCS "Read_MPM/*.h" "Read_MPM/*.cpp")

ADD_LIBRARY(BoundaryConditions SHARED ${BC_SRCS})
ADD_LIBRARY(Custom_Tasks      SHARED ${CUSTOM_TASKS})
ADD_LIBRARY(Elements          SHARED ${Elements_SRCS})
ADD_LIBRARY(Exceptions        SHARED ${EXCEPT_SRCS})
ADD_LIBRARY(GLOBAL_QTS        SHARED ${GLOBAL_QTS})
ADD_LIBRARY(Materials         SHARED ${MATERIALS_SRCS})
ADD_LIBRARY(MPMClasses        SHARED ${MPMCLASSES_SRCS})
ADD_LIBRARY(NairnMPMClass     SHARED ${NAIRNMPMCLASS_SRCS})
ADD_LIBRARY(Nodes             SHARED ${NODES_SRCS})
ADD_LIBRARY(Patches           SHARED ${PATCHES_SRCS})
ADD_LIBRARY(Read_MPM          SHARED ${READ_MPM_SRCS})

TARGET_LINK_LIBRARIES(BoundaryConditions)
TARGET_LINK_LIBRARIES(Custom_Tasks    BoundaryConditions
                                      Nodes    
                                      Read_XML)

TARGET_LINK_LIBRARIES(Read_MPM        Read_XML)
TARGET_LINK_LIBRARIES(Elements        Elements_C)
TARGET_LINK_LIBRARIES(Materials       Materials_C)
TARGET_LINK_LIBRARIES(NairnMPMClass   BoundaryConditions
                                      Custom_Tasks
                                      Elements
                                      GLOBAL_QTS
                                      Materials
                                      # Materials_C
                                      MPMClasses
                                      Nodes
                                      Patches)
