
include_directories("src")
include_directories("..") #NARINMPM
include_directories("../../Common")
include_directories("../..")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


FILE(GLOB_RECURSE BC_SRCS "Boundary_Conditions/*.h" "Boundary_Conditions/*.cpp")
FILE(GLOB_RECURSE Elements_SRCS "Elements/*.h" "Elements/*.cpp")
FILE(GLOB_RECURSE EXCEPT_SRCS "Exceptions/*.h" "Exceptions/*.cpp")

FILE(GLOB_RECURSE NODES_SRCS "Nodes/*.h" "Nodes/*.cpp")
#FILE(GLOB_RECURSE MATERIALS_SRCS "Materials/*.h" "Materials/*.cpp")

FILE(GLOB_RECURSE NAIRNMPMCLASS_SRCS "NairnMPM_Class/*.h" "NairnMPM_Class/*.cpp")

ADD_LIBRARY(BoundaryConditions SHARED ${BC_SRCS})
TARGET_LINK_LIBRARIES(BoundaryConditions)


ADD_LIBRARY(Elements    SHARED ${Elements_SRCS})
ADD_LIBRARY(Exceptions  SHARED ${EXCEPT_SRCS})
#ADD_LIBRARY(Materials   SHARED ${MATERIALS_SRCS})
ADD_LIBRARY(NairnMPMClass   SHARED ${NAIRNMPMCLASS_SRCS})
ADD_LIBRARY(Nodes       SHARED ${NODES_SRCS})