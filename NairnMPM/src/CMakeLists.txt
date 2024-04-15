
include_directories("src")
include_directories("..") #NARINMPM
include_directories("../../Common")
include_directories("../..")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


FILE(GLOB_RECURSE BC_SRCS "Boundary_Conditions/*.h" "Boundary_Conditions/*.cpp")
FILE(GLOB_RECURSE Elements_SRCS "Elements/*.h" "Elements/*.cpp")


ADD_LIBRARY(BondaryConditions SHARED ${BC_SRCS})
TARGET_LINK_LIBRARIES(BondaryConditions)

ADD_LIBRARY(Elements SHARED ${Elements_SRCS})