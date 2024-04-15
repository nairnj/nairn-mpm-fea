
include_directories("src")
include_directories("..") #NARINMPM
include_directories("../../Common")
include_directories("../..")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


FILE(GLOB_RECURSE BC_SRCS "Boundary_Conditions/*.h" "Boundary_Conditions/*.cpp")

ADD_LIBRARY(BondaryConditions SHARED ${BC_SRCS})
TARGET_LINK_LIBRARIES(BondaryConditions)
