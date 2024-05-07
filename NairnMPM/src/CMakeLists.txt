
include_directories("src")
include_directories("..") #NARINMPM
include_directories("../../Common")
include_directories("../..")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


FILE(GLOB_RECURSE BC_SRCS             "Boundary_Conditions/*.h" "Boundary_Conditions/*.cpp")
FILE(GLOB_RECURSE CRACKS_SRCS         "Cracks/*.h" "Cracks/*.cpp")
FILE(GLOB_RECURSE CUSTOM_TASKS        "Custom_Tasks/*.h" "Custom_Tasks/*.cpp")
FILE(GLOB_RECURSE Elements_SRCS       "Elements/*.h" "Elements/*.cpp")
FILE(GLOB_RECURSE EXCEPT_SRCS         "Exceptions/*.h" "Exceptions/*.cpp")
FILE(GLOB_RECURSE GLOBAL_QTS          "Global_Quantities/*.h" "Global_Quantities/*.cpp")

FILE(GLOB_RECURSE NODES_SRCS          "Nodes/*.h" "Nodes/*.cpp")
FILE(GLOB_RECURSE MATERIALS_SRCS      "Materials/*.h"  "Materials/*.cpp" )

FILE(GLOB_RECURSE MPMCLASSES_SRCS     "MPM_Classes/*.h" "MPM_Classes/*.cpp")
FILE(GLOB_RECURSE NAIRNMPMCLASS_SRCS  "NairnMPM_Class/*.h" "NairnMPM_Class/*.cpp")
FILE(GLOB_RECURSE PATCHES_SRCS        "Patches/*.h" "Patches/*.cpp")
FILE(GLOB_RECURSE READ_MPM_SRCS       "Read_MPM/*.h" "Read_MPM/*.cpp")
FILE(GLOB_RECURSE SYSTEM_SRCS         "System/*.h" "System/*.cpp")

FILE(GLOB_RECURSE RXML_SRCS           "../../Common/Read_XML/*.h" "../../Common/Read_XML/*.cpp")

# ADD_LIBRARY(BoundaryConditions SHARED ${BC_SRCS})
# ADD_LIBRARY(Custom_Tasks      STATIC  ${CUSTOM_TASKS})
# ADD_LIBRARY(Elements          STATIC  ${Elements_SRCS})
# ADD_LIBRARY(Exceptions        STATIC  ${EXCEPT_SRCS})
# ADD_LIBRARY(Global_QTs        STATIC  ${GLOBAL_QTS})
# ADD_LIBRARY(Materials         STATIC  ${MATERIALS_SRCS})
# ADD_LIBRARY(MPMClasses        STATIC  ${MPMCLASSES_SRCS})
# ADD_LIBRARY(NairnMPMClass     STATIC  ${NAIRNMPMCLASS_SRCS})
# ADD_LIBRARY(Nodes             STATIC  ${NODES_SRCS})
# ADD_LIBRARY(Patches           STATIC  ${PATCHES_SRCS})
# ADD_LIBRARY(Read_MPM          STATIC  ${READ_MPM_SRCS})
# ADD_LIBRARY(System            STATIC  ${SYSTEM_SRCS})

# TARGET_LINK_LIBRARIES(BoundaryConditions Custom_Tasks
                                         # Elements
                                         # Exceptions
                                         # Nodes
                                         # Patches
                                         # MPMClasses
                                         # Elements
                                         # Read_XML
                                         # System)
# TARGET_LINK_LIBRARIES(Custom_Tasks    
                                      
                                      # Nodes    
                                      # Read_XML)

# TARGET_LINK_LIBRARIES(Global_QTs      BoundaryConditions)

# TARGET_LINK_LIBRARIES(Read_MPM        BoundaryConditions
                                      
                                      # Read_XML
                                      # )
                                     
# TARGET_LINK_LIBRARIES(Elements        Elements_C)
# TARGET_LINK_LIBRARIES(Materials       BoundaryConditions
                                      # Materials_C)
# TARGET_LINK_LIBRARIES(Nodes           Nodes_C
                                      # )

# TARGET_LINK_LIBRARIES(System   System_C)
# TARGET_LINK_LIBRARIES(NairnMPMClass   
                                      
                                      
                                      # Global_QTs
                                      
                                      # # Materials_C
                                      # MPMClasses
                                      
                                      # Patches
                                      # Read_MPM
                                      # BoundaryConditions
                                      # System)

ADD_LIBRARY(NairnMPMClass     STATIC  ${NAIRNMPMCLASS_SRCS}
                                      ${CRACKS_SRCS}
                                      ${BC_SRCS}
                                      ${CUSTOM_TASKS}
                                      ${Elements_SRCS}
                                      ${EXCEPT_SRCS}
                                      ${GLOBAL_QTS}
                                      ${MATERIALS_SRCS}
                                      ${MPMCLASSES_SRCS}
                                      ${NODES_SRCS}
                                      ${PATCHES_SRCS}
                                      ${READ_MPM_SRCS}
                                      ${SYSTEM_SRCS}
                                      ${RXML_SRCS}
                                      )

 TARGET_LINK_LIBRARIES(NairnMPMClass  Materials_C
                                      System_C
                                      Elements_C
                                      Exceptions_C
                                      Nodes_C
                                      System_C)