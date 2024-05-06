

include_directories("Read_XML")


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB_RECURSE EXCEPTIONS_C "Exceptions/*.h" "Exceptions/*.cpp")
FILE(GLOB_RECURSE ELEMENTS_C "Elements/*.h" "Elements/*.cpp")
FILE(GLOB_RECURSE MATERIALS "Materials/*.h" "Materials/*.cpp")
FILE(GLOB_RECURSE NODES_C "Nodes/*.h" "Nodes/*.cpp")

FILE(GLOB_RECURSE BC_SRCS "Read_XML/*.h" "Read_XML/*.cpp")
FILE(GLOB_RECURSE SYSTEM_C "System/*.h" "System/CommonAnalysis.cpp"
                                      "System/CommonArchiveData.cpp"
                                      "System/CommonTask.cpp"
                                      "System/CommonUtilities.cpp"
                                      "System/Gauss.cpp"
                                      "System/LinkedObject.cpp"
                                      "System/Matrix3.cpp"
                                      "System/Matrix4.cpp"
                                      "System/Spline2D.cpp"
                                      "System/svninfo.cpp"
                                      "System/UnitsController.cpp"
                                      )
ADD_LIBRARY(Elements_C SHARED ${ELEMENTS_C})
ADD_LIBRARY(Exceptions_C SHARED ${EXCEPTIONS_C})
ADD_LIBRARY(Materials_C SHARED ${MATERIALS})
ADD_LIBRARY(Nodes_C SHARED ${NODES_C})
ADD_LIBRARY(Read_XML SHARED ${BC_SRCS})
ADD_LIBRARY(System_C SHARED ${SYSTEM_C})

TARGET_LINK_LIBRARIES(Read_XML Materials)
