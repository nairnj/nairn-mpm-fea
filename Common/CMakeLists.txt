

include_directories("Read_XML")


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB_RECURSE EXCEPTIONS_C "Exceptions/*.h" "Exceptions/*.cpp")
FILE(GLOB_RECURSE ELEMENTS_C "Elements/*.h" "Elements/*.cpp")
FILE(GLOB_RECURSE MATERIALS "Materials/*.h" "Materials/*.cpp")
FILE(GLOB_RECURSE NODES_C "Nodes/*.h" "Nodes/*.cpp")

FILE(GLOB_RECURSE READXML_SRCS "Read_XML/*.h" "Read_XML/*.cpp")
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
ADD_LIBRARY(Elements_C STATIC ${ELEMENTS_C})
ADD_LIBRARY(Exceptions_C STATIC ${EXCEPTIONS_C})
ADD_LIBRARY(Materials_C STATIC ${MATERIALS})
ADD_LIBRARY(Nodes_C STATIC ${NODES_C})
# ADD_LIBRARY(Read_XML STATIC ${READXML_SRCS})
ADD_LIBRARY(System_C STATIC ${SYSTEM_C})

# TARGET_LINK_LIBRARIES(Read_XML Materials)
