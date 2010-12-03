# -------------------------------------------------------------------------
# To compile NairnFEA use the following commands
#
# cd to directory with this 'makefile'
# give the command 'make SYSTEM=system' where the lower case system is where
#    you are compiling. Supported systems are:
#       mac - MacOS X
#       carbon - OSU carbon/hedera cluster (the default)
#       cygwin - Standard cygwin install in Windows
#       cluster2 - Grant Smith's cluster at U of Utah
#       inferno - C-SAFE cluster at U of Utah
#       zia - A Debian linux machine
# Other options are
#    make clean - to remove all compiled objects
#    make install - to copy compiled code to desired installation destination
#
# If you need other options, you can edit the makefile or override them in the make command
# The following are the most important variables
#
# 1. $(xercesSo) is path to directory with the compiled Xerces library
# 2. $(headersPath) topath to directory with Xerces include files
# 3. $(CC) is the gcc compiler, $(RM) is the rm command, $(CP) is copy for install
# 4. $(CFLAGS) is the gcc compiler options
#    $(LFLAGS) is the gcc linking options
# 5. $(output) is path to save executable
# 6. $(ioutput) is path to install executable
# 7. $(src) and $(com) are paths to parts of the source code

# -------------------------------------------------------------------------
# 1. Xerces location for Mac, Carbon/Hedera, Windows/cygwin, cluster2, and inferno
SYSTEM=carbon
ifeq ($(SYSTEM),mac)
	xercesSo = /usr/local/lib/libxerces-c.dylib
endif
ifeq ($(SYSTEM),carbon)
	xercesSo = /usr/local/xerces3/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),cygwin)
	xercesSo = -L/usr/local/lib -lxerces-c
endif
ifeq ($(SYSTEM),cluster2)
	xercesSo = /app/software/Xerces/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),inferno)
	xercesSo = /usr/sci/projects/SCIRun/Thirdparty/1.22/Linux/gcc-3.2-32bit/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),zia)
	xercesSo = /usr/lib/libxerces-c.so
endif

# 2. Path to outside headers for (Mac & Carbon/Hedera & Windows/cygwin), cluster2, and inferno 
ifeq ($(SYSTEM),mac)
	headersPath = /usr/local/include
endif
ifeq ($(SYSTEM),carbon)
	headersPath = /usr/local/xerces3/include
endif
ifeq ($(SYSTEM),cygwin)
	headersPath = /usr/local/include
endif
ifeq ($(SYSTEM),cluster2)
	headersPath = /app/software/Xerces/includ
endif
ifeq ($(SYSTEM),inferno)
	headersPath = /usr/sci/projects/SCIRun/Thirdparty/1.22/Linux/gcc-3.2-32bit/include
endif
ifeq ($(SYSTEM),zia)
	headersPath = /usr/include/xercesc
endif

# 3. Define compiler to use
CC= g++
RM= rm
CP= cp

# 4. Compiler options -c=compile only, -O3=optimize most -O0=not optimized -g=debugging info
CFLAGS= -c -O3
LFLAGS=

# 5. Define executable destination
output = ../input/NairnFEA
dtdpath = ../input/NairnFEA.dtd

# 6. Define executable 'make install' location
ioutput = ~/bin

# 7. Define paths to all files (relative to 'makefile')
src = ../src
com = ../../Common
#com = ../../NairnMPM/Common

# Macros for headers and all files
headers = -I$(src) -I$(com) -I$(headersPath)
ArcController = $(com)/Read_XML/ArcController
Area = $(src)/Read_FEA/Area
BitMapFilesCommon = $(com)/Read_XML/BitMapFilesCommon
BitMapFilesFEA = $(src)/Read_FEA/BitMapFilesFEA
BMPLevel = $(com)/Read_XML/BMPLevel
BoxController = $(com)/Read_XML/BoxController
CommonAnalysis = $(com)/System/CommonAnalysis
CommonArchiveData = $(com)/System/CommonArchiveData
CommonException = $(com)/Exceptions/CommonException
CommonReadHandler = $(com)/Read_XML/CommonReadHandler
CommonUtilities = $(com)/System/CommonUtilities
Constraint = $(src)/Boundary_Conditions/Constraint
ConstraintController = $(src)/Read_FEA/ConstraintController
CSTriangle = $(src)/Elements/CSTriangle
EdgeBC = $(src)/Boundary_Conditions/EdgeBC
EdgeBCController = $(src)/Read_FEA/EdgeBCController
EightNodeIsoparam = $(src)/Elements/EightNodeIsoparam
Elastic = $(com)/Materials/Elastic
ElementBase = $(com)/Elements/ElementBase
ElementsController = $(com)/Read_XML/ElementsController
FEAArchiveData = $(src)/System/FEAArchiveData
FEABoundaryCondition = $(src)/Boundary_Conditions/FEABoundaryCondition
FEAReadHandler = $(src)/Read_FEA/FEAReadHandler
FEAStartResults = $(src)/NairnFEA_Class/FEAStartResults
FourNodeIsoparam = $(com)/Elements/FourNodeIsoparam
Gauss = $(src)/Numerical/Gauss
ImperfectInterface = $(src)/Materials/ImperfectInterface
Interface2D = $(src)/Elements/Interface2D
IsotropicMat = $(com)/Materials/IsotropicMat
Keypoint = $(src)/Read_FEA/Keypoint
KeypointsController = $(src)/Read_FEA/KeypointsController
Linear2D = $(com)/Elements/Linear2D
LinearInterface = $(src)/Elements/LinearInterface
LineController = $(com)/Read_XML/LineController
LinkedObject = $(com)/System/LinkedObject
main = $(com)/System/main
MaterialBase = $(com)/Materials/MaterialBase
MaterialBaseFEA = $(src)/Materials/MaterialBaseFEA
MaterialController = $(com)/Read_XML/MaterialController
mathexpr = $(com)/Read_XML/mathexpr
MoreElementBase = $(src)/Elements/MoreElementBase
MoreNodalPoint = $(src)/Nodes/MoreNodalPoint
NairnFEA = $(src)/NairnFEA_Class/NairnFEA
NodalDispBC = $(src)/Boundary_Conditions/NodalDispBC
NodalDispBCController = $(src)/Read_FEA/NodalDispBCController
NodalLoad = $(src)/Boundary_Conditions/NodalLoad
NodalLoadController = $(src)/Read_FEA/NodalLoadController
NodalPoint = $(com)/Nodes/NodalPoint
NodalPoint2D = $(com)/Nodes/NodalPoint2D
NodalPoint3D = $(com)/Nodes/NodalPoint3D
NodesController = $(com)/Read_XML/NodesController
Orthotropic = $(com)/Materials/Orthotropic
ParseController = $(com)/Read_XML/ParseController
Path = $(src)/Read_FEA/Path
PathBCController = $(src)/Read_FEA/PathBCController
PathsController = $(src)/Read_FEA/PathsController
PointController = $(com)/Read_XML/PointController
Quad2D = $(src)/Elements/Quad2D
QuadInterface = $(src)/Elements/QuadInterface
RectController = $(com)/Read_XML/RectController
ShapeController = $(com)/Read_XML/ShapeController
SixNodeTriangle = $(src)/Elements/SixNodeTriangle
StrX = $(com)/Exceptions/StrX
TransIsotropic = $(com)/Materials/TransIsotropic
Utilities = $(src)/System/Utilities

# path from 'makefile' to prefix header file
prefix = $(src)/System/FEAPrefix.hpp

# path from 'makefile' to default headers (this included by $(prefix)
dprefix = $(src)/System/FEAPrefix.hpp $(com)/System/DataTypes.hpp \
		$(com)/System/CommonAnalysis.hpp $(com)/System/LinkedObject.hpp

# header files that include other header files (other than $(dprefix) files)
# ArcController.hpp : LineController.hpp ShapeController.hpp
# BoxController.hpp : ShapeController.hpp
# ConstraintController.hpp : ParseController.hpp
# CSTriangle.hpp : Linear2D.hpp (ElementBase.hpp)
# EdgeBC.hpp : FEABoundaryCondition.hpp
# EdgeBCController.hpp : ParseController.hpp
# EightNodeIsoparam.hpp : Quad2D.hpp (ElementBase.hpp)
# Elastic.hpp : MaterialBase.hpp
# ElementsController.hpp : ParseController.hpp
# FEAArchiveData.hpp : CommonArchiveData.hpp
# FEAReadHandler.hpp : CommonReadHandler.hpp
# FourNodeIsoparam.hpp : Linear2D.hpp (ElementBase.hpp)
# KeypointsController.hpp : ParseController.hpp
# ImperfectInterface.hpp : ElementBase.hpp
# Interface2D.hpp : ElementBase.hpp
# IsotropicMat.hpp : Elastic.hpp (MaterialBase.hpp)
# Linear2D.hpp : ElementBase.hpp
# LinearInterface : Interface2D.hpp (ElementBase.hpp)
# LineController.hpp : ShapeController.hpp
# MaterialController.hpp : ParseController.hpp
# NodalDispBC.hpp : FEABoundaryCondition.hpp
# NodalDispBCController.hpp : ParseController.hpp
# NodalLoad.hpp : FEABoundaryCondition.hpp
# NodalLoadController.hpp : ParseController.hpp
# NodalPoint2D.hpp : NodalPoint.hpp
# NodalPoint3D.hpp : NodalPoint.hpp
# NodesController.hpp : ParseController.hpp
# Orthotropic.hpp : TransIsotropic.hpp (Elastic.hpp MaterialBase.hpp)
# PathBCController.hpp : ShapeController.hpp
# PathsController.hpp : ParseController.hpp
# PointController.hpp : ShapeController.hpp
# Quad2D.hpp : ElementBase.hpp
# QuadInterface : Interface2D.hpp (ElementBase.hpp)
# RectController.hpp : ShapeController.hpp
# SixNodeTriangle.hpp : Quad2D.hpp (ElementBase.hpp)
# TransIsotropic.hpp : Elastic.hpp (MaterialBase.hpp)

# Header files with no dependencies
# Area.hpp CommonArchiveData.hpp CommonException.hpp CommonReadHandler.hpp FEABoundaryCondition.hpp
# ElementBase.hpp Keypoint.hpp MaterialBase.hpp NairnFEA.hpp
# NodalPoint.hpp ParseController.hpp Path.hpp ShapeController.hpp StrX.hpp

# -------------------------------------------------------------------------
# all compiled objects
objects = main.o CommonAnalysis.o CommonArchiveData.o LinkedObject.o CommonUtilities.o \
		CommonException.o ElementBase.o Linear2D.o FourNodeIsoparam.o MaterialBase.o \
		CommonReadHandler.o ParseController.o ElementsController.o NodesController.o \
		LineController.o MaterialController.o IsotropicMat.o TransIsotropic.o BoxController.o \
		Orthotropic.o NodalPoint.o FEAArchiveData.o Utilities.o NairnFEA.o FEAStartResults.o \
		FEAReadHandler.o KeypointsController.o Keypoint.o PathsController.o Path.o Area.o \
		EdgeBCController.o NodalDispBCController.o NodalLoadController.o MaterialBaseFEA.o \
		ImperfectInterface.o MoreElementBase.o Quad2D.o EightNodeIsoparam.o SixNodeTriangle.o \
		CSTriangle.o Interface2D.o LinearInterface.o QuadInterface.o MoreNodalPoint.o \
		NodalDispBC.o NodalLoad.o EdgeBC.o Gauss.o Elastic.o Constraint.o ConstraintController.o \
		mathexpr.o FEABoundaryCondition.o ShapeController.o PointController.o PathBCController.o \
		RectController.o ArcController.o NodalPoint2D.o NodalPoint3D.o BitMapFilesCommon.o \
		BMPLevel.o BitMapFilesFEA.o

# -------------------------------------------------------------------------
# Link all objects
NairnFEA : $(objects)
	$(CC) $(LFLAGS) -o $(output) $(objects) $(xercesSo)

# -------------------------------------------------------------------------
# All source code files
#
# Common: System
main.o : $(main).cpp $(dprefix) $(NairnFEA).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(main).cpp
CommonAnalysis.o : $(CommonAnalysis).cpp $(dprefix) $(FEAArchiveData).hpp $(FEAReadHandler).hpp $(CommonException).hpp \
			$(StrX).hpp $(MaterialBase).hpp $(NodalPoint).hpp $(ElementBase).hpp $(CommonReadHandler).hpp \
			$(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonAnalysis).cpp
CommonArchiveData.o : $(CommonArchiveData).cpp $(dprefix) $(CommonArchiveData).hpp $(NodalPoint).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonArchiveData).cpp
LinkedObject.o : $(LinkedObject).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinkedObject).cpp
CommonUtilities.o : $(CommonUtilities).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonUtilities).cpp

# Common: Read_XML
CommonReadHandler.o : $(CommonReadHandler).cpp $(dprefix) $(CommonReadHandler).hpp $(NairnFEA).hpp \
			$(MaterialController).hpp $(NodesController).hpp $(ElementsController).hpp $(StrX).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonReadHandler).cpp
BitMapFilesCommon.o : $(BitMapFilesCommon).cpp $(dprefix) $(CommonReadHandler).hpp $(BMPLevel).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFilesCommon).cpp
BMPLevel.o : $(BMPLevel).cpp $(dprefix) $(BMPLevel).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BMPLevel).cpp
ParseController.o : $(ParseController).cpp $(dprefix) $(ParseController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ParseController).cpp
ElementsController.o : $(ElementsController).cpp $(dprefix) $(ElementsController).hpp $(FourNodeIsoparam).hpp \
			$(NodesController).hpp $(ParseController).hpp $(Linear2D).hpp $(ElementBase).hpp $(EightNodeIsoparam).hpp \
			$(SixNodeTriangle).hpp $(CSTriangle).hpp $(LinearInterface).hpp $(QuadInterface).hpp $(Interface2D).hpp \
			$(Quad2D).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementsController).cpp
NodesController.o : $(NodesController).cpp $(dprefix) $(NodesController).hpp $(NodalPoint).hpp $(ParseController).hpp \
			$(NodalPoint2D).hpp $(NodalPoint3D.hpp) $(NairnFEA).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodesController).cpp
ShapeController.o : $(ShapeController).cpp $(dprefix) $(ShapeController).hpp $(NodalPoint).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ShapeController).cpp
LineController.o : $(LineController).cpp $(dprefix) $(LineController).hpp $(ShapeController).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LineController).cpp
BoxController.o : $(BoxController).cpp $(dprefix) $(BoxController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BoxController).cpp
ArcController.o : $(ArcController).cpp $(dprefix) $(ArcController).hpp $(LineController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ArcController).cpp
PointController.o : $(PointController).cpp $(dprefix) $(PointController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PointController).cpp
RectController.o : $(RectController).cpp $(dprefix) $(RectController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(RectController).cpp
MaterialController.o : $(MaterialController).cpp  $(dprefix) $(MaterialController).hpp $(ParseController).hpp \
		$(IsotropicMat).hpp $(TransIsotropic).hpp $(Orthotropic).hpp $(ImperfectInterface).hpp \
		$(Elastic).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialController).cpp
mathexpr.o : $(mathexpr).cpp $(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(mathexpr).cpp

# Common: Exceptions
CommonException.o : $(CommonException).cpp $(dprefix) $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonException).cpp

# Common: Materials
MaterialBase.o : $(MaterialBase).cpp $(dprefix) $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialBase).cpp
Elastic.o : $(Elastic).cpp $(dprefix) $(Elastic).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Elastic).cpp
IsotropicMat.o : $(IsotropicMat).cpp $(dprefix) $(IsotropicMat).hpp $(MaterialBase).hpp $(Elastic).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(IsotropicMat).cpp
TransIsotropic.o : $(TransIsotropic).cpp $(dprefix) $(TransIsotropic).hpp $(MaterialBase).hpp $(Elastic).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TransIsotropic).cpp
Orthotropic.o : $(Orthotropic).cpp $(dprefix) $(Orthotropic).hpp $(TransIsotropic).hpp $(MaterialBase).hpp $(Elastic).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Orthotropic).cpp

# Common: Elements
ElementBase.o : $(ElementBase).cpp $(dprefix) $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementBase).cpp
Linear2D.o : $(Linear2D).cpp $(dprefix) $(Linear2D).hpp $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Linear2D).cpp
FourNodeIsoparam.o : $(FourNodeIsoparam).cpp $(dprefix) $(FourNodeIsoparam).hpp $(Linear2D).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FourNodeIsoparam).cpp

# Common: Nodes
NodalPoint.o : $(NodalPoint).cpp $(dprefix) $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint).cpp
NodalPoint2D.o : $(NodalPoint2D).cpp $(dprefix) $(NodalPoint2D).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint2D).cpp
NodalPoint3D.o : $(NodalPoint3D).cpp $(dprefix) $(NodalPoint3D).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint3D).cpp
	
# FEA: System
FEAArchiveData.o : $(FEAArchiveData).cpp $(dprefix) $(FEAArchiveData).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEAArchiveData).cpp
Utilities.o : $(Utilities).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Utilities).cpp

# FEA: NairnFEA_Class
NairnFEA.o : $(NairnFEA).cpp $(dprefix) $(NairnFEA).hpp $(ElementBase).hpp $(CommonException).hpp \
			 $(NodalDispBC).hpp $(NodalLoad).hpp $(EdgeBC).hpp $(NodalPoint).hpp $(FEAArchiveData).hpp \
			 $(CommonArchiveData).hpp $(Constraint).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NairnFEA).cpp
FEAStartResults.o : $(FEAStartResults).cpp $(dprefix) $(NairnFEA).hpp $(ElementBase).hpp $(CommonException).hpp \
			 $(NodalDispBC).hpp $(NodalLoad).hpp $(EdgeBC).hpp $(NodalPoint).hpp $(MaterialBase).hpp \
			 $(mathexpr).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEAStartResults).cpp

# FEA: Read_FEA
FEAReadHandler.o : $(FEAReadHandler).cpp $(dprefix) $(FEAReadHandler).hpp $(CommonReadHandler).hpp \
			$(NairnFEA).hpp $(ElementBase).hpp $(NodalDispBC).hpp $(NodalLoad).hpp $(EdgeBC).hpp $(NodalPoint).hpp \
			$(Keypoint).hpp $(KeypointsController).hpp $(ParseController).hpp $(Path).hpp $(PathsController).hpp \
			$(Area).hpp $(NodesController).hpp $(ElementsController).hpp $(EdgeBCController).hpp \
			$(LineController).hpp $(MaterialController).hpp $(NodalDispBCController).hpp $(NodalLoadController).hpp \
			$(Constraint).hpp $(ConstraintController).hpp $(mathexpr).hpp $(FEABoundaryCondition).hpp \
			$(ShapeController).hpp $(PathBCController).hpp $(PointController).hpp $(ArcController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEAReadHandler).cpp
BitMapFilesFEA.o : $(BitMapFilesFEA).cpp $(dprefix) $(FEAReadHandler).hpp $(CommonArchiveData).hpp $(RectController).hpp \
			$(ShapeController).hpp $(ElementBase).hpp $(FEAArchiveData).hpp $(CommonArchiveData).hpp $(BMPLevel).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFilesFEA).cpp
KeypointsController.o : $(KeypointsController).cpp $(dprefix) $(KeypointsController).hpp $(ParseController).hpp \
			$(Keypoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(KeypointsController).cpp
Keypoint.o : $(Keypoint).cpp $(dprefix) $(Keypoint).hpp $(NodesController).hpp $(ParseController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Keypoint).cpp
PathsController.o : $(PathsController).cpp $(dprefix) $(PathsController).hpp $(ParseController).hpp \
			$(Path).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PathsController).cpp
Path.o : $(Path).cpp $(dprefix) $(Path).hpp $(KeypointsController).hpp $(Keypoint).hpp $(EdgeBC).hpp \
			$(EdgeBCController).hpp $(ParseController).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Path).cpp
Area.o : $(Area).cpp $(dprefix) $(Area).hpp $(KeypointsController).hpp $(Keypoint).hpp $(ParseController).hpp \
			$(EdgeBCController).hpp $(PathsController).hpp $(ElementsController).hpp $(NodesController).hpp \
			$(Path).hpp $(ElementBase).hpp $(EightNodeIsoparam).hpp $(Quad2D).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Area).cpp 
EdgeBCController.o : $(EdgeBCController).cpp $(dprefix) $(EdgeBCController).hpp $(ParseController).hpp \
			$(ElementBase).hpp $(EdgeBC).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EdgeBCController).cpp
NodalDispBCController.o : $(NodalDispBCController).cpp $(dprefix) $(NodalDispBCController).hpp $(ParseController).hpp \
			$(NodalDispBC).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalDispBCController).cpp
NodalLoadController.o : $(NodalLoadController).cpp $(dprefix) $(NodalLoadController).hpp $(ParseController).hpp \
			$(NodalLoad).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalLoadController).cpp
ConstraintController.o : $(ConstraintController).cpp $(dprefix) $(Constraint).hpp $(ParseController).hpp \
			$(ConstraintController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ConstraintController).cpp
PathBCController.o : $(PathBCController).cpp $(dprefix) $(PathBCController).hpp $(ShapeController).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PathBCController).cpp

# FEA: Materials
MaterialBaseFEA.o : $(MaterialBaseFEA).cpp $(dprefix) $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialBaseFEA).cpp
ImperfectInterface.o : $(ImperfectInterface).cpp $(dprefix) $(ImperfectInterface).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ImperfectInterface).cpp

# FEA: Elements
MoreElementBase.o : $(MoreElementBase).cpp $(dprefix) $(NairnFEA).hpp $(ElementBase).hpp $(MaterialBase).hpp \
			$(NodalPoint).hpp $(CommonException).hpp $(FourNodeIsoparam).hpp $(EightNodeIsoparam).hpp \
			$(SixNodeTriangle).hpp $(Quad2D).hpp $(Linear2D).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoreElementBase).cpp
Quad2D.o : $(Quad2D).cpp $(dprefix) $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Quad2D).cpp
EightNodeIsoparam.o : $(EightNodeIsoparam).cpp $(dprefix) $(EightNodeIsoparam).hpp $(Quad2D).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EightNodeIsoparam).cpp
SixNodeTriangle.o : $(SixNodeTriangle).cpp $(dprefix) $(SixNodeTriangle).hpp $(Quad2D).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SixNodeTriangle).cpp
CSTriangle.o : $(CSTriangle).cpp $(dprefix) $(CSTriangle).hpp $(Linear2D).hpp $(ElementBase).hpp \
			$(NairnFEA).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CSTriangle).cpp
Interface2D.o : $(Interface2D).cpp $(dprefix) $(Interface2D).hpp $(ElementBase).hpp \
			$(NairnFEA).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Interface2D).cpp
LinearInterface.o : $(LinearInterface).cpp $(dprefix) $(LinearInterface).hpp $(Interface2D).hpp $(ElementBase).hpp \
			$(NairnFEA).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinearInterface).cpp
QuadInterface.o : $(QuadInterface).cpp $(dprefix) $(QuadInterface).hpp $(Interface2D).hpp $(ElementBase).hpp \
			$(NairnFEA).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(QuadInterface).cpp

# FEA: Nodes
MoreNodalPoint.o : $(MoreNodalPoint).cpp $(dprefix) $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoreNodalPoint).cpp

# FEA: Boundary Conditions
FEABoundaryCondition.o : $(FEABoundaryCondition).cpp $(dprefix) $(FEABoundaryCondition).hpp $(NodalPoint).hpp $(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEABoundaryCondition).cpp
NodalDispBC.o : $(NodalDispBC).cpp $(dprefix) $(NodalDispBC).hpp $(NodalPoint).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalDispBC).cpp
NodalLoad.o : $(NodalLoad).cpp $(dprefix) $(NodalLoad).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalLoad).cpp
EdgeBC.o : $(EdgeBC).cpp $(dprefix) $(EdgeBC).hpp $(ElementBase).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EdgeBC).cpp
Constraint.o : $(Constraint).cpp $(dprefix) $(Constraint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Constraint).cpp

# FEA: Numerical
Gauss.o : $(Gauss).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Gauss).cpp

# -------------------------------------------------------------------------
# To make executable
all : NairnFEA
        
# -------------------------------------------------------------------------
# To clean compiled objects and executable      
clean : 
	$(RM) $(objects) $(output)

# -------------------------------------------------------------------------
# To clean compiled objects only      
mostlyclean :
	$(RM) $(objects)

# -------------------------------------------------------------------------
# To install        
install : 
	$(CP) $(output) $(dtdpath) $(ioutput)
