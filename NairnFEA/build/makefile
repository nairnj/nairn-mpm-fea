# -------------------------------------------------------------------------
# To compile NairnFEA use the following commands
#
# cd to directory with this 'makefile'
# give the command 'make SYSTEM=system' where the lower case system is where
#    you are compiling. Supported systems are:
#       mac -  use Apple's compiler assuming you have install libomp.dylib
#		mac-clang - this is changed in step 1 to your preferred clang version among
#                   (can also specify version to override the default)
#         mac-clang17 - MacOS X using clang-mp 17
#		  mac-clang11 - MacOS X using clang-mp 11
#         mac-clang9 - MacOS X using clang-mp 9.0
#         mac-clang7 - MacOS X using clang-mp 7.0
#       carbon - OSU cluster with installed gcc7 (the default)
#       mpm - OSU machine
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
# 3. $(CC) is the gcc compiler and $(CCLINK) is linker,
#	 $(RM) is the rm command, $(CP) is copy for install,
# 4. $(CFLAGS) is the c++ compiler options
#    $(LFLAGS) is the c++ linking options
# 5. $(output) is path to save executable
#    $(dtdpath) is path to save DTD file
# 6. $(ioutput) is path to install executable
# 7. $(src) and $(com) are paths to parts of the source code
#
# Each of these can changed by editing or can be overriden at make time
# and are documented more in the numbered section below
#
# -------------------------------------------------------------------------
# 1. Path and/or command to link to xerces library
#    The default is standard Linux dynamic library install
SYSTEM ?= carbon
# You can pick a preferred clang-mp when user uses mac-clang without a number
ifeq ($(SYSTEM),mac-clang)
    override SYSTEM=mac-clang17
endif
xercesSo= /usr/local/lib/libxerces-c.so
ifeq ($(findstring mac-clang,$(SYSTEM)),mac-clang)
    xercesSo = /usr/local/lib/libxerces-c.dylib
endif
ifeq ($(SYSTEM),carbon)
	xercesSo = /usr/local/xerces3/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),mpm)
	xercesSo = /usr/lib/x86_64-linux-gnu/libxerces-c.so
endif
ifeq ($(SYSTEM),zia)
	xercesSo = /usr/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),mac)
    xercesSo = /usr/local/lib/libxerces-c.dylib
endif

$(info *** Compiling for system <$(SYSTEM)> using xercesSo <$(xercesSo)>)

# 2. Path to xerces headers (if not at default /usr/local/include)
headersPath= /usr/local/include
ifeq ($(SYSTEM),carbon)
	headersPath = /usr/local/xerces3/include
endif
ifeq ($(SYSTEM),mpm)
	headersPath = /usr/include
endif
ifeq ($(SYSTEM),zia)
	headersPath = /usr/include/xercesc
endif

$(info *** Path to xerces headers is <$(headersPath)>)

# 3. Define compiler to use and other commands
#     CC is path to the c++ compiler, must support OpenMP
#     CCLINK is path to c++ linker
#     RM is path/options to the file remove command (used only in make clean)
#     CP is path/options to copy command (used only in make install
CC= g++
CCLINK= g++
ifeq ($(SYSTEM),carbon)
    CC= /usr/bin/g++
    CCLINK= /usr/bin/g++
endif
ifeq ($(SYSTEM),mac)
    CC= clang
    CCLINK= clang++
endif
ifeq ($(SYSTEM),mac-clang17)
    CC= /opt/local/bin/clang-mp-17
    CCLINK= /opt/local/bin/clang++-mp-17
endif
ifeq ($(SYSTEM),mac-clang11)
    CC= /opt/local/bin/clang-mp-11
    CCLINK= /opt/local/bin/clang++-mp-11
endif
ifeq ($(SYSTEM),mac-clang9)
    CC= /opt/local/bin/clang-mp-9.0
    CCLINK= /opt/local/bin/clang++-mp-9.0
endif
ifeq ($(SYSTEM),mac-clang7)
    CC= /opt/local/bin/clang-mp-7.0
    CCLINK= /opt/local/bin/clang++-mp-7.0
endif
RM= rm
CP= cp

$(info *** Comiler is <$(CC)> and linker is <$(CCLINK)>)

# 4. Compiler options -c=compile only is required, others are optional
#     CFLAGS are flags passed to compile each source file. If must contain -c and
#	  can add any other needed options. Some flags found useful in this project are
#            -O3 = optimize most -O0 = not optimized
#            -fopenmp = enable OpenMP compile (link too)
#            -g = debugging info  -pg = profiling
#            -arch i386 = compile architecture (MacOS X only, other archs are x86_64, ppc, and ppc64)
#     LFLAGS are flags passed at link time. The default is -fopenmp. Some flags found useful
#       in this project are
#            -fopenmp = enable OpenMP link (compile too)
#            -arch i386 = link architecture (MacOS X only, other archs are x86_64, ppc, and ppc64)
#            -pg = profiling
#			 -Wl,-rpath,<path_to_libraries> = add path to libraries and this appears needed
#                   in clang-14 or newer (or maybe caused by MacOS Sonoma command line tools)
#     If using -fopenmp, -arch, or -pg, must have in both CFLAGS and LFLAGS. For MacOS X, the specified arch
#       must match the xerces library used for linking
CFLAGS= -c -O3 -fopenmp -std=c++11
LFLAGS= -fopenmp
ifeq ($(SYSTEM),mac-clang17)
    LFLAGS= -fopenmp -lc++ -Wl,-rpath,/opt/local/libexec/llvm-17/lib
else
    ifeq ($(findstring mac-clang,$(SYSTEM)),mac-clang)
		LFLAGS= -fopenmp -lc++
	endif
endif
ifeq ($(SYSTEM),mac)
    CFLAGS= -x c++ -c -O3 -Xclang -fopenmp -std=c++11
    LFLAGS= -lomp
endif
ifeq ($(SYSTEM),carbon)
    LFLAGS = -fopenmp -Wl,-rpath=/usr/lib64 -L/usr/lib64
endif

$(info *** Linking with flags <$(LFLAGS)>)

# 5. Define executable destination
#     output = relative or full path to save executable or to find it for a subsequent make install
#     dtdpath = relative of full path to DTD file, only used during make install
output = ../input/NairnFEA
dtdpath = ../input/NairnFEA.dtd

# 6. Provide relative of full path to install location, only used during make install
ioutput = ~/bin

# 7. Define paths to all files (relative to 'makefile')
#     The commented com option is needed if you follow the suggested Eclipse procedure
src = ../src
com = ../../Common
std = ../../Common/Headers

headers = -I$(src) -I$(std) -I$(com) -I$(headersPath)
ifeq ($(findstring mac-clang,$(SYSTEM)),mac-clang)
    headers = -I$(src) -I$(std) -I$(com) -I$(headersPath) -I/opt/local/include/libomp
endif
ArcController = $(com)/Read_XML/ArcController
Area = $(src)/Read_FEA/Area
Atomic = $(com)/Read_XML/Atomic
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
CSTriangle = $(com)/Elements/CSTriangle
EdgeBC = $(src)/Boundary_Conditions/EdgeBC
EdgeBCController = $(src)/Read_FEA/EdgeBCController
EightNodeIsoparam = $(src)/Elements/EightNodeIsoparam
Elastic = $(com)/Materials/Elastic
ElementBase = $(com)/Elements/ElementBase
ElementsController = $(com)/Read_XML/ElementsController
Expression = $(com)/Read_XML/Expression
FEAArchiveData = $(src)/System/FEAArchiveData
FEABoundaryCondition = $(src)/Boundary_Conditions/FEABoundaryCondition
FEAReadHandler = $(src)/Read_FEA/FEAReadHandler
FEAStartResults = $(src)/NairnFEA_Class/FEAStartResults
FourNodeIsoparam = $(com)/Elements/FourNodeIsoparam
Gauss = $(com)/System/Gauss
ImperfectInterface = $(src)/Materials/ImperfectInterface
Interface2D = $(src)/Elements/Interface2D
IsotropicMat = $(com)/Materials/IsotropicMat
Keypoint = $(src)/Read_FEA/Keypoint
KeypointsController = $(src)/Read_FEA/KeypointsController
Lagrange2D = $(com)/Elements/Lagrange2D
Linear2D = $(com)/Elements/Linear2D
LinearInterface = $(src)/Elements/LinearInterface
LineController = $(com)/Read_XML/LineController
LinkedObject = $(com)/System/LinkedObject
main = $(com)/System/main
MaterialBase = $(com)/Materials/MaterialBase
MaterialBaseFEA = $(src)/Materials/MaterialBaseFEA
MaterialController = $(com)/Read_XML/MaterialController
MatRegionFEA = $(src)/Read_FEA/MatRegionFEA
MoreElementBase = $(src)/Elements/MoreElementBase
MoreNodalPoint = $(src)/Nodes/MoreNodalPoint
NairnFEA = $(src)/NairnFEA_Class/NairnFEA
NodalDispBC = $(src)/Boundary_Conditions/NodalDispBC
NodalDispBCController = $(src)/Read_FEA/NodalDispBCController
NodalLoad = $(src)/Boundary_Conditions/NodalLoad
NodalLoadController = $(src)/Read_FEA/NodalLoadController
NodalPoint = $(com)/Nodes/NodalPoint
NodesController = $(com)/Read_XML/NodesController
Orthotropic = $(com)/Materials/Orthotropic
OvalController = $(com)/Read_XML/OvalController
ParseController = $(com)/Read_XML/ParseController
Path = $(src)/Read_FEA/Path
PathBCController = $(src)/Read_FEA/PathBCController
PathsController = $(src)/Read_FEA/PathsController
PointController = $(com)/Read_XML/PointController
PolygonController = $(com)/Read_XML/PolygonController
Quad2D = $(com)/Elements/Quad2D
QuadInterface = $(src)/Elements/QuadInterface
RectController = $(com)/Read_XML/RectController
ShapeController = $(com)/Read_XML/ShapeController
SixNodeTriangle = $(src)/Elements/SixNodeTriangle
StrX = $(com)/Exceptions/StrX
svninfo = $(com)/System/svninfo
TransIsotropic = $(com)/Materials/TransIsotropic
UnitsController = $(com)/System/UnitsController
Utilities = $(src)/System/Utilities
XYFileImporter = $(com)/Read_XML/XYFileImporter
XYBMPImporter = $(com)/Read_XML/XYBMPImporter
XYTXTImporter = $(com)/Read_XML/XYTXTImporter

# path from 'makefile' to prefix header file
prefix = $(src)/System/FEAPrefix.hpp

# path from 'makefile' to default headers (this included by $(prefix)
dprefix = $(src)/System/FEAPrefix.hpp $(com)/System/DataTypes.hpp \
		$(com)/System/CommonAnalysis.hpp $(com)/System/LinkedObject.hpp $(com)/System/UnitsController.hpp

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
# Lagrange2D.hpp : Quad2D.hpp (ElementBase.hpp)
# Linear2D.hpp : ElementBase.hpp
# LinearInterface : Interface2D.hpp (ElementBase.hpp)
# LineController.hpp : ShapeController.hpp
# MaterialController.hpp : ParseController.hpp
# NodalDispBC.hpp : FEABoundaryCondition.hpp
# NodalDispBCController.hpp : ParseController.hpp
# NodalLoad.hpp : FEABoundaryCondition.hpp
# NodalLoadController.hpp : ParseController.hpp
# NodesController.hpp : ParseController.hpp
# Orthotropic.hpp : TransIsotropic.hpp (Elastic.hpp MaterialBase.hpp)
# OvalController.hpp : ShapeController.hpp
# PathBCController.hpp : ShapeController.hpp
# PathsController.hpp : ParseController.hpp
# PointController.hpp : ShapeController.hpp
# PolygonController.hpp : ShapeController.hpp
# Quad2D.hpp : ElementBase.hpp
# QuadInterface : Interface2D.hpp (ElementBase.hpp)
# RectController.hpp : ShapeController.hpp
# SixNodeTriangle.hpp : Quad2D.hpp (ElementBase.hpp)
# TransIsotropic.hpp : Elastic.hpp (MaterialBase.hpp)
# XYBMPImporter.hpp : XYFileImporter.hpp, CommonReadHandler.hpp
# XYFileImporter.hpp : CommonReadHandler.hpp
# XYTXTImporter.hpp : XYFileImporter.hpp, CommonReadHandler.hpp

# Header files with no dependencies
# Area.hpp CommonArchiveData.hpp CommonException.hpp CommonReadHandler.hpp FEABoundaryCondition.hpp
# ElementBase.hpp Keypoint.hpp MaterialBase.hpp NairnFEA.hpp
# NodalPoint.hpp ParseController.hpp Path.hpp ShapeController.hpp StrX.hpp
# Expression.hpp Atomic.hpp

# -------------------------------------------------------------------------
# all compiled objects
objects = main.o CommonAnalysis.o CommonArchiveData.o LinkedObject.o CommonUtilities.o svninfo.o \
		CommonException.o ElementBase.o Linear2D.o FourNodeIsoparam.o MaterialBase.o \
		CommonReadHandler.o ParseController.o ElementsController.o NodesController.o \
		LineController.o MaterialController.o IsotropicMat.o TransIsotropic.o BoxController.o \
		Orthotropic.o NodalPoint.o FEAArchiveData.o Utilities.o NairnFEA.o FEAStartResults.o \
		FEAReadHandler.o KeypointsController.o Keypoint.o PathsController.o Path.o Area.o \
		XYFileImporter.o XYBMPImporter.o XYTXTImporter.o \
		EdgeBCController.o NodalDispBCController.o NodalLoadController.o MaterialBaseFEA.o \
		ImperfectInterface.o MoreElementBase.o Quad2D.o EightNodeIsoparam.o SixNodeTriangle.o \
		CSTriangle.o Interface2D.o LinearInterface.o QuadInterface.o MoreNodalPoint.o \
		NodalDispBC.o NodalLoad.o EdgeBC.o Gauss.o Elastic.o Constraint.o ConstraintController.o \
		FEABoundaryCondition.o ShapeController.o PointController.o PathBCController.o \
		RectController.o ArcController.o BitMapFilesCommon.o \
		OvalController.o BMPLevel.o BitMapFilesFEA.o Lagrange2D.o MatRegionFEA.o \
		PolygonController.o UnitsController.o Expression.o Atomic.o

# -------------------------------------------------------------------------
# Link all objects
NairnFEA : $(objects)
	$(CCLINK) $(LFLAGS) -o $(output) $(objects) $(xercesSo)

# -------------------------------------------------------------------------
# All source code files
#
# Common: System
main.o : $(main).cpp $(dprefix) $(NairnFEA).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(main).cpp
CommonAnalysis.o : $(CommonAnalysis).cpp $(dprefix) $(CommonException).hpp \
			$(StrX).hpp $(MaterialBase).hpp $(NodalPoint).hpp $(ElementBase).hpp $(CommonReadHandler).hpp \
			$(CommonArchiveData).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonAnalysis).cpp
CommonArchiveData.o : $(CommonArchiveData).cpp $(dprefix) $(CommonArchiveData).hpp $(NodalPoint).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonArchiveData).cpp
LinkedObject.o : $(LinkedObject).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinkedObject).cpp
CommonUtilities.o : $(CommonUtilities).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonUtilities).cpp
svninfo.o : $(svninfo).cpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(svninfo).cpp
UnitsController.o : $(UnitsController).cpp $(UnitsController).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UnitsController).cpp
Gauss.o : $(Gauss).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Gauss).cpp

# Common: Read_XML
CommonReadHandler.o : $(CommonReadHandler).cpp $(dprefix) $(CommonReadHandler).hpp \
			$(MaterialController).hpp $(NodesController).hpp $(ElementsController).hpp $(StrX).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonReadHandler).cpp
XYFileImporter.o : $(XYFileImporter).cpp $(dprefix) $(XYFileImporter).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(XYFileImporter).cpp
XYBMPImporter.o : $(XYBMPImporter).cpp $(dprefix) $(XYBMPImporter).hpp $(XYFileImporter).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(XYBMPImporter).cpp
XYTXTImporter.o : $(XYTXTImporter).cpp $(dprefix) $(XYTXTImporter).hpp $(XYFileImporter).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(XYTXTImporter).cpp
BitMapFilesCommon.o : $(BitMapFilesCommon).cpp $(dprefix) $(CommonReadHandler).hpp $(BMPLevel).hpp $(MaterialController).hpp \
			$(XYFileImporter).hpp $(XYBMPImporter).hpp $(XYTXTImporter).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFilesCommon).cpp
BMPLevel.o : $(BMPLevel).cpp $(dprefix) $(BMPLevel).hpp $(NairnFEA).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BMPLevel).cpp
ParseController.o : $(ParseController).cpp $(dprefix) $(ParseController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ParseController).cpp
ElementsController.o : $(ElementsController).cpp $(dprefix) $(ElementsController).hpp $(FourNodeIsoparam).hpp \
			$(NodesController).hpp $(ParseController).hpp $(Linear2D).hpp $(ElementBase).hpp $(EightNodeIsoparam).hpp \
			$(SixNodeTriangle).hpp $(CSTriangle).hpp $(LinearInterface).hpp $(QuadInterface).hpp $(Interface2D).hpp \
			$(Quad2D).hpp  $(Lagrange2D).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementsController).cpp
NodesController.o : $(NodesController).cpp $(dprefix) $(NodesController).hpp $(NodalPoint).hpp $(ParseController).hpp \
			$(NairnFEA).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodesController).cpp
ShapeController.o : $(ShapeController).cpp $(dprefix) $(ShapeController).hpp $(NodalPoint).hpp $(CommonReadHandler).hpp \
            $(ElementBase).hpp
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
OvalController.o : $(OvalController).cpp $(dprefix) $(OvalController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(OvalController).cpp
PolygonController.o : $(PolygonController).cpp $(dprefix) $(PolygonController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PolygonController).cpp
MaterialController.o : $(MaterialController).cpp  $(dprefix) $(MaterialController).hpp $(ParseController).hpp \
		$(IsotropicMat).hpp $(TransIsotropic).hpp $(Orthotropic).hpp $(ImperfectInterface).hpp \
		$(Elastic).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialController).cpp
Atomic.o : $(Atomic).cpp $(Atomic).hpp  $(Expression).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Atomic).cpp
Expression.o : $(Expression).cpp $(Expression).hpp $(Atomic).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Expression).cpp

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
Quad2D.o : $(Quad2D).cpp $(dprefix) $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Quad2D).cpp
Lagrange2D.o : $(Lagrange2D).cpp $(dprefix) $(Lagrange2D).hpp $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Lagrange2D).cpp
CSTriangle.o : $(CSTriangle).cpp $(dprefix) $(CSTriangle).hpp $(Linear2D).hpp $(ElementBase).hpp \
			$(NairnFEA).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CSTriangle).cpp

# Common: Nodes
NodalPoint.o : $(NodalPoint).cpp $(dprefix) $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint).cpp
	
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
			 $(FEABoundaryCondition).hpp $(Expression).hpp $(FEAReadHandler).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEAStartResults).cpp

# FEA: Read_FEA
FEAReadHandler.o : $(FEAReadHandler).cpp $(dprefix) $(FEAReadHandler).hpp $(CommonReadHandler).hpp \
			$(NairnFEA).hpp $(ElementBase).hpp $(NodalDispBC).hpp $(NodalLoad).hpp $(EdgeBC).hpp $(NodalPoint).hpp \
			$(Keypoint).hpp $(KeypointsController).hpp $(ParseController).hpp $(Path).hpp $(PathsController).hpp \
			$(Area).hpp $(NodesController).hpp $(ElementsController).hpp $(EdgeBCController).hpp \
			$(LineController).hpp $(MaterialController).hpp $(NodalDispBCController).hpp $(NodalLoadController).hpp \
			$(Constraint).hpp $(ConstraintController).hpp $ $(FEABoundaryCondition).hpp \
			$(ShapeController).hpp $(PathBCController).hpp $(PointController).hpp $(ArcController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEAReadHandler).cpp
BitMapFilesFEA.o : $(BitMapFilesFEA).cpp $(dprefix) $(FEAReadHandler).hpp $(CommonReadHandler).hpp $(RectController).hpp \
			$(ShapeController).hpp $(ElementBase).hpp $(FEAArchiveData).hpp $(CommonArchiveData).hpp $(BMPLevel).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFilesFEA).cpp
MatRegionFEA.o : $(MatRegionFEA).cpp $(dprefix) $(FEAReadHandler).hpp $(CommonReadHandler).hpp $(RectController).hpp \
			$(ShapeController).hpp $(ElementBase).hpp $(MaterialController).hpp $(OvalController).hpp $(ParseController).hpp \
			$(PolygonController).hpp $(Expression).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatRegionFEA).cpp
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
			$(Path).hpp $(ElementBase).hpp $(EightNodeIsoparam).hpp $(Quad2D).hpp $(Expression).hpp
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
EightNodeIsoparam.o : $(EightNodeIsoparam).cpp $(dprefix) $(EightNodeIsoparam).hpp $(Quad2D).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EightNodeIsoparam).cpp
SixNodeTriangle.o : $(SixNodeTriangle).cpp $(dprefix) $(SixNodeTriangle).hpp $(Quad2D).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SixNodeTriangle).cpp
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
FEABoundaryCondition.o : $(FEABoundaryCondition).cpp $(dprefix) $(FEABoundaryCondition).hpp $(NodalPoint).hpp \
			$(Expression).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FEABoundaryCondition).cpp
NodalDispBC.o : $(NodalDispBC).cpp $(dprefix) $(NodalDispBC).hpp $(NodalPoint).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalDispBC).cpp
NodalLoad.o : $(NodalLoad).cpp $(dprefix) $(NodalLoad).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalLoad).cpp
EdgeBC.o : $(EdgeBC).cpp $(dprefix) $(EdgeBC).hpp $(ElementBase).hpp $(FEABoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EdgeBC).cpp
Constraint.o : $(Constraint).cpp $(dprefix) $(Constraint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Constraint).cpp

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
