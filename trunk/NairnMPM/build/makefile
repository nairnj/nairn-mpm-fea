# -------------------------------------------------------------------------
# To compile NairnMPM use the following commands
#
# cd to directory with this 'makefile'
# give the command 'make SYSTEM=system' where the lower case system is where
#    you are compiling. Supported systems are:
#       mac - MacOS X
#       carbon - OSU carbon/hedera cluster (the default)
#       mpm - OSU machine
#       cygwin - Standard cygwin install in Windows
#       cluster2 - Grant Smith's cluster at U of Utah
#       inferno - C-SAFE cluster at U of Utah
#       zia - A Debian linux machine
#       artemis - for Edward Le on FP Innovations cluster
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
#
# Each of these can changed by editing or can be overriden at make time with syntax such as
#
# make CFLAGS="-c -arch x86_64" LFLAGS="-arch x84_64" output="~/Desktop/NairnFEA" SYSTEM=mac
#
# -------------------------------------------------------------------------
# 1. Path and/or command to link to xerces library
#    The default is standard Linux dynamic library install
SYSTEM=carbon
xercesSo= /usr/local/lib/libxerces-c.so
ifeq ($(SYSTEM),mac)
	xercesSo = /usr/local/lib/libxerces-c.dylib
endif
ifeq ($(SYSTEM),carbon)
	xercesSo = /usr/local/xerces3/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),mpm)
	xercesSo = /usr/lib/libxerces-c.so
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
ifeq ($(SYSTEM),artemis)
	xercesSo = /home/edward/xerces/lib/libxerces-c.so
endif

# 2. Path to xerces headers (if not at default /usr/local/include) 
headersPath= /usr/local/include
ifeq ($(SYSTEM),carbon)
	headersPath = /usr/local/xerces3/include
endif
ifeq ($(SYSTEM),mpm)
	headersPath = /usr/include
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
ifeq ($(SYSTEM),artemis)
	headersPath = /home/edward/xerces/include
endif

# 3. Define compiler to use and other commands
#     CC is path to the c++ compiler
#     RM is path/options to the file remove command (used only in make clean)
#     CP is path/options to copy command (used only in make install
CC= g++
RM= rm
CP= cp

# 4. Compiler options -c=compile only, -O3=optimize most -O0=not optimized -g=debugging info
#     CFLAGS are flags passed to compiled for each source file. If must contain -c and
#		can add any other needed options. Some flags found useful in this project are
#            -O3 = optimize most -O0 = not optimized
#            -g = debuggin info  -pg = profiling
#            -arch i386 = compile architecture (MacOS X only, other archs are x86_64, ppc, and ppc64)
#     LFLAGS are flags passed at link time. The default is empty. Some flags found useful
#       in this project are
#            -arch i386 = link architecture (MacOS X only, other archs are x86_64, ppc, and ppc64)
#            -pg = profiling
#     If using -arch or -pg, must have in both CFLAGS and LFLAGS. For MacOS X, the specified arch
#       must match the xerces library used for linking
CFLAGS= -c -O3
LFLAGS=

# 5. Define executable destination
#     output = relative or full path to save executable or to find it for a subsequent make install
#     dtdpath = relative of full path to DTD file, only used during make install
output = ../input/NairnMPM
dtdpath = ../input/NairnMPM.dtd

# 7. Define paths to all files (relative to 'makefile')
#     These are usually not changed unless you reorganize the code
ioutput = ~/bin

# 7. Define paths to all files (relative to 'makefile')
#     The comment com option is needed if you follow the suggested Eclipse procedure
src = ../src
com = ../../Common
#com = ../Common

headers = -I$(src) -I$(com) -I$(headersPath) 
AnisoPlasticity = $(src)/Materials/AnisoPlasticity
ArcController = $(com)/Read_XML/ArcController
ArchiveData = $(src)/System/ArchiveData
BitMapFiles = $(src)/Read_MPM/BitMapFiles
BitMapFilesCommon = $(com)/Read_XML/BitMapFilesCommon
BistableIsotropic = $(src)/Materials/BistableIsotropic
BMPLevel = $(com)/Read_XML/BMPLevel
BodyBoxController = $(src)/Read_MPM/BodyBoxController
BodyForce = $(src)/Global_Quantities/BodyForce
BodyCylinderController = $(src)/Read_MPM/BodyCylinderController
BodyObjectController = $(src)/Read_MPM/BodyObjectController
BodyOvalController = $(src)/Read_MPM/BodyOvalController
BodyPolygonController = $(src)/Read_MPM/BodyPolygonController
BodyPolyhedronController = $(src)/Read_MPM/BodyPolyhedronController
BodyRectController = $(src)/Read_MPM/BodyRectController
BodySphereController = $(src)/Read_MPM/BodySphereController
BoundaryCondition = $(src)/Boundary_Conditions/BoundaryCondition
BoxController = $(com)/Read_XML/BoxController
CalcJKTask = $(src)/Custom_Tasks/CalcJKTask
CohesiveZone = $(src)/Materials/CohesiveZone
CommonAnalysis = $(com)/System/CommonAnalysis
CommonArchiveData = $(com)/System/CommonArchiveData
CommonException = $(com)/Exceptions/CommonException
CommonReadHandler = $(com)/Read_XML/CommonReadHandler
CommonTask = $(com)/System/CommonTask
CommonUtilities = $(com)/System/CommonUtilities
ConductionTask = $(src)/Custom_Tasks/ConductionTask
ContourPoint = $(src)/Cracks/ContourPoint
CrackController = $(src)/Read_MPM/CrackController
CrackHeader = $(src)/Cracks/CrackHeader
CrackNode = $(src)/Cracks/CrackNode
CrackSegment = $(src)/Cracks/CrackSegment
CrackSurfaceContact = $(src)/Cracks/CrackSurfaceContact
CrackVelocityField = $(src)/Nodes/CrackVelocityField
CrackVelocityFieldMulti = $(src)/Nodes/CrackVelocityFieldMulti
CrackVelocityFieldSingle = $(src)/Nodes/CrackVelocityFieldSingle
CubicTraction = $(src)/Materials/CubicTraction
CustomTask = $(src)/Custom_Tasks/CustomTask
DiffusionTask = $(src)/Custom_Tasks/DiffusionTask
EightNodeIsoparamBrick = $(src)/Elements/EightNodeIsoparamBrick
Elastic = $(com)/Materials/Elastic
ElasticMPM = $(src)/Materials/ElasticMPM
ElementBase = $(com)/Elements/ElementBase
ElementBase3D = $(src)/Elements/ElementBase3D
ElementsController = $(com)/Read_XML/ElementsController
FourNodeIsoparam = $(com)/Elements/FourNodeIsoparam
Generators = $(src)/Read_MPM/Generators
GlobalQuantity = $(src)/Global_Quantities/GlobalQuantity
GridForcesTask = $(src)/NairnMPM_Class/GridForcesTask
HEAnisotropic = $(src)/Materials/HEAnisotropic
HillPlastic = $(src)/Materials/HillPlastic
HyperElastic = $(src)/Materials/HyperElastic
InitializationTask = $(src)/NairnMPM_Class/InitializationTask
IsoPlasticity = $(src)/Materials/IsoPlasticity
IsotropicMat = $(com)/Materials/IsotropicMat
JohnsonCook = $(src)/Materials/JohnsonCook
Lagrange2D = $(com)/Elements/Lagrange2D
Linear2D = $(com)/Elements/Linear2D
LinearTraction = $(src)/Materials/LinearTraction
LineController = $(com)/Read_XML/LineController
LinkedObject =$(com)/System/LinkedObject
main = $(com)/System/main
MassAndMomentumTask = $(src)/NairnMPM_Class/MassAndMomentumTask
MaterialBase = $(com)/Materials/MaterialBase
MaterialBaseMPM = $(src)/Materials/MaterialBaseMPM
MaterialController = $(com)/Read_XML/MaterialController
mathexpr = $(com)/Read_XML/mathexpr
MatPoint2D = $(src)/MPM_Classes/MatPoint2D
MatPoint3D = $(src)/MPM_Classes/MatPoint3D
MatPtFluxBC = $(src)/Boundary_Conditions/MatPtFluxBC
MatPtLoadBC = $(src)/Boundary_Conditions/MatPtLoadBC
MatVelocityField = $(src)/Nodes/MatVelocityField
MeshInfo = $(src)/NairnMPM_Class/MeshInfo
MGSCGLMaterial = $(src)/Materials/MGSCGLMaterial
Mooney = $(src)/Materials/Mooney
MoreMPMElementBase = $(src)/Elements/MoreMPMElementBase
MoveCracksTask = $(src)/NairnMPM_Class/MoveCracksTask
MPMBase = $(src)/MPM_Classes/MPMBase
MPMReadHandler = $(src)/Read_MPM/MPMReadHandler
MPMTask = $(src)/NairnMPM_Class/MPMTask
MPMTermination = $(src)/Exceptions/MPMTermination
MPMWarnings = $(src)/Exceptions/MPMWarnings
MpsController = $(src)/Read_MPM/MpsController
NairnMPM = $(src)/NairnMPM_Class/NairnMPM
NodalConcBC = $(src)/Boundary_Conditions/NodalConcBC
NodalPoint = $(com)/Nodes/NodalPoint
NodalPoint2D = $(com)/Nodes/NodalPoint2D
NodalPoint3D = $(com)/Nodes/NodalPoint3D
NodalPointMPM = $(src)/Nodes/NodalPointMPM
NodalTempBC = $(src)/Boundary_Conditions/NodalTempBC
NodalVelBC = $(src)/Boundary_Conditions/NodalVelBC
NodesController = $(com)/Read_XML/NodesController
Orthotropic = $(com)/Materials/Orthotropic
ParseController = $(com)/Read_XML/ParseController
PointController = $(com)/Read_XML/PointController
PolyTriangle = $(src)/Read_MPM/PolyTriangle
PropagateTask = $(src)/Custom_Tasks/PropagateTask
Quad2D = $(com)/Elements/Quad2D
RectController = $(com)/Read_XML/RectController
ResetElementsTask = $(src)/NairnMPM_Class/ResetElementsTask
ReverseLoad = $(src)/Custom_Tasks/ReverseLoad
RigidMaterial = $(src)/Materials/RigidMaterial
RunCustomTasksTask = $(src)/NairnMPM_Class/RunCustomTasksTask
ShapeController = $(com)/Read_XML/ShapeController
SLMaterial = $(src)/Materials/SLMaterial
StartOutput = $(src)/NairnMPM_Class/StartOutput
StrX = $(com)/Exceptions/StrX
ThermalRamp = $(src)/Global_Quantities/ThermalRamp
TractionLaw = $(src)/Materials/TractionLaw
TransIsotropic = $(com)/Materials/TransIsotropic
TransportTask = $(src)/Custom_Tasks/TransportTask
TrilinearTraction = $(src)/Materials/TrilinearTraction
UpdateMomentaTask = $(src)/NairnMPM_Class/UpdateMomentaTask
UpdateParticlesTask = $(src)/NairnMPM_Class/UpdateParticlesTask
UpdateStrainsFirstTask = $(src)/NairnMPM_Class/UpdateStrainsFirstTask
UpdateStrainsLastTask = $(src)/NairnMPM_Class/UpdateStrainsLastTask
Viscoelastic = $(src)/Materials/Viscoelastic
VonMisesHardening = $(src)/Materials/VonMisesHardening
VTKArchive = $(src)/Custom_Tasks/VTKArchive
WoodMaterial = $(src)/Materials/WoodMaterial

# path from 'makefile' to prefix header file
prefix = $(src)/System/MPMPrefix.hpp

# path from 'makefile' to default headers (this included by $(prefix)
dprefix = $(src)/System/MPMPrefix.hpp $(com)/System/DataTypes.hpp \
		$(com)/System/CommonAnalysis.hpp $(com)/System/LinkedObject.hpp

# header files that include other header files (other than $(dprefix) files)
# AnisoPlasticity : Orthotropic.hpp (TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)
# ArcController.hpp : LineController.hpp ShapeController.hpp
# ArchiveData.hpp : CommonArchiveData.hpp
# BistableIsotropic.hpp : IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
# BodyBoxController.hpp : BodyObjectController.hpp
# BodyCylinderController.hpp : BodySphereController.hpp BodyObjectController.hpp
# BodyOvalController.hpp : BodyObjectController.hpp
# BodyPolygonController.hpp : BodyObjectController.hpp
# BodyPolyhedronController.hpp : BodyObjectController.hpp
# BodyRectController.hpp : BodyObjectController.hpp
# BodySphereController.hpp : BodyObjectController.hpp
# BoxController.hpp : ShapeController.hpp
# CalcJKTask.hpp : CustomTask.hpp
# CohesiveZone.hpp : MaterialBase.hpp TractionLaw.hpp
# ConductionTask.hpp : TransportTask.hpp
# CrackController.hpp : ParseController.hpp
# CrackNode.hpp : NodalPoint.hpp
# CrackSurfaceContact.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# CrackVelocityField.hpp : MatVelocityField.hpp
# CrackVelocityFieldMulti.hpp : CrackVelocityField.hpp MatVelocityField.hpp
# CrackVelocityFieldSingle.hpp : CrackVelocityField.hpp MatVelocityField.hpp
# CubicTraction.hpp : MaterialBase.hpp TractionLaw.hpp
# DiffusionTask.hpp : TransportTask.hpp
# EightNodeIsoparamBrick.hpp : ElementBase3D.hpp (ElementBase.hpp)
# Elastic.hpp : MaterialBase.hpp
# ElementsController.hpp : ParseController.hpp
# ElementBase3D : ElementBase
# FourNodeIsoparam.hpp : Linear2D.hpp (ElementBase.hpp)
# GridForcesTask.hpp : MPMTask.hpp (CommonTask.hpp)
# HEAnisotropic.hpp : HyperElastic.hpp (MaterialBase.hpp)
# HillPlastic : AnisoPlasticity.hpp (Orthotropic.hpp TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)
# InitializationTask : MPMTask.hpp (CommonTask.hpp)
# IsoPlasticity : IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
# IsotropicMat.hpp : Elastic.hpp (MaterialBase.hpp)
# JohnsonCook.hpp : IsoPlasticity (IsotropicMat.hpp Elastic.hpp MaterialBase.hpp)
# Lagrange2D.hpp : Quad2D.hpp (ElementBase.hpp)
# Linear2D.hpp : ElementBase.hpp
# LinearTraction.hpp : CohesiveZone.hpp MaterialBase.hpp TractionLaw.hpp
# LineController.hpp : ShapeController.hpp
# MassAndMomentumTask.hpp : MPMTask.hpp (CommonTask.hpp)
# MaterialController.hpp : ParseController.hpp
# MatPtLoadBC.hpp : BoundaryCondition.hpp
# MatPtFluxBC.hpp : BoundaryCondition.hpp
# MatPoint2D.hpp : MPMBase.hpp
# MatPoint3D.hpp : MPMBase.hpp
# Mooney.hpp : HyperElastic.hpp (MaterialBase.hpp)
# MoveCracksTask.hpp : MPMTask.hpp (CommonTask.hpp)
# MPMReadHandler.hpp : CommonReadHandler.hpp
# MPMTask.hpp : CommonTask.hpp
# MPMTermination.hpp : CommonException.hpp
# MpsController.hpp : ParseController.hpp
# NodalConcBC.hpp : BoundaryCondition.hpp
# NodalPoint.hpp : CrackVelocityField.hpp MatVelocityField.hpp
# NodalPoint2D.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# NodalPoint3D.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# NodalTempBC.hpp : BoundaryCondition.hpp
# NodalVelBC.hpp : BoundaryCondition.hpp
# NodesController.hpp : ParseController.hpp
# Orthotropic.hpp : TransIsotropic.hpp (Elastic.hpp MaterialBase.hpp)
# PointController.hpp : ShapeController.hpp
# PropagateTask.hpp : CustomTask.hpp
# Quad2D.hpp : ElementBase.hpp
# RectController.hpp : ShapeController.hpp
# ResetElementsTask.hpp : MPMTask.hpp (CommonTask.hpp)
# ReverseLoad.hpp : CustomTask.hpp
# HyperElastic.hpp : MaterialBase.hpp
# RunCustomTasksTask.hpp : MPMTask.hpp (CommonTask.hpp)
# RigidMaterial.hpp : MaterialBase.hpp
# TractionLaw.hpp : MaterialBase.hpp
# TransIsotropic.hpp : Elastic.hpp (MaterialBase.hpp)
# TrilinearTraction : MaterialBase.hpp TractionLaw.hpp CohesiveZone.hpp
# UpdateMomentaTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateParticlesTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateStrainsFirstTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateStrainsLastTask.hpp : MPMTask.hpp (CommonTask.hpp)
# Viscoelastic.hpp : MaterialBase.hpp
# VonMisesHardening.hpp : IsoPlasticity (IsotropicMat.hpp Elastic.hpp MaterialBase.hpp)
# WoodMaterial : HillPlastic.hpp (AnisoPlasticity.hpp, Orthotropic.hpp TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)

# Header files with no dependencies
# BMPLevel.hpp BodyForce.hpp BoundaryCondition.hpp CommonArchiveData.hpp CommonException.hpp
# CommonReadHandler.hpp ContourPoint.hpp CrackHeader.hpp CrackSegment.hpp
# CustomTask.hpp ElementBase.hpp GlobalQuantity.hpp
# MaterialBase.hpp MPMBase.hpp MPMWarnings.hpp NairnMPM.hpp
# ParseController.hpp ShapeController.hpp StrX.hpp ThermalRamp.hpp TransportTask.hpp
# PolyTriangle.hpp

# -------------------------------------------------------------------------
# all compiled objects
objects = main.o CommonAnalysis.o CommonArchiveData.o LinkedObject.o CommonUtilities.o CommonTask.o \
		CommonException.o ElementBase.o Linear2D.o FourNodeIsoparam.o MaterialBase.o \
		CommonReadHandler.o ParseController.o ElementsController.o NodesController.o \
		LineController.o MaterialController.o IsotropicMat.o TransIsotropic.o BoxController.o \
		Orthotropic.o NodalPoint.o NodalPoint2D.o NodalPoint3D.o ArchiveData.o NairnMPM.o MPMTask.o \
		InitializationTask.o MassAndMomentumTask.o UpdateStrainsFirstTask.o GridForcesTask.o UpdateMomentaTask.o \
		UpdateParticlesTask.o UpdateStrainsLastTask.o RunCustomTasksTask.o MoveCracksTask.o ResetElementsTask.o \
		StartOutput.o MPMReadHandler.o Generators.o BitMapFiles.o BitMapFilesCommon.o BMPLevel.o MpsController.o \
		CrackController.o MPMTermination.o MPMWarnings.o BodyForce.o GlobalQuantity.o \
		ThermalRamp.o MPMBase.o MatPoint2D.o CrackHeader.o CrackSegment.o ContourPoint.o \
		CrackSurfaceContact.o CrackNode.o MaterialBaseMPM.o Viscoelastic.o JohnsonCook.o \
		Mooney.o HEAnisotropic.o VonMisesHardening.o BistableIsotropic.o RigidMaterial.o BoundaryCondition.o \
		MatPtLoadBC.o MatPtFluxBC.o NodalVelBC.o NodalConcBC.o NodalTempBC.o MoreMPMElementBase.o \
		NodalPointMPM.o CustomTask.o CalcJKTask.o PropagateTask.o ReverseLoad.o \
		DiffusionTask.o ConductionTask.o MeshInfo.o TransportTask.o ElementBase3D.o \
		EightNodeIsoparamBrick.o MatPoint3D.o Elastic.o ElasticMPM.o HyperElastic.o \
		mathexpr.o ShapeController.o PointController.o RectController.o BodyObjectController.o \
		BodyBoxController.o BodyOvalController.o BodyPolygonController.o BodyRectController.o \
		BodySphereController.o ArcController.o MatVelocityField.o CrackVelocityField.o BodyCylinderController.o \
		TractionLaw.o CohesiveZone.o LinearTraction.o CubicTraction.o IsoPlasticity.o AnisoPlasticity.o \
		HillPlastic.o MGSCGLMaterial.o VTKArchive.o SLMaterial.o CrackVelocityFieldMulti.o CrackVelocityFieldSingle.o \
		WoodMaterial.o TrilinearTraction.o BodyPolyhedronController.o PolyTriangle.o Quad2D.o Lagrange2D.o

# -------------------------------------------------------------------------
# Link all objects
NairnMPM : $(objects)
	$(CC) $(LFLAGS) -o $(output) $(objects) $(xercesSo)

# -------------------------------------------------------------------------
# Common: System
main.o : $(main).cpp $(dprefix) $(NairnMPM).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(main).cpp
CommonAnalysis.o : $(CommonAnalysis).cpp $(dprefix) $(CommonAnalysis).hpp $(ArchiveData).hpp $(MPMReadHandler).hpp $(CommonException).hpp \
			$(StrX).hpp $(MaterialBase).hpp $(NodalPoint).hpp $(ElementBase).hpp $(CommonReadHandler).hpp \
			$(CommonArchiveData).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonAnalysis).cpp
CommonArchiveData.o : $(CommonArchiveData).cpp $(dprefix) $(CommonArchiveData).hpp $(NodalPoint).hpp $(ElementBase).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonArchiveData).cpp
LinkedObject.o : $(LinkedObject).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinkedObject).cpp
CommonUtilities.o : $(CommonUtilities).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonUtilities).cpp
CommonTask.o : $(CommonTask).cpp $(dprefix) $(CommonTask).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonTask).cpp

# Common: Read_XML
CommonReadHandler.o : $(CommonReadHandler).cpp $(dprefix) $(CommonReadHandler).hpp $(NairnMPM).hpp $(MaterialController).hpp \
		$(NodesController).hpp $(ElementsController).hpp $(StrX).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonReadHandler).cpp
BitMapFilesCommon.o : $(BitMapFilesCommon).cpp $(dprefix) $(CommonReadHandler).hpp $(BMPLevel).hpp $(MaterialController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFilesCommon).cpp
BMPLevel.o : $(BMPLevel).cpp $(dprefix) $(BMPLevel).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BMPLevel).cpp
ParseController.o : $(ParseController).cpp $(dprefix) $(ParseController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ParseController).cpp
ElementsController.o : $(ElementsController).cpp $(dprefix) $(ElementsController).hpp $(FourNodeIsoparam).hpp \
			$(NodesController).hpp $(ParseController).hpp $(Linear2D).hpp $(ElementBase).hpp $(EightNodeIsoparamBrick).hpp \
			$(CommonReadHandler).hpp $(Quad2D).hpp  $(Lagrange2D).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementsController).cpp
NodesController.o : $(NodesController).cpp $(dprefix) $(NodesController).hpp $(NodalPoint).hpp $(ParseController).hpp \
			$(NodalPoint2D).hpp $(NodalPoint3D).hpp $(NairnMPM).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodesController).cpp
ShapeController.o : $(ShapeController).cpp $(dprefix) $(ShapeController).hpp $(NodalPoint).hpp $(MPMBase).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonReadHandler).hpp
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
		$(IsotropicMat).hpp $(TransIsotropic).hpp $(Orthotropic).hpp $(Viscoelastic).hpp $(VonMisesHardening).hpp \
		$(Mooney).hpp $(HEAnisotropic).hpp $(BistableIsotropic).hpp $(RigidMaterial).hpp $(MaterialBase).hpp $(JohnsonCook).hpp \
		$(Elastic).hpp $(HyperElastic).hpp $(NairnMPM).hpp $(TractionLaw).hpp $(CohesiveZone).hpp $(LinearTraction).hpp \
		$(CubicTraction).hpp $(AnisoPlasticity).hpp $(IsoPlasticity.hpp) $(HillPlastic).hpp $(MGSCGLMaterial).hpp \
		$(SLMaterial).hpp $(WoodMaterial).hpp $(TrilinearTraction).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialController).cpp
mathexpr.o : $(mathexpr).cpp $(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(mathexpr).cpp

# Common: Exceptions
CommonException.o : $(CommonException).cpp $(dprefix) $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CommonException).cpp

# Common: Materials
MaterialBase.o : $(MaterialBase).cpp $(dprefix) $(MaterialBase).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialBase).cpp
Elastic.o : $(Elastic).cpp $(dprefix) $(Elastic).hpp $(MaterialBase).hpp $(ThermalRamp).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Elastic).cpp
IsotropicMat.o : $(IsotropicMat).cpp $(dprefix) $(IsotropicMat).hpp $(MaterialBase).hpp \
		$(Elastic).hpp $(ConductionTask).hpp $(DiffusionTask).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(IsotropicMat).cpp
TransIsotropic.o : $(TransIsotropic).cpp $(dprefix) $(TransIsotropic).hpp $(MaterialBase).hpp \
		$(Elastic).hpp $(ConductionTask).hpp $(DiffusionTask).hpp $(TransportTask).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TransIsotropic).cpp
Orthotropic.o : $(Orthotropic).cpp $(dprefix) $(Orthotropic).hpp $(TransIsotropic).hpp $(MaterialBase).hpp \
		 $(Elastic).hpp $(ConductionTask).hpp $(DiffusionTask).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Orthotropic).cpp

# Common: Elements
ElementBase.o : $(ElementBase).cpp $(dprefix) $(ElementBase).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementBase).cpp
Linear2D.o : $(Linear2D).cpp $(dprefix) $(Linear2D).hpp $(ElementBase).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Linear2D).cpp
FourNodeIsoparam.o : $(FourNodeIsoparam).cpp $(dprefix) $(FourNodeIsoparam).hpp $(Linear2D).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(MeshInfo).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FourNodeIsoparam).cpp
Quad2D.o : $(Quad2D).cpp $(dprefix) $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Quad2D).cpp
Lagrange2D.o : $(Lagrange2D).cpp $(dprefix) $(Lagrange2D).hpp $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp \
			$(MeshInfo).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Lagrange2D).cpp

# Common: Nodes
NodalPoint.o : $(NodalPoint).cpp $(dprefix) $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint).cpp
NodalPoint2D.o : $(NodalPoint2D).cpp $(dprefix) $(NodalPoint2D).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint2D).cpp
NodalPoint3D.o : $(NodalPoint3D).cpp $(dprefix) $(NodalPoint3D).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint3D).cpp
MatVelocityField.o : $(MatVelocityField).cpp $(dprefix) $(MatVelocityField).hpp $(NairnMPM).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatVelocityField).cpp
CrackVelocityField.o : $(CrackVelocityField).cpp $(dprefix) $(CrackVelocityField).hpp $(MatVelocityField).hpp $(NairnMPM).hpp \
			$(NodalPoint).hpp $(CommonException).hpp $(MeshInfo).hpp $(BoundaryCondition).hpp $(CrackVelocityFieldMulti).hpp \
			$(CrackVelocityFieldSingle).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackVelocityField).cpp
CrackVelocityFieldMulti.o : $(CrackVelocityFieldMulti).cpp $(dprefix) $(CrackVelocityFieldMulti).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(NairnMPM).hpp $(NodalPoint).hpp $(CommonException).hpp $(MeshInfo).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackVelocityFieldMulti).cpp
CrackVelocityFieldSingle.o : $(CrackVelocityFieldSingle).cpp $(dprefix) $(CrackVelocityFieldSingle).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(NairnMPM).hpp $(NodalPoint).hpp $(CommonException).hpp $(MeshInfo).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackVelocityFieldSingle).cpp
	
# MPM: System
ArchiveData.o : $(ArchiveData).cpp $(dprefix) $(NairnMPM).hpp $(ArchiveData).hpp $(MaterialBase).hpp \
			$(CommonArchiveData).hpp $(CommonException).hpp $(GlobalQuantity).hpp $(ElementBase).hpp $(ThermalRamp).hpp \
			$(CrackHeader).hpp $(MPMBase).hpp $(NodalPoint).hpp $(BoundaryCondition).hpp $(MeshInfo).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ArchiveData).cpp

# MPM: NairnMPM_Class
NairnMPM.o : $(NairnMPM).cpp $(dprefix) $(NairnMPM).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp $(ThermalRamp).hpp \
			$(MaterialBase).hpp $(CustomTask).hpp $(CalcJKTask).hpp $(PropagateTask).hpp $(CommonException).hpp \
			$(MPMTask).hpp $(MPMTermination).hpp $(ElementBase).hpp $(InitializationTask).hpp $(MassAndMomentumTask).hpp \
			$(GridForcesTask).hpp $(UpdateMomentaTask).hpp $(UpdateParticlesTask).hpp $(UpdateStrainsLastTask).hpp \
			$(MoveCracksTask).hpp $(ResetElementsTask).hpp $(CommonTask).hpp $(RunCustomTasksTask).hpp \
			$(UpdateStrainsFirstTask).hpp $(MPMWarnings).hpp $(BoundaryCondition).hpp $(MPMBase).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(CustomTask).hpp $(MPMWarnings).hpp $(CrackHeader).hpp \
			$(CrackSurfaceContact).hpp $(MeshInfo).hpp $(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NairnMPM).cpp
StartOutput.o : $(StartOutput).cpp $(dprefix) $(NairnMPM).hpp $(MaterialBase).hpp $(ThermalRamp).hpp $(ArchiveData).hpp \
			$(CommonArchiveData).hpp $(BodyForce).hpp $(CrackSurfaceContact).hpp $(CommonException).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(NodalConcBC).hpp $(BoundaryCondition).hpp \
			$(NodalTempBC).hpp $(NodalVelBC).hpp $(MatPtLoadBC).hpp $(MatPtFluxBC).hpp $(CrackHeader).hpp \
			$(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(StartOutput).cpp
MeshInfo.o : $(MeshInfo).cpp $(dprefix) $(MeshInfo).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MeshInfo).cpp
MPMTask.o : $(MPMTask).cpp $(dprefix) $(MPMTask).hpp $(CommonTask).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMTask).cpp
InitializationTask.o : $(InitializationTask).cpp $(dprefix) $(InitializationTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(NodalPoint).hpp $(MPMWarnings).hpp $(MatPtLoadBC).hpp $(CrackNode).hpp $(ThermalRamp).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(InitializationTask).cpp
MassAndMomentumTask.o : $(MassAndMomentumTask).cpp $(dprefix) $(MassAndMomentumTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(NodalPoint).hpp $(MatPtLoadBC).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp \
			$(RigidMaterial).hpp $(CrackHeader).hpp $(CrackSurfaceContact).hpp $(TransportTask).hpp $(BoundaryCondition).hpp \
			$(NodalVelBC).hpp $(CommonException).hpp $(NodalTempBC).hpp $(NodalConcBC).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MassAndMomentumTask).cpp
UpdateStrainsFirstTask.o : $(UpdateStrainsFirstTask).cpp $(dprefix) $(UpdateStrainsFirstTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateStrainsFirstTask).cpp
GridForcesTask.o : $(GridForcesTask).cpp $(dprefix) $(GridForcesTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(NodalPoint).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp $(CrackHeader).hpp \
			$(TransportTask).hpp $(NodalVelBC).hpp $(ConductionTask).hpp $(BodyForce).hpp $(CrackNode).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GridForcesTask).cpp
UpdateMomentaTask.o : $(UpdateMomentaTask).cpp $(dprefix) $(UpdateMomentaTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(TransportTask).hpp $(NodalPoint).hpp $(CrackNode).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateMomentaTask).cpp
UpdateParticlesTask.o : $(UpdateParticlesTask).cpp $(dprefix) $(UpdateParticlesTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(TransportTask).hpp $(NodalPoint).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp \
			$(BodyForce).hpp $(ThermalRamp).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateParticlesTask).cpp
UpdateStrainsLastTask.o : $(UpdateStrainsLastTask).cpp $(dprefix) $(UpdateStrainsLastTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(CrackHeader).hpp $(TransportTask).hpp \
			$(NodalPoint).hpp $(CrackSurfaceContact).hpp $(NodalVelBC).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateStrainsLastTask).cpp
RunCustomTasksTask.o : $(RunCustomTasksTask).cpp $(dprefix) $(RunCustomTasksTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(NodalPoint).hpp $(CustomTask).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(RunCustomTasksTask).cpp
MoveCracksTask.o : $(MoveCracksTask).cpp $(dprefix) $(MoveCracksTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(CrackHeader).hpp $(MPMTermination).hpp $(CrackSurfaceContact).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoveCracksTask).cpp
ResetElementsTask.o : $(ResetElementsTask).cpp $(dprefix) $(ResetElementsTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(ElementBase).hpp $(MPMTermination).hpp $(MPMWarnings).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ResetElementsTask).cpp

# MPM: Read_MPM
MPMReadHandler.o : $(MPMReadHandler).cpp $(dprefix) $(NairnMPM).hpp $(ArchiveData).hpp $(MPMReadHandler).hpp \
			$(CommonReadHandler).hpp $(CommonArchiveData).hpp $(CustomTask).hpp $(ReverseLoad).hpp $(ThermalRamp).hpp \
			$(BodyForce).hpp $(CrackSurfaceContact).hpp $(GlobalQuantity).hpp $(MatPoint2D).hpp $(MPMBase).hpp $(NodalConcBC).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(NodalTempBC).hpp $(NodalVelBC).hpp $(MatPtLoadBC).hpp $(MatPtFluxBC).hpp \
			$(BoundaryCondition).hpp $(NodesController).hpp $(ElementsController).hpp $(ParseController).hpp $(MaterialController).hpp \
			$(MpsController).hpp $(MaterialController).hpp $(CrackController).hpp $(ElementBase).hpp $(MaterialBase).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(TransportTask).hpp $(MatPoint3D).hpp $(BodyObjectController).hpp \
			$(BodyBoxController).hpp $(BodyOvalController).hpp $(BodyPolygonController).hpp $(BodyRectController).hpp \
			$(BodySphereController).hpp $(BodyCylinderController).hpp $(RigidMaterial).hpp $(VTKArchive).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(MeshInfo).hpp $(PropagateTask).hpp $(BodyPolyhedronController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMReadHandler).cpp
Generators.o : $(Generators).cpp $(dprefix) $(NairnMPM).hpp $(MPMReadHandler).hpp $(CommonReadHandler).hpp \
			$(MPMBase).hpp $(ElementBase).hpp $(MatPoint2D).hpp $(NodalConcBC).hpp $(NodalTempBC).hpp $(NodalVelBC).hpp \
			$(MatPtLoadBC).hpp $(MatPtFluxBC).hpp $(BoundaryCondition).hpp $(ThermalRamp).hpp $(StrX).hpp $(LineController).hpp \
			$(MpsController).hpp $(FourNodeIsoparam).hpp $(ParseController).hpp $(Linear2D).hpp $(NodalPoint).hpp $(BoxController).hpp \
			$(MeshInfo).hpp $(EightNodeIsoparamBrick).hpp $(ElementBase3D).hpp $(MatPoint3D).hpp $(ShapeController).hpp \
			$(RectController).hpp $(ArcController).hpp $(NodalPoint2D).hpp $(NodalPoint3D).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(BodyObjectController).hpp $(ElementsController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Generators).cpp
BitMapFiles.o : $(BitMapFiles).cpp $(dprefix) $(NairnMPM).hpp $(MPMReadHandler).hpp $(CommonReadHandler).hpp \
			$(BMPLevel).hpp $(MatPoint2D).hpp $(ElementBase).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp \
			$(MpsController).hpp $(MPMBase).hpp $(ParseController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFiles).cpp
MpsController.o : $(MpsController).cpp $(dprefix) $(MpsController).hpp $(MPMBase).hpp $(ParseController).hpp \
			 $(DiffusionTask).hpp $(ThermalRamp).hpp $(NairnMPM).hpp $(ResetElementsTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MpsController).cpp
CrackController.o : $(CrackController).cpp $(dprefix) $(CrackController).hpp $(ParseController).hpp $(CrackHeader).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackController).cpp
BodyObjectController.o : $(BodyObjectController).cpp $(dprefix) $(BodyObjectController).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyObjectController).cpp
BodyRectController.o : $(BodyRectController).cpp $(dprefix) $(BodyRectController).hpp $(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyRectController).cpp
BodyOvalController.o : $(BodyOvalController).cpp $(dprefix) $(BodyOvalController).hpp $(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyOvalController).cpp
BodyBoxController.o : $(BodyBoxController).cpp $(dprefix) $(BodyBoxController).hpp $(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyBoxController).cpp
BodySphereController.o : $(BodySphereController).cpp $(dprefix) $(BodySphereController).hpp $(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodySphereController).cpp
BodyCylinderController.o : $(BodyCylinderController).cpp $(dprefix) $(BodyCylinderController).hpp $(BodySphereController).hpp \
			$(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyCylinderController).cpp
BodyPolygonController.o : $(BodyPolygonController).cpp $(dprefix) $(BodyPolygonController).hpp $(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyPolygonController).cpp
BodyPolyhedronController.o : $(BodyPolyhedronController).cpp $(dprefix) $(BodyPolyhedronController).hpp $(BodyObjectController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyPolyhedronController).cpp
PolyTriangle.o : $(PolyTriangle).cpp $(dprefix) $(PolyTriangle).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PolyTriangle).cpp

# MPM: Exceptions
MPMTermination.o : $(MPMTermination).cpp $(dprefix) $(MPMTermination).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMTermination).cpp
MPMWarnings.o : $(MPMWarnings).cpp $(dprefix) $(MPMWarnings).hpp $(NairnMPM).hpp $(ArchiveData).hpp \
			$(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMWarnings).cpp

# MPM: Global Quantities
GlobalQuantity.o : $(GlobalQuantity).cpp $(dprefix) $(GlobalQuantity).hpp $(NairnMPM).hpp $(MaterialBase).hpp \
			$(ThermalRamp).hpp $(MPMBase).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(ArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GlobalQuantity).cpp
ThermalRamp.o : $(ThermalRamp).cpp $(dprefix) $(ThermalRamp).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ThermalRamp).cpp
BodyForce.o : $(BodyForce).cpp $(dprefix) $(BodyForce).hpp $(mathexpr).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyForce).cpp

# MPM: MPM_Classes
MPMBase.o : $(MPMBase).cpp $(dprefix) $(MPMBase).hpp $(CrackHeader).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMBase).cpp
MatPoint2D.o : $(MatPoint2D).cpp $(dprefix) $(MatPoint2D).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPoint2D).cpp
MatPoint3D.o : $(MatPoint3D).cpp $(dprefix) $(MatPoint3D).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPoint3D).cpp

# MPM: Cracks
CrackHeader.o : $(CrackHeader).cpp $(dprefix) $(MaterialBase).hpp $(NairnMPM).hpp $(MPMTermination).hpp $(MeshInfo).hpp \
			$(ArchiveData).hpp $(CommonException).hpp $(CommonArchiveData).hpp $(CrackHeader).hpp \
			$(ElementBase).hpp $(CrackSurfaceContact).hpp $(ContourPoint).hpp $(CrackSegment).hpp \
			$(NodalPoint).hpp $(MPMBase).hpp $(CrackNode).hpp $(NodalPoint2D).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackHeader).cpp
CrackSegment.o : $(CrackSegment).cpp $(dprefix) $(CrackSegment).hpp $(MaterialBase).hpp $(ArchiveData).hpp \
			$(ElementBase).hpp $(NairnMPM).hpp $(CommonArchiveData).hpp $(CrackHeader).hpp $(NodalPoint).hpp \
			$(TractionLaw).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackSegment).cpp
ContourPoint.o : $(ContourPoint).cpp $(dprefix) $(ContourPoint).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ContourPoint).cpp
CrackSurfaceContact.o : $(CrackSurfaceContact).cpp $(dprefix) $(CrackSurfaceContact).hpp $(NodalPoint).hpp $(NairnMPM).hpp\
			$(MPMBase).hpp $(CrackHeader).hpp $(MaterialBase).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackSurfaceContact).cpp
CrackNode.o : $(CrackNode).cpp $(dprefix) $(CrackNode).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackNode).cpp

# MPM: Materials
MaterialBaseMPM.o : $(MaterialBaseMPM).cpp $(dprefix) $(MaterialBase).hpp $(NairnMPM).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(DiffusionTask).hpp $(TransportTask).hpp $(CommonException).hpp \
			$(CrackSurfaceContact).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialBaseMPM).cpp
ElasticMPM.o : $(ElasticMPM).cpp $(dprefix) $(Elastic).hpp $(MaterialBase).hpp $(ThermalRamp).hpp $(ConductionTask).hpp \
			 $(MPMBase).hpp $(DiffusionTask).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElasticMPM).cpp
IsoPlasticity.o : $(IsoPlasticity).cpp $(dprefix) $(IsoPlasticity).hpp $(IsotropicMat).hpp \
			$(MaterialBase).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp \
			$(TransportTask).hpp $(Elastic).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(IsoPlasticity).cpp
AnisoPlasticity.o : $(AnisoPlasticity).cpp $(dprefix) $(AnisoPlasticity).hpp $(Orthotropic).hpp $(TransIsotropic).hpp \
			$(MaterialBase).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp \
			$(TransportTask).hpp $(Elastic).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(AnisoPlasticity).cpp
Viscoelastic.o : $(Viscoelastic).cpp $(dprefix) $(Viscoelastic).hpp $(MaterialBase).hpp $(MPMBase).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(TransportTask).hpp $(CommonException).hpp $(NairnMPM).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Viscoelastic).cpp
HyperElastic.o : $(HyperElastic).cpp $(dprefix) $(HyperElastic).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HyperElastic).cpp
Mooney.o : $(Mooney).cpp $(dprefix) $(Mooney).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(TransportTask).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Mooney).cpp
HEAnisotropic.o : $(HEAnisotropic).cpp $(dprefix) $(HEAnisotropic).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(TransportTask).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HEAnisotropic).cpp
VonMisesHardening.o : $(VonMisesHardening).cpp $(dprefix) $(VonMisesHardening).hpp $(IsoPlasticity).hpp $(IsotropicMat).hpp \
			$(Elastic).hpp $(MaterialBase).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(VonMisesHardening).cpp
HillPlastic.o : $(HillPlastic).cpp $(dprefix) $(HillPlastic).hpp $(AnisoPlasticity).hpp $(Orthotropic).hpp \
			$(TransIsotropic).hpp $(MaterialBase).hpp $(Elastic).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HillPlastic).cpp
WoodMaterial.o : $(WoodMaterial).cpp $(dprefix) $(WoodMaterial).hpp $(HillPlastic).hpp $(AnisoPlasticity).hpp $(Orthotropic).hpp \
			$(TransIsotropic).hpp $(MaterialBase).hpp $(Elastic).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(WoodMaterial).cpp
JohnsonCook.o : $(JohnsonCook).cpp $(dprefix) $(JohnsonCook).hpp $(IsoPlasticity).hpp $(IsotropicMat).hpp \
			$(Elastic).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ThermalRamp).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(JohnsonCook).cpp
MGSCGLMaterial.o : $(MGSCGLMaterial).cpp $(dprefix) $(MGSCGLMaterial).hpp $(IsoPlasticity).hpp $(IsotropicMat).hpp \
			$(Elastic).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ConductionTask).hpp \
			$(DiffusionTask).hpp $(TransportTask).hpp $(ThermalRamp).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MGSCGLMaterial).cpp
SLMaterial.o : $(SLMaterial).cpp $(dprefix) $(SLMaterial).hpp $(MGSCGLMaterial).hpp $(IsoPlasticity).hpp $(IsotropicMat).hpp \
			$(Elastic).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ConductionTask).hpp \
			$(DiffusionTask).hpp $(TransportTask).hpp $(ThermalRamp).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SLMaterial).cpp
BistableIsotropic.o : $(BistableIsotropic).cpp $(dprefix) $(BistableIsotropic).hpp $(IsotropicMat).hpp \
			$(MaterialBase).hpp $(Elastic).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BistableIsotropic).cpp
RigidMaterial.o : $(RigidMaterial).cpp $(dprefix) $(RigidMaterial).hpp $(MaterialBase).hpp $(CommonException).hpp \
			$(mathexpr).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(RigidMaterial).cpp
TractionLaw.o : $(TractionLaw).cpp $(dprefix) $(TractionLaw).hpp $(MaterialBase).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TractionLaw).cpp
CohesiveZone.o : $(CohesiveZone).cpp $(dprefix) $(CohesiveZone).hpp $(TractionLaw).hpp $(MaterialBase).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CohesiveZone).cpp
TrilinearTraction.o : $(TrilinearTraction).cpp $(dprefix) $(TrilinearTraction).hpp $(CohesiveZone).hpp $(TractionLaw).hpp \
			 $(MaterialBase).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TrilinearTraction).cpp
LinearTraction.o : $(LinearTraction).cpp $(dprefix) $(LinearTraction).hpp $(CohesiveZone).hpp $(TractionLaw).hpp \
			$(MaterialBase).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinearTraction).cpp
CubicTraction.o : $(CubicTraction).cpp $(dprefix) $(CubicTraction).hpp $(TractionLaw).hpp $(MaterialBase).hpp \
			$(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CubicTraction).cpp

# MPM: Boundary_Conditions
BoundaryCondition.o : $(BoundaryCondition).cpp $(dprefix) $(BoundaryCondition).hpp \
			$(mathexpr).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BoundaryCondition).cpp
MatPtLoadBC.o : $(MatPtLoadBC).cpp $(dprefix) $(MatPtLoadBC).hpp $(BoundaryCondition).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(MeshInfo).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPtLoadBC).cpp
MatPtFluxBC.o : $(MatPtFluxBC).cpp $(dprefix) $(MatPtFluxBC).hpp $(MatPtLoadBC).hpp $(BoundaryCondition).hpp \
			$(MPMBase).hpp $(MaterialBase).hpp $(MeshInfo).hpp $(mathexpr).hpp $(NairnMPM).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPtFluxBC).cpp
NodalVelBC.o : $(NodalVelBC).cpp $(dprefix) $(NodalVelBC).hpp $(BoundaryCondition).hpp $(NodalPoint).hpp \
			$(NairnMPM).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalVelBC).cpp
NodalConcBC.o : $(NodalConcBC).cpp $(dprefix) $(NodalConcBC).hpp $(BoundaryCondition).hpp $(NodalPoint).hpp \
			 $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalConcBC).cpp 
NodalTempBC.o : $(NodalTempBC).cpp $(dprefix) $(NodalTempBC).hpp $(BoundaryCondition).hpp $(NodalPoint).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalTempBC).cpp

# MPM: Custom_Tasks
CustomTask.o : $(CustomTask).cpp $(dprefix) $(CustomTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CustomTask).cpp
CalcJKTask.o : $(CalcJKTask).cpp $(dprefix) $(CalcJKTask).hpp $(MaterialBase).hpp $(CustomTask).hpp \
			$(ArchiveData).cpp $(CommonArchiveData).hpp $(NodalPoint).hpp $(NairnMPM).hpp $(MPMBase).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CalcJKTask).cpp
PropagateTask.o : $(PropagateTask).cpp $(dprefix) $(PropagateTask).hpp $(CalcJKTask).hpp $(NairnMPM).hpp \
			$(MaterialBase).hpp $(ElementBase).hpp $(CustomTask).hpp $(MPMBase).hpp $(ConductionTask).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PropagateTask).cpp
ReverseLoad.o : $(ReverseLoad).cpp $(dprefix) $(ReverseLoad).hpp $(PropagateTask).hpp $(NairnMPM).hpp \
			$(MatPtLoadBC).hpp $(MPMTermination).hpp $(MaterialBase).hpp $(CustomTask).hpp $(BoundaryCondition).hpp \
			$(CommonException).hpp $(CrackHeader).hpp $(MatPtLoadBC).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ReverseLoad).cpp
TransportTask.o : $(TransportTask).cpp $(dprefix) $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TransportTask).cpp
DiffusionTask.o : $(DiffusionTask).cpp $(dprefix) $(DiffusionTask).hpp $(NairnMPM).hpp $(ElementBase).hpp \
			$(MaterialBase).hpp $(NodalConcBC).hpp $(MatPtFluxBC).hpp $(BoundaryCondition).hpp $(MatPtLoadBC).hpp \
			$(MPMBase).hpp $(NodalPoint).hpp $(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(DiffusionTask).cpp
ConductionTask.o : $(ConductionTask).cpp $(dprefix) $(NairnMPM).hpp $(ElementBase).hpp $(MaterialBase).hpp \
			$(NodalTempBC).hpp $(ConductionTask).hpp $(BoundaryCondition).hpp $(MPMBase).hpp \
			$(CrackHeader).hpp $(NodalPoint).hpp $(CrackSegment).hpp $(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ConductionTask).cpp
VTKArchive.o : $(VTKArchive).cpp $(dprefix) $(VTKArchive).hpp $(CustomTask).hpp $(NairnMPM).hpp $(NodalPoint).hpp \
			$(CommonException).hpp $(ArchiveData).hpp $(MeshInfo).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(VTKArchive).cpp

# MPM: Elements
MoreMPMElementBase.o : $(MoreMPMElementBase).cpp $(dprefix) $(ElementBase).hpp $(NodalPoint).hpp \
			$(MeshInfo).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoreMPMElementBase).cpp
ElementBase3D.o : $(ElementBase3D).cpp $(dprefix) $(ElementBase).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementBase3D).cpp
EightNodeIsoparamBrick.o : $(EightNodeIsoparamBrick).cpp $(dprefix) $(ElementBase3D).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EightNodeIsoparamBrick).cpp
	
# MPM: Nodes
NodalPointMPM.o : $(NodalPointMPM).cpp $(dprefix) $(NodalPoint).hpp $(NairnMPM).hpp $(ArchiveData).hpp \
			$(CommonArchiveData).hpp $(CrackSurfaceContact).hpp $(CommonException).hpp $(MPMWarnings).hpp \
			$(CrackNode).hpp $(MeshInfo).hpp $(CrackSegment).hpp $(CrackHeader).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(BoundaryCondition).hpp $(CrackVelocityFieldMulti).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPointMPM).cpp

# -------------------------------------------------------------------------
# To make executable
all : NairnMPM
        
# -------------------------------------------------------------------------
# To clean compiled objects        
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
