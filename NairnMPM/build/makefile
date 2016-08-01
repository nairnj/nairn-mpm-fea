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
#		cygwinbin - Standard cygwin install using binary xerces
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
#    $(dtdpath) is path to save DTD file
# 6. $(ioutput) is path to install executable
# 7. $(src) and $(com) are paths to parts of the source code
#
# Each of these can changed by editing or can be overriden at make time with syntax such as
# and are documented more in the numbered section below
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
ifeq ($(SYSTEM),cygwinbin)
	xercesSo = -L/lib -lxerces-c
endif
ifeq ($(SYSTEM),inferno)
	xercesSo = /usr/sci/projects/SCIRun/Thirdparty/1.22/Linux/gcc-3.2-32bit/lib/libxerces-c.so
endif
ifeq ($(SYSTEM),zia)
	xercesSo = /usr/lib/libxerces-c.so
endif

# 2. Path to xerces headers (if not at default /usr/local/include) 
headersPath= /usr/local/include
ifeq ($(SYSTEM),carbon)
	headersPath = /usr/local/xerces3/include
endif
ifeq ($(SYSTEM),mpm)
	headersPath = /usr/include
endif
ifeq ($(SYSTEM),cygwinbin)
	headersPath= /usr/include
endif
ifeq ($(SYSTEM),inferno)
	headersPath = /usr/sci/projects/SCIRun/Thirdparty/1.22/Linux/gcc-3.2-32bit/include
endif
ifeq ($(SYSTEM),zia)
	headersPath = /usr/include/xercesc
endif

# 3. Define compiler to use and other commands
#     CC is path to the c++ compiler, must support OpenMP is used
#     RM is path/options to the file remove command (used only in make clean)
#     CP is path/options to copy command (used only in make install
CC= g++
ifeq ($(SYSTEM),carbon)
    CC= /usr/local/gcc48/bin/g++
endif
ifeq ($(SYSTEM),mac)
    CC= /opt/local/bin/g++-mp-4.9
endif
RM= rm
CP= cp

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
#     If using -fopenmp, -arch, or -pg, must have in both CFLAGS and LFLAGS. For MacOS X, the specified arch
#       must match the xerces library used for linking
CFLAGS= -c -O3 -fopenmp
LFLAGS= -fopenmp

# 5. Define executable destination
#     output = relative or full path to save executable or to find it for a subsequent make install
#     dtdpath = relative of full path to DTD file, only used during make install
output = ../input/NairnMPM
dtdpath = ../input/NairnMPM.dtd

# 6. Define paths to all files (relative to 'makefile')
#     These are usually not changed unless you reorganize the code
ioutput = ~/bin

# 7. Define paths to all files (relative to 'makefile')
#     The commented com option is needed if you follow the suggested Eclipse procedure
src = ../src
com = ../../Common
#com = ../Common

headers = -I$(src) -I$(com) -I$(headersPath)

# NairnMPM file abbreviations
AdjustTimeStepTask = $(src)/Custom_Tasks/AdjustTimeStepTask
AnisoPlasticity = $(src)/Materials/AnisoPlasticity
ArcController = $(com)/Read_XML/ArcController
ArchiveData = $(src)/System/ArchiveData
BitMapFiles = $(src)/Read_MPM/BitMapFiles
BitMapFilesCommon = $(com)/Read_XML/BitMapFilesCommon
BistableIsotropic = $(src)/Materials/BistableIsotropic
BMPLevel = $(com)/Read_XML/BMPLevel
BodyForce = $(src)/Global_Quantities/BodyForce
BoundaryCondition = $(src)/Boundary_Conditions/BoundaryCondition
BoxController = $(com)/Read_XML/BoxController
CalcJKTask = $(src)/Custom_Tasks/CalcJKTask
CarnotCycle = $(src)/Custom_Tasks/CarnotCycle
ClampedNeohookean = $(src)/Materials/ClampedNeohookean
CohesiveZone = $(src)/Materials/CohesiveZone
CommonAnalysis = $(com)/System/CommonAnalysis
CommonArchiveData = $(com)/System/CommonArchiveData
CommonException = $(com)/Exceptions/CommonException
CommonReadHandler = $(com)/Read_XML/CommonReadHandler
CommonTask = $(com)/System/CommonTask
CommonUtilities = $(com)/System/CommonUtilities
ConductionTask = $(src)/Custom_Tasks/ConductionTask
ContactLaw = $(src)/Materials/ContactLaw
ContourPoint = $(src)/Cracks/ContourPoint
CoulombFriction = $(src)/Materials/CoulombFriction
CoupledSawTooth = $(src)/Materials/CoupledSawTooth
CrackController = $(src)/Read_MPM/CrackController
CrackHeader = $(src)/Cracks/CrackHeader
CrackLeaf = $(src)/Cracks/CrackLeaf
CrackNode = $(src)/Cracks/CrackNode
CrackSegment = $(src)/Cracks/CrackSegment
CrackSurfaceContact = $(src)/Cracks/CrackSurfaceContact
CrackVelocityField = $(src)/Nodes/CrackVelocityField
CrackVelocityFieldMulti = $(src)/Nodes/CrackVelocityFieldMulti
CrackVelocityFieldSingle = $(src)/Nodes/CrackVelocityFieldSingle
CrossedCrack = $(src)/Cracks/CrossedCrack
CubicTraction = $(src)/Materials/CubicTraction
CustomTask = $(src)/Custom_Tasks/CustomTask
DDBHardening = $(src)/Materials/DDBHardening
DiffusionTask = $(src)/Custom_Tasks/DiffusionTask
EightNodeIsoparamBrick = $(src)/Elements/EightNodeIsoparamBrick
Elastic = $(com)/Materials/Elastic
ElasticMPM = $(src)/Materials/ElasticMPM
ElementBase = $(com)/Elements/ElementBase
ElementBase3D = $(src)/Elements/ElementBase3D
ElementsController = $(com)/Read_XML/ElementsController
ExtrapolateRigidBCsTask = $(src)/NairnMPM_Class/ExtrapolateRigidBCsTask
FourNodeIsoparam = $(com)/Elements/FourNodeIsoparam
Generators = $(src)/Read_MPM/Generators
GhostNode = $(src)/Patches/GhostNode
GlobalQuantity = $(src)/Global_Quantities/GlobalQuantity
GridArchive = $(src)/Custom_Tasks/GridArchive
GridForcesTask = $(src)/NairnMPM_Class/GridForcesTask
GridPatch = $(src)/Patches/GridPatch
HardeningLawBase = $(src)/Materials/HardeningLawBase
HEIsotropic = $(src)/Materials/HEIsotropic
HEMGEOSMaterial = $(src)/Materials/HEMGEOSMaterial
HillPlastic = $(src)/Materials/HillPlastic
HistoryArchive = $(src)/Custom_Tasks/HistoryArchive
HyperElastic = $(src)/Materials/HyperElastic
IdealGas = $(src)/Materials/IdealGas
InitializationTask = $(src)/NairnMPM_Class/InitializationTask
InitVelocityFieldsTask = $(src)/NairnMPM_Class/InitVelocityFieldsTask
IsoPlasticity = $(src)/Materials/IsoPlasticity
IsotropicMat = $(com)/Materials/IsotropicMat
JohnsonCook = $(src)/Materials/JohnsonCook
Lagrange2D = $(com)/Elements/Lagrange2D
Linear2D = $(com)/Elements/Linear2D
LinearHardening = $(src)/Materials/LinearHardening
LinearInterface = $(src)/Materials/LinearInterface
LinearTraction = $(src)/Materials/LinearTraction
LineController = $(com)/Read_XML/LineController
LinkedObject =$(com)/System/LinkedObject
main = $(com)/System/main
MassAndMomentumTask = $(src)/NairnMPM_Class/MassAndMomentumTask
MaterialBase = $(com)/Materials/MaterialBase
MaterialBaseMPM = $(src)/Materials/MaterialBaseMPM
MaterialController = $(com)/Read_XML/MaterialController
MaterialInterfaceNode = $(src)/Nodes/MaterialInterfaceNode
mathexpr = $(com)/Read_XML/mathexpr
MatPointAS = $(src)/MPM_Classes/MatPointAS
MatPoint2D = $(src)/MPM_Classes/MatPoint2D
MatPoint3D = $(src)/MPM_Classes/MatPoint3D
MatPtFluxBC = $(src)/Boundary_Conditions/MatPtFluxBC
MatPtHeatFluxBC = $(src)/Boundary_Conditions/MatPtHeatFluxBC
MatPtLoadBC = $(src)/Boundary_Conditions/MatPtLoadBC
MatPtTractionBC = $(src)/Boundary_Conditions/MatPtTractionBC
Matrix3 = $(com)/System/Matrix3
MatVelocityField = $(src)/Nodes/MatVelocityField
MeshInfo = $(src)/NairnMPM_Class/MeshInfo
Mooney = $(src)/Materials/Mooney
MoreIsotropicMat = $(src)/Materials/MoreIsotropicMat
MoreMPMElementBase = $(src)/Elements/MoreMPMElementBase
MoveCracksTask = $(src)/NairnMPM_Class/MoveCracksTask
MPMBase = $(src)/MPM_Classes/MPMBase
MPMReadHandler = $(src)/Read_MPM/MPMReadHandler
MPMTask = $(src)/NairnMPM_Class/MPMTask
MPMWarnings = $(src)/Exceptions/MPMWarnings
MpsController = $(src)/Read_MPM/MpsController
NairnMPM = $(src)/NairnMPM_Class/NairnMPM
Neohookean = $(src)/Materials/Neohookean
NodalConcBC = $(src)/Boundary_Conditions/NodalConcBC
NodalPoint = $(com)/Nodes/NodalPoint
NodalPoint2D = $(com)/Nodes/NodalPoint2D
NodalPoint3D = $(com)/Nodes/NodalPoint3D
NodalPointMPM = $(src)/Nodes/NodalPointMPM
NodalTempBC = $(src)/Boundary_Conditions/NodalTempBC
NodalVelBC = $(src)/Boundary_Conditions/NodalVelBC
NodesController = $(com)/Read_XML/NodesController
NonlinearHardening = $(src)/Materials/NonlinearHardening
Nonlinear2Hardening = $(src)/Materials/Nonlinear2Hardening
Orthotropic = $(com)/Materials/Orthotropic
OvalController = $(com)/Read_XML/OvalController
ParseController = $(com)/Read_XML/ParseController
PointController = $(com)/Read_XML/PointController
PolygonController = $(com)/Read_XML/PolygonController
PolyhedronController = $(src)/Read_MPM/PolyhedronController
PolyTriangle = $(src)/Read_MPM/PolyTriangle
PostExtrapolationTask = $(src)/NairnMPM_Class/PostExtrapolationTask
PostForcesTask = $(src)/NairnMPM_Class/PostForcesTask
PressureLaw = $(src)/Materials/PressureLaw
ProjectRigidBCsTask = $(src)/NairnMPM_Class/ProjectRigidBCsTask
PropagateTask = $(src)/Custom_Tasks/PropagateTask
Quad2D = $(com)/Elements/Quad2D
RectController = $(com)/Read_XML/RectController
ResetElementsTask = $(src)/NairnMPM_Class/ResetElementsTask
ReverseLoad = $(src)/Custom_Tasks/ReverseLoad
RigidMaterial = $(src)/Materials/RigidMaterial
RunCustomTasksTask = $(src)/NairnMPM_Class/RunCustomTasksTask
SCGLHardening = $(src)/Materials/SCGLHardening
SetRigidContactVelTask = $(src)/NairnMPM_Class/SetRigidContactVelTask
ShapeController = $(com)/Read_XML/ShapeController
ShellController = $(src)/Read_MPM/ShellController
SLMaterial = $(src)/Materials/SLMaterial
SphereController = $(src)/Read_MPM/SphereController
StartOutput = $(src)/NairnMPM_Class/StartOutput
StrX = $(com)/Exceptions/StrX
svninfo = $(com)/System/svninfo
TaitLiquid = $(src)/Materials/TaitLiquid
ThermalRamp = $(src)/Global_Quantities/ThermalRamp
TorusController = $(src)/Read_MPM/TorusController
TractionLaw = $(src)/Materials/TractionLaw
TransIsotropic = $(com)/Materials/TransIsotropic
TransportTask = $(src)/Custom_Tasks/TransportTask
TrilinearTraction = $(src)/Materials/TrilinearTraction
UnitsController = $(com)/System/UnitsController
UpdateMomentaTask = $(src)/NairnMPM_Class/UpdateMomentaTask
UpdateParticlesTask = $(src)/NairnMPM_Class/UpdateParticlesTask
UpdateStrainsFirstTask = $(src)/NairnMPM_Class/UpdateStrainsFirstTask
UpdateStrainsLastTask = $(src)/NairnMPM_Class/UpdateStrainsLastTask
UpdateStrainsLastContactTask = $(src)/NairnMPM_Class/UpdateStrainsLastContactTask
Viscoelastic = $(src)/Materials/Viscoelastic
VTKArchive = $(src)/Custom_Tasks/VTKArchive
WoodMaterial = $(src)/Materials/WoodMaterial

# path from 'makefile' to prefix header file
prefix = $(src)/System/MPMPrefix.hpp

# path from 'makefile' to default headers (this included by $(prefix)
dprefix = $(src)/System/MPMPrefix.hpp $(com)/System/DataTypes.hpp $(com)/System/Matrix3.hpp\
		$(com)/System/CommonAnalysis.hpp $(com)/System/LinkedObject.hpp $(com)/System/UnitsController.hpp

# NairnMPM header files that include other header files (other than $(dprefix) files)
# AdjustTimeStepTask.hpp : CustomTask.hpp
# AnisoPlasticity.hpp : Orthotropic.hpp (TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)
# ArcController.hpp : LineController.hpp ShapeController.hpp
# ArchiveData.hpp : CommonArchiveData.hpp
# BistableIsotropic.hpp : IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
# PolygonController.hpp : ShapeController.hpp
# PolyhedronController.hpp : ShapeController.hpp
# BoxController.hpp : ShapeController.hpp
# CalcJKTask.hpp : CustomTask.hpp
# CarnotCycle.hpp : CustomTask.hpp
# ClampedNeohookean.hpp : Neohookean.hpp HyperElastic.hpp (MaterialBase.hpp)
# CohesiveZone.hpp : MaterialBase.hpp TractionLaw.hpp
# ConductionTask.hpp : TransportTask.hpp
# ContactLaw.hpp : MaterialBase.hpp
# CoulombFriction.hpp : ContactLaw.hpp MaterialBase.hpp
# CoupledSawTooth.hpp : CohesiveZone.hpp MaterialBase.hpp TractionLaw.hpp
# CrackController.hpp : ParseController.hpp
# CrackNode.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# CrackSegment.hpp : CrackHeader.hpp
# CrackSurfaceContact.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# CrackVelocityField.hpp : MatVelocityField.hpp
# CrackVelocityFieldMulti.hpp : CrackVelocityField.hpp MatVelocityField.hpp
# CrackVelocityFieldSingle.hpp : CrackVelocityField.hpp MatVelocityField.hpp
# CubicTraction.hpp : MaterialBase.hpp TractionLaw.hpp
# DDDHardening.hpp: HardeningLawBase.hpp, MaterialBase.hpp
# DiffusionTask.hpp : TransportTask.hpp
# EightNodeIsoparamBrick.hpp : ElementBase3D.hpp (ElementBase.hpp)
# Elastic.hpp : MaterialBase.hpp
# ElementsController.hpp : ParseController.hpp
# ElementBase3D : ElementBase
# ExtrapolateRigidBCsTask.hpp : MPMTask.hpp (CommonTask.hpp)
# FourNodeIsoparam.hpp : Linear2D.hpp (ElementBase.hpp)
# GridArchive.hpp : CustomTask.hpp
# GridForcesTask.hpp : MPMTask.hpp (CommonTask.hpp)
# HardeningLawBase.hpp: MaterialBase.hpp
# HEIsotropic.hpp : HyperElastic.hpp (MaterialBase.hpp)
# HEMGEOSMaterial.hpp : HEIsotropic.hpp (HyperElastic.hpp, MaterialBase.hpp)
# HillPlastic : AnisoPlasticity.hpp (Orthotropic.hpp TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)
# HistoryArchive.hpp : CustomTask.hpp
# IdealGas.hpp : HyperElastic.hpp (MaterialBase.hpp)
# InitializationTask : MPMTask.hpp (CommonTask.hpp)
# InitVelocityFieldsTask : MPMTask.hpp (CommonTask.hpp)
# IsoPlasticity : IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
# IsotropicMat.hpp : Elastic.hpp (MaterialBase.hpp)
# JohnsonCook.hpp : HardeningLawBase.hpp, MaterialBase.hpp
# Lagrange2D.hpp : Quad2D.hpp (ElementBase.hpp)
# Linear2D.hpp : ElementBase.hpp
# LinearHardening.hpp: HardeningLawBase.hpp, MaterialBase.hpp
# LinearInterface.hpp : MaterialBase.hpp
# LinearTraction.hpp : CohesiveZone.hpp MaterialBase.hpp TractionLaw.hpp
# LineController.hpp : ShapeController.hpp
# MassAndMomentumTask.hpp : MPMTask.hpp (CommonTask.hpp)
# MaterialController.hpp : ParseController.hpp
# MaterialIntefaceNode.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# MatPtLoadBC.hpp : BoundaryCondition.hpp
# MatPtTractionBC.hpp : MatPtLoadBC.hpp BoundaryCondition.hpp
# MatPtFluxBC.hpp : BoundaryCondition.hpp
# MatPtHeatFluxBC.hpp : BoundaryCondition.hpp
# MatPointAS.hpp : MatPoint2D.hpp MPMBase.hpp
# MatPoint2D.hpp : MPMBase.hpp
# MatPoint3D.hpp : MPMBase.hpp
# Mooney.hpp : HyperElastic.hpp (MaterialBase.hpp)
# MoveCracksTask.hpp : MPMTask.hpp (CommonTask.hpp)
# MPMReadHandler.hpp : CommonReadHandler.hpp
# MPMTask.hpp : CommonTask.hpp
# MpsController.hpp : ParseController.hpp
# Neohookean.hpp : HyperElastic.hpp (MaterialBase.hpp)
# NodalConcBC.hpp : BoundaryCondition.hpp
# NodalPoint.hpp : CrackVelocityField.hpp MatVelocityField.hpp
# NodalPoint2D.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# NodalPoint3D.hpp : NodalPoint.hpp CrackVelocityField.hpp MatVelocityField.hpp
# NodalTempBC.hpp : BoundaryCondition.hpp
# NodalVelBC.hpp : BoundaryCondition.hpp
# NodesController.hpp : ParseController.hpp
# NonlinearHardening.hpp: HardeningLawBase.hpp, MaterialBase.hpp
# Nonlinear2Hardening.hpp: NonlinearHardening.hpp, HardeningLawBase.hpp, MaterialBase.hpp
# Orthotropic.hpp : TransIsotropic.hpp (Elastic.hpp MaterialBase.hpp)
# OvalController.hpp : ShapeController.hpp
# PointController.hpp : ShapeController.hpp
# PolygonController.hpp : ShapeController.hpp
# PostExtrapolationTask.hpp : MPMTask.hpp (CommonTask.hpp)
# PostForcesTask.hpp : MPMTask.hpp (CommonTask.hpp)
# PresureLaw.hpp : MaterialBase.hpp TractionLaw.hpp
# ProjectRigidBCsTask.hpp : MPMTask.hpp (CommonTask.hpp)
# PropagateTask.hpp : CustomTask.hpp
# Quad2D.hpp : ElementBase.hpp
# RectController.hpp : ShapeController.hpp
# ResetElementsTask.hpp : MPMTask.hpp (CommonTask.hpp)
# ReverseLoad.hpp : CustomTask.hpp
# HyperElastic.hpp : MaterialBase.hpp
# RunCustomTasksTask.hpp : MPMTask.hpp (CommonTask.hpp)
# RigidMaterial.hpp : MaterialBase.hpp
# SCGLHardening.hpp : HardeningLawBase.hpp, MaterialBase.hpp
# SetRigidContactVelTask.hpp : MPMTask.hpp (CommonTask.hpp)
# ShellController.hpp : ShapeController.hpp
# SLMaterial.hpp :  SCGLHardening.hpp, HardeningLawBase.hpp, MaterialBase.hpp
# SphereController.hpp : ShapeController.hpp
# TaitLiquid.hpp : HyperElastic.hpp (MaterialBase.hpp)
# TorusController.hpp : ShapeController.hpp
# TractionLaw.hpp : MaterialBase.hpp
# TransIsotropic.hpp : Elastic.hpp (MaterialBase.hpp)
# TrilinearTraction : MaterialBase.hpp TractionLaw.hpp CohesiveZone.hpp
# UpdateMomentaTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateParticlesTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateStrainsFirstTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateStrainsLastTask.hpp : MPMTask.hpp (CommonTask.hpp)
# UpdateStrainsLastContactTask.hpp : MPMTask.hpp (CommonTask.hpp)
# Viscoelastic.hpp : MaterialBase.hpp
# WoodMaterial : HillPlastic.hpp (AnisoPlasticity.hpp, Orthotropic.hpp TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)
# VTKArchive.hpp : CustomTask.hpp, GridArchive.hpp

# NairnMPM header files with no dependencies
# BMPLevel.hpp
# BodyForce.hpp
# BoundaryCondition.hpp
# CommonArchiveData.hpp
# CommonException.hpp
# CommonReadHandler.hpp
# ContourPoint.hpp
# CrackHeader.hpp
# CrackLeaf.hpp
# CrossedCrack.hpp
# CustomTask.hpp
# ElementBase.hpp
# GlobalQuantity.hpp
# GridPatch.hpp
# MaterialBase.hpp
# MatVelocityField.hpp
# MPMBase.hpp
# MPMWarnings.hpp
# NairnMPM.hpp
# ParseController.hpp
# ShapeController.hpp
# StrX.hpp
# ThermalRamp.hpp
# TransportTask.hpp
# PolyTriangle.hpp

# -------------------------------------------------------------------------
# all compiled objects
objects = main.o CommonAnalysis.o CommonArchiveData.o LinkedObject.o CommonUtilities.o CommonTask.o svninfo.o \
		CommonException.o ElementBase.o Linear2D.o FourNodeIsoparam.o MaterialBase.o Matrix3.o \
		CommonReadHandler.o ParseController.o ElementsController.o NodesController.o TorusController.o \
		LineController.o MaterialController.o IsotropicMat.o TransIsotropic.o BoxController.o \
		Orthotropic.o NodalPoint.o NodalPoint2D.o NodalPoint3D.o ArchiveData.o NairnMPM.o MPMTask.o \
		InitializationTask.o MassAndMomentumTask.o UpdateStrainsFirstTask.o GridForcesTask.o UpdateMomentaTask.o SetRigidContactVelTask.o \
		UpdateParticlesTask.o UpdateStrainsLastTask.o UpdateStrainsLastContactTask.o RunCustomTasksTask.o MoveCracksTask.o ResetElementsTask.o \
		StartOutput.o MPMReadHandler.o Generators.o BitMapFiles.o BitMapFilesCommon.o BMPLevel.o MpsController.o \
		CrackController.o MPMWarnings.o BodyForce.o GlobalQuantity.o BoundaryCondition.o DDBHardening.o CrossedCrack.o \
		ThermalRamp.o MPMBase.o MatPoint2D.o MatPointAS.o CrackHeader.o CrackSegment.o CrackLeaf.o ContourPoint.o \
		CrackSurfaceContact.o CrackNode.o MaterialInterfaceNode.o MaterialBaseMPM.o Viscoelastic.o JohnsonCook.o \
		Mooney.o HEIsotropic.o BistableIsotropic.o RigidMaterial.o HardeningLawBase.o LinearHardening.o \
		MatPtLoadBC.o MatPtTractionBC.o MatPtFluxBC.o NodalVelBC.o NodalConcBC.o NodalTempBC.o MoreMPMElementBase.o \
		NodalPointMPM.o CustomTask.o CalcJKTask.o PropagateTask.o ReverseLoad.o IdealGas.o NonlinearHardening.o \
		DiffusionTask.o ConductionTask.o MeshInfo.o TransportTask.o ElementBase3D.o SCGLHardening.o TaitLiquid.o \
		EightNodeIsoparamBrick.o MatPoint3D.o Elastic.o ElasticMPM.o HyperElastic.o Nonlinear2Hardening.o \
		mathexpr.o ShapeController.o PointController.o RectController.o OvalController.o PolygonController.o \
		SphereController.o ArcController.o MatVelocityField.o CrackVelocityField.o HistoryArchive.o CarnotCycle.o \
		TractionLaw.o CohesiveZone.o CoupledSawTooth.o LinearTraction.o CubicTraction.o IsoPlasticity.o AnisoPlasticity.o \
		HillPlastic.o VTKArchive.o AdjustTimeStepTask.o SLMaterial.o CrackVelocityFieldMulti.o LinearInterface.o \
		CrackVelocityFieldSingle.o WoodMaterial.o TrilinearTraction.o PolyhedronController.o PolyTriangle.o Quad2D.o Lagrange2D.o \
        HEMGEOSMaterial.o GridPatch.o GhostNode.o MatPtHeatFluxBC.o PressureLaw.o ShellController.o UnitsController.o  \
		Neohookean.o ClampedNeohookean.o GridArchive.o InitVelocityFieldsTask.o MoreIsotropicMat.o PostForcesTask.o \
		CoulombFriction.o ContactLaw.o PostExtrapolationTask.o ProjectRigidBCsTask.o ExtrapolateRigidBCsTask.o

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
			$(CommonArchiveData).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(MPMBase).hpp
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
Matrix3.o : $(Matrix3).cpp $(dprefix)
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Matrix3).cpp
svninfo.o : $(svninfo).cpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(svninfo).cpp
UnitsController.o : $(UnitsController).cpp $(UnitsController).hpp $(CommonReadHandler).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UnitsController).cpp

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
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonReadHandler).hpp $(ElementBase).hpp
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
		$(IsotropicMat).hpp $(TransIsotropic).hpp $(Orthotropic).hpp $(Viscoelastic).hpp $(PressureLaw).hpp \
		$(Mooney).hpp $(BistableIsotropic).hpp $(RigidMaterial).hpp $(MaterialBase).hpp \
		$(Elastic).hpp $(HyperElastic).hpp $(NairnMPM).hpp $(TractionLaw).hpp $(CohesiveZone).hpp $(LinearTraction).hpp \
		$(CubicTraction).hpp $(AnisoPlasticity).hpp $(IsoPlasticity.hpp) $(HillPlastic).hpp \
		$(WoodMaterial).hpp $(TrilinearTraction).hpp $(IdealGas).hpp $(CoupledSawTooth).hpp $(TaitLiquid).hpp \
		$(HEIsotropic).hpp $(HEMGEOSMaterial).hpp $(Neohookean).hpp $(ClampedNeohookean).hpp $(ContactLaw).hpp \
		$(CoulombFriction).hpp $(LinearInterface).hpp
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
TransIsotropic.o : $(TransIsotropic).cpp $(dprefix) $(TransIsotropic).hpp $(MaterialBase).hpp $(CommonException).hpp \
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
			$(NodalPoint).hpp $(MeshInfo).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(NairnMPM).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(FourNodeIsoparam).cpp
Quad2D.o : $(Quad2D).cpp $(dprefix) $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Quad2D).cpp
Lagrange2D.o : $(Lagrange2D).cpp $(dprefix) $(Lagrange2D).hpp $(Quad2D).hpp $(ElementBase).hpp $(NodalPoint).hpp \
			$(MeshInfo).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Lagrange2D).cpp

# Common: Nodes
NodalPoint.o : $(NodalPoint).cpp $(dprefix) $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(MaterialInterfaceNode).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint).cpp
NodalPoint2D.o : $(NodalPoint2D).cpp $(dprefix) $(NodalPoint2D).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint2D).cpp
NodalPoint3D.o : $(NodalPoint3D).cpp $(dprefix) $(NodalPoint3D).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPoint3D).cpp
MatVelocityField.o : $(MatVelocityField).cpp $(dprefix) $(MatVelocityField).hpp $(NairnMPM).hpp $(CommonException).hpp \
			$(NodalPoint).hpp $(CrackHeader).hpp $(CrackVelocityField).hpp $(TransportTask).hpp $(ConductionTask).hpp $(BodyForce).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatVelocityField).cpp
CrackVelocityField.o : $(CrackVelocityField).cpp $(dprefix) $(CrackVelocityField).hpp $(MatVelocityField).hpp $(NairnMPM).hpp \
			$(NodalPoint).hpp $(CommonException).hpp $(MeshInfo).hpp $(BoundaryCondition).hpp $(CrackVelocityFieldMulti).hpp \
			$(CrackVelocityFieldSingle).hpp $(MaterialBase).hpp $(ConductionTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackVelocityField).cpp
CrackVelocityFieldMulti.o : $(CrackVelocityFieldMulti).cpp $(dprefix) $(CrackVelocityFieldMulti).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(NairnMPM).hpp $(NodalPoint).hpp $(CommonException).hpp $(MeshInfo).hpp $(BoundaryCondition).hpp \
			$(MaterialInterfaceNode).hpp $(MPMBase).hpp $(ConductionTask).hpp $(MaterialBase).hpp $(ContactLaw).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackVelocityFieldMulti).cpp
CrackVelocityFieldSingle.o : $(CrackVelocityFieldSingle).cpp $(dprefix) $(CrackVelocityFieldSingle).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(NairnMPM).hpp $(NodalPoint).hpp $(CommonException).hpp $(MeshInfo).hpp $(BoundaryCondition).hpp \
            $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackVelocityFieldSingle).cpp
MaterialInterfaceNode.o : $(MaterialInterfaceNode).cpp $(dprefix) $(MaterialInterfaceNode).hpp $(NodalPoint).hpp \
            $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialInterfaceNode).cpp
	
# MPM: System
ArchiveData.o : $(ArchiveData).cpp $(dprefix) $(NairnMPM).hpp $(ArchiveData).hpp $(MaterialBase).hpp \
			$(CommonArchiveData).hpp $(CommonException).hpp $(GlobalQuantity).hpp $(ElementBase).hpp $(ThermalRamp).hpp \
			$(CrackHeader).hpp $(MPMBase).hpp $(NodalPoint).hpp $(BoundaryCondition).hpp $(MeshInfo).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ArchiveData).cpp

# MPM: NairnMPM_Class
NairnMPM.o : $(NairnMPM).cpp $(dprefix) $(NairnMPM).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp $(ThermalRamp).hpp \
			$(MaterialBase).hpp $(CustomTask).hpp $(CalcJKTask).hpp $(PropagateTask).hpp $(CommonException).hpp $(SetRigidContactVelTask).hpp \
			$(MPMTask).hpp $(ElementBase).hpp $(InitializationTask).hpp $(MassAndMomentumTask).hpp $(GridPatch).hpp $(ExtrapolateRigidBCsTask).hpp \
			$(GridForcesTask).hpp $(UpdateMomentaTask).hpp $(UpdateParticlesTask).hpp $(UpdateStrainsLastTask).hpp $(UpdateStrainsLastContactTask).hpp \
			$(MoveCracksTask).hpp $(ResetElementsTask).hpp $(CommonTask).hpp $(RunCustomTasksTask).hpp $(MatPtTractionBC).hpp \
			$(UpdateStrainsFirstTask).hpp $(MPMWarnings).hpp $(BoundaryCondition).hpp $(MPMBase).hpp $(MatPtLoadBC).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(CustomTask).hpp $(MPMWarnings).hpp $(CrackHeader).hpp $(MatPtFluxBC).hpp \
			$(CrackSurfaceContact).hpp $(MeshInfo).hpp $(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(MatPtHeatFluxBC).hpp $(InitVelocityFieldsTask).hpp $(ProjectRigidBCsTask).hpp $(PostExtrapolationTask).hpp $(PostForcesTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NairnMPM).cpp
StartOutput.o : $(StartOutput).cpp $(dprefix) $(NairnMPM).hpp $(MaterialBase).hpp $(ThermalRamp).hpp $(ArchiveData).hpp \
			$(CommonArchiveData).hpp $(BodyForce).hpp $(CrackSurfaceContact).hpp $(CommonException).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(NodalConcBC).hpp $(BoundaryCondition).hpp \
			$(NodalTempBC).hpp $(NodalVelBC).hpp $(MatPtLoadBC).hpp $(MatPtFluxBC).hpp $(CrackHeader).hpp $(MatPtHeatFluxBC).hpp \
			$(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(MatPtTractionBC).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(StartOutput).cpp
MeshInfo.o : $(MeshInfo).cpp $(dprefix) $(MeshInfo).hpp $(GridPatch).hpp $(MPMBase).hpp $(CommonException).hpp $(ElementBase).hpp \
			$(BoundaryCondition).hpp $(NairnMPM).hpp $(NodalPoint).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MeshInfo).cpp
MPMTask.o : $(MPMTask).cpp $(dprefix) $(MPMTask).hpp $(CommonTask).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp $(GridPatch).hpp \
            $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMTask).cpp
InitializationTask.o : $(InitializationTask).cpp $(dprefix) $(InitializationTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(NodalPoint).hpp $(MPMWarnings).hpp $(MatPtLoadBC).hpp $(CrackNode).hpp $(ThermalRamp).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(BoundaryCondition).hpp $(MaterialInterfaceNode).hpp \
            $(GridPatch).hpp $(MPMBase).hpp $(ElementBase).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(InitializationTask).cpp
InitVelocityFieldsTask.o : $(InitVelocityFieldsTask).cpp $(dprefix) $(InitVelocityFieldsTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(GridPatch).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp $(CrackHeader).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp  $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(InitVelocityFieldsTask).cpp
MassAndMomentumTask.o : $(MassAndMomentumTask).cpp $(dprefix) $(MassAndMomentumTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(NodalPoint).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp $(RigidMaterial).hpp \
			$(CommonException).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(GridPatch).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MassAndMomentumTask).cpp
SetRigidContactVelTask.o : $(SetRigidContactVelTask).cpp $(dprefix) $(SetRigidContactVelTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(RigidMaterial).hpp $(MPMBase).hpp \
			$(BoundaryCondition).hpp $(NodalVelBC).hpp $(NodalTempBC).hpp $(NodalConcBC).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SetRigidContactVelTask).cpp
ExtrapolateRigidBCsTask.o : $(ExtrapolateRigidBCsTask).cpp $(dprefix) $(ExtrapolateRigidBCsTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(RigidMaterial).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp \
			$(BoundaryCondition).hpp $(NodalVelBC).hpp $(NodalTempBC).hpp $(NodalConcBC).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ExtrapolateRigidBCsTask).cpp
ProjectRigidBCsTask.o : $(ProjectRigidBCsTask).cpp $(dprefix) $(ProjectRigidBCsTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(RigidMaterial).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp $(CommonException).hpp \
			$(BoundaryCondition).hpp $(NodalVelBC).hpp $(NodalTempBC).hpp $(NodalConcBC).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ProjectRigidBCsTask).cpp
PostExtrapolationTask.o : $(PostExtrapolationTask).cpp $(dprefix) $(PostExtrapolationTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(CommonException).hpp $(CrackHeader).hpp $(CrackSurfaceContact).hpp $(TransportTask).hpp \
			$(CrackNode).hpp $(MaterialInterfaceNode).hpp $(NodalVelBC).hpp $(BoundaryCondition).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PostExtrapolationTask).cpp
UpdateStrainsFirstTask.o : $(UpdateStrainsFirstTask).cpp $(dprefix) $(UpdateStrainsFirstTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(NodalPoint).hpp $(MaterialBase).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateStrainsFirstTask).cpp
GridForcesTask.o : $(GridForcesTask).cpp $(dprefix) $(GridForcesTask).hpp $(MPMTask).hpp $(CommonTask).hpp $(GridPatch).hpp \
			$(NairnMPM).hpp $(NodalPoint).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp $(TransportTask).hpp $(ConductionTask).hpp \
			$(CommonException).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GridForcesTask).cpp
PostForcesTask.o : $(PostForcesTask).cpp $(dprefix) $(PostForcesTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(CrackHeader).hpp $(ConductionTask).hpp $(TransportTask).hpp $(BodyForce).hpp \
			$(CrackNode).hpp $(MaterialInterfaceNode).hpp $(NodalVelBC).hpp $(MatPtTractionBC).hpp $(BoundaryCondition).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PostForcesTask).cpp
UpdateMomentaTask.o : $(UpdateMomentaTask).cpp $(dprefix) $(UpdateMomentaTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(TransportTask).hpp $(NodalPoint).hpp $(CrackNode).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateMomentaTask).cpp
UpdateParticlesTask.o : $(UpdateParticlesTask).cpp $(dprefix) $(UpdateParticlesTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(TransportTask).hpp $(NodalPoint).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ElementBase).hpp \
			$(BodyForce).hpp $(ThermalRamp).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(ConductionTask).hpp \
			$(CommonException).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateParticlesTask).cpp
UpdateStrainsLastTask.o : $(UpdateStrainsLastTask).cpp $(dprefix) $(UpdateStrainsLastTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(TransportTask).hpp $(UpdateStrainsFirstTask).cpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateStrainsLastTask).cpp
UpdateStrainsLastContactTask.o : $(UpdateStrainsLastContactTask).cpp $(dprefix) $(UpdateStrainsLastContactTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(TransportTask).hpp $(MassAndMomentumTask).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(NodalVelBC).hpp $(BoundaryCondition).hpp \
			$(UpdateStrainsFirstTask).cpp $(CrackNode).hpp $(GridPatch).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(UpdateStrainsLastContactTask).cpp
RunCustomTasksTask.o : $(RunCustomTasksTask).cpp $(dprefix) $(RunCustomTasksTask).hpp $(MPMTask).hpp $(CommonTask).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(NodalPoint).hpp $(CustomTask).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(RunCustomTasksTask).cpp
MoveCracksTask.o : $(MoveCracksTask).cpp $(dprefix) $(MoveCracksTask).hpp $(MPMTask).hpp $(CommonTask).hpp $(TractionLaw).hpp \
			$(MaterialBase).hpp $(NairnMPM).hpp $(CrackHeader).hpp $(CrackSurfaceContact).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoveCracksTask).cpp
ResetElementsTask.o : $(ResetElementsTask).cpp $(dprefix) $(ResetElementsTask).hpp $(MPMTask).hpp $(CommonTask).hpp $(BodyForce).hpp \
			$(NairnMPM).hpp $(MPMBase).hpp $(ElementBase).hpp $(MPMWarnings).hpp $(CommonException).hpp $(GridPatch).hpp $(MeshInfo).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ResetElementsTask).cpp

# MPM: Read_MPM
MPMReadHandler.o : $(MPMReadHandler).cpp $(dprefix) $(NairnMPM).hpp $(ArchiveData).hpp $(MPMReadHandler).hpp \
			$(CommonReadHandler).hpp $(CommonArchiveData).hpp $(CustomTask).hpp $(ReverseLoad).hpp $(ThermalRamp).hpp \
			$(BodyForce).hpp $(CrackSurfaceContact).hpp $(GlobalQuantity).hpp $(MatPoint2D).hpp $(MPMBase).hpp $(NodalConcBC).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(NodalTempBC).hpp $(NodalVelBC).hpp $(MatPtLoadBC).hpp $(MatPtFluxBC).hpp \
			$(BoundaryCondition).hpp $(NodesController).hpp $(ElementsController).hpp $(ParseController).hpp $(MaterialController).hpp \
			$(MpsController).hpp $(MaterialController).hpp $(CrackController).hpp $(ElementBase).hpp $(MaterialBase).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(TransportTask).hpp $(MatPoint3D).hpp $(MatPtTractionBC).hpp $(AdjustTimeStepTask).hpp \
			$(PolygonController).hpp $(ShapeController).hpp $(SphereController).hpp $(ShellController).hpp $(RigidMaterial).hpp \
			$(VTKArchive).hpp $(HistoryArchive).hpp $(CarnotCycle).hpp $(MatPtHeatFluxBC).hpp $(GridArchive).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(MeshInfo).hpp $(PropagateTask).hpp $(PolyhedronController).hpp \
            $(MatPointAS).hpp $(PressureLaw).hpp $(TaitLiquid).hpp $(ContactLaw).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMReadHandler).cpp
Generators.o : $(Generators).cpp $(dprefix) $(NairnMPM).hpp $(MPMReadHandler).hpp $(CommonReadHandler).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(ElementBase).hpp $(MatPoint2D).hpp $(NodalConcBC).hpp $(NodalTempBC).hpp $(NodalVelBC).hpp \
			$(MatPtLoadBC).hpp $(MatPtFluxBC).hpp $(BoundaryCondition).hpp $(ThermalRamp).hpp $(StrX).hpp $(LineController).hpp \
			$(MpsController).hpp $(FourNodeIsoparam).hpp $(ParseController).hpp $(Linear2D).hpp $(NodalPoint).hpp $(BoxController).hpp \
			$(MeshInfo).hpp $(EightNodeIsoparamBrick).hpp $(ElementBase3D).hpp $(MatPoint3D).hpp $(ShapeController).hpp \
			$(RectController).hpp $(ArcController).hpp $(NodalPoint2D).hpp $(NodalPoint3D).hpp $(CrackVelocityField).hpp \
			$(OvalController).hpp $(MatVelocityField).hpp $(ElementsController).hpp $(MatPtTractionBC).hpp $(MatPointAS).hpp \
            $(TorusController).hpp $(MatPtHeatFluxBC).hpp $(ResetElementsTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Generators).cpp
BitMapFiles.o : $(BitMapFiles).cpp $(dprefix) $(NairnMPM).hpp $(MPMReadHandler).hpp $(CommonReadHandler).hpp \
			$(BMPLevel).hpp $(MatPoint2D).hpp $(MatPointAS).hpp $(ElementBase).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp \
			$(MpsController).hpp $(MPMBase).hpp $(ParseController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BitMapFiles).cpp
MpsController.o : $(MpsController).cpp $(dprefix) $(MpsController).hpp $(MPMBase).hpp $(ParseController).hpp \
			 $(DiffusionTask).hpp $(ThermalRamp).hpp $(NairnMPM).hpp $(ResetElementsTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MpsController).cpp
CrackController.o : $(CrackController).cpp $(dprefix) $(CrackController).hpp $(ParseController).hpp $(CrackHeader).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackController).cpp
SphereController.o : $(SphereController).cpp $(dprefix) $(SphereController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SphereController).cpp
ShellController.o : $(ShellController).cpp $(dprefix) $(ShellController).hpp $(ShapeController).hpp $(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ShellController).cpp
TorusController.o : $(TorusController).cpp $(dprefix) $(TorusController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TorusController).cpp
PolyhedronController.o : $(PolyhedronController).cpp $(dprefix) $(PolyhedronController).hpp $(ShapeController).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PolyhedronController).cpp
PolyTriangle.o : $(PolyTriangle).cpp $(dprefix) $(PolyTriangle).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PolyTriangle).cpp

# MPM: Exceptions
MPMWarnings.o : $(MPMWarnings).cpp $(dprefix) $(MPMWarnings).hpp $(NairnMPM).hpp $(ArchiveData).hpp \
			$(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMWarnings).cpp

# MPM: Global Quantities
GlobalQuantity.o : $(GlobalQuantity).cpp $(dprefix) $(GlobalQuantity).hpp $(NairnMPM).hpp $(MaterialBase).hpp \
			$(ThermalRamp).hpp $(MPMBase).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(ArchiveData).hpp $(BoundaryCondition).hpp $(NodalVelBC).hpp $(CommonArchiveData).hpp $(NodalTempBC).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GlobalQuantity).cpp
ThermalRamp.o : $(ThermalRamp).cpp $(dprefix) $(ThermalRamp).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ThermalRamp).cpp
BodyForce.o : $(BodyForce).cpp $(dprefix) $(BodyForce).hpp $(mathexpr).hpp $(MPMBase).hpp $(NodalPoint).hpp \
            $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BodyForce).cpp

# MPM: MPM_Classes
MPMBase.o : $(MPMBase).cpp $(dprefix) $(MPMBase).hpp $(CrackHeader).hpp $(MaterialBase).hpp $(MatPtFluxBC).hpp $(MatPtHeatFluxBC).hpp \
            $(MatPtTractionBC).hpp $(MatPtLoadBC).hpp $(BoundaryCondition).hpp $(MeshInfo).hpp $(ElementBase).hpp $(BodyForce).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MPMBase).cpp
MatPoint2D.o : $(MatPoint2D).cpp $(dprefix) $(MatPoint2D).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(MeshInfo).hpp \
			$(NodalPoint).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(CommonException).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPoint2D).cpp
MatPointAS.o : $(MatPointAS).cpp $(dprefix) $(MatPointAS).hpp $(MatPoint2D).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(MeshInfo).hpp \
			$(NodalPoint).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(CommonException).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPointAS).cpp
MatPoint3D.o : $(MatPoint3D).cpp $(dprefix) $(MatPoint3D).hpp $(MPMBase).hpp $(MaterialBase).hpp $(ElementBase).hpp $(MeshInfo).hpp \
			$(NodalPoint).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(CommonException).hpp $(BoundaryCondition).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPoint3D).cpp

# MPM: Cracks
CrackHeader.o : $(CrackHeader).cpp $(dprefix) $(MaterialBase).hpp $(NairnMPM).hpp $(MeshInfo).hpp \
			$(ArchiveData).hpp $(CommonException).hpp $(CommonArchiveData).hpp $(CrackHeader).hpp $(CommonException).hpp \
			$(ElementBase).hpp $(CrackSurfaceContact).hpp $(ContourPoint).hpp $(CrackSegment).hpp $(CrackLeaf).hpp \
			$(NodalPoint).hpp $(MPMBase).hpp $(CrackNode).hpp $(NodalPoint2D).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(CrossedCrack).hpp $(ParseController).hpp $(PropagateTask).hpp $(CustomTask).hpp $(ConductionTask).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackHeader).cpp
CrackLeaf.o : $(CrackLeaf).cpp $(dprefix) $(CrackLeaf).hpp $(CrackHeader).hpp $(CrackSegment).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackLeaf).cpp
CrackSegment.o : $(CrackSegment).cpp $(dprefix) $(CrackSegment).hpp $(MaterialBase).hpp $(ArchiveData).hpp \
			$(ElementBase).hpp $(NairnMPM).hpp $(CommonArchiveData).hpp $(CrackHeader).hpp $(NodalPoint).hpp \
			$(TractionLaw).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackSegment).cpp
ContourPoint.o : $(ContourPoint).cpp $(dprefix) $(ContourPoint).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ContourPoint).cpp
CrackSurfaceContact.o : $(CrackSurfaceContact).cpp $(dprefix) $(CrackSurfaceContact).hpp $(NodalPoint).hpp $(NairnMPM).hpp \
			$(MPMBase).hpp $(CrackHeader).hpp $(MaterialBase).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(MeshInfo).hpp \
            $(ConductionTask).hpp $(ContactLaw).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackSurfaceContact).cpp
CrackNode.o : $(CrackNode).cpp $(dprefix) $(CrackNode).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrackNode).cpp
CrossedCrack.o : $(CrossedCrack).cpp $(dprefix) $(CrossedCrack).hpp $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CrossedCrack).cpp

# MPM: Materials
MaterialBaseMPM.o : $(MaterialBaseMPM).cpp $(dprefix) $(MaterialBase).hpp $(NairnMPM).hpp $(ThermalRamp).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(DiffusionTask).hpp $(TransportTask).hpp $(CommonException).hpp \
			$(CrackSurfaceContact).hpp $(HardeningLawBase).hpp $(LinearHardening).hpp $(DDBHardening).hpp \
			$(NonlinearHardening).hpp $(JohnsonCook).hpp $(SCGLHardening).hpp $(SLMaterial).hpp $(Nonlinear2Hardening).hpp \
			$(ContactLaw).hpp $(BodyForce).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MaterialBaseMPM).cpp
ElasticMPM.o : $(ElasticMPM).cpp $(dprefix) $(Elastic).hpp $(MaterialBase).hpp $(ThermalRamp).hpp $(ConductionTask).hpp \
			 $(MPMBase).hpp $(DiffusionTask).hpp $(TransportTask).hpp $(MeshInfo).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElasticMPM).cpp
MoreIsotropicMat.o : $(MoreIsotropicMat).cpp $(dprefix) $(Elastic).hpp $(MaterialBase).hpp $(ConductionTask).hpp \
			 $(MPMBase).hpp $(DiffusionTask).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoreIsotropicMat).cpp
IsoPlasticity.o : $(IsoPlasticity).cpp $(dprefix) $(IsoPlasticity).hpp $(IsotropicMat).hpp $(HardeningLawBase).hpp \
			$(MaterialBase).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp $(LinearHardening).hpp \
			$(TransportTask).hpp $(Elastic).hpp $(CommonException).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(IsoPlasticity).cpp
AnisoPlasticity.o : $(AnisoPlasticity).cpp $(dprefix) $(AnisoPlasticity).hpp $(Orthotropic).hpp $(TransIsotropic).hpp \
			$(MaterialBase).hpp $(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp $(MPMWarnings).hpp \
			$(TransportTask).hpp $(Elastic).hpp $(CommonException).hpp $(NairnMPM).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(AnisoPlasticity).cpp
Viscoelastic.o : $(Viscoelastic).cpp $(dprefix) $(Viscoelastic).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ThermalRamp).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(TransportTask).hpp $(CommonException).hpp $(NairnMPM).hpp $(MPMWarnings).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Viscoelastic).cpp
HyperElastic.o : $(HyperElastic).cpp $(dprefix) $(HyperElastic).hpp $(MaterialBase).hpp $(TransportTask).hpp \
			$(ConductionTask).hpp $(DiffusionTask).hpp $(ThermalRamp).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HyperElastic).cpp
Mooney.o : $(Mooney).cpp $(dprefix) $(Mooney).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(TransportTask).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Mooney).cpp
Neohookean.o : $(Neohookean).cpp $(dprefix) $(Neohookean).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(TransportTask).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Neohookean).cpp
ClampedNeohookean.o : $(ClampedNeohookean).cpp $(dprefix) $(ClampedNeohookean).hpp $(Neohookean).hpp $(MaterialBase).hpp $(HyperElastic).hpp  \
			$(NairnMPM).hpp $(MPMBase).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ClampedNeohookean).cpp
HEIsotropic.o : $(HEIsotropic).cpp $(dprefix) $(HEIsotropic).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(TransportTask).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp $(HardeningLawBase).hpp $(LinearHardening).hpp \
			$(CommonException).hpp $(NairnMPM).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HEIsotropic).cpp
HEMGEOSMaterial.o : $(HEMGEOSMaterial).cpp $(dprefix) $(HEMGEOSMaterial).hpp $(HEIsotropic).hpp $(NairnMPM).hpp \
			$(HyperElastic).hpp $(MaterialBase).hpp $(MPMBase).hpp $(ConductionTask).hpp $(HardeningLawBase).hpp \
			$(DiffusionTask).hpp $(TransportTask).hpp $(ThermalRamp).hpp $(CommonException).hpp $(MPMWarnings).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HEMGEOSMaterial).cpp
IdealGas.o : $(IdealGas).cpp $(dprefix) $(IdealGas).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(TransportTask).hpp \
			$(DiffusionTask).hpp $(ConductionTask).hpp $(MPMBase).hpp $(ThermalRamp).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(IdealGas).cpp
TaitLiquid.o : $(TaitLiquid).cpp $(dprefix) $(TaitLiquid).hpp $(MaterialBase).hpp $(HyperElastic).hpp $(MPMBase).hpp \
			$(CommonException).hpp $(DiffusionTask).hpp $(ThermalRamp).hpp $(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TaitLiquid).cpp
HillPlastic.o : $(HillPlastic).cpp $(dprefix) $(HillPlastic).hpp $(AnisoPlasticity).hpp $(Orthotropic).hpp \
			$(TransIsotropic).hpp $(MaterialBase).hpp $(Elastic).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HillPlastic).cpp
WoodMaterial.o : $(WoodMaterial).cpp $(dprefix) $(WoodMaterial).hpp $(HillPlastic).hpp $(AnisoPlasticity).hpp $(Orthotropic).hpp \
			$(TransIsotropic).hpp $(MaterialBase).hpp $(Elastic).hpp $(MPMBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(WoodMaterial).cpp
BistableIsotropic.o : $(BistableIsotropic).cpp $(dprefix) $(BistableIsotropic).hpp $(IsotropicMat).hpp \
			$(MaterialBase).hpp $(Elastic).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BistableIsotropic).cpp
RigidMaterial.o : $(RigidMaterial).cpp $(dprefix) $(RigidMaterial).hpp $(MaterialBase).hpp $(CommonException).hpp \
			$(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(RigidMaterial).cpp
    
# MPM: Hardening Laws
HardeningLawBase.o : $(HardeningLawBase).cpp $(dprefix) $(HardeningLawBase).hpp $(MaterialBase).hpp $(MPMBase).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HardeningLawBase).cpp
LinearHardening.o : $(LinearHardening).cpp $(dprefix) $(LinearHardening).hpp $(HardeningLawBase).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinearHardening).cpp
NonlinearHardening.o : $(NonlinearHardening).cpp $(dprefix) $(NonlinearHardening).hpp $(HardeningLawBase).hpp $(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NonlinearHardening).cpp
Nonlinear2Hardening.o : $(Nonlinear2Hardening).cpp $(dprefix) $(Nonlinear2Hardening).hpp $(NonlinearHardening).hpp $(HardeningLawBase).hpp \
			$(MaterialBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(Nonlinear2Hardening).cpp
JohnsonCook.o : $(JohnsonCook).cpp $(dprefix) $(JohnsonCook).hpp $(HardeningLawBase).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(ThermalRamp).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(JohnsonCook).cpp
SCGLHardening.o : $(SCGLHardening).cpp $(dprefix) $(SCGLHardening).hpp $(HardeningLawBase).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(ThermalRamp).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SCGLHardening).cpp
SLMaterial.o : $(SLMaterial).cpp $(dprefix) $(SLMaterial).hpp $(SCGLHardening).hpp $(HardeningLawBase).hpp \
			$(MaterialBase).hpp $(MPMBase).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(SLMaterial).cpp
DDBHardening.o : $(DDBHardening).cpp $(dprefix) $(DDBHardening).hpp $(HardeningLawBase).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(ConductionTask).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(DDBHardening).cpp

# MPM: Traction Laws
TractionLaw.o : $(TractionLaw).cpp $(dprefix) $(TractionLaw).hpp $(MaterialBase).hpp $(CrackSegment).hpp $(CrackHeader).hpp \
			$(CommonArchiveData).hpp $(ArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TractionLaw).cpp
CohesiveZone.o : $(CohesiveZone).cpp $(dprefix) $(CohesiveZone).hpp $(TractionLaw).hpp $(MaterialBase).hpp $(CrackSegment).hpp $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CohesiveZone).cpp
PressureLaw.o : $(PressureLaw).cpp $(dprefix) $(PressureLaw).hpp $(TractionLaw).hpp $(MaterialBase).hpp $(CrackSegment).hpp \
			$(CrackHeader).hpp $(mathexpr).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PressureLaw).cpp
CoupledSawTooth.o : $(CoupledSawTooth).cpp $(dprefix) $(CoupledSawTooth).hpp $(CohesiveZone).hpp $(TractionLaw).hpp $(MaterialBase).hpp \
			$(CrackSegment).hpp $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CoupledSawTooth).cpp
TrilinearTraction.o : $(TrilinearTraction).cpp $(dprefix) $(TrilinearTraction).hpp $(CohesiveZone).hpp $(TractionLaw).hpp \
            $(MaterialBase).hpp $(CrackSegment).hpp $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TrilinearTraction).cpp
LinearTraction.o : $(LinearTraction).cpp $(dprefix) $(LinearTraction).hpp $(CohesiveZone).hpp $(TractionLaw).hpp \
			$(MaterialBase).hpp $(CrackSegment).hpp $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinearTraction).cpp
CubicTraction.o : $(CubicTraction).cpp $(dprefix) $(CubicTraction).hpp $(TractionLaw).hpp $(MaterialBase).hpp \
			$(CrackSegment).hpp $(CrackHeader).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CubicTraction).cpp

# MPM: Contact Laws
ContactLaw.o : $(ContactLaw).cpp $(dprefix) $(ContactLaw).hpp $(MaterialBase).hpp $(MaterialController).hpp \
			$(CoulombFriction).hpp $(LinearInterface).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ContactLaw).cpp
CoulombFriction.o : $(CoulombFriction).cpp $(dprefix) $(CoulombFriction).hpp $(ContactLaw).hpp  $(MaterialBase).hpp \
			$(CrackSurfaceContact).hpp $(ConductionTask).hpp $(NairnMPM).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CoulombFriction).cpp
LinearInterface.o : $(LinearInterface).cpp $(dprefix) $(LinearInterface).hpp $(ContactLaw).hpp  $(MaterialBase).hpp \
			$(UnitsController).hpp $(MeshInfo).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(LinearInterface).cpp
	
# MPM: Boundary_Conditions
BoundaryCondition.o : $(BoundaryCondition).cpp $(dprefix) $(BoundaryCondition).hpp \
			$(mathexpr).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(BoundaryCondition).cpp
MatPtLoadBC.o : $(MatPtLoadBC).cpp $(dprefix) $(MatPtLoadBC).hpp $(BoundaryCondition).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(MeshInfo).hpp $(ElementBase).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPtLoadBC).cpp
MatPtTractionBC.o : $(MatPtTractionBC).cpp $(dprefix) $(MatPtTractionBC).hpp $(MatPtLoadBC).hpp $(BoundaryCondition).hpp \
			$(MaterialBase).hpp $(MPMBase).hpp $(MeshInfo).hpp $(ElementBase).hpp $(MaterialBase).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPtTractionBC).cpp
MatPtFluxBC.o : $(MatPtFluxBC).cpp $(dprefix) $(MatPtFluxBC).hpp $(MatPtLoadBC).hpp $(BoundaryCondition).hpp \
			$(MPMBase).hpp $(MaterialBase).hpp $(MeshInfo).hpp $(mathexpr).hpp $(NairnMPM).hpp $(NodalPoint).hpp \
            $(CrackVelocityField).hpp $(MatVelocityField).hpp $(TransportTask).hpp $(DiffusionTask).cpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPtFluxBC).cpp
MatPtHeatFluxBC.o : $(MatPtHeatFluxBC).cpp $(dprefix) $(MatPtHeatFluxBC).hpp $(MatPtLoadBC).hpp $(BoundaryCondition).hpp \
			$(MPMBase).hpp $(MaterialBase).hpp $(MeshInfo).hpp $(mathexpr).hpp $(NairnMPM).hpp $(NodalPoint).hpp \
            $(CrackVelocityField).hpp $(MatVelocityField).hpp $(TransportTask).hpp $(ConductionTask).cpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MatPtHeatFluxBC).cpp
NodalVelBC.o : $(NodalVelBC).cpp $(dprefix) $(NodalVelBC).hpp $(BoundaryCondition).hpp $(NodalPoint).hpp $(MeshInfo).hpp \
			$(NairnMPM).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalVelBC).cpp
NodalConcBC.o : $(NodalConcBC).cpp $(dprefix) $(NodalConcBC).hpp $(BoundaryCondition).hpp $(NodalPoint).hpp \
			 $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalConcBC).cpp 
NodalTempBC.o : $(NodalTempBC).cpp $(dprefix) $(NodalTempBC).hpp $(BoundaryCondition).hpp $(NodalPoint).hpp $(NairnMPM).hpp \
			$(CrackVelocityField).hpp $(MatVelocityField).hpp $(TransportTask).hpp $(ConductionTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalTempBC).cpp

# MPM: Custom_Tasks
CustomTask.o : $(CustomTask).cpp $(dprefix) $(CustomTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CustomTask).cpp
CalcJKTask.o : $(CalcJKTask).cpp $(dprefix) $(CalcJKTask).hpp $(MaterialBase).hpp $(CustomTask).hpp $(GridPatch).hpp \
			$(ArchiveData).cpp $(CommonArchiveData).hpp $(NodalPoint).hpp $(NairnMPM).hpp $(MPMBase).hpp $(ElementBase).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(MPMTask).hpp $(CommonTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CalcJKTask).cpp
PropagateTask.o : $(PropagateTask).cpp $(dprefix) $(PropagateTask).hpp $(CalcJKTask).hpp $(NairnMPM).hpp $(ArchiveData).hpp \
			$(MaterialBase).hpp $(ElementBase).hpp $(CustomTask).hpp $(MPMBase).hpp $(ConductionTask).hpp $(CommonArchiveData).hpp \
			$(CrackHeader).hpp $(CrackSegment).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(PropagateTask).cpp
ReverseLoad.o : $(ReverseLoad).cpp $(dprefix) $(ReverseLoad).hpp $(PropagateTask).hpp $(NairnMPM).hpp $(ArchiveData).hpp \
			$(MatPtLoadBC).hpp $(MaterialBase).hpp $(CustomTask).hpp $(BoundaryCondition).hpp $(GlobalQuantity).hpp \
			$(CommonException).hpp $(CrackHeader).hpp $(MatPtLoadBC).hpp $(MPMBase).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ReverseLoad).cpp
TransportTask.o : $(TransportTask).cpp $(dprefix) $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(TransportTask).cpp
DiffusionTask.o : $(DiffusionTask).cpp $(dprefix) $(DiffusionTask).hpp $(NairnMPM).hpp $(ElementBase).hpp $(CommonException).hpp \
			$(MaterialBase).hpp $(NodalConcBC).hpp $(MatPtFluxBC).hpp $(BoundaryCondition).hpp $(MatPtLoadBC).hpp \
			$(MPMBase).hpp $(NodalPoint).hpp $(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(DiffusionTask).cpp
ConductionTask.o : $(ConductionTask).cpp $(dprefix) $(NairnMPM).hpp $(ElementBase).hpp $(MaterialBase).hpp $(RigidMaterial).hpp \
			$(NodalTempBC).hpp $(ConductionTask).hpp $(BoundaryCondition).hpp $(MPMBase).hpp $(ThermalRamp).hpp $(CommonException).hpp \
			$(CrackHeader).hpp $(NodalPoint).hpp $(CrackSegment).hpp $(TransportTask).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp \
			$(MatPtHeatFluxBC).hpp $(CrackSurfaceContact).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ConductionTask).cpp
GridArchive.o : $(GridArchive).cpp $(dprefix) $(GridArchive).hpp $(CustomTask).hpp $(NairnMPM).hpp $(ArchiveData).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GridArchive).cpp
VTKArchive.o : $(VTKArchive).cpp $(dprefix) $(VTKArchive).hpp $(CustomTask).hpp $(GridArchive).hpp $(NairnMPM).hpp $(NodalPoint).hpp \
			$(CommonException).hpp $(ArchiveData).hpp $(MeshInfo).hpp $(MaterialBase).hpp \
			$(MPMBase).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(VTKArchive).cpp
HistoryArchive.o : $(HistoryArchive).cpp $(dprefix) $(HistoryArchive).hpp $(CustomTask).hpp $(NairnMPM).hpp \
			$(ArchiveData).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(HistoryArchive).cpp
AdjustTimeStepTask.o : $(AdjustTimeStepTask).cpp $(dprefix) $(AdjustTimeStepTask).hpp $(CustomTask).hpp $(NairnMPM).hpp \
			$(ArchiveData).hpp $(MeshInfo).hpp $(MaterialBase).hpp $(MPMBase).hpp $(CommonArchiveData).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(AdjustTimeStepTask).cpp
CarnotCycle.o : $(CarnotCycle).cpp $(dprefix) $(CarnotCycle).hpp $(CustomTask).hpp $(NairnMPM).hpp $(ThermalRamp).hpp \
			$(CommonException).hpp $(MaterialBase).hpp $(MPMBase).hpp $(CustomTask).hpp \
			$(ConductionTask).hpp $(MaterialBase).hpp $(BodyForce).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(CarnotCycle).cpp

# MPM: Elements
MoreMPMElementBase.o : $(MoreMPMElementBase).cpp $(dprefix) $(ElementBase).hpp $(NodalPoint).hpp $(MPMBase).hpp $(NairnMPM).hpp \
			$(MeshInfo).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(MoreMPMElementBase).cpp
ElementBase3D.o : $(ElementBase3D).cpp $(dprefix) $(ElementBase).hpp $(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(ElementBase3D).cpp
EightNodeIsoparamBrick.o : $(EightNodeIsoparamBrick).cpp $(dprefix) $(ElementBase3D).hpp $(ElementBase).hpp \
			$(NodalPoint).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp 
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(EightNodeIsoparamBrick).cpp
	
# MPM: Nodes
NodalPointMPM.o : $(NodalPointMPM).cpp $(dprefix) $(NodalPoint).hpp $(NairnMPM).hpp $(ArchiveData).hpp $(MaterialBase).hpp \
			$(CommonArchiveData).hpp $(CrackSurfaceContact).hpp $(CommonException).hpp $(MPMWarnings).hpp \
			$(CrackNode).hpp $(MeshInfo).hpp $(CrackSegment).hpp $(CrackHeader).hpp $(CrackVelocityField).hpp \
			$(MatVelocityField).hpp $(BoundaryCondition).hpp $(CrackVelocityFieldMulti).hpp $(TransportTask).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(NodalPointMPM).cpp

# MPM: Patches
GridPatch.o : $(GridPatch).cpp $(dprefix) $(GridPatch).hpp $(GhostNode).hpp $(MeshInfo).hpp $(NodalPoint).hpp \
			$(MPMBase).hpp $(MaterialBase).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp $(CommonException).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GridPatch).cpp
GhostNode.o : $(GhostNode).cpp $(dprefix) $(GhostNode).hpp $(MeshInfo).hpp $(NodalPoint).hpp $(CommonException).hpp \
			$(NodalPoint2D).hpp $(NodalPoint3D).hpp $(CrackVelocityField).hpp $(MatVelocityField).hpp
	$(CC) $(CFLAGS) $(headers) -include $(prefix) $(GhostNode).cpp


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
