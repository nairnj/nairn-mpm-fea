/*********************************************************************
    MPMPrefix.hpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on feb 3, 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		DataTypes.hpp
*********************************************************************/

#define MPM_CODE

// To compile with OpenMP
#define USE_OPENMP

// include this to write tasks to a .log file during the calculation
// If a crash occurs, the log file can be consulted to see where it
// crashed in the time step
// For quicker check, activate COUT_PROGRESS in NairnMPM.cpp
//#define LOG_PROGRESS

//  compile directives in OSParticulas and in NairnMPM

// New transport analysis methods using FMPM approach
//        FMPM getting gradients in post M&M extrapolation using lumped values (seems best of FMPM)
//        Other options to get gradient (in particle update or in M&M will full mass matrix) cause
//            time step need to be too small
#define TRANSPORT_FMPM

// Activate option for Poroelasticy calculations
#define POROELASTICITY

// option to restart time step if acceleration is too high
#define RESTART_OPTION

// Grid boundary conditions
// When activated, the copied momenta are adjusted for BCs too (regular ones always are)
// Set to 0 (do not adjust), 1, to adjust only those on symmetry planes, 2 to adjust all with BCs
// If 2 is best, probably better to copy after contact and BCs are done (or not copy at all?)
#define ADJUST_COPIED_PK 1

// This option adds new terms in multimaterial mode
// 0: no changes for contact, 1: Correction in paper applied to delta p^\alpha only, 2: future
#define MM_XPIC 1

// MOVECRACKS_EARLY: move cracks right after particle update instead of after custom tasks
//		Issues: does it affect J calculation?
#define MOVECRACKS_EARLY

// Define this to use ASCII map to speed up expressions. This approach limited to small set
// of single character variables, but only a small set is used here
#define USE_ASCII_MAP

// C includes
#ifdef WINDOWS_EXE
#include <stdio.h>
#endif
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <unordered_map>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <ctime>

// For c++89 compliance
using namespace std;

// This includes all data types and some common headers
#include "System/DataTypes.hpp"
#ifdef LOG_PROGRESS
#include "System/ArchiveData.hpp"
#endif

// math utilities
#define fmin(a,b) (((a)>(b))?(b):(a))
#define fmax(a,b) (((a)>(b))?(a):(b))

// general constants
#define TRUE 1
#define FALSE 0
#define PI_CONSTANT 3.141592653589793
#define SQRT_TWOTHIRDS 0.8164965809277260
#define TWOTHIRDS 0.6666666666666667
#define ONETHIRD 0.3333333333333333
#define SQRT_EIGHT27THS 0.5443310539518174

// maxShapeNodes - max number non-zero shape functions +1
extern int maxShapeNodes;
#define MAX_SHAPE_NODES 40

// Max number of nodes in an element + 1
#define MaxElNd 10

// bits for J and K or set both bits
#define NEED_J 1
#define NEED_K 2
#define NEED_JANDK 3

// MPM method
enum { USF_METHOD=0,UNUSED_METHOD,USAVG_METHOD,USL_METHOD};
        
// analysis type
// also defined in FEAPrefix.hpp - keep same
enum { PLANE_STRAIN=0,PLANE_STRESS,AXI_SYM,THREE_DIM,END_FEA_TYPES,
		BEGIN_MPM_TYPES=9,PLANE_STRAIN_MPM,PLANE_STRESS_MPM,THREED_MPM,AXISYMMETRIC_MPM,END_MPM_TYPES };

// propagation status
enum { NOGROWTH=0,GROWNOW };

// crack velocity fields
enum { NO_CRACK=0,ABOVE_CRACK,BELOW_CRACK};

// grid value calculations
#define VELOCITY_FOR_STRAIN_UPDATE 0

// GIMP methods
#define POINT_GIMP 0
#define UNIFORM_GIMP 1
#define LINEAR_CPDI 2
#define QUADRATIC_CPDI 3
#define UNIFORM_GIMP_AS 4
#define LINEAR_CPDI_AS 5
#define FINITE_GIMP 6
#define BSPLINE_GIMP 7
#define BSPLINE 8
#define BSPLINE_CPDI 9
#define BSPLINE_GIMP_AS 10
