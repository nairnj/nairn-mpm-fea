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
//#define LOG_PROGRESS

// Define to test multiple patches, but not be in parallel
// This overrides incompatiible omp calls but does not deactivate #pragma omp,
// which therefore should be commented out. Also assume omp is enabled to
// _OPENMP is defined, even though not being used
//#define _OMPTEST

// Activate debugging messages that look for nan in results that should always be valid numbers
//#define CHECK_NAN

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
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <ctime>

// For c++89 compliance
using namespace std;

// This includes all data types and some common headers
#include "System/DataTypes.hpp"

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

// Max number of particles in an element
#define MaxElParticles 27

// bits for J and K or set both bits
#define NEED_J 1
#define NEED_K 2
#define NEED_JANDK 3

// MPM method
enum { USF_METHOD=0,USL_METHOD,USAVG_METHOD,SZS_METHOD};
        
// analysis type
// also defined in FEAPrefix.hpp - keep same
enum { PLANE_STRAIN=0,PLANE_STRESS,AXI_SYM,THREE_DIM,END_FEA_TYPES,
		BEGIN_MPM_TYPES=9,PLANE_STRAIN_MPM,PLANE_STRESS_MPM,THREED_MPM,AXISYMMETRIC_MPM,END_MPM_TYPES };

// propagation status
enum { NOGROWTH=0,GROWNOW };

// crack velocity fields
enum { NO_CRACK=0,ABOVE_CRACK,BELOW_CRACK};

// GIMP methods
#define POINT_GIMP 0
#define UNIFORM_GIMP 1
#define LINEAR_CPDI 2
#define QUADRATIC_CPDI 3
#define UNIFORM_GIMP_AS 4
#define LINEAR_CPDI_AS 5




