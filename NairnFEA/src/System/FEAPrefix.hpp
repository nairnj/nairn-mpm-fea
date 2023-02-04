/*********************************************************************
    FEAPrefix.hpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on feb 3, 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		DataTypes.hpp
*********************************************************************/

#define FEA_CODE

#define OSPARTICULAS

// Note that when OpenMP is available, it defines _OPENMP
// This define just says to include omp.h
#define USE_OPENMP

// C includes
#ifdef WINDOWS_EXE
#include <stdio.h>
#endif
#include <cstdio>
#include <string.h>
#include <stdlib.h>
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

// math utilities
#define fmin(a,b) (((a)>(b))?(b):(a))
#define fmax(a,b) (((a)>(b))?(a):(b))

// general constants
#define TRUE 1
#define FALSE 0
#define PI_CONSTANT 3.141592653589793
// max nodes in an element + 1
#define MaxElNd 10

// analysis type
// also defined in MPMPrefix.hpp - keep same
enum { PLANE_STRAIN=0,PLANE_STRESS,AXI_SYM,THREE_DIM,END_FEA_TYPES,
        BEGIN_MPM_TYPES=9,PLANE_STRAIN_MPM,PLANE_STRESS_MPM,THREED_MPM,AXISYMMETRIC_MPM,END_MPM_TYPES };

// FEA utilities
int IsEven(int);
void GetLineParameters(int,double,double *,double *);
double GetRatio(int,double);
