/*********************************************************************
    DataTypes.hpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on feb 3, 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	This is included in prefix files and it includes
		CommonAnalysis.hpp
		LinkedObject.hpp
*********************************************************************/

#ifndef _DATA_TYPES_

#define _DATA_TYPES_

// vector 3D
typedef struct {
    double x;
    double y;
	double z;
} Vector;

// index to the sig components
enum { XX=0,YY,ZZ,YZ,XZ,XY};

#ifdef MPM_CODE
	// symmetric tensor in contracted notation (3D)
	typedef struct {
		double xx;
		double yy;
		double zz;
		double yz;
		double xz;
		double xy;
	} Tensor;
	
	// antisymmetric tensor in contracted notation
	typedef struct {
		double yz;
		double xz;
		double xy;
	} TensorAntisym;
	
	// full tensor in 3D
	typedef struct {
		double comp[3][3];
	} TensorFull;

	// crack field data
	typedef struct {
		short loc;
		int crackNum;
		Vector norm;
	} CrackField;
	
	// grid properties structure
	typedef struct {
		Vector du;
		Vector dv;
		double kinetic;
		double work;
		Tensor stress;
		double *matWeight;
	} DispField;
	
	// contact law properties, nextFriction and matID only used for material contact
	typedef struct {
		int law;
		double friction;
		double Dn;
		double Dnc;
		double Dt;
		char *nextFriction;
		int matID;
	} ContactDetails;

	// for each node in CPDI domain give its element,
	// its naturual coordinates, weighting factor
	// for finding the gradient, and weighting factor
	// for shape function
	typedef struct {
		int inElem;
		Vector ncpos;
		Vector wg;
		double ws;
	} CPDIDomain;

	// Transport Properties
	// conductivity is divided by rho (in g/mm^3)
	typedef struct {
		Tensor diffusionTensor;
		Tensor kCondTensor;
	} TransportProperties;

#else
	// tensor (2D and axisymmetric and plane strain)
	typedef struct {
		double xx;
		double yy;
		double zz;
		double xy;
	} Tensor;

	// nodal forces and stresses
	typedef struct {
		Tensor stress;
		Vector force;
		int numElems;
	} ForceField;
	
	// periodic properties
	typedef struct {
		int dof;
		bool fixDu;
		double du;
		bool fixDudy;
		double dudy;
		double xmin;
		double xmax;
		bool fixDv;
		double dv;
		bool fixDvdx;
		double dvdx;
		double ymin;
		double ymax;
	} PeriodicInfo;

#endif

// common utility methods
void PrintSection(const char *text);
int DbleEqual(double,double);
int Reverse(char *,int);
void PrimeFactors(int,vector<int> &);

Vector MakeVector(double,double,double);
Vector *ZeroVector(Vector *);
Vector *CopyVector(Vector *,const Vector *);
Vector *ScaleVector(Vector *,const double);
Vector *CopyScaleVector(Vector *,const Vector *,const double);
Vector *AddVector(Vector *,const Vector *);
Vector *SubVector(Vector *,const Vector *);
Vector *AddScaledVector(Vector *,const Vector *,const double);
Vector *CrossProduct(Vector *,const Vector *,const Vector *);
double DotVectors(const Vector *,const Vector *);
double DotVectors2D(const Vector *,const Vector *);
void PrintVector(const char *,const Vector *);

Tensor *ZeroTensor(Tensor *);
Tensor *AddTensor(Tensor *,Tensor *);
Tensor *ScaleTensor(Tensor *,double);
double Tensor_i(Tensor *,int);
double Tensor_ij(Tensor *,int,int);
void PrintTensor(const char *,Tensor *);

#ifdef MPM_CODE
TensorAntisym *ZeroTensorAntisym(TensorAntisym *);
#endif

#include "System/CommonAnalysis.hpp"
#include "System/LinkedObject.hpp"

#ifdef MPM_CODE
#include "System/Matrix3.hpp"
#endif

#endif
