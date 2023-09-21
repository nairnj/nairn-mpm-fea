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

class Matrix3;

// vector 3D
typedef struct {
    double x;
    double y;
	double z;
} Vector;

// rect in 2D
typedef struct {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
} Rect;


// index to the sig components
enum { XX=0,YY,ZZ,YZ,XZ,XY,ZY,ZX,YX};

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
        double mass;
        double vx;
        double vy;
		double *matWeight;
	} DispField;
	
	// contact law for material pair
	typedef struct {
		int lawID;
		int matID;
		char *nextFriction;
	} ContactPair;

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

	// variables for multimaterial conduction calculations
	typedef struct {
		double gTValue;		  // material transport value
							  // remove next line if #ifdef USE_GTVALUEREL is removed
		double gTValueRel;	  // to handle diffusion with variable csat (c/csat actual)
		double gTRelValue;	  // to handle diffusion with variable csat (csat relative)
		double gVCT;		  // transport capacity
		double gQ;			  // transport velocity
		double gTstar;		  // for XPIC
		double gTprev;		  // for XPIC
		double gTnext;		  // for XPIC
		bool gFirstBC;	  // flag for first BC on a node
	} TransportField;

	// For diffusion calculations
	typedef struct {
		double conc;
		double prevConc;
		Vector grad;
		double bufferSource;
	} DiffusionInfo;

	// For residual strains in constitutive laws
	typedef struct {
		double dT;
		double dC;
		double doopse;		// for generalized plane stress or strain
	} ResidualStrains;

	// plastic law properties
	typedef struct {
		double alpint;
		double dalpha;
		double strialmag;
	} HardeningAlpha;

	typedef struct {
		Vector Sacc;			// FLIP and XPIC only
		Vector Svtilde;
		Vector Svlumped;		// XPIC(k>1) only
		double gridAlpha;
		double particleAlpha;
		int m;					// m>1 FMPM,m=0 FLIP,m<0 XPIC
	} GridToParticleExtrap;

	// For Finite Gimp,each particle saves a list of elements
	// also the moments of the intersection
	typedef struct {
		int *InTheseElements;
		double *moment_0;
		double *moment_x;
		double *moment_y;
		double *moment_xy;
	} FiniteGIMPInfo;

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
bool DbleEqual(double,double);
int Reverse(char *,int);
void PrimeFactors(int,vector<int> &);

Vector MakeVector(double,double,double);
Vector MakeScaledVector(double,double,double,double);
Vector SetDiffVectors(const Vector *,const Vector *);
Vector SetSumVectors(const Vector *,const Vector *);
Vector SetScaledVector(const Vector *,const double );
Vector *ZeroVector(Vector *);
Vector *CopyVector(Vector *,const Vector *);
Vector *ScaleVector(Vector *,const double);
Vector *CopyScaleVector(Vector *,const Vector *,const double);
Vector *AddVector(Vector *,const Vector *);
Vector *SubVector(Vector *,const Vector *);
Vector *AddScaledVector(Vector *,const Vector *,const double);
Vector *CrossProduct(Vector *,const Vector *,const Vector *);
double CrossProduct2D(const Vector *,const Vector *);
double DotVectors(const Vector *,const Vector *);
double DotVectors2D(const Vector *,const Vector *);
void PrintVector(const char *,const Vector *);

Tensor MakeTensor(double, double, double, double, double, double);
Tensor MakeTensor2D(double, double, double, double);
Tensor *ZeroTensor(Tensor *);
Tensor *AddTensor(Tensor *,Tensor *);
Tensor *SubTensor(Tensor *,Tensor *);
Tensor *ScaleTensor(Tensor *,double);
double DotTensors2D(const Tensor *, const Tensor *);
double Tensor_i(Tensor *,int);
double Tensor_ij(Tensor *,int,int);
void PrintTensor(const char *,Tensor *);

#ifdef MPM_CODE
double DotTensors(const Tensor *, const Tensor *);
Vector TensorEigenvalues(Tensor *,bool,bool);
Matrix3 TensorToMatrix(Tensor *,bool);
Matrix3 TensorToMatrix2D(Tensor *,bool);
#endif

bool XYInRect(double,double,Rect *);
bool PtInRect(Vector *,Rect *);
double RectArea(Rect *);

unsigned charConvertToLower(unsigned char);
int CIstrcmp(const char *, const char *);
void GetFileExtension(const char *,char *,int);
char *MakeDOSPath(char *);

void dout(const char* fmt...);
void doutCritical(const char* fmt...);
void doutMM(const char* fmt...);
void doutMMCritical(const char* fmt...);

#ifdef MPM_CODE
TensorAntisym *ZeroTensorAntisym(TensorAntisym *);
#endif

void InitRandom(unsigned int);
long RandomLong(long,long);
double RandomRange(double,double);
double Random(void);
double NormalCDFInverse(double p);
double RationalApproximation(double t);
bool RealQuadraticRoots(double,double,double,double &,double &);
double SmoothStepRange(double,double,double,int,int);
double SmoothStep(double,int,int);

#ifdef MPM_CODE
double gsl_sf_lambert_W0(double x);
double gsl_sf_lambert_Wm1(double x);
#endif

double gamma_fxn(double z);

#include "System/CommonAnalysis.hpp"
#include "System/LinkedObject.hpp"

#ifdef MPM_CODE
#include "System/Matrix3.hpp"
#include "System/Matrix4.hpp"
#endif

#endif
