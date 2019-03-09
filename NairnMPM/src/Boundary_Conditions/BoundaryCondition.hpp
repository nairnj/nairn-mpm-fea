/********************************************************************************
    BoundaryCondition.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Mar 29 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _BOUNDARYCONDITION_

#define _BOUNDARYCONDITION_

class Expression;

// BC directions (x and y are used at bit locations too)
#define X_DIRECTION 1
#define Y_DIRECTION 2
#define XY_SKEWED_DIRECTION 3
#define Z_DIRECTION 4
#define XZ_SKEWED_DIRECTION 5
#define YZ_SKEWED_DIRECTION 6
#define XYZ_SKEWED_DIRECTION 7

// Higher bits
#define TEMP_DIRECTION 8
#define CONC_DIRECTION 16
#define XSYMMETRYPLANE_DIRECTION 32
#define YSYMMETRYPLANE_DIRECTION 64
#define ZSYMMETRYPLANE_DIRECTION 128
#define ANYSYMMETRYPLANE_DIRECTION 224

// non-bit based

// for input only
#define Z_DIRECTION_INPUT 3
#define XY_SKEWED_INPUT 12
#define XZ_SKEWED_INPUT 13
#define YZ_SKEWED_INPUT 23
#define XYZ_SKEWED_INPUT 123

// Traction
#define N_DIRECTION 11
#define T1_DIRECTION 12
#define T2_DIRECTION 13


// Flux conditions
#define EXTERNAL_FLUX 1
#define COUPLED_FLUX 2

// loading boundary conditions
enum {	CONSTANT_VALUE=1,LINEAR_VALUE,SINE_VALUE,COSINE_VALUE,SILENT,FUNCTION_VALUE };

class BoundaryCondition : public LinkedObject
{
    public:
 		
		// constructors
		BoundaryCondition(int,double,double);
		virtual ~BoundaryCondition();
        virtual int GetSetDirection(void) const;
        virtual BoundaryCondition *UnsetDirection(void);
		virtual BoundaryCondition *SetRigidProperties(int,int,int,double);
		
		// pure virtual
        virtual BoundaryCondition *PrintBC(ostream &);
		
		// virtual methods
		virtual double BCValue(double);
		int GetNodeNum(double);
		int GetNodeNum(void);
		virtual void SetFunction(char *);
		virtual void GetPosition(unordered_map<string,double> & );
		virtual void PrintFunction(ostream &);
		int GetID(void);
		void SetID(int);
		double *GetBCValuePtr(void);
		virtual void SetBCValue(double);
		void SetBCValueCU(double);
		double GetBCValueOut(void);
		virtual void SetBCFirstTime(double);
		void SetBCFirstTimeCU(double);
		double GetBCValue(void);
		double GetBCFirstTime(void);
		virtual double GetBCFirstTimeOut(void);
		void SetBCOffset(double);
		double GetBCOffset(void);
		int GetBCStyle();
	
	protected:
		Expression *function;
		double scale;
		int nodeNum;
		int bcID;
		int style;
	
	private:
		double value,ftime,offset;
};

#endif

