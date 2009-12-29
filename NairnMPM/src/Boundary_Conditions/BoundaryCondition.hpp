/********************************************************************************
    BoundaryCondition.hpp
    NairnMPM
    
    Created by John Nairn on Mon Mar 29 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _BOUNDARYCONDITION_

#define _BOUNDARYCONDITION_

class ROperation;

// BC directions (x and y are used at bit locations too)
#define SKEW_DIRECTION 0
#define X_DIRECTION 1
#define Y_DIRECTION 2
#define Z_DIRECTION 4
#define TEMP_DIRECTION 8
#define CONC_DIRECTION 16
#define EXTERNAL_FLUX 1
#define COUPLED_FLUX 2

// loading boundary conditions
enum {	CONSTANT_VALUE=1,LINEAR_VALUE,SINE_VALUE,COSINE_VALUE,SILENT,FUNCTION_VALUE };

class BoundaryCondition : public LinkedObject
{
    public:
        long nodeNum;
        int style;
        double value,ftime,offset;
		
		// constructors
		BoundaryCondition(int,double,double);
		virtual ~BoundaryCondition();
		virtual BoundaryCondition *UnsetDirection(void);
		virtual BoundaryCondition *SetRigidProperties(long,int,int,double);
		
		// pure virtual
        virtual BoundaryCondition *PrintBC(ostream &) = 0;
		
		// virtual methods
		virtual double BCValue(double);
		long GetNodeNum(double);
		long GetNodeNum(void);
		virtual void SetFunction(char *);
		virtual char *GetFunctionString(void);
		virtual void PrintFunction(ostream &);
		virtual void GetPosition(double *,double *,double *,double *);
	
	protected:
		ROperation *function;
		static double varTime,varXValue,varYValue,varZValue,varRotValue;
		double scale;
};

#endif

