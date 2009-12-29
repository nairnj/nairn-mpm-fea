/********************************************************************************
    FEABoundaryCondition.hpp
    NairnFEA
    
    Created by John Nairn on July 24, 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _FEABOUNDARYCONDITION_

#define _FEABOUNDARYCONDITION_

class FEABoundaryCondition : public LinkedObject
{
    public:
		
		// constructors
		FEABoundaryCondition();
		virtual ~FEABoundaryCondition();
		
		// methods
		virtual LinkedObject *DecrementNodeNum(long,LinkedObject *,bool *);
		virtual void SetValue(double,char *);
		virtual double *GetValuePtr(void);
		
	protected:
        long nodeNum;
		int direction;
		double bcValue;
		static double varXValue,varYValue,varRValue,varZValue;
};

#endif

