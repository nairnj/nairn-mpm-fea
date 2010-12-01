/********************************************************************************
    Constraint.hpp
    NairnFEA
    
    Created by John Nairn on 2/6/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _CONSTRAINT_

#define _CONSTRAINT_

class Constraint : public LinkedObject
{
    public:
		// constructors and destructors
		Constraint(int,int);
        
        // methods
		void AddNode(int,double);
		void AddConstant(double);
		int GetBandWidth(int,int);
		Constraint *MapNodes(int *);
		
		// Accessors
		int NumberNodes(void);
		int NodalDof(int,int);
		double GetQ(void);
		double GetCoeff(int);
		void SetLambdaNum(int);
		int GetLambdaNum(void);
		void Describe(void);
	
	private:
		int lambda,dof;
		vector< int > nodes;
		vector< double > coeff;
		double qvalue;
};

extern Constraint *firstConstraint;
extern int numConstraints;

#endif
