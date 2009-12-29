/********************************************************************************
    NodalDispBC.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _NODALDISPBC_

#define _NODALDISPBC_

#include "Boundary_Conditions/FEABoundaryCondition.hpp"

class NodalDispBC : public FEABoundaryCondition
{
    public:
		double angle;
        
        // constructors and destructors
        NodalDispBC(long,int);
		NodalDispBC(long,int,double);
		void DefaultSetup(long,int,double);
		
        // methods
        NodalDispBC *PrintBC(ostream &);
        NodalDispBC *FixOrRotate(double **,double *,long,long,int);
		NodalDispBC *Unrotate(double *,int);
        void PrintReaction(void);
		NodalDispBC *MapNodes(int *);
		bool SameDofSetting(NodalDispBC *);
	
	private:
		int axis;
};

extern NodalDispBC *firstDispBC;

#endif

