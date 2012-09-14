/********************************************************************************
	MatPtTractionBC.hpp
	NairnMPM

	Created by John Nairn on 9/13/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.

	Dependencies
	BoundaryCondition.hpp
********************************************************************************/

#ifndef _MATPTTRACTIONBC_

#define _MATPTTRACTIONBC_

#include "Boundary_Conditions/MatPtLoadBC.hpp"

class MatPtTractionBC: public MatPtLoadBC
{
	public:
		int face;
	
		// constructors and destructors
		MatPtTractionBC(int,int,int,int);
	
		// virtual methods
		virtual BoundaryCondition *PrintBC(ostream &);
		MatPtTractionBC *AddMPTraction(double);
	
		// class methods
		static void SetParticleSurfaceTractions(double);
};

extern MatPtTractionBC *firstTractionPt;

#endif
