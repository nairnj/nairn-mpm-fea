/********************************************************************************
    MaterialController.hpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParseController.hpp
********************************************************************************/

#ifndef _MATERIALCONTROLLER_

#define _MATERIALCONTROLLER_

#include "Read_XML/ParseController.hpp"

class MaterialController : public ParseController
{
    public:
#ifdef MPM_CODE
		double friction,Dn,Dnc,Dt;			// temporary varoables when reading friction
		int otherMatID;
#endif
	
		// methods
		int AddMaterial(int,char *);
		int SetMaterialArray(void);
		char *InputPointer(char *,int &);
		void SetMatColor(float,float,float);
#ifdef MPM_CODE
		void SetCriterion(int,int);
		void SetDirection(int,int);
		void SetTractionMat(int,int);
		void SetMaterialFriction(void);
#endif
};

extern MaterialController *matCtrl;

#endif
