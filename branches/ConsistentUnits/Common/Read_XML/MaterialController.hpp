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
		double friction,Dn,Dnc,Dt;			// temporary variables when reading friction
		int otherMatID;
#endif
		ParseController *nameCtrl;
	
		MaterialController(void);
		~MaterialController(void);
	
		// methods
		int AddMaterial(int,char *);
		const char *SetMaterialArray(void);
		char *InputPointer(char *,int &,double &);
		const char *MaterialType(void);
		void SetMatColor(float,float,float,float);
#ifdef MPM_CODE
		void SetCriterion(int,int);
		void SetDirection(int,int);
		void SetTractionMat(int,int);
		void SetMaterialFriction(void);
#endif
		int GetIDFromName(char *);
		int GetIDFromNewName(char *);
};

extern MaterialController *matCtrl;

#endif
