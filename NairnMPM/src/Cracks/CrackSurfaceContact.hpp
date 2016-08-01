/********************************************************************************
    CrackSurfaceContact.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Header for dealing with crack surface contact laws
		Global contact

	Dependencies
		NodalPoint.hpp, CrackVelocityField.hpp, MatVelocityField.hpp
********************************************************************************/

#ifndef _CRACKSURFACECONTACT_

#define _CRACKSURFACECONTACT_

#include "Nodes/NodalPoint.hpp"

class MPMBase;
class NodalPoint;
class CrackVelocityField;
class ContactLaw;

// crack contact law
enum { NOCONTACT=0,STICK,FRICTIONLESS,FRICTIONAL,IMPERFECT_INTERFACE};

enum { IN_CONTACT=0,SEPARATED };

// method to find the normal (0|1|2|3)
enum { MAXIMUM_VOLUME_GRADIENT=0,MAXIMUM_VOLUME,AVERAGE_MAT_VOLUME_GRADIENTS,EACH_MATERIALS_MASS_GRADIENT,
		SPECIFIED_NORMAL };

class CrackSurfaceContact
{
    public:
		int materialContactLawID,crackContactLawID;
		ContactLaw *materialContactLaw;
		bool hasImperfectInterface,displacementCheck;
		double materialContactVmin,rigidGradientBias;
		int materialNormalMethod;
		Vector contactNormal;
		
        // constructors and destructors
        CrackSurfaceContact();
    
        // methods
		short HasContact(int);
		short IsImperfect(int);
		bool GetDeltaMomentum(NodalPoint *np,Vector *,CrackVelocityField *,CrackVelocityField *,Vector *,int,bool,double,int *);
		short MaterialContact(Vector *,Vector *,double,bool,double,NodalPoint *);
		void Output(void);
		void CustomCrackContactOutput(int &,int);
		void MaterialOutput(void);
		bool GetInterfaceForceOnCrack(NodalPoint *,Vector *,CrackVelocityField *,CrackVelocityField *,Vector *,int,double *,double);
		bool GetMoveOnlySurfaces(void) const;
		void SetMoveOnlySurfaces(bool);
		bool GetPreventPlaneCrosses(void) const;
		void SetPreventPlaneCrosses(bool);
		ContactLaw *GetMaterialContactLaw(int,int);
		void SetContactNormal(double,double);
	
		void MaterialContactPairs(int);
	
	private:
		ContactLaw **crackContactLaw;
		ContactLaw ***mmContactLaw;
		bool moveOnlySurfaces;
        bool preventPlaneCrosses;
};

extern CrackSurfaceContact contact;
    
#endif
