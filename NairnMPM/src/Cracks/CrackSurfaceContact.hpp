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

class CrackSurfaceContact
{
    public:
		int crackContactLawID;
		bool crackContactByDisplacements;	// true if crack contact using displacements, false if need to adjust normal COD
		double crackPositionCutoff;			// cut off for normal contact when using position for cracks
	
        // constructors and destructors
        CrackSurfaceContact();
    
        // methods
		short HasContact(int);
		bool GetDeltaMomentum(NodalPoint *np,Vector *,CrackVelocityField *,CrackVelocityField *,Vector *,int,int,double,int *);
		double MaterialSeparation(double,double,Vector *,NodalPoint *,bool,double);
		void Output(void);
		void CustomCrackContactOutput(int &,int);
		bool GetMoveOnlySurfaces(void) const;
		void SetMoveOnlySurfaces(bool);
		bool GetPreventPlaneCrosses(void) const;
		void SetPreventPlaneCrosses(bool);
	
	private:
		ContactLaw **crackContactLaw;
		bool moveOnlySurfaces;
        bool preventPlaneCrosses;
};

extern CrackSurfaceContact contact;
    
#endif
