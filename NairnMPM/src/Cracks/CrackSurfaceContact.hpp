/********************************************************************************
    CrackSurfaceContact.hpp
    NairnMPM
    
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

// crack contact law
enum { NOCONTACT=0,STICK,FRICTIONLESS,FRICTIONAL,IMPERFECT_INTERFACE};

enum { IN_CONTACT=0,SEPARATED };

// method to find the normal
enum { MAXIMUM_VOLUME_GRADIENT=0,MAXIMUM_VOLUME,EACH_MATERIALS_MASS_GRADIENT,AVERAGE_MAT_VOLUME_GRADIENTS};

class CrackSurfaceContact
{
    public:
		double friction,Dn,Dnc,Dt;
		double materialFriction,materialDn,materialDnc,materialDt;
		bool hasImperfectInterface,displacementCheck;
		double positionCutoff,materialContactVmin,rigidGradientBias;
		short materialContactLaw;
		int materialNormalMethod;
		
        // constructors and destructors
        CrackSurfaceContact();
    
        // methods
		void AddDisplacementVolumeTask1(short,int,NodalPoint *,MPMBase *,double);
		void AddDisplacementTask6(short,int,NodalPoint *,MPMBase *,double);
		short HasContact(int);
		short IsImperfect(int);
#ifndef _BC_CRACK_SIDE_ONLY_
		short GetDeltaMomentum(NodalPoint *np,Vector *,CrackVelocityField *,CrackVelocityField *,Vector *,int,bool,double);
#else
		short GetDeltaMomentum(NodalPoint *np,Vector *,Vector *,CrackVelocityField *,CrackVelocityField *,Vector *,int,bool,double);
#endif
		short MaterialContact(Vector *,Vector *,Vector *,double,bool,double);
		void Output(int);
		void CrackOutput(bool,double,double,double,double,int);
		void MaterialOutput(void);
		short GetInterfaceForce(NodalPoint *,Vector *,CrackVelocityField *,CrackVelocityField *,Vector *,int,double *);
		bool GetMoveOnlySurfaces(void);
		void SetMoveOnlySurfaces(bool);
		bool GetPreventPlaneCrosses(void);
		void SetPreventPlaneCrosses(bool);
		bool GetContactByDisplacements(void);
		void SetContactByDisplacements(bool);
		void SetNormalCODCutoff(double);
		double GetNormalCODCutoff(void);
		int GetMaterialContactLaw(int,int);
		double GetMaterialFriction(int,int);
        void GetMaterialInterface(int,int,double *,double *,double *);
	
		void MaterialContactPairs(int);
	
	private:
		short ContactLaw;
		ContactDetails *CrackContactLaw;
		ContactDetails **mmContact;
		bool moveOnlySurfaces;
		bool contactByDisplacements;
        bool preventPlaneCrosses;
		double normalCODAdjust;

		void TangentialSlipDeltaP(Vector *,Vector *);
		void NormalSlipDeltaP(Vector *,Vector *);
		void FrictionalDeltaP(Vector *,Vector *,int);
};

extern CrackSurfaceContact contact;
    
#endif
