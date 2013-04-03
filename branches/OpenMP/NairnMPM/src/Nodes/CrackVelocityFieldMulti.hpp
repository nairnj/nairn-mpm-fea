/********************************************************************************
	CrackVelocityFieldMulti.hpp
	NairnMPM
 
	Created by John Nairn on 11 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
 
	Dependencies
		CrackVelocityField.hpp, MatVelocityField.hpp
********************************************************************************/

#ifndef _CRACKVELOCITYFIELDMULTI_

#define _CRACKVELOCITYFIELDMULTI_

#include "Nodes/CrackVelocityField.hpp"

class CrackVelocityFieldMulti : public CrackVelocityField
{
	public:
		
        // constructors and destructors
		CrackVelocityFieldMulti(short,int);
		virtual void ZeroMatFields(void);
	
		// specific task methods
		virtual void AddMassTask1(int,double);
		virtual double GetTotalMassAndCount(void);
		virtual void AddVolumeGradient(int,MPMBase *,double,double,double);
		virtual void CombineRigidFrom(CrackVelocityFieldMulti *,int);
		virtual void CopyRigidFrom(CrackVelocityFieldMulti *,int);
	
		virtual void AddFintSpreadTask3(Vector *);
		virtual void AddFextSpreadTask3(Vector *);
		virtual void CalcFtotTask3(double);
	
		virtual void UpdateMomentaOnField(double);
	
		virtual void RezeroNodeTask6(double);
	
		virtual void MaterialContact(int,int,bool,double);
        virtual void GetFrictionalDeltaMomentum(Vector *,Vector *,double,double);
		virtual void GetVolumeGradient(int,NodalPoint *,Vector *,double);
		virtual void RigidMaterialContact(int,int,int,bool,double);
		virtual bool GetInterfaceForcesForNode(Vector *,Vector *,double,double,
										   double,Vector *,double *,double,Vector *,double,bool,bool,double);
		virtual void CalcVelocityForStrainUpdate(void);
	
		// boundary conditions
        virtual void SetMomVel(Vector *);
        virtual void AddMomVel(Vector *,double);
        virtual void SetFtot(Vector *,double);
        virtual void AddFtot(Vector *,double,double);
	
		// accessors
		virtual int GetNumberPointsNonrigid(void);
		virtual void SumAndClearRigidContactForces(Vector *,bool);
		virtual double GetTotalMass(void);
		virtual void AddKineticEnergyAndMass(double &,double &);
		virtual double GetVolumeNonrigid(void);
		virtual double GetVolumeTotal(double);
		virtual Vector GetCMatMomentum(void);
		virtual Vector GetCMDisplacement(void);
		virtual Vector GetCMatFtot(void);
		virtual MatVelocityField *GetRigidMaterialField(int *);
		virtual void ChangeMomentum(Vector *,bool,double);
		virtual int CopyFieldMomenta(Vector *,int);
		virtual int PasteFieldMomenta(Vector *,int);
		virtual void Describe(void);
	
	private:
		// variables (changed in MPM time step)
		int numberMaterials;		// number of materials in this crack velocity field
		int numberRigidPoints;		// number of rigid particles in this velocity field
	
		// constants (not changed in MPM time step)

};

#endif
