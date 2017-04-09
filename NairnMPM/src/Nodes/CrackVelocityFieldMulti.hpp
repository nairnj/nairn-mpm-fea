/********************************************************************************
	CrackVelocityFieldMulti.hpp
	nairn-mpm-fea
 
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
		CrackVelocityFieldMulti(int,short,int);
		virtual ~CrackVelocityFieldMulti();
		virtual void ZeroMatFields(void);
		virtual void AddMatVelocityField(int);
        virtual bool NeedsMatVelocityField(int) const;
		virtual void MatchMatVelocityFields(MatVelocityField **);
	
		// specific task methods
		virtual void AddMassTask1(int,double,int);
		virtual double GetTotalMassAndCount(void);
		virtual void AddVolumeGradient(int,MPMBase *,double,double,double);
		virtual void CopyVolumeGradient(int,Vector *);
		virtual void CopyMassAndMomentum(NodalPoint *);
        virtual void CopyMassAndMomentumLast(NodalPoint *);
	
		virtual void AddFtotSpreadTask3(Vector *);
		virtual void CopyGridForces(NodalPoint *);
		virtual void AddGravityAndBodyForceTask3(Vector *);
		virtual void RestoreMomenta(void);
	
		virtual void UpdateMomentaOnField(double);

		virtual void RezeroNodeTask6(double);

		virtual void MaterialContactOnCVF(NodalPoint *,double,int);
		virtual bool HasVolumeGradient(int) const;
		virtual void GetVolumeGradient(int,const NodalPoint *,Vector *,double) const;
		virtual void RigidMaterialContactOnCVF(int,bool,NodalPoint *,double,int);
		virtual double GetContactArea(NodalPoint *,double,double,Vector *,Vector *,Vector *,double *,double *) const;
		virtual void CalcVelocityForStrainUpdate(void);
	
		// boundary conditions
        virtual void SetMomVel(Vector *);
        virtual void AddMomVel(Vector *,double);
		virtual void ReflectMomVel(Vector *,CrackVelocityField *,double);
        virtual void SetFtotDirection(Vector *,double,Vector *);
        virtual void AddFtotDirection(Vector *,double,double,Vector *);
		virtual void ReflectFtotDirection(Vector *,double,CrackVelocityField *,double,Vector *);
	
		// accessors
		virtual bool HasPointsNonrigid(void) const;
		virtual void SumAndClearRigidContactForces(Vector *,bool,double,Vector *);
		virtual double GetTotalMass(bool) const;
		virtual void AddKineticEnergyAndMass(double &,double &);
		virtual double GetVolumeNonrigid(bool) const;
		virtual double GetVolumeTotal(NodalPoint *) const;
		virtual Vector GetCMatMomentum(bool &,double *) const;
		virtual Vector GetCMDisplacement(NodalPoint *,bool) const;
		virtual Vector GetCMatFtot(void);
		virtual void ChangeCrackMomentum(Vector *,bool,double);
		virtual int CopyFieldMomenta(Vector *,int);
#ifdef ADJUST_EXTRAPOLATED_PK_FOR_SYMMETRY
		virtual void AdjustForSymmetryBC(NodalPoint *);
#endif
		virtual int PasteFieldMomenta(Vector *,int);
		virtual void Describe(void) const;
	
		// NFM rigid materials ignore cracks
		virtual void CopyRigidFrom(MatVelocityField *,int);
		virtual MatVelocityField *GetRigidMaterialField(int *);
	
	private:
		// variables (changed in MPM time step)
		int numberMaterials;		// number of materials in this crack velocity field
		int numberRigidPoints;		// number of rigid particles in this velocity field
	
		// constants (not changed in MPM time step)

};

#endif
