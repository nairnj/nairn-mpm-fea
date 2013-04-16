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
		virtual void AddMatVelocityField(int);
		virtual void MatchMatVelocityFields(MatVelocityField **);
	
		// specific task methods
		virtual void AddMassTask1(int,double,int);
		virtual double GetTotalMassAndCount(void);
		virtual void AddVolumeGradient(int,MPMBase *,double,double,double);
		virtual void CopyVolumeGradient(int,Vector *);
		virtual void CombineRigidFrom(CrackVelocityFieldMulti *,int);
		virtual void CopyRigidFrom(CrackVelocityFieldMulti *,int);
		virtual void CopyMassAndMomentum(NodalPoint *,int);
	
		virtual void AddFtotSpreadTask3(Vector *);
		virtual void CopyGridForces(NodalPoint *,int);
		virtual void AddGridDampingTask3(double);
	
		virtual void UpdateMomentaOnField(double);
	
		virtual void RezeroNodeTask6(double);
	
		virtual void MaterialContactOnCVF(NodalPoint *,int,double,bool,MaterialInterfaceNode **,MaterialInterfaceNode **);
        virtual void GetFrictionalDeltaMomentum(Vector *,Vector *,double,double);
		virtual void GetVolumeGradient(int,const NodalPoint *,Vector *,double) const;
		virtual void RigidMaterialContactOnCVF(int,NodalPoint *,int,double,bool,MaterialInterfaceNode **,MaterialInterfaceNode **);
		virtual bool GetInterfaceForcesForNode(Vector *,Vector *,double,double,
										   double,Vector *,double *,double,Vector *,double,bool,bool,double);
		virtual void CalcVelocityForStrainUpdate(void);
	
		// boundary conditions
        virtual void SetMomVel(Vector *);
        virtual void AddMomVel(Vector *,double);
        virtual void SetFtotDirection(Vector *,double);
        virtual void AddFtotDirection(Vector *,double,double);
	
		// accessors
		virtual int GetNumberPointsNonrigid(void);
		virtual void SumAndClearRigidContactForces(Vector *,bool);
		virtual double GetTotalMass(void);
		virtual void AddKineticEnergyAndMass(double &,double &);
		virtual double GetVolumeNonrigid(void);
		virtual double GetVolumeTotal(double) const;
		virtual Vector GetCMatMomentum(void);
		virtual Vector GetCMDisplacement(void) const;
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
