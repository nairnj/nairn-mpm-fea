/********************************************************************************
	CrackVelocityField.hpp
	NairnMPM
 
	Created by John Nairn on 11 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
 
	Dependencies
		MatVelocityField.hpp
********************************************************************************/

#ifndef _CRACKVELOCITYFIELD_

#define _CRACKVELOCITYFIELD_

#define FIRST_CRACK 0
#define SECOND_CRACK 1

#include "Nodes/MatVelocityField.hpp"

class MPMBase;
class NodalPoint;

class CrackVelocityField
{
	public:
		// variables (changed in MPM time step)
		short loc[2];				// crack location's
		int crackNum[2];			// crack number's/crackNum[0] used for total particles when moving crack planes by cm
		Vector norm[2];				// crack normal's/norm[0] used for velocity when moving crack planes by cm
		DispField *df;				// For J and K calculations
		
		// constants (not changed in MPM time step)
	
		// constructors and destructors
        CrackVelocityField(short,int);
        virtual ~CrackVelocityField();
		virtual void Zero(short,int,bool);
		virtual void ZeroMatFields() = 0;
		
		// specific task methods
		void AddMomentumTask1(int,Vector *,Vector *);
		virtual void AddMass(int,double);
		virtual void AddMassTask1(int,double);
		virtual double GetTotalMassAndCount(void) = 0;
		virtual void AddVolumeGradient(int,MPMBase *,double,double,double);
	
		void AddFintTask3(int,Vector *);
		virtual void AddFintSpreadTask3(Vector *) = 0;
		void AddFextTask3(int,Vector *);
		virtual void AddFextSpreadTask3(Vector *) = 0;
		virtual void CalcFtotTask3(double) = 0;
	
		virtual void UpdateMomentaOnField(double) = 0;
	
		void IncrementDelvaTask5(int,double,Vector *,Vector *);
	
		virtual void RezeroNodeTask6(double) = 0;
		void AddMomentumTask6(int,double,Vector *);
	
		void CreateStrainField(void);
		void DeleteStrainField(void);
		
		short IncrementDelvTask8(double,Vector *,double *);
		int CollectMomentaTask8(Vector *);
		void SetCMVelocityTask8(Vector *,int);
		bool GetCMVelocityTask8(Vector *);
	
		void AddNormals(Vector *,int);
		void AddDisplacement(int,double,Vector *);
		void AddVolume(int,double);
	
		// methods
		virtual void MaterialContact(int,int,bool,double);
		virtual void GetVolumeGradient(int,NodalPoint *,Vector *,double);
		virtual void CalcVelocityForStrainUpdate(void) = 0;
	
		// boundary conditions
        virtual void SetMomVel(int) = 0;
        virtual void AddMomVel(int,double) = 0;
        virtual void SetFtot(int,double) = 0;
        virtual void AddFtot(int,double,double) = 0;
	
		// accessors
		short location(int);
		int crackNumber(int);
		int OppositeCrackTo(int,int);
		void SetLocationAndCrack(short,int,int);
		virtual double GetTotalMass(void) = 0;
		virtual void AddKineticEnergyAndMass(double &,double &) = 0;
		virtual double GetVolumeNonrigid(void) = 0;
		virtual double GetVolumeTotal(double) = 0;
		virtual Vector GetCMatMomentum(void) = 0;
		virtual Vector GetCMDisplacement(void) = 0;
		virtual Vector GetCMatFtot(void) = 0;
		virtual void ChangeMomentum(Vector *,bool,double) = 0;
		virtual int CopyFieldMomenta(Vector *,int) = 0;
		virtual int PasteFieldMomenta(Vector *,int) = 0;
		Vector GetVelocity(int);
		Vector GetContactForce(int);
		virtual int GetNumberPoints(void);
		virtual int GetNumberPointsNonrigid(void);
		virtual void Describe(void);
		virtual void SumAndClearRigidContactForces(Vector *,bool);
	
		// class methods
		static bool ActiveField(CrackVelocityField *);
		static bool ActiveNonrigidField(CrackVelocityField *cvf);
		static CrackVelocityField *CreateCrackVelocityField(short,int);
	
	protected:
		// variables (changed in MPM time step)
		int numberPoints;			// total number of materials points in this field/field [0] changed to sum of all in task 8
		MatVelocityField **mvf;		// material velocity fields
		// unscaled nonrigid volume (ignores dilation) only used for imperfect interface forces and material contact
		// unscaleRigidVolume is due to rigid contaft materials (type 8) (always zero unless multimaterial mode)
};

#endif
