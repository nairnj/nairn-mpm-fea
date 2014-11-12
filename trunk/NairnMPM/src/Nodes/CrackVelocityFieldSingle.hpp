/********************************************************************************
	CrackVelocityFieldSingle.hpp
	nairn-mpm-fea
 
	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
 
	Dependencies
		CrackVelocityField.hpp,MatVelocityField.hpp
********************************************************************************/

#ifndef _CRACKVELOCITYFIELDSINGLE_

#define _CRACKVELOCITYFIELDSINGLE_

#include "Nodes/CrackVelocityField.hpp"

class CrackVelocityFieldSingle : public CrackVelocityField
{
	public:
		
        // constructors and destructors
		CrackVelocityFieldSingle(short,int);
		virtual void ZeroMatFields(void);
	
		// specific task methods
		virtual double GetTotalMassAndCount(void);
	
		virtual void AddFtotSpreadTask3(Vector *);
		virtual void UpdateMomentaOnField(double);
	
		virtual void RezeroNodeTask6(double);
	
		virtual void CalcVelocityForStrainUpdate(void);
	
		// boundary conditions
        virtual void SetMomVel(Vector *);
        virtual void AddMomVel(Vector *,double);
		virtual void ReflectMomVel(Vector *,CrackVelocityField *);
        virtual void SetFtotDirection(Vector *,double,Vector *);
        virtual void AddFtotDirection(Vector *,double,double,Vector *);
		virtual void ReflectFtotDirection(Vector *,double,CrackVelocityField *,Vector *);
	
		// accessors
		virtual double GetTotalMass(void) const;
		virtual void AddKineticEnergyAndMass(double &,double &);
		virtual double GetVolumeNonrigid(void);
		virtual double GetVolumeTotal(NodalPoint *) const;
		virtual Vector GetCMatMomentum(void) const;
		virtual Vector GetCMDisplacement(NodalPoint *) const;
		virtual Vector GetCMatFtot(void);
		virtual void ChangeMomentum(Vector *,bool,double);
		virtual int CopyFieldMomenta(Vector *,int);
		virtual int PasteFieldMomenta(Vector *,int);
		virtual void Describe(void) const;
	
};

#endif
