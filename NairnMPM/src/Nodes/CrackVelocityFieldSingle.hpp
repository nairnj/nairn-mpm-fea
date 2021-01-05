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
		CrackVelocityFieldSingle(int,short,int);
        virtual ~CrackVelocityFieldSingle();
		virtual void ZeroMatFields(void);
	
		// specific task methods
		virtual double GetTotalMassAndCount(bool &);
	
		virtual void AddFtotSpreadTask3(Vector *);
		virtual void AddGravityAndBodyForceTask3(Vector *,double,double);
#ifdef RESTART_OPTION
        virtual bool IsTravelTooMuch(double,double) const;
#endif
		virtual void RestoreMomenta(void);
	
		virtual void UpdateMomentum(double);
		virtual void RezeroNodeTask6(double);
	
		virtual void GridValueCalculation(int);
	
		// boundary conditions
        virtual void ZeroVelocityBC(Vector *,int,double,Vector *);
        virtual void AddVelocityBC(Vector *,double,int,double,Vector *);
		virtual void ReflectVelocityBC(Vector *,CrackVelocityField *,double,double,int,double,Vector *);
	
		// contact accessors
		virtual double GetContactVolumeNonrigid(bool) const;
		virtual Vector GetCMDisplacement(NodalPoint *,bool,bool) const;
		virtual Vector GetCMatFtot(void);
	
		// accessors
		virtual double GetTotalMass(bool) const;
		virtual void AddKineticEnergyAndMass(double &,double &);
		virtual Vector GetCMatMomentum(bool &,double *,Vector *,bool) const;
		virtual void ChangeCrackMomentum(Vector *,int,double);
		virtual int CopyFieldMomenta(Vector *,int);
#if ADJUST_COPIED_PK == 1
		virtual void AdjustForSymmetryBC(NodalPoint *);
#endif
		virtual int PasteFieldMomenta(Vector *,int);
		virtual void Describe(void) const;
	
		// XPIC
		virtual void XPICSupport(int,int,NodalPoint *,double,int,int,double);
};

#endif
