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

class MaterialContactNode;
class ContactLaw;

#include "Nodes/CrackVelocityField.hpp"

#define MASS_MIN 1.e-5
#define MASS_MAX 0.99999

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
		virtual double GetTotalMassAndCount(bool &);
		virtual void AddVolumeGradient(int,MPMBase *,double,double,double);
		virtual void AddVolumeGradient(int,Vector *);
		virtual void CopyMassAndMomentum(NodalPoint *);
        virtual void CopyMassAndMomentumLast(NodalPoint *);
	
		virtual void AddFtotSpreadTask3(Vector *);
		virtual void CopyGridForces(NodalPoint *);
		virtual void AddGravityAndBodyForceTask3(Vector *,double,double);
#ifdef RESTART_OPTION
        virtual bool IsTravelTooMuch(double,double) const;
#endif
		virtual void RestoreMomenta(void);
	
		virtual void UpdateMomentum(double);
		virtual void RezeroNodeTask6(double);

		// standard mathod
		virtual void MaterialContactOnCVF(MaterialContactNode *,double,int);
		void MaterialContactOnCVFLumped(MaterialContactNode *,double,int,int *,int,Vector,double);
		virtual void RigidMaterialContactOnCVF(int,bool,MaterialContactNode *mcn,double,int);
	
		// contact support methods
		Vector GetNormalVector(MaterialContactNode *,int,int,double,Vector *,double,double,bool &,Vector *);
		bool NonRigidCustomNormal(NodalPoint *,int,int,Vector &);
		bool RigidCustomNormal(NodalPoint *,int,int,Vector &);
		virtual Vector GetDisplacementVector(NodalPoint *,Vector *,double,Vector *,double,bool,double,Vector *,Vector *,Vector *,bool &);
		virtual void GetVolumeGradient(int,const NodalPoint *,Vector *,double) const;
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
		virtual bool HasPointsNonrigid(void) const;
		virtual int HasPointsThatSeeCracks(void);
		virtual void SumAndClearRigidContactForces(Vector *,bool,double,Vector *);
		virtual double GetTotalMass(bool) const;
		virtual void AddKineticEnergyAndMass(double &,double &);
		virtual Vector GetCMatMomentum(bool &,double *,Vector *,bool) const;
		virtual void ChangeCrackMomentum(Vector *,int,double);
		virtual int CopyFieldMomenta(Vector *,int);
#if ADJUST_COPIED_PK == 1
		virtual void AdjustForSymmetryBC(NodalPoint *);
#endif
		virtual int PasteFieldMomenta(Vector *,int);
		virtual int GetNumberMaterials(void);
		virtual int GetNumberNonrigidMaterials(void);
		virtual void Describe(void) const;

		// XPIC methods
		virtual void XPICSupport(int,int,NodalPoint *,double,int,int,double);
	
		// class methods
		static double GetTangentCOD(Vector *,Vector *,Vector *);
		static double GetContactArea(NodalPoint *,double,double,Vector *,double *,Vector *);
	
		virtual void MirrorFieldsThatIgnoreCracks(CrackVelocityFieldMulti *);
		virtual bool HasFieldsThatIgnoreCracks(void) const;
		virtual bool MVFInMemory(int mi) const;
    
        // machine learning
        virtual Vector LinearRegressionNormal(MaterialContactNode *,int,int,int,bool &,Vector *);
        virtual Vector LogisticRegressionNormal(MaterialContactNode *,int,int,int,bool &,Vector *);
        virtual void FindSepFromNormalAndPointCloud(Vector *,MaterialContactNode *,int,int,int,Vector *);
    
	private:
		// variables (changed in MPM time step)
		int numberMaterials;		// number of materials in this crack velocity field
		int numberRigidPoints;		// number of rigid particles in this velocity field

};

#endif
