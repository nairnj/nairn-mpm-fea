/********************************************************************************
    NodalPoint.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _NODALPOINT_

#define _NODALPOINT_

// To do imperfect interface force after mass and momentum (OK to first order) define this term
// Undefine to do after momentum update, which should be second order
// Make sure all code includes NodalPoint.hpp
//#define MANDMIMPINT

#ifdef MPM_CODE

class CrackSegment;
class CrackNode;
class TransportTask;
#include "Nodes/CrackVelocityField.hpp"
#include "Cracks/CrackHeader.hpp"

#endif

class NodalPoint : public LinkedObject
{
    public:
		// variables (changed in MPM time step)
#ifdef MPM_CODE
		CrackVelocityField **cvf;	// crack velocity fields
		TransportField gCond;		// conduction
		//double gTemperature,gMpCp,fcond;		// conduction
		TransportField gDiff;		// diffusion
		//double gConcentration,gVolume,fdiff;	// diffusion
		unsigned short fixedDirection;
	
		static double interfaceEnergy;		// total tracked each time step
		static double frictionWork;			// cumulative friction work
#else
		double gTemperature;		// delta T in FEA
#endif
	
		// constants (not changed in MPM time step)
        double x,y,z;
        int	num;

#ifdef FEA_CODE
        ForceField *fs;
#endif
       
        // constructors and destructors
        NodalPoint(int);
#ifdef MPM_CODE
		NodalPoint(NodalPoint *);
#endif
		virtual ~NodalPoint();
        
        // methods - general and abstract
        virtual void PrintNodalPoint(ostream &) = 0;

#ifdef MPM_CODE
        // methods - MPM only
		void InitializeForTimeStep();
		void CopyFieldInitialization(NodalPoint *);
		void UseTheseFields(CrackVelocityField **);
		short AddCrackVelocityField(int,CrackField *);
		void AddMatVelocityField(short,int);
        bool NeedsMatVelocityField(short,int) const;
	
		void AddMassMomentum(MPMBase *,short,int,double,double,double,double,int,bool);
        void AddMassMomentumLast(MPMBase *,short,int,double,double,double,double);
		void AddMomentumTask1(short,int,double,Vector *,int);
		void AddMass(short,int,double);
		void AddMassTask1(short,int,double,int);
		void CopyVolumeGradient(short,int,Vector *);
		void CopyMassAndMomentum(NodalPoint *);
        void CopyMassAndMomentumLast(NodalPoint *);
        void RezeroNodeTask6(double);
        void AddMomentumTask6(short,int,double,Vector *);
	
		void AddFtotTask3(short,int,Vector *);
		void AddFtotSpreadTask3(short,Vector);
		void AddTractionTask3(MPMBase *,short,int,Vector *);
		void AddGravityAndBodyForceTask3(Vector *);
		void CopyGridForces(NodalPoint *);
		void UpdateMomentaOnNode(double);
		void IncrementDelvaTask5(short,int,double,GridToParticleExtrap *) const;
	
		bool IncrementDelvSideTask8(short,int,double,Vector *,double *,CrackSegment *) const;
		short GetFieldForSurfaceParticle(short,int,CrackSegment *) const;
		void SurfaceCrossesCracks(double,double,double,double,CrackField *) const;
		int SurfaceCrossesOneCrack(double,double,double,double,int) const;
		int SurfaceCrossesOtherCrack(double,double,double,double,int) const;
		void CalcCMVelocityTask8(void);
		bool GetCMVelocityTask8(Vector *) const;
	
		// specific task methods
		void PrepareForFields(void);
        void ZeroDisp(void);
		int GetFieldForCrack(bool,bool,DispField **,int);
        void ZeroDisp(NodalPoint *);
        void CopyUGradientStressEnergy(NodalPoint *);
        void DeleteDisp(void);
        void DeleteDisp(NodalPoint *);
		int NumberParticles(void);
		int NumberNonrigidCracks(void);
		bool NodeHasNonrigidParticles(void) const;
		double GetNodalMass(bool) const;
		void Describe(void) const;
		void AddDisplacement(short,int,double,Vector *);
		void AddVolume(short,int,double);
        void AddUGradient(short,double,double,double,double,double,int,double);
		void AddMatWeights(double,double *);
		// GRID_JTERMS
        void AddGridVelocity(short,double,double,double);
        void AddEnergy(short,double,double,double,double);
        void AddStress(short,double,Tensor *);
        Vector GetVelocity(short,int);
		void AddKineticEnergyAndMass(double &,double &);
		void CalcVelocityForStrainUpdate(void);
        short GetCMVelocity(Vector *);
        void CalcStrainField(void);
		void Interpolate(NodalPoint *,NodalPoint *,double,int);
        void CrackContact(int,double,CrackNode **,CrackNode **);
		void AdjustDelPiForBCs(Vector *) const;
		void CrackContactThree(int,int,double);
		void MaterialContactOnNode(double,int);
        void GetMatVolumeGradient(int,Vector *) const;
        void SetMomVel(Vector *);
        void AddMomVel(Vector *,double);
		void ReflectMomVel(Vector *,NodalPoint *,double);
        void SetFtotDirection(Vector *,double,Vector *);
        void AddFtotDirection(Vector *,double,double,Vector *);
		void ReflectFtotDirection(Vector *,double,NodalPoint *,double,Vector *);
		void SetFixedDirection(int);
		void UnsetFixedDirection(int);
		void CalcTotalMassAndCount(void);
		void RestoreMomenta(void);
		void AddGetContactForce(bool,Vector *,double,Vector *);
		void AddRigidBCInfo(MPMBase *,double,int,Vector *);
		int ReadAndZeroRigidBCInfo(Vector *,double *,double *);
#else
        void InitForceField(void);
        void PrintAvgStress(void);
#endif
	
		// NairnMPM methods
#ifdef MPM_CODE
		void CopyRigidParticleField(void);
#endif

		// class methods
#ifdef MPM_CODE
		static void PreliminaryCalcs(void);
		static void GetGridVelocitiesForStrainUpdate(void);
		static void GetGridCMVelocitiesTask8(void);
	
	protected:
		double nodalMass;				// total mass
#endif
    
    private:
		
#ifdef MPM_CODE
        //methods - MPM only
		void AverageStrain(DispField *,DispField *,DispField *,double);
        void AdjustContact(short,short,Vector *,int,int,double);
#endif

};

// List of nodes stored as nd[1] to nd[nnodes]
extern NodalPoint **nd;
extern int nnodes;

#endif
