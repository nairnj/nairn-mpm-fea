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

// Define to have BC apply only to velocity field on same side as crack, otherwise apply to
//	all velocity fields
//#define _BC_CRACK_SIDE_ONLY_

#ifdef MPM_CODE
class CrackSegment;
class MaterialInterfaceNode;
class CrackNode;
class TransportTask;
#include "Nodes/CrackVelocityField.hpp"
#endif

class NodalPoint : public LinkedObject
{
    public:
		// variables (changed in MPM time step)
		double gTemperature;		// absolute T in MPM, delta T in FEA
#ifdef MPM_CODE
		double mass;				// total mass
		CrackVelocityField **cvf;	// material velocity fields
		double gVolume;
		double gConcentration;
		double fdiff;				// diffusion
		double gMpCp;
		double fcond;				// conduction
		unsigned short fixedDirection;
	
		static double interfaceEnergy;		// total tracked each time step
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
#ifdef USE_FEXT
		void AddFextSpreadTask3(short,Vector);
#endif
		void AddTractionTask3(MPMBase *,int,Vector *);
		void CopyGridForces(NodalPoint *);
	
		void UpdateMomentaOnNode(double);
	
		void IncrementDelvaTask5(short,int,double,Vector *,Vector *) const;
	
	
		short IncrementDelvSideTask8(short,int,double,Vector *,double *,CrackSegment *);
		void SurfaceCrossesCracks(double,double,double,double,CrackField *);
		int SurfaceCrossesOneCrack(double,double,double,double,int);
		bool SurfaceCrossesOtherCrack(double,double,double,double,int);
		void CalcCMVelocityTask8(void);
		bool GetCMVelocityTask8(Vector *);
	
		// specific task methods
		void PrepareForFields(void);
        void ZeroDisp(void);
		int GetFieldForCrack(int,int,DispField **,DispField *);
        void ZeroDisp(NodalPoint *);
        void CopyUGradientStressEnergy(NodalPoint *);
        void DeleteDisp(void);
        void DeleteDisp(NodalPoint *);
		int NumberParticles(void);
		int NumberNonrigidParticles(void);
		void Describe(void) const;
		void AddDisplacement(short,int,double,Vector *);
		void AddVolume(short,int,double);
        void AddUGradient(short,double,double,double,double,double,int,double);
		void AddMatWeights(double,double *);
        void AddEnergy(short,double,double,double,double);
        void AddStress(short,double,Tensor *);
        Vector GetVelocity(short,int);
		void AddKineticEnergyAndMass(double &,double &);
		Vector GetContactForce(short,int);
		void CalcVelocityForStrainUpdate(void);
        short GetCMVelocity(Vector *);
        void CalcStrainField(void);
        void Interpolate(NodalPoint *,NodalPoint *,double,bool,int);
        void CrackContact(bool,double,CrackNode **,CrackNode **);
		void CrackContactThree(int,bool,double);
		void CrackInterfaceForce(void);
		void InterfaceForceThree(int);
		void MaterialContactOnNode(double,int,MaterialInterfaceNode **,MaterialInterfaceNode **);
        void MaterialInterfaceForce(MaterialInterfaceNode *);
		void GetVolumeGradient(short,int,Vector *,double) const;
		void GetMatVolumeGradient(int,Vector *) const;
        void SetMomVel(Vector *);
        void AddMomVel(Vector *,double);
		void ReflectMomVel(Vector *,NodalPoint *);
        void SetFtotDirection(Vector *,double,Vector *);
        void AddFtotDirection(Vector *,double,double,Vector *);
		void ReflectFtotDirection(Vector *,double,NodalPoint *,Vector *);
        Vector GetCMatMomentum(void);
		void SetFixedDirection(int);
		void UnsetFixedDirection(int);
		void CalcTotalMassAndCount(void);
#ifdef COMBINE_RIGID_MATERIALS
		void CopyRigidParticleField(void);
#endif
		Vector GetTotalContactForce(bool);
#else
        void InitForceField(void);
        void PrintAvgStress(void);
#endif

		// class methods
#ifdef MPM_CODE
		static void PreliminaryCalcs(void);
		static void GetGridVelocitiesForStrainUpdate(void);
		static void GetGridCMVelocitiesTask8(void);
#endif
    
    private:
		
#ifdef MPM_CODE
        //methods - MPM only
		void AverageStrain(DispField *,DispField *,DispField *,double);
		int WeightAverageStrain(int,int,DispField *);
    int WeightAverageStrain(int fld1,int fld2,int fld3,DispField *dest);
        void AdjustContact(short,short,Vector *,int,bool,double);
		void AddInterfaceForce(short,short,Vector *,int);
#endif

};

// List of nodes stored as nd[1] to nd[nnodes]
extern NodalPoint **nd;
extern int nnodes;

#endif
