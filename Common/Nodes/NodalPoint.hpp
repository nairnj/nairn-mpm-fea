/********************************************************************************
    NodalPoint.hpp
    NairnMPM
    
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
		double gRhoVCp;
		double fcond;				// conduction
		unsigned char fixedDirection;
	
		static double interfaceEnergy;		// total tracked each time step
#endif
	
		// constants (not changed in MPM time step)
        double x,y,z;
        long num;
#ifdef MPM_CODE
        char above,below;
#endif

#ifdef FEA_CODE
        ForceField *fs;
#endif
       
        // constructors and destructors
        NodalPoint(long);
        virtual ~NodalPoint();
        
        // methods - general and abstract
        virtual void PrintNodalPoint(ostream &) = 0;

#ifdef MPM_CODE
        // methods - MPM only
		void InitializeForTimeStep();
	
		short AddMomentumTask1(int,CrackField *,double,Vector *);
		void AddMass(short,int,double);
		void AddMassTask1(short,int);
		void AddMassGradient(short,int,double,double,double,double);
	
		void AddFintTask3(short,int,Vector *);
		void AddFintSpreadTask3(short,Vector);
		void AddFextTask3(short,int,Vector *);
		void AddFextSpreadTask3(short,Vector);
		void CalcFtotTask3(double);
	
		void UpdateMomentaTask4(double);
	
		void IncrementDelvaTask5(short,int,double,Vector *,Vector *);
	
		void RezeroNodeTask6(double);
		void AddMomentumTask6(short,int,double,Vector *);
	
		short IncrementDelvSideTask8(short,int,double,Vector *,double *,CrackSegment *);
		void SurfaceCrossesCracks(double,double,double,double,CrackField *);
		int SurfaceCrossesOneCrack(double,double,double,double,int);
		bool SurfaceCrossesOtherCrack(double,double,double,double,int);
		void CalcCMVelocityTask8(void);
		bool GetCMVelocityTask8(Vector *);
	
		// specific task methods
		void PrepareForFields(void);
        void ZeroDisp(void);
        void DeleteDisp(void);
		int NumberParticles(void);
		int NumberNonrigidParticles(void);
		void Describe(void);
		void AddDisplacement(short,int,double,Vector *);
		void AddUnscaledVolume(short,double);
        void AddUGradient(short,double,double,double,double,double);
        void AddEnergy(short,double,double,double,double);
        void AddStress(short,double,Tensor *);
        Vector GetVelocity(short,int);
		void CalcVelocityForStrainUpdate(void);
        short GetCMVelocity(Vector *);
        void CalcStrainField(void);
        void Interpolate(NodalPoint *,NodalPoint *,double,bool);
        void CrackContact(int,bool,double);
		void CrackContactThree(int,bool,double);
		void InterfaceForce(void);
		void InterfaceForceThree(int);
		void MaterialContactOnNode(bool,double);
		void GetMassGradient(short,int,Vector *,double);
		double GetMass(short,int);
		void SetXMomVel(void);
		void SetYMomVel(void);
		void SetZMomVel(void);
		void SetSkewMomVel(double);
        void AddXMomVel(double);
        void AddYMomVel(double);
        void AddZMomVel(double);
        void AddSkewMomVel(double,double);
        void SetXFtot(double);
        void SetYFtot(double);
        void SetZFtot(double);
        void SetSkewFtot(double,double);
        void AddXFtot(double,double);
        void AddYFtot(double,double);
        void AddZFtot(double,double);
        void AddSkewFtot(double,double,double);
		void SetFixedDirection(int);
		void UnsetFixedDirection(int);
		void CalcTotalMassAndCount(void);
		void CombineRigidParticles(void);
#else
        void InitForceField(void);
        void PrintAvgStress(void);
#endif

		// class methods
#ifdef MPM_CODE
		static void PreliminaryCalcs(void);
		static void CombineRigidMaterials(void);
		static void MaterialContact(bool,bool,double);
		static void GetGridVelocitiesForStrainUpdate(void);
		static void GetGridCMVelocitiesTask8(void);
#endif
    
    private:
		
#ifdef MPM_CODE
        //methods - MPM only
		void AverageStrain(DispField *,DispField *,DispField *,double);
        void AdjustContact(short,short,Vector *,int,bool,double);
		void AddInterfaceForce(short,short,Vector *,int);
#endif

};

// List of nodes stored as nd[1] to nd[nnodes]
extern NodalPoint **nd;
extern long nnodes;

#endif
