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

#ifdef MPM_CODE

#include "System/DataTypes.hpp"

class CrackSegment;
class CrackNode;
class TransportTask;
class MaterialContactNode;
#include "Nodes/CrackVelocityField.hpp"
#include "Nodes/MatVelocityField.hpp"
#include "Cracks/CrackHeader.hpp"

#endif

class NodalPoint : public LinkedObject
{
    public:
		// variables
#ifdef MPM_CODE
		CrackVelocityField **cvf;	// crack velocity fields
		TransportField gCond;		// conduction
		//double gTemperature,gMpCp,fcond;		// conduction
		TransportField *gDiff;		// diffusion(s)
		//double gConcentration,gVolume,fdiff;	// diffusion
		unsigned short fixedDirection;
		MaterialContactNode *contactData;
	
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
		NodalPoint(int,double,double);
		NodalPoint(int,double,double,double);
		void CreateNodalPoint(int,double,double,double);
#ifdef MPM_CODE
		NodalPoint(NodalPoint *);
#endif
		virtual ~NodalPoint();
        
        // methods - general and abstract
        virtual void PrintNodalPoint(ostream &);

#ifdef MPM_CODE
        // methods - MPM only
		void InitializeForTimeStep();
		void CopyFieldInitialization(NodalPoint *);
		void UseTheseFields(CrackVelocityField **);
		short AddCrackVelocityField(int,CrackField *);
		void AddMatVelocityField(short,int);
        bool NeedsMatVelocityField(short,int) const;
		void CreateDiffusionVariables(void);

		void AddMassMomentum(MPMBase *,short,int,double,double,double,double,int,bool);
        void AddMassMomentumLast(MPMBase *,short,int,double,double,double,double);
		void AddMomentumTask1(short,int,double,Vector *,int);
		void AddMass(short,int,double);
		void AddMassTask1(short,int,double,int);
		void CopyMassAndMomentum(NodalPoint *);
        void CopyMassAndMomentumLast(NodalPoint *);
        void RezeroNodeTask6(double);
        void AddMomentumTask6(short,int,double,Vector *);
	
		void AddFtotTask3(short,int,Vector *);
		void AddFtotSpreadTask3(short,Vector);
		bool AddTractionTask3(MPMBase *,short,int,Vector *);
		void AddGravityAndBodyForceTask3(Vector *);
 #ifdef RESTART_OPTION
        bool IsTravelTooMuch(double,double) const;
 #endif
		void CopyGridForces(NodalPoint *);
		void UpdateMomentum(double);
		void IncrementDelvaTask5(short,int,double,GridToParticleExtrap *) const;

		bool IncrementDelvSideTask8(short,int,double,Vector *,Vector *,double *,CrackSegment *) const;
        bool IncrementDispField(short,int,double,DispField *,CrackSegment *) const;
		bool GetCMVelocityTask8(Vector *,Vector *) const;
		short GetFieldForSurfaceParticle(short,int,CrackSegment *,bool) const;
		void SurfaceCrossesCracks(Vector *,Vector *,CrackField *) const;
		int SurfaceCrossesOneCrack(Vector *,Vector *,int) const;
		int SurfaceCrossesOtherCrack(Vector *,Vector *,int) const;
		bool GetCenterOfMassVelocity(Vector *,bool);
	
		// specific task methods
		void PrepareForFields(void);
        void ZeroDisp(void);
		int GetFieldForCrack(bool,bool,DispField **,int,int &);
        void ZeroDisp(NodalPoint *);
        void CopyUGradientStressEnergy(NodalPoint *);
        void DeleteDisp(void);
        void DeleteDisp(NodalPoint *);
		int NumberParticles(void);
		bool NeedsParticleListed(short vfld);
		int NumberNonrigidCracks(void);
		bool NodeHasNonrigidParticles(void) const;
		bool NodeHasParticles(void) const;
		void Describe(bool) const;
		void AddContactTerms(short,int,ContactTerms *);
        void AddUGradient(short,double,double,double,double,double,int,double);
		void AddMatWeights(double,double *);
		// GRID_JTERMS
        void AddGridVelocity(short,double,double,double);
        void AddEnergy(short,double,double,double,double);
        void AddStress(short,double,Tensor *);
        Vector GetVelocity(short,int);
		void AddKineticEnergyAndMass(double &,double &);
		void GridValueCalculation(int);
        short GetCMVelocity(Vector *);
        void CalcStrainField(void);
		void Interpolate(NodalPoint *,NodalPoint *,double,int,int *);
		int HasCrackContact(void);
		void CrackContact(int,double,int);
		void CrackContactThree(int,int,double);
		void MaterialContactOnNode(double,int,MaterialContactNode *);
        void GetMatVolumeGradient(int,Vector *) const;
        void ZeroVelocityBC(Vector *,int,double,Vector *);
        void AddVelocityBC(Vector *,double,int,double,Vector *);
		void ReflectVelocityBC(Vector *,NodalPoint *,double,double,int,double,Vector *);
		void SetFixedDirection(int);
		void UnsetFixedDirection(int);
		bool CalcTotalMassAndCount(void);
		void RestoreMomenta(void);
		void AddGetContactForce(bool,Vector *,double,Vector *);
		void AddRigidBCInfo(MPMBase *,double,int,Vector *);
		int ReadAndZeroRigidBCInfo(Vector *,double *,double *);
		void MirrorIgnoredCrackFields(void);
	
		// XPIC
		void XPICSupport(int,int,NodalPoint *,double,int,int,double);
		void AddVStarNext(short,int,Vector *,Vector *,Vector *,Matrix3 *,double,double);
		virtual Vector *GetVStarPrev(short,int) const;
		double GetMaterialMass(short,int) const;
	
#else // not MPM_CODE
		// FEA code
        void InitForceField(void);
        void PrintAvgStress(void);
#endif // end MPM_CODE
	
		// class methods
#ifdef MPM_CODE
		static void PrepareNodeCrackFields(void);
		static NodalPoint *CreateGhostFromReal(NodalPoint *);
#endif
		static NodalPoint *Create2DNode(int,double,double);
		static NodalPoint *Create3DNode(int,double,double,double);
	
	protected:
#ifdef MPM_CODE
		double nodalMass;				// total mass
		bool hasParticles;				// true if node has particles
#endif
    
    private:
#ifdef MPM_CODE
        // methods - MPM only
		void AverageStrain(DispField *,DispField *,DispField *,double);
        void AdjustContact(short,short,Vector *,int,int,double);
#endif

};

// List of nodes stored as nd[1] to nd[nnodes]
extern NodalPoint **nd;
extern int nnodes;
extern int *nda;

#endif
