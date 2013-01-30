/********************************************************************************
    MaterialBase.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		none
********************************************************************************/

#ifndef _MATERIAL_BASE_

#define _MATERIAL_BASE_

// prototypes
#ifdef MPM_CODE
	class CrackHeader;
	class CrackSegment;
	class DiffusionTask;
	class ConductionTask;
	class MPMBase;
#endif

#define UNSPECIFIED -1
enum { NO_PROPAGATION=0,MAXHOOPSTRESS,STEADYSTATEGROWTH,TOTALENERGYBALANCE,
       STRAINENERGYDENSITY,EMPIRICALCRITERION,MAXCTODCRITERION,CRITICALERR };
enum { DEFAULT_DIRECTION=0,SELF_SIMILAR,NORMAL_TO_COD,HOOP_FROM_COD,INITIAL_DIRECTION };

class MaterialBase : public LinkedObject
{
    public:
		// variables (changed in MPM time step)
	
		// constants (not changed in MPM time step)
        char *name;
        double rho,concSaturation,betaI;
		float red,green,blue;
#ifdef FEA_CODE
        double mdm[5][5];
        double me0[5];
#else
		int criterion[2];
        double KIc,KIIc,KIexp,KIIexp,JIc,JIIc,gamma,delIc,delIIc,nmix;
		Tensor diffusionTensor;
		Tensor kCondTensor;                 // divided by rho (in g/mm^3)
		double pCrit3,gain,initSpeed;
		double initTime,maxLength;
		int constantDirection;
		Vector growDir;
		int matPropagateDirection[2];
		int tractionMat[2];
        double diffusionCon,kCond;			// for isotropic properties
		ContactDetails *lastFriction;
		static vector<int> fieldMatIDs;
		static int incrementalDefGradTerms;
#endif
        
        // constructors and destructors
        MaterialBase();
        MaterialBase(char *);
		virtual ~MaterialBase();
		
		// initialization and verification
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		static void PrintProperty(const char *,double,const char *);
		static void PrintProperty(const char *,bool);
#ifdef MPM_CODE
		virtual void ValidateForUse(int);
        virtual char *InitHistoryData(void);
		virtual void InitialLoadTransProps(void);
		virtual void PrintTransportProperties(void);
        virtual void SetInitialParticleState(MPMBase *,int);
        virtual void SetHardeningLaw(char *);
		double *CreateAndZeroDoubles(int);
#endif

		// initialization (base class only)
        void PrintMaterial(int);
		void PrintCommonProperties(void);
#ifdef MPM_CODE
		void PrintCriterion(int,int);
#endif
         
        // Methods
        virtual void LoadMechProps(int,double,int);
#ifdef MPM_CODE
        virtual void LoadMechanicalProps(MPMBase *,int);
		virtual void LoadTransportProps(MPMBase *,int);
		virtual double GetHeatCapacity(MPMBase *);
		virtual double GetHeatCapacityVol(MPMBase *);
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
#endif

		// Methods (base class only)
#ifdef MPM_CODE
		void Hypo2DCalculations(MPMBase *,double,double,double,double);
		void Hypo3DCalculations(MPMBase *,double,double,double,double *);
#endif

        // base class fracture methods
#ifdef MPM_CODE
		virtual MaterialBase *SetFinalPropagate(void);
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
		virtual int CriterionNeeds(int);
        virtual int ShouldPropagate(CrackSegment *,Vector &,CrackHeader *,int,int);
        virtual bool ControlCrackSpeed(CrackSegment *,double &);
		virtual bool SelectDirection(CrackSegment *,Vector &,CrackHeader *,int);
		virtual void HoopDirection(double,double,Vector *);
		void RotateDirection(Vector *,double,double);
        virtual double CrackPropagationAngleFromStrainEnergyDensityCriterion(double,double,double);
		static const char *PreferredDirection(int);
#endif

		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag(void) = 0;
#ifdef MPM_CODE
        virtual double WaveSpeed(bool,MPMBase *) = 0;
		virtual double ShearWaveSpeed(bool,MPMBase *);
        virtual double CurrentWaveSpeed(bool,MPMBase *);
		virtual double MaximumDiffusion(void);
        virtual double MaximumDiffusivity(void);
        virtual double GetHistory(int,char *);
		virtual bool Rigid(void);
		virtual short RigidBC(void);
		virtual short RigidContact(void);
		virtual bool isTractionLaw(void);
		virtual void SetFriction(double,int,double,double,double);
		virtual ContactDetails *GetContactToMaterial(int);
		virtual void ContactOutput(int);
        virtual bool HasPlasticStrainForGradient(void);
#else
        virtual double GetStressStrainZZ(double,double,double,double,double,int);
#endif

		// accessors (base class only)
#ifdef MPM_CODE
		virtual int SetField(int,bool,int);
		static short GetMVFIsRigid(int matfld);
		static int GetFieldMatID(int matfld);
		int GetField(void);
 		Tensor *GetkCondTensor(void);
		Tensor *GetDiffusionTensor(void);
        virtual double GetCurrentRelativeVolume(MPMBase *);
        virtual Tensor GetStress(Tensor *sp,double pressure);
#endif
		
	protected:
		// variables (changed in MPM time step)
#ifdef MPM_CODE
		double heatCapacity;			// changed if depends on particle state
		double heatCapacityVol;			// changed if depends on particle state
#endif
	
		// constants (changed in MPM time step)
        double lastMatAngle;
		double C11,C12,C13,C22,C23,C33,C44,C55,C66;
		double CTE1,CTE2,CTE3;
#ifdef MPM_CODE
		double CME1,CME2,CME3;
		int field;
#endif
};

// List of Materials from theMaterials[0] to theMaterials[nmat-1]
// Material IDs are particles, etc., are stored from 1 to nmat
extern MaterialBase **theMaterials;
extern int nmat;

#endif

