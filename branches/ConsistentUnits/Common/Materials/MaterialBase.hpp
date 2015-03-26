/********************************************************************************
    MaterialBase.hpp
    nairn-mpm-fea
    
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
class HardeningLawBase;

#else

// C is stiffness matrix and some other things
// alpha is thermal expansion coefficients and other things
typedef struct {
	double C[5][5];
	double alpha[5];
} ElasticProperties;

#endif

#define UNSPECIFIED -1
enum { NO_PROPAGATION=0,MAXHOOPSTRESS,STEADYSTATEGROWTH,DELETED_TOTALENERGYBALANCE_USECUSTOMTASKIFNEEDED,
       STRAINENERGYDENSITY,EMPIRICALCRITERION,MAXCTODCRITERION,CRITICALERR };
enum { DEFAULT_DIRECTION=0,SELF_SIMILAR,NORMAL_TO_COD,HOOP_FROM_COD,INITIAL_DIRECTION };

class MaterialBase : public LinkedObject
{
    public:
		// variables (changed in MPM time step)
	
		// constants (not changed in MPM time step)
        char *name;
        double rho,concSaturation,betaI;
		float red,green,blue,alpha;
#ifdef FEA_CODE
		ElasticProperties pr;
#else
		int criterion[2];
        double KIc,KIIc,KIexp,KIIexp,JIc,JIIc,delIc,delIIc,nmix;
		double initSpeed;
		double initTime,maxLength;
		bool constantDirection;
		int constantTip;
		Vector growDir;
		int matPropagateDirection[2];
		int tractionMat[2];
        double diffusionCon,kCond;			// for isotropic properties
		ContactDetails *lastFriction;
	
		// class varibles
		static vector<int> fieldMatIDs;		// 0 to # material in multimatrial mode
		static vector<int> activeMatIDs;	// 0 to # active, non-rigid fields any mode
		static int incrementalDefGradTerms;
        static bool isolatedSystemAndParticles;
        static int maxPropertyBufferSize;
        static int maxAltBufferSize;
#endif
        
        // constructors and destructors
        MaterialBase();
        MaterialBase(char *);
		virtual ~MaterialBase();
		
		// initialization and verification
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
#ifdef MPM_CODE
		virtual char *InitHistoryData(void);
		virtual void FillTransportProperties(TransportProperties *);
		virtual void SetHardeningLaw(char *);
        virtual bool AcceptHardeningLaw(HardeningLawBase *,int );
		virtual void ValidateForUse(int) const;
		virtual void PrintTransportProperties(void) const;
        virtual void SetInitialParticleState(MPMBase *,int) const;
		double *CreateAndZeroDoubles(int) const;
#endif

		// initialization (base class only)
		void PrintMaterial(int) const;
		void PrintCommonProperties(void) const;
#ifdef MPM_CODE
		void PrintCriterion(int,int) const;
#endif
         
        // Methods
#ifdef MPM_CODE
		virtual void GetTransportProps(MPMBase *,int,TransportProperties *) const;
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *) const;
		virtual double GetHeatCapacity(MPMBase *) const;
		virtual double GetCpHeatCapacity(MPMBase *) const;
		virtual double GetCpMinusCv(MPMBase *) const;
        virtual void IncrementHeatEnergy(MPMBase *,double,double,double) const;
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
#else
		virtual void LoadMechanicalPropertiesFEA(int,double,int);
#endif

		// Methods (base class only)
#ifdef MPM_CODE
		void Hypo2DCalculations(MPMBase *,double,double,double,double) const;
		void Hypo3DCalculations(MPMBase *,double,double,double,double *) const;
#endif

        // base class fracture methods
#ifdef MPM_CODE
		virtual MaterialBase *SetFinalPropagate(void);
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
		virtual Vector IsotropicJToK(Vector,Vector,Vector,int,double,double);
		virtual int CriterionNeeds(int);
        virtual int ShouldPropagate(CrackSegment *,Vector &,CrackHeader *,int,int);
        virtual bool ControlCrackSpeed(CrackSegment *,double &);
		virtual bool SelectDirection(CrackSegment *,Vector &,CrackHeader *,int);
		virtual void HoopDirection(double,double,Vector *);
		void RotateDirection(Vector *,double,double);
        virtual double CrackPropagationAngleFromStrainEnergyDensityCriterion(double,double,double);
#endif

		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag(void) const = 0;
#ifdef MPM_CODE
        virtual double WaveSpeed(bool,MPMBase *) const = 0;
		virtual double ShearWaveSpeed(bool,MPMBase *) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *) const;
		virtual double MaximumDiffusion(void) const;
        virtual double MaximumDiffusivity(void) const;
        virtual double GetHistory(int,char *) const;
		virtual bool Rigid(void) const;
		virtual short RigidBC(void) const;
		virtual short RigidContact(void) const;
		virtual bool isTractionLaw(void) const;
		virtual bool isMembrane(void) const;
		virtual int KeepsCrackTip(void) const;
		virtual void SetFriction(double,int,double,double,double);
		virtual ContactDetails *GetContactToMaterial(int);
		virtual void ContactOutput(int);
        virtual bool PartitionsElasticAndPlasticStrain(void) const;
        virtual double GetArtificalViscosity(double,double) const;
		virtual bool SupportsArtificialViscosity(void) const;
		virtual int GetShareMatField(void) const;
#else
        virtual double GetStressStrainZZ(double,double,double,double,double,int);
#endif

		// accessors (base class only)
#ifdef MPM_CODE
		virtual int SetField(int,bool,int,int &);
		int GetField(void) const;
		int GetActiveField(void) const;
        virtual double GetCurrentRelativeVolume(MPMBase *) const;
        virtual Tensor GetStress(Tensor *,double) const;
#endif
	
		// class methods
		static void PrintProperty(const char *,double,const char *);
		static void PrintProperty(const char *,bool);
#ifdef MPM_CODE
		static const char *PreferredDirection(int);
		static short GetMVFIsRigid(int matfld);
		static int GetFieldMatID(int matfld);
		static int GetActiveMatID(int matfld);
#endif
		
	protected:
		// variables (changed in MPM time step)
#ifdef MPM_CODE
		double heatCapacity;			// changed if depends on particle state
        bool artificialViscosity;       // true to false for artifical viscosity
        double avA1,avA2;               // artificial viscosity coefficients
		TransportProperties tr;			// transport tensors
#endif
	
		// constants (changed in MPM time step)
		double C11,C12,C13,C22,C23,C33,C44,C55,C66;
		double CTE1,CTE2,CTE3;
#ifdef MPM_CODE
		double CME1,CME2,CME3;
		int field,activeField;
		int shareMatField;
#endif
};

// List of Materials from theMaterials[0] to theMaterials[nmat-1]
// Material IDs are particles, etc., are stored from 1 to nmat
extern MaterialBase **theMaterials;
extern int nmat;

#endif

