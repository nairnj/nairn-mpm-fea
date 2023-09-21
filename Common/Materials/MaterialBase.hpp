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

// softening
class FailureSurface;
class SofteningLaw;
class InitialCondition;
enum { VARY_STRENGTH=1,VARY_TOUGHNESS,VARY_STRENGTH_AND_TOUGHNESS };
enum { SOFTDIST_NONE=0,SOFTDIST_NORMAL,SOFTDIST_WEIBULL };

// softening history variables
enum { SOFT_DAMAGE_STATE=0,DELTANORMAL,DELTASHEAR,DELTASHEAR2,DAMAGENORMAL,DAMAGESHEAR,DAMAGESHEAR2,
	NORMALDIR1,NORMALDIR2,NORMALDIR3,GCSCALING,RELATIVE_STRENGTH,RELATIVE_TOUGHNESS,
	SOFT_NUMBER_HISTORY };
enum { CUBOID_SURFACE=0,CYLINDER_SURFACE,OVOID_SURFACE,COUPLED_CUBOID_SURFACE};

enum { NO_DIFFUSION=0,MOISTURE_DIFFUSION,POROELASTICITY_DIFFUSION,
    FRACTURE_PHASE_FIELD,BATTERY_PHASE_FIELD,CONDUCTION_PHASE_FIELD,NUM_DUFFUSION_OPTIONS };

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
enum { SOLID_MAT=0,MEMBRANE_MAT,TRACTION_MAT,CONTACT_MAT };
enum { CURRENT_CONFIGURATION=0,INITIAL_CONFIGURATION,INITIALMATERIAL};

// NOTHING:							alt strain is not used
// ENG_BIOT_PLASTIC_STRAIN:			small strain plasticity it biot plastic strain in alt strain to total strain
//									in strain
// LEFT_CAUCHY_ELASTIC_B_STRAIN:	hyperelastic plastic, Alt strain has elastic B. Need to calculate
//									plastic strain
// LEFT_CAUCHY_TOTAL_B_STRAIN:		hyperelastic with B = FF^T (and no plastic strain
// MEMBRANE_DEFORMATION:			membrane deformation
enum { NOTHING,ENG_BIOT_PLASTIC_STRAIN,LEFT_CAUCHY_ELASTIC_B_STRAIN,LEFT_CAUCHY_TOTAL_B_STRAIN,MEMBRANE_DEFORMATION};

class MaterialBase : public LinkedObject
{
    public:
		// variables
        char *name;
 		float red,green,blue,alpha;
#ifdef FEA_CODE
		ElasticProperties pr;
#else
		// crack settings
		int criterion[2];
        double KIc,KIIc,KIexp,KIIexp,JIc,JIIc,delIc,delIIc,nmix;
		double initSpeed;
		double initTime,maxLength;
		bool constantDirection;
		int constantTip;
		Vector growDir;
		int matPropagateDirection[2];
		int tractionMat[2];
		
		// class varibles
		static vector<int> fieldMatIDs;		// 0 to # material in multimatrial mode
		static vector<int> activeMatIDs;	// 0 to # active, non-rigid fields any mode
		static int incrementalDefGradTerms;
        static bool isolatedSystemAndParticles;
        static int maxPropertyBufferSize;
        static int maxAltBufferSize;
		static bool extrapolateRigidBCs;
#endif
        
        // constructors and destructors
        MaterialBase();
        MaterialBase(char *,int);
		virtual ~MaterialBase();
		
		// initialization and verification
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        void SwapProperties(double &,double &,double &,double &);
        void SwapProperties(double &,double &);
        static void pswap(double &,double &);
#ifdef MPM_CODE
		// history data
		virtual int SizeOfHistoryData(void) const;
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
		virtual double GetHistory(int,char *) const;
		double *CreateAndZeroDoubles(char *,int) const;
		void ZeroDoubles(char *,int) const;
        virtual void SetInitialConditions(InitialCondition *,int,bool);
		virtual Vector GetDamageNormal(MPMBase *,bool) const;
	
		virtual void FillTransportProperties(TransportProperties *);
		virtual void SetHardeningLaw(char *);
		virtual bool AcceptHardeningLaw(HardeningLawBase *,int);
		virtual HardeningLawBase *GetPlasticLaw(void) const;

		virtual void ValidateForUse(int) const;
		virtual void PrintTransportProperties(void) const;
        virtual void SetInitialParticleState(MPMBase *,int,int) const;
#endif

		// initialization (base class only)
		virtual void PrintMaterial(int) const;
		virtual void PrintCommonProperties(void) const;
#ifdef MPM_CODE
		void PrintCriterion(int,int) const;
#endif
         
        // Methods
#ifdef MPM_CODE
		virtual void GetTransportProps(MPMBase *,int,TransportProperties *) const;
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual double GetHeatCapacity(MPMBase *) const;
		//virtual double GetCpHeatCapacity(MPMBase *) const;
		virtual double GetCpMinusCv(MPMBase *) const;
		virtual double GetDiffusionCT(void) const;
		virtual bool NeedsCrelExtrap(void) const;
		virtual double GetCsatRelative(MPMBase *) const;
		virtual double GetMatDiffusionSource(int,MPMBase *,double,double,double,double,double,double) const;
		virtual bool GetParticleDiffusionSource(DiffusionTask *,MPMBase *,double,double *) const;
    
        virtual void IncrementHeatEnergy(MPMBase *,double,double) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int,Tensor *) const;
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual double GetIncrementalResJ(MPMBase *,ResidualStrains *,double) const;
    	virtual Matrix3 LRGetStrainIncrement(int,MPMBase *,Matrix3,Matrix3 *,Matrix3 *,Matrix3 *,Matrix3 *) const;
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
		virtual const int MaterialID(void) const;
#ifdef MPM_CODE
		virtual double GetRho(MPMBase *) const;
        virtual double WaveSpeed(bool,MPMBase *) const = 0;
		virtual double ShearWaveSpeed(bool,MPMBase *,int) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
		virtual double MaximumDiffusion(void) const;
        virtual double MaximumDiffusivity(void) const;
		virtual bool IsRigid(void) const;
		virtual bool IsRigidBC(void) const;
		virtual bool IsRigidContact(void) const;
		virtual bool IsRigidBlock(void) const;
		virtual bool IsRigidSpringBlock(void) const;
		virtual int MaterialStyle(void) const;
		virtual int KeepsCrackTip(void) const;
		virtual void SetFriction(int,int);
		virtual int GetContactToMaterial(int);
		virtual void ContactOutput(int);
        virtual double GetArtificialViscosity(double,double,MPMBase *) const;
		virtual bool SupportsArtificialViscosity(void) const;
		virtual bool SupportsDiffusion(void) const;
		virtual int GetShareMatField(void) const;
		virtual int AltStrainContains(void) const;
		virtual bool GetCrackingStrain(MPMBase *,Tensor *,bool,Matrix3 *) const;
		virtual Vector GetCrackingCOD(MPMBase *,bool) const;
		virtual double GetMaterialConcentrationSaturation(MPMBase *) const;
		virtual double GetMGEOSXmax(double,double,double,double,double &);
		virtual double GetMGEOSKRatio(double,double,double,double,double) const;
        virtual void SetRelativeStrength(MPMBase *,double);
        virtual void SetRelativeToughness(MPMBase *,double);
#else
        virtual double GetStressStrainZZ(double,double,double,double,double,int);
#endif

		// accessors (base class only)
#ifdef MPM_CODE
		virtual int SetField(int,bool,int,int &);
		int GetField(void) const;
		int GetActiveField(void) const;
		virtual double GetCurrentRelativeVolume(MPMBase *,int) const;
        virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
		virtual void SetStress(Tensor *,MPMBase *) const;
		virtual void IncrementThicknessStress(double,MPMBase *) const;
		virtual Tensor GetStressPandDev(Tensor *,double,MPMBase *) const;
		virtual void SetStressPandDev(Tensor *,MPMBase *) const;
		virtual void IncrementThicknessStressPandDev(double,MPMBase *) const;
	
		// material damping
		virtual void SetDamping(double);
		virtual double GetMaterialDamping(double) const;
    
		// for liquid contact
		virtual double GetViscosity(double) const;
		virtual double BracketContactLawShearRate(double,double,double &,double &,double &,double &) const;
	
        // material in cracks
        virtual int AllowsCracks(void) const;

		// for softening materials
		virtual void SetInitiationLaw(char *);
		virtual bool AcceptInitiationLaw(FailureSurface *,int);
		virtual void SetSofteningLaw(char *,int);
		virtual bool AcceptSofteningLaw(SofteningLaw *,int);
		virtual double *GetSoftHistoryPtr(MPMBase *) const;
		virtual int GetTractionFailureSurface(void) const;
	
#ifdef POROELASTICITY
		virtual void UndrainedPressIncrement(MPMBase *,double) const;
		virtual void UndrainedPressIncrement(MPMBase *,double,double,double) const;
#endif
	
#ifdef SUPPORT_MEMBRANES
		// membrane materials
		virtual void MPMMembrane2DLaw(MPMBase *,double,double &,double &,double &,double &,bool &,double,int,void *,ResidualStrains *) const;
		virtual void MPMMembrane3DLaw(MPMBase *,double,double,double &,double,double,
								  double &,double &,double &,bool &,double,int,void *,ResidualStrains *) const;
		static short GetMVFIsMembrane(int);
#endif // end SUPPORT_MEMBRANES
	
#endif	// MPM_CODE
	
		// class methods
		static void PrintProperty(const char *,double,const char *);
		static void PrintProperty(const char *,bool);
#ifdef MPM_CODE
		static const char *PreferredDirection(int);
		static int GetMVFFlags(int);
		static int GetFieldMatID(int);
		static MaterialBase *GetFieldMaterial(int);
		static int GetActiveMatID(int);
		static int GetContactLawNum(int);
#endif
		
	protected:
		int materialID;
		double rho;
		double concSaturation,betaI;	// not used in FEA, but needed for compiling
#ifdef MPM_CODE
		double heatCapacity;			// changed if depends on particle state
        bool artificialViscosity;       // true to false for artifical viscosity
        double avA1,avA2;               // artificial viscosity coefficients
		TransportProperties tr;			// transport tensors
		double diffusionCon,kCond;		// for isotropic properties
		double diffusionCT;				// eta*alpha and eta/Q for poroelasticity, 1 for diffusion
		ContactPair *lastFriction;
		double matPdamping;				// particle damping
		bool matUsePDamping;            // true it particle damping or PIC damping changed in this material
        int allowsCracks;				// false to ignore cracks (option ignored in multimaterial mode)
#ifdef POROELASTICITY
		double Darcy;					// isotropic coefficient or permeability
		double alphaPE,Qalpha;			// Poroelasticity Biot coefficent
		double Ku;						// Poroelasticity undrained bulk modulus
#endif
	
#endif	// MPM_CODE
	
		// constants (changed in MPM time step)
		double C11,C12,C13,C22,C23,C33,C44,C55,C66;
		double S13,S23,S33;				// for generalized plane stress and strain
		double CTE1,CTE2,CTE3;
#ifdef MPM_CODE
		double CME1,CME2,CME3;
		int field,activeField;
		int shareMatField;
#endif
        int swapz;
};

// List of Materials from theMaterials[0] to theMaterials[nmat-1]
// Material IDs are particles, etc., are stored from 1 to nmat
extern MaterialBase **theMaterials;
extern int nmat;

#endif

