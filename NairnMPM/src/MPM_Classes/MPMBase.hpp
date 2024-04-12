/********************************************************************************
    MPMBase.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _MPM_BASE_

#define _MPM_BASE_

// If define, traction BC always integrate on deformed edge, regardless of shape function
// being used. If not defined, uGIMP and and non-GIMP use undeformed edges
//#define TRACTION_ALWAYS_DEFORMED

#define DEFORMED_VOLUME 0
#define DEFORMED_AREA 1
#define DEFORMED_AREA_FOR_GRADIENT 2

#define GRAD_GLOBAL 0
#define GRAD_SECOND 3
#define GRAD_THIRD 6
enum { gGRADx=0,gGRADy,gGRADz };

// maxElementIntersections - max number of elements for finite GIMP
extern int maxElementIntersections;

class MaterialBase;

class MPMBase : public LinkedObject
{
    public:
		//  variables (changed in MPM time step)
		Vector pos,vel;
		char *vfld;
	
		// for conduction, pTemp have 3, 6, or 9 depending on gradient needs
		double pTemperature,pPreviousTemperature,*pTemp;
	
		// dT and dC for for particle residual strains
		ResidualStrains dTrans;

		// conc potential (0 to 1) (archived * concSaturation)
		DiffusionInfo **pDiff;
	
		// for generalized plane stress or strain
		double oopIncrement;

		// constants (not changed in MPM time step)
        double mp;
		Vector origpos;
               
        // constructors and destructors and other intializers
        MPMBase();
        MPMBase(int,int,double);
		virtual ~MPMBase();
		void AllocateJStructures(void);
        bool AllocateCPDIorGIMPStructures(int,bool);
		bool AllocateGIMPStructures(int,bool);
		void ResetMaterialPoint(void);

        // virtual methods
		virtual ResidualStrains ScaledResidualStrains(int);
        virtual double thickness(void) = 0;
        virtual void SetOrigin(Vector *) = 0;
        virtual void SetPosition(Vector *) = 0;
        virtual void SetVelocity(Vector *) = 0;
		virtual void UpdateStrain(double,int,int,void *,int,bool) = 0;
		virtual void GetFintPlusFext(Vector *,double,double,double,double) = 0;
		virtual void MoveParticle(GridToParticleExtrap *) = 0;
		virtual void MovePosition(double) = 0;
		virtual void SetVelocitySpeed(double) = 0;
		virtual void AddTemperatureGradient(int);
		virtual void AddTemperatureGradient(int,Vector *) = 0;
		virtual void AddConcentrationGradient(int);
		virtual void AddConcentrationGradient(int,Vector *) = 0;
		virtual double FCond(int,double,double,double,TransportProperties *) = 0;
		virtual double FDiff(double,double,double,TransportProperties *,int) = 0;
		virtual double KineticEnergy(void) = 0;
		virtual Matrix3 GetDeformationGradientMatrix(void) const = 0;
		virtual Matrix3 GetElasticBiotStrain(void) = 0;
		virtual void SetDeformationGradientMatrix(Matrix3) = 0;
		virtual Matrix3 GetDisplacementGradientMatrix(void) const = 0;
		virtual Matrix3 GetDisplacementGradientForJ(const MaterialBase *) = 0;
        virtual Matrix3 GetElasticLeftCauchyMatrix(void) = 0;
        virtual void GetDeformationGradient(double F[][3]) const = 0;
        virtual double GetRelativeVolume(void) = 0;
		virtual double GetVolume(int) = 0;
        virtual void GetSemiSideVectors(Vector *,Vector *,Vector *) const = 0;
		virtual void ScaleSemiSideVectorsForCPDI(Vector *,Vector *,Vector *) const = 0;
		virtual double GetDeformedRadius(Vector *) const = 0;
		virtual void GetUndeformedSemiSides(double *,double *,double *) const = 0;
		virtual void GetCPDINodesAndWeights(int) = 0;
		virtual Vector GetSurfaceInfo(int,int,int *,Vector *,Vector *,int *,double *) = 0;
		virtual void SetDimensionlessSize(Vector *);
		virtual void SetDimensionlessByPts(int);
        virtual void GetDimensionlessSize(Vector &) const;
		virtual void GetInitialSize(Vector &) const;
		virtual Vector GetParticleSize(void) const;
		virtual double GetParticleXSize(void) const;
		virtual double GetParticleYSize(void) const;
		virtual double GetParticleZSize(void) const;
		virtual double GetMinParticleLength(void) const;
 
        // defined virtual methods
		virtual double GetUnscaledVolume(void);
		virtual void IncrementDeformationGradientZZ(double dezz);
 
		// base only methods (make virtual if need to override)
		int MatID(void) const;
		int ArchiveMatID(void) const;
		bool InReservoir(void) const;
		int ElemID(void) const;
		void ChangeElemID(int,bool);
		int ArchiveElemID(void);
		void ReverseParticle(bool,bool);
		void StopParticle(void);
		int GetElementCrossings(void);
		void SetElementCrossings(int);
		void IncrementElementCrossings(void);
		bool HasLeftTheGridBefore(void);
		void SetHasLeftTheGridBefore(bool);
		void Get2DSinCos(double *,double *);
		Matrix3 *GetRtotPtr(void);
		void SetRtot(Matrix3);
		void InitRtot(Matrix3);
		Matrix3 GetRtot(void);
		double GetRotationZ(void);
		double GetRotationY(void);
		double GetRotationX(void);
		double GetParticleRotationZ(void);
		double GetParticleRotationY(void);
		double GetParticleRotationX(void);
		double GetRotationZInDegrees(void);
		double GetRotationYInDegrees(void);
		double GetRotationXInDegrees(void);
		void SetAnglez0(double);
		void SetAngley0(double);
		void SetAnglex0(double);
		void SetAnglez0InDegrees(double);
		void SetAngley0InDegrees(double);
		void SetAnglex0InDegrees(double);
		double GetAnglez0InDegrees(void);
		double GetAngley0InDegrees(void);
		double GetAnglex0InDegrees(void);
        double GetAnglez0InRadians(void);
        double GetAngley0InRadians(void);
        double GetAnglex0InRadians(void);
		Matrix3 GetBiotStrain(void) const;
		virtual Matrix3 GetInitialRotation(void);
		void IncrementRotationStrain(double);
		void IncrementRotationStrain(double,double,double);
		virtual void InitializeMass(double,double,bool);
		void SetConcentration(double,bool);
		void SetTemperature(double,double);
        void SetVelocityGradient(double,double,double,double,int);
		Vector *GetPFext(void);
		Vector *GetNcpos(void);
		CPDIDomain **GetCPDIInfo(void);
		Vector *GetAcc(void);
		Tensor *GetVelGrad(void);
		double GetPlastEnergy(void);
		void AddPlastEnergy(double);
    	double GetClear_dTad(void);
    	void Add_dTad(double);
		double GetWorkEnergy(void);
		void SetWorkEnergy(double);
		void AddWorkEnergy(double);
		double GetResidualEnergy(void);
		void SetResidualEnergy(double);
		void AddResidualEnergy(double);
		void AddWorkEnergyAndResidualEnergy(double,double);
		double GetStrainEnergy(void);
        double GetHeatEnergy(void);
        void SetHeatEnergy(double);
        void AddHeatEnergy(double);
        double GetEntropy(void);
        void SetEntropy(double);
        void AddEntropy(double,double);
		void AddEntropy(double,double,double);
        double GetInternalEnergy(void);
        void IncrementPressure(double);
        void SetPressure(double);
        double GetPressure(void);
		Tensor ReadStressTensor(void);
		void StoreStressTensor(Tensor *);
		void StoreThicknessStressIncrement(double);
		void StoreThicknessStrainIncrement(double);
		Tensor *GetStressTensor(void);
		Tensor *GetStrainTensor(void);
 		Tensor *GetAltStrainTensor(void);
		TensorAntisym *GetRotationStrainTensor(void);
		char *GetHistoryPtr(int);
		void SetHistoryPtr(char *);
		double GetHistoryDble(int,int);
		void SetHistoryDble(int,double,int);
		char *GetCustomHistoryPtr();
		void SetCustomHistoryPtr(char *);
		double GetCustomHistoryDble(int);
		void SetCustomHistoryDble(int,double);
        void Describe(void);
		double GetRho(void);
		double GetConcSaturation();
        void AddParticleDiffusionSource(int,double);
        bool GetClearParticleDiffusionSource(int,double &);
        void SetRelativeStrength(double);
        void SetRelativeToughness(double);
		virtual void GetExactTractionInfo(int,int,int *,Vector *,Vector *,int *) const;
        void SetNum(int zeroBasedNumber);
        int GetNum(void) const;
    
        // These always defined, but a reserved for a future feature
		Matrix3 GetParticleGradVp(bool);
		Vector GetParticleAngMomentum(void);
	
		bool AllocateFiniteGIMPStructures(bool);
		virtual void GetFiniteGIMP_Integrals(void);
		FiniteGIMPInfo *GetFiniteGIMPInfo(void);

	protected:
		// variables (changed in MPM time step)
		Vector mpm_lp;				// Dimensionless size relative to current element (radius in -1 to 1 natural coordinates)
		Vector pFext;				// external force
		Vector ncpos;				// natural coordinates position
		char *cpdi_or_gimp;         // Should make pointer and allocate only what is needed
        Vector *faceArea;           // make pointer when needed
		Vector acc;					// acceleration (hold velocity of rigid particle in hold phase)
		Tensor *velGrad;			// used for J Integral only on non-rigid particles only
		Tensor sp;					// stress tensor (init 0)
        double pressure;            // for use if materials wants to, otherwise it is zero
		Tensor ep;					// total strain tensor (init 0)
		Tensor eplast;				// plastic strain tensor (init 0)
		TensorAntisym wrot;			// rotation strain tensor (init 0)
		double plastEnergy;			// total plastic energy
		double prev_dTad;			// adiabatic temperature rise in previous step
    	double buffer_dTad;			// adiabatic temperature rise current step
		double workEnergy;			// total work energy  sigma.de
        double heatEnergy;          // total heat flow on the particle
        double entropy;             // total entropy on the particle
		double resEnergy;			// total residual energy sigma.dres
		char *matData;				// material history if needed (init NULL)
		char *customMatData;		// particle history by custom tasks
		Matrix3 *Rtot;				// only track for large rotation hypo and 3D aniso small rotation
	
		// constants (not changed in MPM time step)
 		double anglez0;				// initial cw x rotation angle (2D or 3D) (stored in radians)
		double angley0;				// initial cw y rotation (3D)
		double anglex0;				// initial cw x rotation (3D)
	
    private:
		// variables (changed in MPM time step)
		int inElem;
		int elementCrossings;		// fabs() is # element crossinsgs, when <0 particle has left the grid
	
		// constants (not changed in MPM time step)
        int matnum;
        int num;            // address in mpm[] array is num-1
};

// Lists of material points from mpm[0] to mpm[nmpms-1]
// ordered such that nonrigid are from mpm[0] to mpm[nmpmsNR-1], rigid block [mpm[nmpmsNR] to mpm[nmpmsRB-1],
//      rigid contact mpm[nmpmsRB] to  mpm[nmpmsRC-1] and rigid BC mpm[nmpmsRC] to  mpm[nmpms-1]
extern MPMBase **mpm;
extern int nmpmsNR,nmpmsRB,nmpmsRC,nmpms;

#endif
