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

// to activate simplifed heat energy and entropy calculations
// When permanent, delete dT from IncrementHeatEnergy() arguments
#define NEW_HEAT_METHOD

#define DEFORMED_VOLUME 0
#define DEFORMED_AREA 1
#define DEFORMED_AREA_FOR_GRADIENT 2

#define GRAD_GLOBAL 0
#define GRAD_SECOND 3
#define GRAD_THIRD 6
enum { gGRADx=0,gGRADy,gGRADz };

// define to load nodes and node IDs in initialization rather than on the fly each
// time they are needed (to remove also remove GIMPNodes in DataTypes)
//#define LOAD_GIMP_INFO

class MaterialBase;

class MPMBase : public LinkedObject
{
    public:
		//  variables (changed in MPM time step)
		Vector pos,vel;
		char *vfld;
	
		// for conduction, pTemp have 3, 6, or 9 depending on gradient needs
		double pTemperature,pPreviousTemperature,*pTemp;
	
		// conc potential (0 to 1) (archived * concSaturation)
		double pConcentration,pPreviousConcentration,*pDiffusion;
	
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
    
        // virtual methods
        virtual double thickness(void) = 0;
        virtual void SetOrigin(Vector *) = 0;
        virtual void SetPosition(Vector *) = 0;
        virtual void SetVelocity(Vector *) = 0;
        virtual void SetStress(Tensor *) = 0;
		virtual void UpdateStrain(double,int,int,void *,int) = 0;
		virtual void GetFintPlusFext(Vector *,double,double,double,double) = 0;
		virtual void MovePosition(double,Vector *,Vector *,double) = 0;
		virtual void MoveVelocity(double,Vector *) = 0;
		virtual void MovePosition(double) = 0;
		virtual void SetVelocitySpeed(double) = 0;
		virtual void AddTemperatureGradient(int);
		virtual void AddTemperatureGradient(int,Vector *) = 0;
		virtual double FCond(int,double,double,double,TransportProperties *) = 0;
		virtual void AddConcentrationGradient(void);
		virtual void AddConcentrationGradient(Vector *) = 0;
		virtual double FDiff(double,double,double,TransportProperties *) = 0;
		virtual double KineticEnergy(void) = 0;
		virtual Matrix3 GetDeformationGradientMatrix(void) const = 0;
		virtual Matrix3 GetElasticBiotStrain(void) = 0;
		virtual void SetDeformationGradientMatrix(Matrix3) = 0;
		virtual Matrix3 GetDisplacementGradientMatrix(void) const = 0;
        virtual Matrix3 GetElasticLeftCauchyMatrix(void) = 0;
        virtual void GetDeformationGradient(double F[][3]) const = 0;
        virtual double GetRelativeVolume(void) = 0;
		virtual double GetVolume(int) = 0;
        virtual void GetSemiSideVectors(Vector *,Vector *,Vector *) const = 0;
		virtual void GetUndeformedSemiSides(double *,double *,double *) const = 0;
		virtual void GetCPDINodesAndWeights(int) = 0;
		virtual double GetTractionInfo(int,int,int *,Vector *,Vector *,int *) = 0;
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
		void SetConcentration(double,double);
		void SetTemperature(double,double);
        void SetVelocityGradient(double,double,double,double,int);
		Vector *GetPFext(void);
		Vector *GetNcpos(void);
		CPDIDomain **GetCPDIInfo(void);
#ifdef LOAD_GIMP_INFO
		GIMPNodes *GetGIMPInfo(void);
#endif
		Vector *GetAcc(void);
		Tensor *GetVelGrad(void);
		double GetPlastEnergy(void);
		void AddPlastEnergy(double);
#ifdef NEW_HEAT_METHOD
    	double GetClear_dTad(void);
#else
    	double GetClearPrevious_dTad(void);
    	double GetBufferClear_dTad(void);
#endif
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
        void AddEntropy(double);
        double GetInternalEnergy(void);
        void IncrementPressure(double);
        void SetPressure(double);
        double GetPressure(void);
		Tensor ReadStressTensor(void);
		Tensor *GetStressTensor(void);
		Tensor *GetStrainTensor(void);
 		Tensor *GetAltStrainTensor(void);
		TensorAntisym *GetRotationStrainTensor(void);
		char *GetHistoryPtr(int);
		void SetHistoryPtr(char *);
		double GetHistoryDble(int,int);
		void SetHistoryDble(int,double,int);
        void Describe(void);
		double GetRho(void);
		double GetConcSaturation(void);
    
	protected:
		// variables (changed in MPM time step)
		Vector mpm_lp;				// Dimensionless size relative to current element (radius in -1 to 1 natural coordinates)
		Vector pFext;				// external force
		Vector ncpos;				// natural coordinates position
		char *cpdi_or_gimp;          // Should make pointer and allocate only what is needed
        Vector *faceArea;           // make pointer then needed
		Vector acc;					// acceleration (hold velocity of rigid particle in hold phase)
		Tensor *velGrad;			// used for J Integral only on non-rigid particles only
		Tensor sp;					// stress tensor (init 0)
        double pressure;            // for use if materials wants to, otherwise it is zero
		Tensor ep;					// total strain tensor (init 0)
		Tensor eplast;				// plastic strain tensor (init 0)
		TensorAntisym wrot;			// rotation strain tensor (init 0)
		double plastEnergy;			// total plastic energy
#ifndef NEW_HEAT_METHOD
    	double prev_dTad;			// adiabatic temperature rise in previous step
#endif
    	double buffer_dTad;			// adiabatic temperature rise current step
		double workEnergy;			// total work energy  sigma.de
        double heatEnergy;          // total heat flow on the particle
        double entropy;             // total entropy on the particle
		double resEnergy;			// total residual energy sigma.dres
		char *matData;				// material history if needed (init NULL)
		Matrix3 *Rtot;				// only track for large rotation hypo and 3D aniso small rotation
	
		// constants (not changed in MPM time step)
 		double anglez0;				// initial cw x rotation angle (2D or 3D) (stored in radians)
		double angley0;				// initial cw y rotation (3D)
		double anglex0;				// initial cw x rotation (3D)
    	
    private:
		// variables (changed in MPM time step)
		int inElem;
		int elementCrossings;		// abs() is # element crossinsgs, when <0 particle has left the grid
	
		// constants (not changed in MPM time step)
        int matnum;
};

// Lists of material points from mpm[0] to mpm[nmpms-1]
// ordered such that nonrigid are from mpm[0] to mpm[nmpmsNR-1]
//      rigid contact mpm[nmpmsNR] to  mpm[nmpmsRC-1] and rigid BC mpm[nmpmsRC] to  mpm[nmpms-1]
extern MPMBase **mpm;
extern int nmpms,nmpmsNR,nmpmsRC;

#endif
