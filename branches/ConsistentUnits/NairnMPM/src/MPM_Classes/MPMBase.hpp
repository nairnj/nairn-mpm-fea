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

#define DEFORMED_VOLUME 0
#define DEFORMED_AREA 1

#define GRAD_GLOBAL 0
#define GRAD_SECOND 3
#define GRAD_THIRD 6
enum { gGRADx=0,gGRADy,gGRADz };

class MaterialBase;

class MPMBase : public LinkedObject
{
    public:
		//  variables (changed in MPM time step)
		Vector pos,vel;
		double pTemperature,pPreviousTemperature;
		double pConcentration,pPreviousConcentration;	// conc potential (0 to 1) (archived * concSaturation)
		char *vfld;
		double *pTemp;
		double *pDiffusion;
	
		// constants (not changed in MPM time step)
        double mp;
		Vector origpos;
               
        // constructors and destructors and other intializers
        MPMBase();
        MPMBase(int,int,double);
		virtual ~MPMBase();
		void AllocateDiffusion(bool);
		void AllocateTemperature(void);
		void AllocateJStructures(void);
        bool AllocateCPDIStructures(int,bool);
    
        // virtual methods
        virtual double thickness(void) = 0;
        virtual void SetOrigin(Vector *) = 0;
        virtual void SetPosition(Vector *) = 0;
        virtual void SetVelocity(Vector *) = 0;
		virtual void UpdateStrain(double,int,int,void *,int) = 0;
		virtual void GetFintPlusFext(Vector *,double,double,double,double) = 0;
        virtual void MovePosition(double,Vector *,double,double) = 0;
        virtual void MoveVelocity(double,double) = 0;
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
		virtual void GetCPDINodesAndWeights(int) = 0;
		virtual double GetTractionInfo(int,int,int *,Vector *,Vector *,int *) = 0;

        // defined virtual methods
        virtual double GetUnscaledVolume(void);
        double GetMassForGradient(void);

		// base only methods (make virtual if need to override)
		int MatID(void) const;
		int ArchiveMatID(void) const;
		int ElemID(void);
		void ChangeElemID(int);
		int ArchiveElemID(void);
		void ReverseParticle(bool,bool);
		void StopParticle(void);
		int GetElementCrossings(void);
		void SetElementCrossings(int);
		void IncrementElementCrossings(void);
		bool HasLeftTheGridBefore(void);
		void SetHasLeftTheGridBefore(bool);
#ifndef USE_PSEUDOHYPERELASTIC
        bool PartitionsElasticAndPlasticStrain(void) const;
#endif
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
		void InitializeMass(double);
		void SetConcentration(double,double);
		void SetTemperature(double,double);
        void SetVelocityGradient(double,double,double,double,int);
		Vector *GetPFext(void);
		Vector *GetNcpos(void);
		CPDIDomain **GetCPDIInfo(void);
		Vector *GetAcc(void);
		Tensor *GetVelGrad(void);
		double GetPlastEnergy(void);
		void AddPlastEnergy(double);
		double GetDispEnergy(void);
		void AddDispEnergy(double);
		void SetDispEnergy(double);
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
		char *GetHistoryPtr(void);
		void SetHistoryPtr(char *);
		double GetHistoryDble(void);
		void SetHistoryDble(double);
		double GetHistoryDble(int);
		void SetHistoryDble(int,double);
        void Describe(void);
    
	protected:
		// variables (changed in MPM time step)
		Vector pFext;				// external force
		Vector ncpos;				// natural coordinates position
		CPDIDomain **cpdi;          // Should makle pointer and allocate only what is needed
        Vector *faceArea;           // make pointer then needed
		Vector acc;					// acceleration (rigid BC particle use to store held velocity)
		Tensor *velGrad;			// used for J Integral only on non-rigid particles only
		Tensor sp;					// stress tensor (init 0)
        double pressure;            // for use if materials wants to, otherwise it is zero
		Tensor ep;					// total strain tensor (init 0)
		Tensor eplast;				// plastic strain tensor (init 0)
		TensorAntisym wrot;			// rotation strain tensor (init 0)
		double plastEnergy;			// total plastic energy
		double dispEnergy;			// dissipated energy in current step
		double workEnergy;			// total work energy  sigma.de
        double heatEnergy;          // total heat flow on the particle
        double entropy;             // total entropy on the particle
		double resEnergy;			// total residual energy sigma.dres
		char *matData;				// material history if needed (init NULL)
	
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
