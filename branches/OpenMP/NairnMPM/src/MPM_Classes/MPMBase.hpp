/********************************************************************************
    MPMBase.hpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _MPM_BASE_

#define _MPM_BASE_

#define DEFORMED_VOLUME 0
#define DEFORMED_AREA 1

class MaterialBase;

// variables for conduction calculations
typedef struct {
	Vector DT;			// conc potential gradient (archived * concSaturation)
} TemperatureField;

// variables for diffusion calculations
typedef struct {
	Vector Dc;			// conc potential gradient (archived * concSaturation)
} DiffusionField;

class MPMBase : public LinkedObject
{
    public:
		//  variables (changed in MPM time step)
		Vector pos,vel;
		double pTemperature,pPreviousTemperature;
		double pConcentration,pPreviousConcentration;	// conc potential (0 to 1) (archived * concSaturation)
		char *vfld;
		TemperatureField *pTemp;
		DiffusionField *pDiffusion;
	
		// constants (not changed in MPM time step)
        double mp;
		Vector origpos;
               
        // constructors and destructors and other intializers
        MPMBase();
        MPMBase(int,int,double);
		virtual ~MPMBase();
		void AllocateDiffusion(void);
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
        virtual void MovePosition(double,Vector *) = 0;
        virtual void MoveVelocity(double,double,Vector *) = 0;
		virtual void SetVelocitySpeed(double) = 0;
		virtual void AddTemperatureGradient(void);
		virtual void AddTemperatureGradient(Vector *) = 0;
		virtual double FCond(double,double,double,TransportProperties *) = 0;
		virtual void AddConcentrationGradient(void);
		virtual void AddConcentrationGradient(Vector *) = 0;
		virtual double FDiff(double,double,double,TransportProperties *) = 0;
		virtual double KineticEnergy(void) = 0;
		virtual Matrix3 GetDeformationGradientMatrix(void) = 0;
        virtual Matrix3 GetElasticLeftCauchyMatrix(void) = 0;
        virtual void GetDeformationGradient(double F[][3]) = 0;
        virtual double GetRelativeVolume(void) = 0;
		virtual double GetVolume(bool) = 0;
		virtual void GetCPDINodesAndWeights(int) = 0;
		virtual double GetTractionInfo(int,int,int *,Vector *,Vector *,int *) = 0;

        // defined virtual methods
        virtual double GetUnscaledVolume(void);
        double GetMassForGradient(void);

		// base only methods (make virtual if need to override)
		int MatID(void);
		int ArchiveMatID(void);
		int ElemID(void);
		void ChangeElemID(int);
		int ArchiveElemID(void);
		void ReverseParticle(void);
		void StopParticle(void);
		int GetElementCrossings(void);
		void SetElementCrossings(int);
		void IncrementElementCrossings(void);
		bool HasLeftTheGrid(void);
		void SetHasLeftTheGrid(bool);
        bool PartitionsElasticAndPlasticStrain(void);
		double GetDuDx(void);					// du/dr in axisym
		double GetDuDy(void);					// du/dz in axisym
		double GetDvDx(void);					// dw/dr in axisym
		double GetDvDy(void);					// dw/dz in axisym
		double GetDwDz(void);					// v/r = etheta in axisym
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
		void AddPlastEnergy(double energyInc);
		double GetDispEnergy(void);
		void AddDispEnergy(double energyInc);
		void SetDispEnergy(double energy);
		double GetStrainEnergy(void);
		void SetStrainEnergy(double energyTot);
		void AddStrainEnergy(double energyInc);
        double GetHeatEnergy(void);
        void SetHeatEnergy(double energyTot);
        void AddHeatEnergy(double energyInc);
        double GetInternalEnergy(void);
		double GetExtWork(void);
        void IncrementPressure(double);
        void SetPressure(double);
        double GetPressure(void);
		Tensor ReadStressTensor(void);
		Tensor *GetStressTensor(void);
		Tensor *GetStrainTensor(void);
 		Tensor *GetPlasticStrainTensor(void);
        Tensor *GetElasticLeftCauchyTensor(void);
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
		Vector acc;					// acceleration
		Tensor *velGrad;			// used for J Integral only on non-rigid particles only
		Tensor sp;					// stress tensor (init 0)
        double pressure;            // for use if materials wants to, otherwise it is zero
		Tensor ep;					// total strain tensor (init 0)
		Tensor eplast;				// plastic strain tensor (init 0)
		TensorAntisym wrot;			// rotation strain tensor (init 0)
		double plastEnergy;			// total plastic energy
		double dispEnergy;			// dissipated energy in current step
		double strainEnergy;		// total strain energy
        double heatEnergy;          // total heat flow on the particle
		double extWork;				// total external work
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
