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

// variables for conduction calculations
typedef struct {
	Vector DT;			// conc potential gradient (archived * concSaturation)
} TemperatureField;

// variables for diffusion calculations
typedef struct {
	double flux;		// potential*volume/sec for external fluxes
	Vector Dc;			// conc potential gradient (archived * concSaturation)
} DiffusionField;

class MPMBase : public LinkedObject
{
    public:
		//  variables (changed in MPM time step)
		Vector pos,vel;
		double pTemperature,pPreviousTemperature;
		double pConcentration,pPreviousConcentration;	// conc potential (0 to 1) (archived * concSaturation)
		char vfld[MaxShapeNds];
		TemperatureField *pTemp;
		DiffusionField *pDiffusion;
		static int currentParticleNum;					// zero based particle number in some loops
	
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
		virtual void UpdateStrain(double,int,int) = 0;
		virtual void Fint(Vector &,double,double,double) = 0;
		virtual void Fext(Vector &,double) = 0;
        virtual void MovePosition(double,Vector *) = 0;
        virtual void MoveVelocity(double,double,Vector *) = 0;
		virtual void SetVelocitySpeed(double) = 0;
		virtual void AddTemperatureGradient(void) = 0;
		virtual void AddTemperatureGradient(Vector *) = 0;
		virtual double FCond(double,double,double) = 0;
		virtual void AddConcentrationGradient(void) = 0;
		virtual void AddConcentrationGradient(Vector *) = 0;
		virtual double FDiff(double,double,double) = 0;
		virtual double KineticEnergy(void) = 0;
        virtual void GetDeformationGradient(double F[][3]) = 0;
        virtual double GetRelativeVolume(void) = 0;
		virtual void GetCPDINodesAndWeights(int) = 0;
		virtual double GetTractionInfo(int,int,int *,Vector *,Vector *,int *) = 0;

       
		// base only methods (make virtual if need to override)
		void SetDilatedVolume(void);
		double GetVolume();
		int MatID(void);
		int ArchiveMatID(void);
		int ElemID(void);
		void ChangeElemID(int);
		int ArchiveElemID(void);
		void ReverseParticle(void);
		void StopParticle(void);
		int GetResetElementCrossings(void);
		double GetDuDy(void);
		double GetDvDx(void);
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
		double GetExtWork(void);
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
    
		// class methods
		static void FullStrainUpdate(double,int,int);
	
	protected:
		// variables (changed in MPM time step)
		Vector pFext;				// external force
		Vector ncpos;				// natural coordinates position
		CPDIDomain **cpdi;          // Should makle pointer and allocate only what is needed
        Vector *faceArea;           // make pointer then needed
		Vector acc;					// acceleration
		Tensor *velGrad;			// used for J Integral only
		Tensor sp;					// stress tensor (init 0)
		Tensor ep;					// total strain tensor (init 0)
		Tensor eplast;				// plastic strain tensor (init 0)
		TensorAntisym wrot;			// rotation strain tensor (init 0)
		double volume;
		double plastEnergy;			// total plastic energy
		double dispEnergy;			// dissipated energy in current step
		double strainEnergy;		// total strain energy
		double extWork;				// total external work
		char *matData;				// material history if needed (init NULL)
	
		// constants (not changed in MPM time step)
 		double anglez0;				// initial cw x rotation angle (2D or 3D) (stored in radians)
		double angley0;				// initial cw y rotation (3D)
		double anglex0;				// initial cw x rotation (3D)
        
    private:
		// variables (changed in MPM time step)
		int inElem;
		int elementCrossings;
	
		// constants (not changed in MPM time step)
        int matnum;
};

// Lists of material points from mpm[0] to mpm[nmpms-1]
extern MPMBase **mpm;
extern int nmpms;

#endif
