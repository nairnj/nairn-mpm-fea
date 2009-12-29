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
		// properties initialized while reading the input file
        double mp;
		Vector pos,vel,origpos;
        double pTemperature,pPreviousTemperature;
		double pConcentration,pPreviousConcentration;	// conc potential (0 to 1) (archived * concSaturation)
		
		// properties set during calculations
		double volume;
        char vfld[MaxShapeNds];
		
		// dynamic transport data
		TemperatureField *pTemp;
		DiffusionField *pDiffusion;
               
        // constructors and destructors and other intializers
        MPMBase();
        MPMBase(int,int,double);
		virtual ~MPMBase();
		void AllocateDiffusion(void);
		void AllocateTemperature(void);
		void AllocateJStructures(void);
        
        // virtual methods
        virtual double thickness(void) = 0;
        virtual void SetOrigin(Vector *) = 0;
        virtual void SetPosition(Vector *) = 0;
        virtual void SetVelocity(Vector *) = 0;
		virtual void SetDilatedVolume(void) = 0;
		virtual void UpdateStrain(double,int,int) = 0;
		virtual Vector Fint(double,double,double) = 0;
		virtual Vector Fext(double fni) = 0;
        virtual void MovePosition(double,Vector *) = 0;
        virtual void MoveVelocity(double,double,Vector *) = 0;
		virtual void SetVelocitySpeed(double) = 0;
		virtual void AddTemperatureGradient(void) = 0;
		virtual void AddTemperatureGradient(Vector *) = 0;
		virtual double FCond(double,double,double) = 0;
		virtual void AddConcentrationGradient(void) = 0;
		virtual void AddConcentrationGradient(Vector *) = 0;
		virtual double FDiff(double,double,double) = 0;
       
		// base only methods (make virtual if need to override)
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
		double GetRotation(void);
		double GetRotationInDegrees(void);
		void SetAngle0InDegrees(double);
		double GetAngle0InDegrees(void);
		void IncrementRotationStrain(double);
		void IncrementRotationStrain(double,double,double);
		void InitializeMass(double);
		void SetConcentration(double,double);
		void SetTemperature(double,double);
        void SetVelocityGradient(double,double,double,double,int);
		Vector *GetPFext(void);
		Vector *GetNcpos(void);
		Vector *GetAcc(void);
		Tensor *GetVelGrad(void);
		double GetPlastEnergy(void);
		void AddPlastEnergy(double energyInc);
		double GetDispEnergy(void);
		void AddDispEnergy(double energyInc);
		void SetDispEnergy(double energy);
		double GetStrainEnergy(void);
		void AddStrainEnergy(double energyInc);
		double GetExtWork(void);
		Tensor *GetStressTensor(void);
		Tensor *GetStrainTensor(void);
 		Tensor *GetPlasticStrainTensor(void);
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
		// used in calculations
		Vector pFext;				// external force
		Vector ncpos;				// natural coordinates position
 		Vector acc;					// acceleration
		
		Tensor *velGrad;			// used for J Integral only
		Tensor sp;					// stress tensor (init 0)
        Tensor ep;					// total strain tensor (init 0)
		Tensor eplast;				// plastic strain tensor (init 0)
        TensorAntisym wrot;			// rotation strain tensor (init 0)
		
        double plastEnergy;			// total plastic energy
		double dispEnergy;			// dissipated energy in current step
        double strainEnergy;		// total strain energy
		double extWork;				// total external work
		double angle0;				// initial 2D rotation angle (stored in radians)
        
		// dynamic data
        char *matData;				// material history if needed (init NULL)
		
    private:
		// used on initialization
        int inElem;
        int matnum;
		
		// used in calculations
		int elementCrossings;
};

// Lists of material points from mpm[0] to mpm[nmpms-1]
extern MPMBase **mpm;
extern long nmpms;

#endif
