/*********************************************************************
    NairnMPM.hpp
    nairn-mpm-fea
    
    Created by jnairn on Mon Nov 19 2001.
    Copyright (c) 2001, All rights reserved.
	
	Header for NairnMPM.hpp, ReadInput.cpp, and OutputSetup.hpp

	Dependencies
		none
*********************************************************************/

#ifndef _NAIRNMPM_

#define _NAIRNMPM_

class MPMBase;
class BoundaryCondition;
class MatPtLoadBC;

// global variables
extern double mtime,propTime,timestep,strainTimestepFirst,strainTimestepLast,fractionUSF;
extern int maxCrackFields,maxMaterialFields,numActiveMaterials;

class NairnMPM : public CommonAnalysis
{
    public:
		int mpmApproach;				// mpm method for updating strain
		int ptsPerElement;				// points per element
		int propagate[2];				// progation method
		int propagateDirection[2];		// optional crack direction setting
		int propagateMat[2];			// optional traction material to create when propagates
		bool hasTractionCracks;			// TRUE if any crack segment has traction law material
		int mstep;						// step number
		double maxtime;					// maximum time for analysis
		int warnParticleLeftGrid;		// warning ID
		bool multiMaterialMode;			// TRUE to use separate velocity fields for each material
		bool hasNoncrackingParticles;	// TRUE is some particles are ignoring the cracks
		bool skipPostExtrapolation;		// Skip post update extrapolation
		double timeStepMinMechanics;	// time step for  mechanics
	
        //  Constructors and Destructor
		NairnMPM();
		
		// methods
		void StartAnalysis(bool);
		virtual void MyStartResultsOutput(void);
		void OutputSetup(void);
		void MPMAnalysis(bool);
		void MPMStep(void);
		void PreliminaryCalcs(void);
		void CreateTasks(void);
		void ReorderPtBCs(MatPtLoadBC *,int,int);
		void SetForceBCs(void);
		void ValidateOptions(void);
		void Usage(void);
    
		// accessors
		virtual void PrintAnalysisTitle(void);
		virtual void PrintAnalysisType(void);
		virtual const char *CodeName(void) const;
		virtual bool ValidAnalysisType(void);
		virtual void SetHasTractionCracks(bool);
        double *GetCFLPtr(void);
		double *GetTransCFLPtr(void);
		double GetCFLCondition(void);
		double GetTransCFLCondition(void);
        double GetPropagationCFLCondition(void);
		virtual const char *MPMAugmentation(void);
		bool HasDiffusion(void);
		bool HasPoroelasticity(void);
		bool HasFluidTransport(void);
	
    protected:
        double FractCellTime;			// fraction of cell crossed at wave speed (<1)
        double PropFractCellTime;       // separate fraction of cell crossed for propagation time steps (currently not settable)
		double TransFractCellTime;		// fraction used for finding transport time step (<1)
};

extern NairnMPM *fmobj;

#endif

