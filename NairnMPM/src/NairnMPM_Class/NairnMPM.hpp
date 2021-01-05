/*********************************************************************
    NairnMPM.hpp
    nairn-mpm-fea
    
    Created by jnairn on Mon Nov 19 2001.
    Copyright (c) 2001, All rights reserved.
	
	Header for NairnMPM.hpp, ReadInput.cpp, and OutputSetup.cpp

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
		int ptsPerElement,ptsPerSide;	// points per element and side
		int customPtsPerElement,customPtsPerSide;		// to overide in current region
		int propagate[2];				// progation method
		int propagateDirection[2];		// optional crack direction setting
		int propagateMat[2];			// optional traction material to create when propagates
		bool hasTractionCracks;			// TRUE if any crack segment has traction law material
		int mstep;						// step number
		double maxtime;					// maximum time for analysis
		int warnParticleLeftGrid;		// warning ID for particle leaving the grid
        bool deleteLeavingParticles;    // to delete (true) or push back (false) leaving particles
        int warnParticleDeleted;        // warning ID for deleting a nan particle
		bool multiMaterialMode;			// TRUE to use separate velocity fields for each material
		bool hasNoncrackingParticles;	// TRUE is some particles are ignoring the cracks
		bool skipPostExtrapolation;		// Skip post update extrapolation
		double timeStepMinMechanics;	// time step for  mechanics
		bool exactTractions;			// implement exact tractions
#ifdef RESTART_OPTION
        double restartScaling;          // if >0, restart time step if acceleration too high with scaled time step
        double restartCFL;              // fraction of cell size travel to trigger a restart
        int warnRestartTimeStep;        // warn the first time time step is restarted
#endif
	
        //  Constructors and Destructor
		NairnMPM();
		
		// output
		virtual void PrintAnalysisTitle(void);
		virtual void CMStartResultsOutput(void);
		virtual void PrintAnalysisMethod(void);
		void OutputBCMassAndGrid(void);
	
		// analysis
		virtual void CMPreparations(void);
		virtual void CMAnalysis(bool);
		void MPMStep(void);
	
		// preparation tasks
		void CreateTransportTasks(void);
		void PreliminaryParticleCalcs(void);
		void GridAndElementCalcs(void);
		void SetupMaterialModeContactXPIC(void);
		void PreliminaryCrackCalcs(void);
		void CreateWarnings(void);
		void CreateTasks(void);
	
		// support
		void ReorderParticles(int,int);
        void SwapMaterialPoints(int,int) const;
		void CFLTimeStep(void);
		void ReorderPtBCs(MatPtLoadBC *,int,int);
		void SetForceBCs(void);
		void ValidateOptions(void);
		void Usage(void);
    
		// accessors
		virtual bool ValidAnalysisType(void);
		virtual void SetHasTractionCracks(bool);
        double *GetCFLPtr(void);
		double *GetTransCFLPtr(void);
		double GetCFLCondition(void);
		double GetTransCFLCondition(void);
        double GetPropagationCFLCondition(void);
		bool HasDiffusion(void);
#ifdef POROELASTICITY
		bool HasPoroelasticity(void);
#endif
		bool HasFluidTransport(void);
	
		// archiver access while reading
		virtual void ArchiveNodalPoints(int);
		virtual void ArchiveElements(int);
		virtual void SetInputDirPath(const char *,bool);
	
		// output accessors
		virtual CommonReadHandler *GetReadHandler(void);
		virtual void GetAnalysisType(int,char *);
		virtual const char *CodeName(void) const;
		virtual const char *NodesAndElementsTitle(void) const;
		virtual const char *MPMAugmentation(void);

	protected:
        double FractCellTime;			// fraction of cell crossed at wave speed (<1)
        double PropFractCellTime;       // separate fraction of cell crossed for propagation time steps (currently not settable)
		double TransFractCellTime;		// fraction used for finding transport time step (<1)
};

extern NairnMPM *fmobj;

#endif

