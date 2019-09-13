/*********************************************************************
    NairnFEA.hpp
    Nairn Research Group FEA Code
    
    Created by John Nairn on Feb 2 2003.
    Copyright (c) 2003, All rights reserved.
	
	Header for NairnFEA.cpp, FEAReadInput.cpp, and BeginResults.cpp
	
	Dependencies
		none
*********************************************************************/

#ifndef _NAIRNFEA_

#define _NAIRNFEA_

class CommonReadHandler;

// FEA output flags
enum { DISPLACEMENT_OUT=0,FORCE_OUT,ELEMSTRESS_OUT,AVGSTRESS_OUT,
        REACT_OUT,ENERGY_OUT,NUMBER_OUT };

class NairnFEA : public CommonAnalysis
{
    public:
		int nsize;						// size of the problem
		int nband;						// bandwidth of the problem
		double **st;					// stiffness matrix st[i][j] - i,j 1 based
		double *stiffnessMemory;		// actually location of stiffness matrix in contiguous memory for speed
		double *rm;						// reaction vector
		char *temperatureExpr;			// temperature expression
		double stressFreeTemperature;	// stress free temperature
		char xax,yax,zax;				// axis names
		char outFlags[NUMBER_OUT+1];	// output flags
		vector< int > selectedNodes;	// selected nodes for conditional output
		PeriodicInfo periodic;			// periodic settings
		
        //  Constructors and Destructor
		NairnFEA();
		
		// start analysis
		virtual void PrintAnalysisTitle(void);
		virtual void CMAnalysis(bool);
		virtual void CMStartResultsOutput(void);
	
		// FEA methods
		void BeginResults(void);
		void Usage();
		void ForcesOnEdges(void);
		void BuildStiffnessMatrix(void);
		int GetBandWidth(void);
		void DisplacementResults(void);
		void ForceStressEnergyResults(void);
		void AvgNodalStresses(void);
		void ReactionResults(void);
		void EnergyResults(void);
	
		// archiver access while reading
		virtual void ArchiveNodalPoints(int);
		virtual void ArchiveElements(int);
		virtual void SetInputDirPath(const char *,bool);
	
		// accessors
		virtual CommonReadHandler *GetReadHandler(void);
		virtual void GetAnalysisType(int,char *);
		virtual const char *CodeName(void) const;
		virtual const char *NodesAndElementsTitle(void) const;
		virtual bool ValidAnalysisType(void);
};

extern NairnFEA *fmobj;

#endif

