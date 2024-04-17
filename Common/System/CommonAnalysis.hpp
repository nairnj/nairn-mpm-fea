/********************************************************************************
    CommonAnalysis.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 13, 2006
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _COMMONANALYSIS_

#define _COMMONANALYSIS_

#define NUMBER_DEVELOPMENT_FLAGS 20

class CommonReadHandler;

class CommonAnalysis
{
    public:
		int np;									// analysis method
		int nfree;								// degrees of freedom (here 2D)
		int dflag[NUMBER_DEVELOPMENT_FLAGS];	// flags that can be set for use in development
		int randseed;							// random seend if > 0
	
        //  Constructors and Destructor
        CommonAnalysis();
		virtual ~CommonAnalysis();
		
		// abstract methods
		virtual CommonReadHandler *GetReadHandler(void) = 0;
		virtual void PrintAnalysisTitle(void) = 0;
		virtual void GetAnalysisType(int,char *) = 0;
		virtual const char *CodeName(void) const = 0;
		virtual const char *NodesAndElementsTitle(void) const = 0;
		virtual bool ValidAnalysisType(void) = 0;
		virtual void CMStartResultsOutput(void) = 0;
		virtual void CMAnalysis(bool) = 0;
	
		// abstract archiver methods
		virtual void ArchiveNodalPoints(int) = 0;
		virtual void ArchiveElements(int) = 0;
		virtual void SetInputDirPath(const char *,bool) = 0;
	
		// start the analysis
		virtual void StartAnalysis(bool);
		virtual void StartResultsOutput(void);
		virtual void PrintAnalysisMethod(void);
		virtual void CMPreparations(void);
	
		// methods
		int ReadFile(const char *,bool);
		void CoutCodeVersion(void);
		double CPUTime(void);
		double ElapsedTime(void);

		// accessors
		void SetValidate(bool);
		bool GetValidate(void);
		void SetReverseBytes(bool);
		bool GetReverseBytes(void);
		void SetDescription(const char *);
		char *GetDescription(void);
		bool IsThreeD(void);
		bool IsAxisymmetric(void);
		virtual const char *MPMAugmentation(void);
		void SetNumberOfProcessors(int);
		int GetNumberOfProcessors(void);
		int GetTotalNumberOfPatches(void);
	
	protected:
		bool validate,reverseBytes;
		int version,subversion,buildnumber,numProcs;
		char *description;
	
#ifdef _OPENMP
		double startTime;
#else
		time_t startTime;			// timers
#endif
		clock_t startCPU;
};

// Input types
// Note: Keep all SOFT(type)_LAW_SELECTION in same order
enum { NO_INPUT=0,TEXT_BLOCK,INT_NUM,DOUBLE_NUM,NODE_BLOCK,
	BC_BLOCK,MPMORDER_BLOCK,CRACKORDER_BLOCK,NOT_NUM,
	STRESS_LIST,OUTFLAGS_BLOCK,ANALYSIS_NUM,FUNCTION_BLOCK,
	SETTING_FUNCTION_BLOCK,SETTING_FUNCTION2_BLOCK,SETTING_FUNCTION3_BLOCK,VALUE_FUNCTION_BLOCK,
	POLYHEDRON_BLOCK,HARDENING_LAW_SELECTION,STRESS_FUNCTION_BLOCK,
	GRID_X_BODY_FORCE_FUNCTION_BLOCK,GRID_Y_BODY_FORCE_FUNCTION_BLOCK,GRID_Z_BODY_FORCE_FUNCTION_BLOCK,
	PRESSURE_FUNCTION_BLOCK,TEXT_PARAMETER,INITIATION_LAW_SELECTION,
	SOFTI_LAW_SELECTION,SOFTII_LAW_SELECTION,SOFTAI_LAW_SELECTION,SOFTTII_LAW_SELECTION,SOFTAII_LAW_SELECTION,
	SOFTZZ_LAW_SELECTION,SOFTXZX_LAW_SELECTION,SOFTXZZ_LAW_SELECTION,
	SOFTYZZ_LAW_SELECTION,PHASE_NUM,PATCH_SIZES,EXPRESSION_STR };

void ThrowSAXException(const char *);
void ThrowSAXException(const char *,const char *);

#endif
