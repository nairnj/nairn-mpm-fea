/********************************************************************************
    CommonAnalysis.hpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on Jan 13, 2006
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _COMMONANALYSIS_

#define _COMMONANALYSIS_

#define NUMBER_DEVELOPMENT_FLAGS 5

class CommonAnalysis
{
    public:
		int np;									// analysis method
		int nfree;								// degrees of freedom (here 2D)
		int dflag[NUMBER_DEVELOPMENT_FLAGS];	// flags that can be set for use in development
	
        //  Constructors and Destructor
        CommonAnalysis();
		virtual ~CommonAnalysis();
		
		// abstract methods
		virtual void PrintAnalysisTitle(void) = 0;
		virtual void PrintAnalysisType(void) = 0;
		virtual char *CodeName(void) = 0;
		virtual bool ValidAnalysisType(void) = 0;
		virtual void MyStartResultsOutput(void) = 0;
	
		// methods
		int ReadFile(const char *,bool);
		void StartResultsOutput(void);
		void CoutCodeVersion(void);

		// accessors
		void SetValidate(bool);
		bool GetValidate(void);
		void SetReverseBytes(bool);
		bool GetReverseBytes(void);
		void SetDescription(const char *);
		char *GetDescription(void);
		bool IsThreeD(void);
		static bool IsThreeD(int);
		bool IsAxisymmetric(void);
		static bool IsAxisymmetric(int);
		
	protected:
		bool validate,reverseBytes;
		int version,subversion,buildnumber;
		char *description;
};

// Input types
enum { NO_INPUT=0,TEXT_BLOCK,INT_NUM,DOUBLE_NUM,NODE_BLOCK,
	BC_BLOCK,MPMORDER_BLOCK,CRACKORDER_BLOCK,NOT_NUM,
	STRESS_LIST,OUTFLAGS_BLOCK,ANALYSIS_NUM,FUNCTION_BLOCK,
	SETTING_FUNCTION_BLOCK,SETTING_FUNCTION2_BLOCK,SETTING_FUNCTION3_BLOCK,VALUE_FUNCTION_BLOCK,
	POLYHEDRON_BLOCK };

void ThrowSAXException(const char *);
void ThrowSAXException(const char *,const char *);

#endif
