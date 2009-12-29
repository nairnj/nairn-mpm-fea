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

class CommonAnalysis
{
    public:
		int np;							// analysis method
		int nfree;						// degrees of freedom (here 2D)
	
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
		int ReadFile(const char *);
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

#endif
