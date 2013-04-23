/********************************************************************************
    VTKArchive.hpp
    NairnMPM
    
    Created by John Nairn on 12/5/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _VTKARCHIVETASK_

#define _VTKARCHIVETASK_

#include "Custom_Tasks/CustomTask.hpp"

class VTKArchive : public CustomTask
{
    public:
        
        // constructors and destructors
        VTKArchive();
        
        // standard methods
		virtual const char *TaskName(void);
		virtual char *InputParam(char *,int &);
        virtual CustomTask *Initialize(void);
	
        virtual CustomTask *PrepareForStep(bool &);
		virtual CustomTask *StepCalculation(void);
        virtual CustomTask *FinishForStep(void);
	
        virtual CustomTask *BeginExtrapolations(void);
        virtual CustomTask *EndExtrapolations(void);
        virtual CustomTask *NodalExtrapolation(NodalPoint *,MPMBase *,short,int,double,short);
        
    private:
		vector< int > quantity;
		vector< int > quantitySize;
		vector< char * > quantityName;
		double customArchiveTime,nextCustomArchiveTime;
		int bufferSize;				// if task has quantity that must be extrapolated
		bool getVTKExtraps;			// flag to do extrapolations this step
		bool doVTKExport;			// flag to export the file this step
		double **vtk;				// buffer when extrapolating to the nodes (1-based array for node extrapolations)
		
};

#endif

