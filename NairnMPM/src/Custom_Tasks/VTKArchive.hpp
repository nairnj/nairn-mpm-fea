/********************************************************************************
    VTKArchive.hpp
    NairnMPM
    
    Created by John Nairn on 12/5/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp, GridArvhive.hpp
********************************************************************************/

#ifndef _VTKARCHIVETASK_

#define _VTKARCHIVETASK_

#include "Custom_Tasks/CustomTask.hpp"
#include "Custom_Tasks/GridArchive.hpp"

class VTKArchive : public GridArchive
{
    public:
        
        // constructors and destructors
        VTKArchive();
        
        // standard methods
		virtual const char *TaskName(void);
		virtual char *InputParam(char *,int &,double &);
        virtual CustomTask *Initialize(void);
	
		// grid archive methods
		virtual bool CheckExportForExtrapolations(void);
		virtual void AllocateExtrapolationBuffers(void);
		virtual CustomTask *NodalExtrapolation(NodalPoint *,MPMBase *,short,int,double,short);
		virtual void FinishExtrapolationCalculations(void);
		virtual void ExportExtrapolationsToFiles(void);
    
    private:
		vector< int > quantity;
		vector< int > quantitySize;
		vector< char * > quantityName;
        vector< int > qparam;
        int intIndex;
        int intArgs[MAX_INTEGER_ARGUMENTS];
		int bufferSize;				// if task has quantity that must be extrapolated
		double **vtk;				// buffer when extrapolating to the nodes (1-based array for node extrapolations)
		
};

#endif

