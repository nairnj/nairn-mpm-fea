/********************************************************************************
	GridArchive.hpp
	NairnMPM - FracGeo branch

	Created by John Nairn on 11/4/14.

	Dependencies
		CustomTask.hpp
 *******************************************************************************/

#ifndef _GRIDARCHIVETASK_

#define _GRIDARCHIVETASK_

#include "Custom_Tasks/CustomTask.hpp"

#define MAX_INTEGER_ARGUMENTS 10

class GridArchive : public CustomTask
{
	public:
	
		// constructors and destructors
		GridArchive();
	
		// standard custom task methods methods
		virtual char *InputParam(char *,int &,double &);
		virtual CustomTask *Initialize(void);
		virtual CustomTask *PrepareForStep(bool &);
		virtual CustomTask *StepCalculation(void);
		virtual CustomTask *BeginExtrapolations(void);
		virtual CustomTask *EndExtrapolations(void);
	
		// special grid archive methods
		virtual bool CheckExportForExtrapolations(void);
		virtual void AllocateExtrapolationBuffers(void);
		//virtual CustomTask *NodalExtrapolation(NodalPoint *,MPMBase *,short,int,double,short);
		virtual void FinishExtrapolationCalculations(void);
		virtual void ExportExtrapolationsToFiles(void);
	
	protected:
		double customArchiveTime,nextCustomArchiveTime;
		bool getExtrapolations;			// flag to do extrapolations this step
		bool doExport;					// flag to export the file this step
};

#endif

