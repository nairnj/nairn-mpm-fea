/********************************************************************************
	Delete Damaged Particles Custon Task
	nairn-mpm-fea
 
	Created by Chad Hammerquist, 2017
	Copyright (c) 2017 John A. Nairn, All rights reserved.
 
	Dependencies
 		CustomTask.hpp
********************************************************************************/

#ifndef _DELETEDAMAGED_

#define _DELETEDAMAGED_

enum { TOTAL_COD=1,TENSILE_COD,SHEAR_COD,SHEARXY_COD,SHEARXZ_COD };

#include "Custom_Tasks/CustomTask.hpp"

class DeleteDamaged : public CustomTask
{
	public:
		// constructors and destructors
		DeleteDamaged();

		// standard methods
		virtual const char *TaskName(void);
		virtual char *InputParam(char *, int &, double &);
		virtual void SetTextParameter(char *,char *);
		virtual CustomTask *Initialize(void);
		virtual CustomTask *StepCalculation(void);
        virtual bool HasReport(void);
        virtual CustomTask *Report(void);

	private:
		// input values
		char *matName;
		int material;
		int damage_direction;
		Vector Store;
		double minRelativeCod;
        int numDeleted;
        double deleteTime;
};

#endif
