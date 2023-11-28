/********************************************************************************
    PhaseFieldDiffusion.hpp
	nairn-mpm-fea
	
	Created by John Nairn on 2/17/2022
	Copyright (c) 2022 John A. Nairn, All rights reserved.
	
	Dependencies
		DiffusionTasks.hpp, TransportTask.hpp
********************************************************************************/

#ifndef _PHASEFIELDDIFFUSION_

#define _PHASEFIELDDIFFUSION_

#include "Custom_Tasks/DiffusionTask.hpp"

class MPMBase;
class NodalPoint;
class MatVelocityField;

class PhaseFieldDiffusion : public DiffusionTask
{
	public:
	
		// constructors and destructors
        PhaseFieldDiffusion(int);
	
		// mass and momentum
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double,short,int);
        virtual double GetVpCTp(MPMBase *);

		// grid forces
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int);
	
		// update particles
		virtual void AdjustRateAndValue(MPMBase *,double &,double &,double &,double) const;
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double,double) const;

		// update particle strains
		virtual void GetDeltaValue(MPMBase *,double,double *) const;
	
		// accessors
		virtual const char *StyleName(void);
		virtual TransportTask *TransportTimeStepFactor(int,double *);
	
};

#endif
