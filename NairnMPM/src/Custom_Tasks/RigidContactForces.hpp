/*
 *  RigidContactForces.h
 *  NairnMPM
 *
 *  Created by John Nairn on 7/29/10.
 *  Copyright 2010 Oregon State University. All rights reserved.
 *
 */

/********************************************************************************
	RigidContactForces.hpp
	NairnMPM

	Created by John Nairn on  7/29/10.
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _RIGIDCONTACTFORCES_

#define _RIGIDCONTACTFORCES_

#import "Custom_Tasks/CustomTask.hpp"

class RigidContactForces : public CustomTask
{
    public:
        
        // constructors and destructors
        RigidContactForces();
        
        // standard methods
		virtual const char *TaskName(void);
        virtual CustomTask *Initialize(void);
        virtual CustomTask *PrepareForStep(bool &);
	virtual CustomTask *BeginExtrapolations(void);
        virtual CustomTask *ParticleCalculation(NodalPoint *,MPMBase *,short,int,double,double,double,double,short);
        virtual CustomTask *ParticleExtrapolation(MPMBase *,short);
        
    private:
		int getForcesThisStep;
};


#endif
