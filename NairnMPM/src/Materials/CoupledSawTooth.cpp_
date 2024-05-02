/********************************************************************************
    CoupledSawTooth.hpp
    nairn-mpm-fea

    Created by John Nairn on 8/30/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
 
    This material is based on mixed-mode damage law developed by Hogberg such
 	that tractions are
 
       Tn = (sigma_c*S(lamdaw)/(lambdaw*un^c) * un
       Tt = (tau_c*S(lamdaw)/(lambdaw*un^t) * ut
 
    where lambdaw is maximum lambda after failure. The lambda term is
 
       lambda = sqrt((un/un^c)^2 + (ut/ut^c)^2)
 
    is a damage variable. The damage initiates with lambda first reaches
    lambda_p(theta) where
 
	    tan(theta) = un*ut^c/(un^c*ut)
 
    Failure occurs when lambdaw reaches 1. Hogberg used a saw tooth law
 
        S(lambdaw) = (1-lambdaw)/(1-lambda_p(theta))
 
 	The model also allows exponential law with
 
        S(lambdaw) = exp(-alpha*(lambdaw-lambda_p(theta))/(1.-lambda_p(theta)))-exp(-alpha)
                     -----------------------------------------------------------------------
                                        1 - exp(-alpha)
 
    Comanho (and others) are a special case of Hogberg that assume un^c=ut_c = uc. Then
 
        lambda = sqrt(un^2 + ut^2)/uc
 
    This model is technically only valid lambda_p(theta) is independent of theta. The
 	later applies if input properties have
 
        un^e/un^c = ut^e/ut^c
 
    Or on input (where umidI and umidII are relative to delIIc), they
 	are equal. If specifying different properties, would need to pick
 
            JIc/(sn^c un^c) = JIIc/(st^c ut^c)
            kIe un^e/sn^c = kIIe ut^e/st^c
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/CoupledSawTooth.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"
#include "Materials/ExponentialTraction.hpp"

extern double mtime;

#pragma mark CoupledSawTooth::Constructors and Destructors

// Constructor 
CoupledSawTooth::CoupledSawTooth(char *matName,int matID) : CohesiveZone(matName,matID)
{
    model = TRIANGULARTRACTIONMATERIAL;
    alpha = -1.;
    
	// rest are parent class properties
}

#pragma mark CoupledSawTooth::Initialization

// Read properties (read read in super classes)
char *CoupledSawTooth::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"alpha")==0)
    {   input=DOUBLE_NUM;
        return((char *)&alpha);
    }
    
    else if(strcmp(xName,"model")==0)
    {   input=INT_NUM;
        return((char *)&model);
    }
    
    return CohesiveZone::InputTractionLawProperty(xName,input,gScaling);
}

// Calculate properties used in analyses - here triangular law
// Do mode I and mode II separately
const char *CoupledSawTooth::VerifyAndLoadProperties(int np)
{
	const char *msg = CohesiveZone::VerifyAndLoadProperties(np);
	if(msg!=NULL) return msg;
	
    if(model==EXPONENTIALTRACTIONMATERIAL)
    {   // adjust uc of each mode to keep Gc as input (or culculated)
        if(alpha<=0.)
            return "Coupled exponential cohesive laws must specify alpha > 0";
        expalpha = exp(-alpha);
        romexpalpha = 1./(1.-expalpha);
        double fealpha = 2.*(1./alpha - expalpha*romexpalpha);
        
        // rescale uc to keep Gc constant (but this complicates setting uc to valid values)
        // delIc = (2.*JIc/stress1 - (1.-fealpha)*umidI)/fealpha;
        // delIIc = (2.*JIIc/stress2 - (1.-fealpha)*umidII)/fealpha;
        
        // recalculated Gc to input k, ue, uc, and alpha
        JIc = 0.5*stress1*(umidI + fealpha*(delIc-umidI));
        JIIc = 0.5*stress2*(umidII + fealpha*(delIIc-umidII));
    }
    else if(model!=TRIANGULARTRACTIONMATERIAL)
        return "Coupled traction law must be sawtooth or exponential";
	
    unpbar = umidI/delIc;
    unpbar2 = unpbar*unpbar;
    utpbar = umidII/delIIc;
    utpbar2 = utpbar*utpbar;
    
	return NULL;
}

// print to output window
void CoupledSawTooth::PrintMechanicalProperties(void) const
{
    PrintSawToothModel("I",JIc,stress1,delIc,kI1,umidI,-1.);
    PrintSawToothModel("II",JIIc,stress2,delIIc,kII1,umidII,-1.);
    switch(model)
    {   case EXPONENTIALTRACTIONMATERIAL:
            cout << "Exponential decay with alpha = " << alpha << endl;
            break;
        default:
            cout << "Linear decay:" << endl;
            break;
    }
    if(!DbleEqual(unpbar,utpbar))
    {   cout << "WARNING: Because ue/uc for normal and tangential directions differ, this cohesive" << endl;
        cout << "         law will be invalid for mixed-mode loading." << endl;
    }
}

#pragma mark CoupledSawTooth::History Data Methods

// history variables:
// h[0=HOG_LAMBDA] is max lambda = lambdaw (it is <0 until initiation)
// h[1=HOG_DW] is cumulative work
// h[2=HOG_UN] is nCod stored to get dun on next step
// h[3=HOG_UT] is tCod stored to get dut on next step
char *CoupledSawTooth::InitHistoryData(char *pchr,MPMBase *mptr)
{
    numTractionHistory = 4;
	double *p = CreateAndZeroDoubles(pchr,numTractionHistory);
    p[HOG_LAMBDA] = -1.;
	return (char *)p;
}

#pragma mark CoupledSawTooth::Traction Law

// Traction law - assume trianglar shape with unloading down slope back to the origin
void CoupledSawTooth::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,Vector *n,Vector *t,double area)
{
    double Tn=0.,Tt=0.;
	double *h = (double *)cs->GetHistoryData();
	
	// get theta in current deformation state
	double lambda;
	double unbar=0.,utbar;
	double sinTheta=0.,cosTheta=1.;
    double dun=0.,dut;
	if(nCod>0.)
	{	unbar = nCod/delIc;
		utbar = tCod/delIIc;
		lambda = sqrt(unbar*unbar + utbar*utbar);
        // This is finding the Hogberg theta
		sinTheta = unbar/lambda;
		cosTheta = utbar/lambda;
		dun = nCod-h[HOG_UN];
		h[HOG_UN] = nCod;
	}
	else
	{   // we assume nCod=0, unbar=0, dun=0, theta=0
        utbar = tCod/delIIc;
		lambda = fabs(utbar);
		h[HOG_UN] = 0.;
	}
    dut = tCod-h[HOG_UT];
    h[HOG_UT] = tCod;
    double sinTheta2 = sinTheta*sinTheta, cosTheta2 = cosTheta*cosTheta;
	
	// Find current lambda_p(theta)
	double lampTheta = 1./sqrt(sinTheta2/unpbar2 + cosTheta2/utpbar2);
    
    // check for damage initiation
    if(h[HOG_LAMBDA]<0.)
    {   if(lambda<=lampTheta)
        {   // the cohesive zone has not initiated damage yet
            // Find tractions and return
            if(nCod>0.) Tn = stress1*unbar/lampTheta;
            Tt = stress2*utbar/lampTheta;
            
            // work increment (loading and unloading)
            // Note energy error here is lampTheta depends on theta and theta changes in this phase
            h[HOG_DW] += Tn*dun + Tt*dut;
            
            // force is traction times area projected onto plane of unit vectors (units F)
            // tract = -area*(Tn*n + Tt*t)
            // In 2D, if t=(dx,dy), then n=(-dy,dx)
            cs->tract.x = -area*(Tn*n->x + Tt*t->x);
            cs->tract.y = -area*(Tn*n->y + Tt*t->y);
            cs->tract.z = -area*(Tn*n->z + Tt*t->z);

            return;
        }
        
        // damage has initiated, set h[HOG_LAMBDA] to initiation lambda and proceed to damage evolution
        h[HOG_LAMBDA] = lampTheta;
    }
	
	// If damage evolution calculate energy dissipation
	if(lambda > h[HOG_LAMBDA])
	{	// truncate to maximum value
        double dlambda;
        if(lambda>=1.)
        {   dlambda = 1.-h[HOG_LAMBDA];
            h[HOG_LAMBDA] = 1.;
        }
        else
        {   dlambda = lambda-h[HOG_LAMBDA];
            h[HOG_LAMBDA] = lambda;
        }
		
		// Energy dissipation
        if(model == TRIANGULARTRACTIONMATERIAL)
        {   dlambda /= (1.-lampTheta);
            cs->czmdG.x += JIc*sinTheta2*dlambda;
            cs->czmdG.y += JIIc*cosTheta2*dlambda;
        }
        else
        {   double arg = exp(-alpha*(h[HOG_LAMBDA]-lampTheta)/(1.-lampTheta));
            dlambda *= romexpalpha*(arg*(1+alpha*h[HOG_LAMBDA]/(1.-lampTheta))-expalpha);
            cs->czmdG.x += 0.5*stress1*delIc*sinTheta2*dlambda;
            cs->czmdG.y += 0.5*stress1*delIc*cosTheta2*dlambda;
        }
        cs->czmdG.z = 1.;
	}
	
	// Did it fail?
	if(lambda >= 1.)
    {   cs->SetMatID(0);                        // now debonded
		// calculate mode mixity
		ReportDebond(mtime,cs,cs->czmdG.x/(cs->czmdG.x+cs->czmdG.y),cs->czmdG.x+cs->czmdG.y);
		cs->tract.x = 0.;
		cs->tract.y = 0.;
		cs->tract.z = 0.;
		return;
	}
    
    // tractions - either elastic or after updating h[HOG_LAMBDA] in damage evolution
    double Hlambda;
    if(model == TRIANGULARTRACTIONMATERIAL)
        Hlambda = (1.-h[HOG_LAMBDA])/(1.-lampTheta);
    else
    {   double arg = exp(-alpha*(h[HOG_LAMBDA]-lampTheta)/(1.-lampTheta));
        Hlambda =romexpalpha*(arg-expalpha);
    }
	
	// Get S(lamw)/lamw (the fmin makes sure never gets stiffer than undamaged zone)
    double SLamLam = fmin(Hlambda/h[HOG_LAMBDA],1./lampTheta);
	
	// get current tractions
    if(nCod>0.) Tn = SLamLam*stress1*unbar;
    Tt = SLamLam*stress2*utbar;
	
	// work increment (loading and unloading)
	h[HOG_DW] += Tn*dun + Tt*dut;
	
	// force is traction times area projected onto plane of unit vectors (units F)
	// tract = -area*(Tn*n + Tt*t)
	// In 2D, if t=(dx,dy), then n=(-dy,dx)
	cs->tract.x = -area*(Tn*n->x + Tt*t->x);
	cs->tract.y = -area*(Tn*n->y + Tt*t->y);
	cs->tract.z = -area*(Tn*n->z + Tt*t->z);	
}

// Return current traction law work energy (Int T.du).
//	This energy is needed for J integral (and only used in J Integral)
// units of F/L
double CoupledSawTooth::CrackWorkEnergy(CrackSegment *cs,double nCod,double tCod)
{	double *h = (double *)cs->GetHistoryData();
	return h[HOG_DW];
}

// Return mode I and II energy that has been released by current segment.
void CoupledSawTooth::CrackDissipatedEnergy(CrackSegment *cs,double &GI,double &GII)
{	GI = cs->czmdG.x;
	GII = cs->czmdG.y;
}

#pragma mark CoupledSawTooth::Accessors

// return material type
const char *CoupledSawTooth::MaterialType(void) const { return "Hogberg Mixed-Mode Cohesive Law"; }

