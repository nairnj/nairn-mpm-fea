/********************************************************************************
    Mooney.cpp
    NairnMPM
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Mooney.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
 
#pragma mark Mooney::Constructors and Destructors

// Constructors
Mooney::Mooney()
{
}
// Constructors with arguments 
Mooney::Mooney(char *matName) : RubberElastic(matName)
{
    int i;
    
	betaI=0.;
    for(i=0;i<MOONEY_PROPS;i++)
        read[i]=0;
	
}

#pragma mark Mooney::Initialization

// print mechanical properties output window
void Mooney::PrintMechanicalProperties(void)
{
	PrintProperty("C1",C1,"");
	PrintProperty("C2",C2,"");
	PrintProperty("a",aI,"");
	cout << endl;
}
	
// Read material properties
char *Mooney::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"C1")==0)
    {	read[C1_PROP]=1;
        return((char *)&C1);
    }
    
    else if(strcmp(xName,"C2")==0)
    {	read[C2_PROP]=1;
        return((char *)&C2);
    }
    
    else if(strcmp(xName,"alpha")==0)
    {	read[CTE_PROP]=1;
        return((char *)&aI);
    }
    
    else if(strcmp(xName,"beta")==0)
        return((char *)&betaI);
    
    return(MaterialBase::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *Mooney::VerifyProperties(int np)
{
	// check input
	if(np==PLANE_STRAIN || np==PLANE_STRAIN_MPM)
		return "Plane strain not valid for Mooney Rivlin Material";
    
    for(int i=0;i<MOONEY_PROPS;i++)
    {	if(!read[i])
			return "A required material property is missing";
    } 

	// call super class
	return MaterialBase::VerifyProperties(np);
}

// Private properties used in constitutive law
void Mooney::InitialLoadMechProps(int makeSpecific,int np)
{
	hasMatProps=TRUE;
	
	// C1 and C2 Specific units
	// for MPM (units N/m^2 cm^3/g)
	C1sp=C1*1e+06/rho;
	C2sp=C2*1e+06/rho;
	
	// nothing needed from superclass
}

#pragma mark Mooney::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
	Does not support thermal or moisture strains
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
								double delTime,int np)
{
	// Strain Increment
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx+=dvxx;
    ep->yy+=dvyy;
    ep->xy+=dvxy+dvyx;
	mptr->IncrementRotationStrain(dvyx-dvxy);

	// save initial stresses
	Tensor *sp=mptr->GetStressTensor();
	Tensor st0=*sp;

	// Simplifying large expressions
	double exx=ep->xx;
	double eyy=ep->yy; 
	double euxy=mptr->GetDuDy();
	double evxy=mptr->GetDvDx();
		
	//ezz as incompressible
	ep->zz=1/((1+exx)*(1+eyy))-1;

	// find stresses for total particle strain
	sp->xx=2*C1sp*(((1+exx)*(1+exx))+(euxy*euxy)-(1/((1+exx)*(1+exx)*(1+eyy)*(1+eyy))))
			+2*C2sp*((1+exx)*(1+exx)*(1+eyy)*(1+eyy)-((evxy*evxy+(1+eyy)*(1+eyy))/(pow((1+exx-euxy*evxy+eyy+exx*eyy),2))));

	sp->yy=2*C1sp*((evxy*evxy)-(1/((1+exx)*(1+exx)*(1+eyy)*(1+eyy)))+(1+eyy)*(1+eyy))
			+2*C2sp*((1+exx)*(1+exx)*(1+eyy)*(1+eyy)-((1+2*exx+exx*exx+euxy*euxy)/pow((1+exx-euxy*evxy+eyy+exx*eyy),2)));

	sp->xy=2*C1sp*(euxy+evxy+exx*evxy+euxy*eyy)
			+2*C2sp*((euxy+evxy+exx*evxy+euxy*eyy)/(pow((1+exx-euxy*evxy+eyy+exx*eyy),2)));


	// strain energy (by midpoint rule)
    double dgam=dvxy+dvyx;
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx
                            + (st0.yy+sp->yy)*dvyy
                            + (st0.xy+sp->xy)*dgam));
}

/* For 3D MPM analysis, take increments in strain and calculate new
    Particle: nothing yet
    dvij are (gradient rates X time increment) to give deformation gradient change
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
}

#pragma mark Mooney::Accessors

// Return the material tag
int Mooney::MaterialTag(void) { return MOONEYRIVLIN; }

/*	calculate wave speed in mm/sec (because C1 and C2 in MPa and rho in g/cm^3)
        For undeformed rubber use sqrt((2C1+2C2)/rho)
*/
double Mooney::WaveSpeed(bool threeD) { return sqrt(2.e9*(C1+C2)/rho); }

// return material type
const char *Mooney::MaterialType(void) { return "Mooney-Rivlin Rubber"; }

// remove when above 3D law written
bool  Mooney::ThreeDMaterial(void) { return false; }



