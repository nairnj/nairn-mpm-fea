//
//  CamClayModified.cpp
//  NairnMPM
//
//  Created by Raydel Lorenzo on 6/11/12.
//  Copyright (c) 2012 __Geotecnia-UNB__. All rights reserved.
//

#include <iostream>
#include <math.h>

/********************************************************************************
 VonMisesHardening.cpp
 NairnMPM
  ********************************************************************************/

#include "Materials/CamClayModified.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/ElasticMatrixCamClay.cpp"

#pragma mark CamClayModified::Constructors and Destructors

// Constructors
CamClayModified::CamClayModified()
{
}

// Constructors
CamClayModified::CamClayModified(char *matName) : GeoMaterials(matName)
{
    // default values
    popCC=0.;
    kapaCC=0.;
    lamdaCC=0.;
    Mcc=0.;
    e0=0.;
}

#pragma mark CamClayModified::Initialization

// print to output window
void CamClayModified::PrintMechanicalProperties(void)
{	
    PrintProperty("po",popCC,"");
    PrintProperty("M",Mcc,"");
    PrintProperty("lamda",lamdaCC,"");
    PrintProperty("kapa",kapaCC,"");
    PrintProperty("e0",e0,"");
    
    IsotropicMat::PrintMechanicalProperties();
	PrintYieldProperties();
}

// print just yield properties to output window
void CamClayModified::PrintYieldProperties(void)
{
}

// Read material properties
char *CamClayModified::InputMat(char *xName,int &input)
{
    // 
    if(strcmp(xName,"phi")==0)
    {   input=DOUBLE_NUM;
        readphi=true;
        return((char *)&phi);
    }
    // Slope of the curve void ratio vs isotropic stress in the elasto-plastic part.
    else if(strcmp(xName,"lamdaCC")==0)
    {   input=DOUBLE_NUM;
        readlamdaCC=true;
        return((char *)&lamdaCC);
    }
    // Slope of the curve void ratio vs isotropic stress in the elastic part.
    else if(strcmp(xName,"kapaCC")==0)
    {   input=DOUBLE_NUM;
        readkapaCC=true;
        return((char *)&kapaCC);
    }
    // Slope of the critical state line
    else if(strcmp(xName,"Mcc")==0)
    {   input=DOUBLE_NUM;
        readMcc=true;
        return((char *)&Mcc);
    }
    // void ratio
    else if(strcmp(xName,"e0")==0)
    {   input=DOUBLE_NUM;
        reade0=true;
        return((char *)&e0);
    }
    return(GeoMaterials::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *CamClayModified::VerifyProperties(int np)
{
	// check properties
	if(!readlamdaCC) return "The lamda is missing";
    
	if(!readkapaCC) return "The kapa is missing";
    
    if(!readphi) return "The phi is missing";
    
    if(!reade0) return "The e0 is missing";
        
	// call super class
	return GeoMaterials::VerifyProperties(np);
}

// Constant reduced properties used in constitutive law
void CamClayModified::InitialLoadMechProps(int makeSpecific,int np)
{
            
    // reduced properties
    Mccred = 6 * sin(phi * PI_CONSTANT / 180.) / (3 - sin(phi * PI_CONSTANT  / 180.));
    // reduced properties
    lamdaCCred = lamdaCC /rho; 
    // reduced properties
    kapaCCred = kapaCC / rho;
    // reduced properties
    popCCred = popCC * 1000. / rho;     // in KPa and here transform to N/m^2 equal to stress. pre-consolidation stress
	
	// beta and npow (if used) are dimensionless
	// if either is entered, Ep and ET are ignored
	
	GeoMaterials::InitialLoadMechProps(makeSpecific,np);
}

#pragma mark CamClayModified::Methods


//Get po of Cam clay for each state, here -dvplas because dvplas (+) is compresion in the expression of hardening and pop is the initial po
//computed from OCR and the initial state
double CamClayModified::GetPoActual(MPMBase *mptr,int np)
{
    double p0act, pop;
    Vector si;

    si.y = -(yref - mptr->origpos.y) * rho * 10.;
    si.x = (si.y * k0);
    si.z = (si.y * k0);
    

    pop = -(si.x + si.y + si.z)/3. * OCR / rho;                                 // in N/m^2, pre-consolidation stress and reduced propertie.

    if (pop<10000.) pop=10000.;
        
    p0act = pop * pow(EULER, (1.+ e0) * (-dvplas)/(lamdaCCred - kapaCCred));
    
    return p0act;
}
 

// Return value of f for current conditions
double CamClayModified::GetFtrial(MPMBase *mptr,int np, Tensor *st)
{
    double f, poact, pmedia, qdesv;

    poact = GetPoActual(mptr,np);
    pmedia = GetPcc(mptr, np, st);
    qdesv = Getqcc(mptr,np,st);

    f = (Mccred * Mccred * pmedia * pmedia) - (Mccred * Mccred * poact * pmedia) + (qdesv * qdesv);

    return f;
}

// Return mean stress
double CamClayModified::GetPcc(MPMBase *mptr,int np,Tensor *st)
{   
    double p;
    
    p=((st->xx + st->yy + st->zz)/3.);
    
    if (p<10000.)
    {
        p=10000.;
    }
    return p;
}

// Return elastic stress tensor for trial update
Tensor CamClayModified::GetElasticIncrement(MPMBase *mptr, int np, Tensor * sti, double dxx, double dyy, double dzz, double dxy, double dyz,double dxz)
{
    double pm;
    Tensor ds;
    
    pm = GetPcc(mptr,np,sti);
    matrix De(kapaCC,pm,nu,e0,rho);      //create and object of matrix in ElasticMatrixCamClay.cpp

    //elastic stress increment
    
    ds.xx = De.get(0,0) * dxx + De.get(0,1) * dyy + De.get(0,2) * dzz;
    ds.yy = De.get(1,0) * dxx + De.get(1,1) * dyy + De.get(1,2) * dzz;
    ds.zz = De.get(2,0) * dxx + De.get(2,1) * dyy + De.get(2,2) * dzz;
    ds.xy = De.get(3,3) * dxy;
    ds.yz = De.get(4,4) * dyz;
    ds.xz = De.get(5,5) * dxz;
        
    return ds;
}


// Return q do cam clay sqrt(3J2)
double CamClayModified::Getqcc(MPMBase *mptr,int np,Tensor *st)
{
    double q;
    q = sqrt(3. * (1./6. * (pow(st->xx-st->yy,2.) + pow(st->yy-st->zz,2.) + pow(st->xx-st->zz,2.)) + st->xy * st->xy + st->yz * st->yz + st->xz * st->xz));
    return q;
}

//Return dpds to use in dfdsigma
Tensor CamClayModified::GetDpDsigma()
{
    Tensor dpds;
    dpds.xx=1./3.;
    dpds.yy=1./3.;
    dpds.zz=1./3.;
    dpds.xy=0.;
    dpds.yz=0.;
    dpds.xz=0.;
    
    return dpds;
}
//Return dqds to use in dfdsigma
Tensor CamClayModified::GetDqdsigma(Tensor *st)
{
    Tensor dqds;
    double aux;
    
    aux=(pow(st->xx-st->yy,2.) + pow(st->yy-st->zz,2.) + pow(st->xx-st->zz,2.)) + 6. * st->xy * st->xy + 6. * st->yz * st->yz + 6. * st->xz * st->xz;
    
    dqds.xx=1./sqrt(2.) * pow(aux,-1./2.) * (2. * st->xx - st->yy - st->zz);
    dqds.yy=1./sqrt(2.) * pow(aux,-1./2.) * (2. * st->yy - st->xx - st->zz);
    dqds.zz=1./sqrt(2.) * pow(aux,-1./2.) * (2. * st->zz - st->xx - st->yy);
    dqds.xy=6./sqrt(2.) * pow(aux,-1./2.) * st->xy;
    dqds.yz=6./sqrt(2.) * pow(aux,-1./2.) * st->yz;
    dqds.xz=6./sqrt(2.) * pow(aux,-1./2.) * st->xz;

    return dqds;
}

//Return derivative of df/dsigma in a vector, declaration tensor becouse have 6 components.
Tensor CamClayModified::GetDfDsigmaGeo(MPMBase *mptr, int np, Tensor *st)
{
    double p,q,poact,dfdp,dfdq;
    Tensor dfds;
    Tensor dpds = GetDpDsigma();
    Tensor dqds = GetDqdsigma(st);
    p = GetPcc(mptr,np,st);
    q = Getqcc(mptr,np,st);
    poact = GetPoActual(mptr,np);
    dfdp = (2. * p - poact) * Mccred * Mccred;
    dfdq = 2. * q;
                       
    dfds.xx = dfdp * dpds.xx + dfdq * dqds.xx;
    dfds.yy = dfdp * dpds.yy + dfdq * dqds.yy;
    dfds.zz = dfdp * dpds.zz + dfdq * dqds.zz;
    dfds.xy = dfdp * dpds.xy + dfdq * dqds.xy;
    dfds.yz = dfdp * dpds.yz + dfdq * dqds.yz;
    dfds.xz = dfdp * dpds.xz + dfdq * dqds.xz;
    
    return dfds;    
}

//Return Derivative F/volumetric plastic deformation
double CamClayModified::GetDfDdefvp(MPMBase * mptr, int np, Tensor * st)
{
    double pm, pact;
    pm = GetPcc(mptr,np,st);
    pact = GetPoActual(mptr,np);

    return -(Mccred * Mccred * pm * pact * (1.+ e0) / (lamdaCCred - kapaCCred));
    
}


//Return Plastic Multiplier
double CamClayModified::GetMultPlast(MPMBase *mptr,int np, Tensor * st,double dxx, double dyy, double dzz, double dxy, double dyz,double dxz)
{
    int i;
    Tensor dfdss;
    double numerador, BDA, dfdevp, BD[6], Aii, multplas, pm;
    pm = GetPcc(mptr,np,st);
    
    matrix De(kapaCCred,pm,nu,e0,rho);
    
    dfdss = GetDfDsigmaGeo(mptr, np, st);
    dfdevp = GetDfDdefvp(mptr, np, st);
    
    for  (i=0;i<6;i++)
    {
        BD[i]=dfdss.xx * De.get(0,i) + dfdss.yy * De.get(1,i) + dfdss.zz * De.get(2,i) + dfdss.xy * De.get(3,i) + dfdss.yz * De.get(4,i) + dfdss.xz * De.get(5,i);
    }
    numerador = BD[0] * dxx + BD[1] * dyy + BD[2] * dzz + BD[3] * dxy + BD[4] * dyz + BD[5] * dxz;
    BDA = BD[0] * dfdss.xx + BD[1] * dfdss.yy + BD[2] * dfdss.zz + BD[3] * dfdss.xy + BD[4] * dfdss.yz + BD[5] * dfdss.xz;
    Aii = dfdss.xx + dfdss.yy + dfdss.zz;
    multplas = numerador / (BDA - dfdevp * Aii);
    
    return multplas;
}


#pragma mark CamClayModified::Accessors

// Return the material tag
int CamClayModified::MaterialTag(void) { return CAMCLAYMODIFIED; }

// return material type
const char *CamClayModified::MaterialType(void) { return "Modified Cam Clay"; }

