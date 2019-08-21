/********************************************************************************
    TransIsotropic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tue Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/TransIsotropic.hpp"
#ifdef MPM_CODE
	#include "Custom_Tasks/ConductionTask.hpp"
	#include "Custom_Tasks/DiffusionTask.hpp"
	#include "MPM_Classes/MPMBase.hpp"
	#include "Exceptions/CommonException.hpp"
#endif
#include "System/UnitsController.hpp"

#pragma mark TransIsotropic::Constructors and Destructors

// Constructor
TransIsotropic::TransIsotropic(char *matName,int matID) : Elastic(matName,matID)
{
    int i;
    
    for(i=0;i<ORTHO_PROPS;i++)
        read[i]=0;
    
    nuA=0.33;
    aA=40;
    aT=40;

#ifdef MPM_CODE
	// Here "A" is "y" when othotropic and "T" is "x"
	diffA=0.;
	diffT=0.;
	kCondA=0.;
	kCondT=0.;
#else
	// For FEA_CODE
    hasMatProps=FALSE;
#endif
	betaA=0.;
	betaT=0.;
}

#pragma mark TransIsotropic::Initialization

#ifdef MPM_CODE
// print mechanical properties to output window
void TransIsotropic::PrintMechanicalProperties(void) const
{
	PrintProperty("Ea",EA*UnitsController::Scaling(1.e-6),"");
	PrintProperty("Et",ET*UnitsController::Scaling(1.e-6),"");
	PrintProperty("va",nuA,"");
	PrintProperty("vt",nuT,"");
	cout << endl;
	
	PrintProperty("Gt",GT*UnitsController::Scaling(1.e-6),"");
	PrintProperty("Ga",GA*UnitsController::Scaling(1.e-6),"");
	cout << endl;
	
	PrintProperty("aa",aA,"");
	PrintProperty("at",aT,"");
	cout << endl;
}

// print transport properties to output window
void TransIsotropic::PrintTransportProperties(void) const
{
	// Diffusion constants
	if(DiffusionTask::HasFluidTransport())
	{	PrintProperty("Da",diffA,"mm^2/s");
		PrintProperty("Dt",diffT,"mm^2/s");
		PrintProperty("csat",concSaturation,"");
		cout << endl;
		
		PrintProperty("ba",betaA,"1/wt fr");
		PrintProperty("bt",betaT,"1/wt fr");
		cout << endl;
	}
	// Conductivity constants
	if(ConductionTask::active)
	{	PrintProperty("ka",rho*kCondA*UnitsController::Scaling(1.e-6),UnitsController::Label(CONDUCTIVITY_UNITS));
		PrintProperty("kt",rho*kCondT*UnitsController::Scaling(1.e-6),UnitsController::Label(CONDUCTIVITY_UNITS));
		PrintProperty("C",heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		cout << endl;
	}
}
#else

// print mechanical properties to output window
void TransIsotropic::PrintMechanicalProperties(void) const
{
	PrintProperty("Ea",EA,"");
	PrintProperty("Et",ET,"");
	PrintProperty("va",nuA,"");
	PrintProperty("vt",nuT,"");
	cout << endl;
	
	PrintProperty("Gt",GT,"");
	PrintProperty("Ga",GA,"");
	cout << endl;
	
	PrintProperty("aa",aA,"");
	PrintProperty("at",aT,"");
	cout << endl;
}

#endif

// Read material properties
char *TransIsotropic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;

#ifdef MPM_CODE
    if(strcmp(xName,"EA")==0)
    {	read[EA_PROP]=1;
        return UnitsController::ScaledPtr((char *)&EA,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"ET")==0)
    {	read[ET_PROP]=1;
        return UnitsController::ScaledPtr((char *)&ET,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"GA")==0)
    {	read[GA_PROP]=1;
        return UnitsController::ScaledPtr((char *)&GA,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"GT")==0)
    {	read[GT_PROP]=1;
        return UnitsController::ScaledPtr((char *)&GT,gScaling,1.e6);
    }
#else
    if(strcmp(xName,"EA")==0)
    {	read[EA_PROP]=1;
		return (char *)&EA;
    }
    
    else if(strcmp(xName,"ET")==0)
    {	read[ET_PROP]=1;
		return (char *)&ET;
    }
    
    else if(strcmp(xName,"GA")==0)
    {	read[GA_PROP]=1;
		return (char *)&GA;
    }
    
    else if(strcmp(xName,"GT")==0)
    {	read[GT_PROP]=1;
		return (char *)&GT;
    }
#endif
    
    else if(strcmp(xName,"nuT")==0)
    {	read[NUT_PROP]=1;
        return((char *)&nuT);
    }
    
    else if(strcmp(xName,"nuA")==0)
        return((char *)&nuA);
    
    else if(strcmp(xName,"alphaA")==0)
        return((char *)&aA);
    
    else if(strcmp(xName,"alphaT")==0)
        return((char *)&aT);

#ifdef MPM_CODE
    else if(strcmp(xName,"betaA")==0)
        return((char *)&betaA);
    
    else if(strcmp(xName,"betaT")==0)
        return((char *)&betaT);

    else if(strcmp(xName,"DA")==0)
        return((char *)&diffA);
    
    else if(strcmp(xName,"DT")==0)
        return((char *)&diffT);
		
	else if(strcmp(xName,"kCondA")==0)
		return UnitsController::ScaledPtr((char *)&kCondA,gScaling,1.e6);
	
	else if(strcmp(xName,"kCondT")==0)
		return UnitsController::ScaledPtr((char *)&kCondT,gScaling,1.e6);
#endif
    
    return Elastic::InputMaterialProperty(xName,input,gScaling);
}

// calculate properties used in analyses
const char *TransIsotropic::VerifyAndLoadProperties(int np)
{
    // finish input and verify all there
    if(!read[GT_PROP])
    {	GT=ET/(2.*(1.+nuT));
        read[GT_PROP]=1;
    }
    else if(!read[ET_PROP])
    {	ET=2.*GT*(1.+nuT);
        read[ET_PROP]=1;
    }
    else if(!read[NUT_PROP])
    {	nuT=ET/(2.*GT)-1.;
        read[NUT_PROP]=1;
    }
    else
		return "ET, nuT, and GT all specified. Only two allowed";

    int i;
    for(i=0;i<TRANS_PROPS;i++)
    {	if(!read[i])
			return "A required material property is missing";
    }
    KT=0.5/((1.-nuT)/ET - 2.*nuA*nuA/EA);
	nuAp = ET*nuA/EA;

    // set properties
	const char *err;
    if(AxialDirection()==AXIAL_Z)
    {	err=SetAnalysisProps(np,ET,ET,EA,nuT,nuAp,nuAp,
					GT,GA,GA,1.e-6*aT,1.e-6*aT,1.e-6*aA,
					betaT*concSaturation,betaT*concSaturation,betaA*concSaturation);
    }
    else
    {	err=SetAnalysisProps(np,ET,EA,ET,nuAp,nuT,nuA,
					GA,GA,GT,1.e-6*aT,1.e-6*aA,1.e-6*aT,
					betaT*concSaturation,betaA*concSaturation,betaT*concSaturation);
    }
	if(err!=NULL) return err;

#ifdef MPM_CODE
    // make conductivity (input as (N/(sec-K)) specific (N mm^3/(sec-K-g))
    kCondA /= rho;
    kCondT /= rho;
#endif
	
	// load elastic properties with constant values
	FillUnrotatedElasticProperties(&pr,np);
	
	// superclass call
	return MaterialBase::VerifyAndLoadProperties(np);
}

#ifdef MPM_CODE
// If needed, a material can initialize particle state
// For subclasses of TransIsotropic, rotation matrix is tracked in large rotation mode
//		and in small rotation is 3D
void TransIsotropic::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
	// always track in 3D so transport properties do not need second decomposition
	if(np==THREED_MPM)
		mptr->InitRtot(mptr->GetInitialRotation());
	
	// call super class
    Elastic::SetInitialParticleState(mptr,np,offset);
}
#endif

#pragma mark TransIsotropic::Methods

#ifdef MPM_CODE

// buffer size for mechanical properties
int TransIsotropic::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = 0;
 	if(useLargeRotation)
		return 0;
	else
		return sizeof(ElasticProperties);
}

// Isotropic material can use read-only initial properties
void *TransIsotropic::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	// if large rotation mode (in material axes) or if isotropic in 2D plane, use initial properties
	if(useLargeRotation || (AxialDirection()==AXIAL_Z && np!=THREED_MPM))
		return (void *)&pr;
	
	// create new elastic properties
	ElasticProperties *p = (ElasticProperties *)matBuffer;
	if(np!=THREED_MPM)
	{	double s,c;
		mptr->Get2DSinCos(&s,&c);
		FillElasticProperties2D(p,TRUE,s,c,np);
	}
	else
		FillElasticProperties3D(mptr,p,np);
	return p;
}

#pragma mark TransIsotropic::Methods (Small Rotation)

// Fill ElasticProperties variable with current particle state
void TransIsotropic::FillElasticProperties3D(MPMBase *mptr,ElasticProperties *p,int np) const
{
	// Full rotation using rotation matrix (save after calculationg for later use
	//    by aniostropic plastic or by transport properties)
	Matrix3 Rtot = mptr->GetRtot();
	mptr->SetRtot(Rtot);
	
	// Get 6X6 rotation matrix
	double R[6][6];
	Rtot.GetRStress(R);
	
	double C[6][6];
	C[0][0] = C11;
	C[0][1] = C12;
	C[0][2] = C13;
	C[1][0] = C12;
	C[1][1] = C22;
	C[1][2] = C23;
	C[2][0] = C13;
	C[2][1] = C23;
	C[2][2] = C33;
	C[3][3] = C44;
	C[4][4] = C55;
	C[5][5] = C66;
	
	int i,j,k,k3,l;
	for(i=0;i<6;i++)
	{	for(j=i;j<6;j++)
		{	double cterm = 0.;
			for(k=0;k<3;k++)
			{	k3 = k+3;
				cterm += R[i][k3]*R[j][k3]*C[k3][k3];
				for(l=0;l<3;l++)
					cterm += R[i][k]*R[j][l]*C[k][l];
			}
			p->C[i][j] = cterm;
		}
	}
	
	// make specific and symmetric
	double rrho=1./rho;
	for(i=0;i<6;i++)
	{	p->C[i][i]*=rrho;
		for(j=i+1;j<6;j++)
		{	p->C[i][j]*=rrho;
			p->C[j][i]=p->C[i][j];
		}
	}
	
	// thermal and moisture expansion
	
	// only these are needed and different from Rstress
	//Rtot.GetRStrain(R);
	R[3][0] = 2.*Rtot(1,0)*Rtot(2,0);
	R[4][0] = 2.*Rtot(0,0)*Rtot(2,0);
	R[5][0] = 2.*Rtot(0,0)*Rtot(1,0);
	R[3][1] = 2.*Rtot(1,1)*Rtot(2,1);
	R[4][1] = 2.*Rtot(0,1)*Rtot(2,1);
	R[5][1] = 2.*Rtot(1,1)*Rtot(0,1);
	R[3][2] = 2.*Rtot(2,2)*Rtot(1,2);
	R[4][2] = 2.*Rtot(2,2)*Rtot(0,2);
	R[5][2] = 2.*Rtot(0,2)*Rtot(1,2);
	
	p->alpha[0] = R[0][0]*CTE1 + R[0][1]*CTE2 + R[0][2]*CTE3;
	p->alpha[1] = R[1][0]*CTE1 + R[1][1]*CTE2 + R[1][2]*CTE3;
	p->alpha[2] = R[2][0]*CTE1 + R[2][1]*CTE2 + R[2][2]*CTE3;
	p->alpha[3] = R[3][0]*CTE1 + R[3][1]*CTE2 + R[3][2]*CTE3;
	p->alpha[4] = R[4][0]*CTE1 + R[4][1]*CTE2 + R[4][2]*CTE3;
	p->alpha[5] = R[5][0]*CTE1 + R[5][1]*CTE2 + R[5][2]*CTE3;
	
	p->beta[0] = R[0][0]*CME1 + R[0][1]*CME2 + R[0][2]*CME3;
	p->beta[1] = R[1][0]*CME1 + R[1][1]*CME2 + R[1][2]*CME3;
	p->beta[2] = R[2][0]*CME1 + R[2][1]*CME2 + R[2][2]*CME3;
	p->beta[3] = R[3][0]*CME1 + R[3][1]*CME2 + R[3][2]*CME3;
	p->beta[4] = R[4][0]*CME1 + R[4][1]*CME2 + R[4][2]*CME3;
	p->beta[5] = R[5][0]*CME1 + R[5][1]*CME2 + R[5][2]*CME3;

	/*
	p->alpha[0] = CTE1*cy2*cz2 + CTE2*(cx*s2z*sx*sy + cz2*sx2*sy2 + cx2*sz2) + CTE3*(-(cx*s2z*sx*sy) + cx2*cz2*sy2 + sx2*sz2);
	p->alpha[1] = CTE1*cy2*sz2 + CTE3*(cz2*sx2 + cx*s2z*sx*sy + cx2*sy2*sz2) + CTE2*(cx2*cz2 - cx*s2z*sx*sy + sx2*sy2*sz2);
	p->alpha[2] = CTE3*cx2*cy2 + CTE2*cy2*sx2 + CTE1*sy2;
	p->alpha[3] = -(CTE1*s2y*sz)/2. + CTE3*(cx*cy*cz*sx + (cx2*s2y*sz)/2.) + CTE2*(-(cx*cy*cz*sx) + (s2y*sx2*sz)/2.);
	p->alpha[4] = (CTE1*cz*s2y)/2. + CTE2*(-(cz*s2y*sx2)/2. - cx*cy*sx*sz) + CTE3*(-(cx2*cz*s2y)/2. + cx*cy*sx*sz);
	p->alpha[5] = -(CTE1*cy2*s2z)/2. + CTE3*((s2z*sx2)/2. - c2z*cx*sx*sy - (cx2*s2z*sy2)/2.) + CTE2*((cx2*s2z)/2. + c2z*cx*sx*sy - (s2z*sx2*sy2)/2.);
	
	p->beta[0] = CME1*cy2*cz2 + CME2*(cx*s2z*sx*sy + cz2*sx2*sy2 + cx2*sz2) + CME3*(-(cx*s2z*sx*sy) + cx2*cz2*sy2 + sx2*sz2);
	p->beta[1] = CME1*cy2*sz2 + CME3*(cz2*sx2 + cx*s2z*sx*sy + cx2*sy2*sz2) + CME2*(cx2*cz2 - cx*s2z*sx*sy + sx2*sy2*sz2);
	p->beta[2] = CME3*cx2*cy2 + CME2*cy2*sx2 + CME1*sy2;
	p->beta[3] = -(CME1*s2y*sz)/2. + CME3*(cx*cy*cz*sx + (cx2*s2y*sz)/2.) + CME2*(-(cx*cy*cz*sx) + (s2y*sx2*sz)/2.);
	p->beta[4] = (CME1*cz*s2y)/2. + CME2*(-(cz*s2y*sx2)/2. - cx*cy*sx*sz) + CME3*(-(cx2*cz*s2y)/2. + cx*cy*sx*sz);
	p->beta[5] =-(CME1*cy2*s2z)/2. + CME3*((s2z*sx2)/2. - c2z*cx*sx*sy - (cx2*s2z*sy2)/2.) + CME2*((cx2*s2z)/2. + c2z*cx*sx*sy - (s2z*sx2*sy2)/2.);
	*/
}

#else

// Fill the object in the material class
void TransIsotropic::LoadMechanicalPropertiesFEA(int makeSpecific,double angle,int np)
{	// if filled with same angle then no need to fill again
	if(hasMatProps && DbleEqual(angle,lastMatAngle)) return;
	lastMatAngle=angle;
	hasMatProps=TRUE;
	FillElasticProperties2D(&pr,FALSE,sin(angle),cos(angle),np);
}

#endif

// fill in stiffness matrix if necessary
// makeSpecific divides by density, but only used in MPM code
// Used by TranIsoptropic 2 and by Orthotropic
// sn = sin(angle) and cn = cos(angle) where angle is cw angle from current axes to the
//     material axes or ccw angle from material axes to the current axes
void TransIsotropic::FillElasticProperties2D(ElasticProperties *p,int makeSpecific,double sn,double cs,int np) const
{
    // If angle not zero do rotation
    if(!DbleEqual(sn,0.))
	{	// analysis axes are ccw from material axes (note the -sin(ang))
        double c2=cs*cs;
        double c3=c2*cs;
        double c4=c2*c2;
        double s2=sn*sn;
        double s3=s2*sn;
        double s4=s2*s2;
        double c2s2=c2*s2;
        double cs3=cs*s3;
        double c3s=c3*sn;
        double c2ms2=c2-s2;
		double cssn=cs*sn;

        // Rotated stiffness matrix ccw to return to analysis axes
        p->C[1][1]=c4*C11+2.*c2s2*(C12+2.*C66)+s4*C22;
        p->C[1][2]=c2s2*(C11+C22-4*C66)+(c4+s4)*C12;
        p->C[1][3]=c3s*(C11-C12-2*C66)+cs3*(C12-C22+2.*C66);
        p->C[2][2]=s4*C11+2.*c2s2*(C12+2.*C66)+c4*C22;
        p->C[2][3]=cs3*(C11-C12-2.*C66)+c3s*(C12-C22+2.*C66);
        p->C[3][3]=c2s2*(C11+C22-2.*C12)+c2ms2*c2ms2*C66;

#ifdef MPM_CODE
		p->C[4][1]=c2*C13+s2*C23;
		p->C[4][2]=s2*C13+c2*C23;
		p->C[4][3]=cssn*(C13-C23);
		p->C[4][4]=C33;
		
		// rotated S13,S23,and S36 for generalized plane stress or strain
		p->C[5][1]=S13*c2+S23*s2;
		p->C[5][2]=S13*s2+S23*c2;
		p->C[5][3]=2.*(S13-S23)*cssn;
		
		p->alpha[5]=c2*prop1+s2*prop2;
		p->alpha[6]=s2*prop1+c2*prop2;
		p->alpha[7]=2.*cssn*(prop1-prop2);
		
		// concentration strains
		p->beta[1]=c2*CME1+s2*CME2;
		p->beta[2]=s2*CME1+c2*CME2;
		p->beta[3]=2.*cssn*(CME1-CME2);
		p->beta[4]=CME3;
		
#endif
		
		// initial strains - all thermal and strain per temperature change
		p->alpha[1]=c2*CTE1+s2*CTE2;
		p->alpha[2]=s2*CTE1+c2*CTE2;
		p->alpha[3]=2.*cssn*(CTE1-CTE2);
		p->alpha[4]=CTE3;
		
    }
    
    else
    {	// Stiffness matrix
        p->C[1][1]=C11;
        p->C[1][2]=C12;
        p->C[1][3]=0.;
        p->C[2][2]=C22;
        p->C[2][3]=0.;
        p->C[3][3]=C66;
		
#ifdef MPM_CODE
		p->C[4][1]=C13;
		p->C[4][2]=C23;
		p->C[4][3]=0.;
		p->C[4][4]=C33;
		
		// rotated S13,S23,and S36 for generalized plane stress
		p->C[5][1]=S13;
		p->C[5][2]=S23;
		p->C[5][3]=0.;
		
		p->alpha[5]=prop1;
		p->alpha[6]=prop2;
		p->alpha[7]=0.;
		
		// concentration strains 
		p->beta[1]=CME1;
		p->beta[2]=CME2;
		p->beta[3]=0.;
		p->beta[4]=CME3;
#endif

		// initial strains - all thermal and strain per temperature change
		p->alpha[1]=CTE1;
		p->alpha[2]=CTE2;
		p->alpha[3]=0.;
		p->alpha[4]=CTE3;
    }

#ifdef MPM_CODE
    // for MPM (units N/m^2 mm^3/g)
    if(makeSpecific)
    {	double rrho=1./rho;
    	p->C[1][1]*=rrho;
        p->C[1][2]*=rrho;
        p->C[1][3]*=rrho;
        p->C[2][2]*=rrho;
        p->C[2][3]*=rrho;
        p->C[3][3]*=rrho;
		if(np==PLANE_STRAIN_MPM || np==AXISYMMETRIC_MPM)
    	{	p->C[4][1]*=rrho;
			p->C[4][2]*=rrho;
			p->C[4][3]*=rrho;
			p->C[4][4]*=rrho;
			// for generalized plane strain
			p->C[5][1]/=S33;
			p->C[5][2]/=S33;
			p->C[5][3]/=S33;
		}
		else if(np==PLANE_STRESS_MPM)
		{	// for generalized plane stress
			p->C[4][4]*=rrho;
			p->C[5][1]*=rho;
			p->C[5][2]*=rho;
			p->C[5][3]*=rho;
		}
    }
#endif

    // Fill bottom half of stiffness matrix
    p->C[2][1]=p->C[1][2];
    p->C[3][1]=p->C[1][3];
    p->C[3][2]=p->C[2][3];
	
#ifdef FEA_CODE
    /* For Plane Strain analysis, save term for strain energy (4 = E3*a3^2*T^2))
          or energy +=0.5*p->C[4][4]*volume*(delta T)^2 */
    if(np==PLANE_STRAIN)
        p->C[4][4]=prop3*CTE3;

    /* For axisymmetric, put theta direction properties in 4th
			column and row */
    else if(np==AXI_SYM)
    {	p->C[4][1]=p->C[1][4]=C13;
        p->C[4][2]=p->C[2][4]=C23;
        p->C[4][3]=p->C[3][4]=0.;
        p->C[4][4]=C33;
    }
#endif
}

#ifdef MPM_CODE

// Called before analysis, material can fill in things that never change during the analysis
// Note: no angle, because cannot depend on material angle
void TransIsotropic::FillTransportProperties(TransportProperties *t)
{
	if(AxialDirection()==AXIAL_Z)
	{	t->diffusionTensor.xx = diffT;
		t->diffusionTensor.yy = diffT;
		t->kCondTensor.xx = kCondT;
		t->kCondTensor.yy = kCondT;
	}
	else
	{	t->diffusionTensor.xx=diffT;
		t->diffusionTensor.yy=diffA;
		t->kCondTensor.xx = kCondT;
		t->kCondTensor.yy = kCondA;
	}
	t->diffusionTensor.zz = GetDiffZ();
	t->kCondTensor.zz = GetKcondZ();
	t->diffusionTensor.xy = 0.;
	t->diffusionTensor.xz = 0.;
	t->diffusionTensor.yz = 0.;
	t->kCondTensor.xy = 0.;
	t->kCondTensor.xz = 0.;
	t->kCondTensor.yz = 0.;
}

// fill in specific transport tensor if necessary
// Used by TranIsoptropic 1 and 2 and by Orthotropic
void TransIsotropic::GetTransportProps(MPMBase *mptr,int np,TransportProperties *t) const
{	
	if(np!=THREED_MPM)
	{	// if isotropic in 2D plane, use initial properties
		if(AxialDirection()==AXIAL_Z)
		{	*t = tr;
			return;
		}
		
		// get rotation matrix info
		double cs,sn;
		mptr->Get2DSinCos(&sn,&cs);

		// analysis axes are ccw from material axes
		double c2=cs*cs;
		double s2=sn*sn;
		double cssn=cs*sn;
		
		// diffusion and conductivity tensors = R.Tens.RT
		if(DiffusionTask::HasFluidTransport())
		{	t->diffusionTensor.xx = diffA*s2 + diffT*c2;
			t->diffusionTensor.yy = diffA*c2 + diffT*s2;
			t->diffusionTensor.xy = (diffT-diffA)*cssn;
		}
		if(ConductionTask::active)
		{	t->kCondTensor.xx = kCondA*s2 + kCondT*c2;
			t->kCondTensor.yy = kCondA*c2 + kCondT*s2;
			t->kCondTensor.xy = (kCondT-kCondA)*cssn;
		}
	}
	
	else
	{
		double R[6][6];
		Matrix3 *Rtot = mptr->GetRtotPtr();
		Rtot->GetRStress(R);
		
		if(DiffusionTask::HasFluidTransport())
		{	double diffz = GetDiffZ();
			t->diffusionTensor.xx = R[0][0]*diffT + R[0][1]*diffA + R[0][2]*diffz;
			t->diffusionTensor.yy = R[1][0]*diffT + R[1][1]*diffA + R[1][2]*diffz;
			t->diffusionTensor.zz = R[2][0]*diffT + R[2][1]*diffA + R[2][2]*diffz;
			t->diffusionTensor.yz = R[3][0]*diffT + R[3][1]*diffA + R[3][2]*diffz;
			t->diffusionTensor.xz = R[4][0]*diffT + R[4][1]*diffA + R[4][2]*diffz;
			t->diffusionTensor.xy = R[5][0]*diffT + R[5][1]*diffA + R[5][2]*diffz;
		}
		
		if(ConductionTask::active)
		{	double kz = GetKcondZ();
			t->kCondTensor.xx = R[0][0]*kCondT + R[0][1]*kCondA + R[0][2]*kz;
			t->kCondTensor.yy = R[1][0]*kCondT + R[1][1]*kCondA + R[1][2]*kz;
			t->kCondTensor.zz = R[2][0]*kCondT + R[2][1]*kCondA + R[2][2]*kz;
			t->kCondTensor.yz = R[3][0]*kCondT + R[3][1]*kCondA + R[3][2]*kz;
			t->kCondTensor.xz = R[4][0]*kCondT + R[4][1]*kCondA + R[4][2]*kz;
			t->kCondTensor.xy = R[5][0]*kCondT + R[5][1]*kCondA + R[5][2]*kz;
		}
	}
}

#endif

#pragma mark TransIsotropic::Accessors

// Return base axial direction
int TransIsotropic::AxialDirection(void) const { return materialID==TRANSISO1 ? AXIAL_Z : AXIAL_Y; }

// return material type
const char *TransIsotropic::MaterialType(void) const
{	if(AxialDirection()==AXIAL_Z)
        return "Tranversely isotropic (unrotated A axis along z axis)";
    else
        return "Tranversely isotropic (unrotated A axis along y axis)";
}

#ifdef MPM_CODE

/* Calculate maximum wave speed in mm/sec (moduli in MPa, rho in g/mm^3)
	TRANSISO1
		wave speeds are GT/rho and (KT+GT)/rho - return larger one
	TRANSISO2
		wave speeds are GT/rho, GA/rho, (KT+GT)/rho, and (EA + 4KT nuA^2)/rho
		return largest (assumes shear ones are not largest)
*/
double TransIsotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{
    if(AxialDirection()==AXIAL_Z && !threeD)
        return sqrt((KT+GT)/rho);
    else
        return sqrt(fmax(GA,fmax(KT+GT,EA+4.*KT*nuA*nuA))/rho);
}

// maximum diffusion coefficient in mm^2/sec (diff in mm^2/sec)
double TransIsotropic::MaximumDiffusion(void) const { return fmax(diffA,diffT); }

// maximum diffusivity in mm^2/sec
// specific k is nJ mm^2/(sec-K-g) and Cp is nJ/(g-K) so k/Cp = mm^2 /sec
double TransIsotropic::MaximumDiffusivity(void) const { return fmax(kCondA,kCondT)/heatCapacity; }

// diffusion and conductivity in the z direction
double TransIsotropic::GetDiffZ(void) const { return AxialDirection()==AXIAL_Z ? diffA : diffT; }
double TransIsotropic::GetKcondZ(void) const { return AxialDirection()==AXIAL_Z ? kCondA : kCondT; }

// not supported yet, need to deal with aniostropi properties
bool TransIsotropic::SupportsDiffusion(void) const
{
    return DiffusionTask::HasPoroelasticity() ? false : true;
}

#endif	// MPM_CODE


