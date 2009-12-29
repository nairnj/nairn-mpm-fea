/********************************************************************************
    TransIsotropic.cpp
    NairnMPM
    
    Created by John Nairn on Tue Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/TransIsotropic.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#ifdef MPM_CODE
	#include "Custom_Tasks/ConductionTask.hpp"
	#include "Custom_Tasks/DiffusionTask.hpp"
	#include "MPM_Classes/MPMBase.hpp"
#endif

#pragma mark TransIsotropic::Constructors and Destructors

// Constructors
TransIsotropic::TransIsotropic() {}

// Constructors
TransIsotropic::TransIsotropic(char *matName,int matID) : Elastic(matName)
{
    int i;
    
	tiType=matID;
    for(i=0;i<ORTHO_PROPS;i++)
        read[i]=0;

#ifdef MPM_CODE
	diffA=0.;
	diffT=0.;
	kcondA=0.;
	kcondT=0.;
#endif
	betaA=0.;
	betaT=0.;
}

#pragma mark TransIsotropic::Initialization

// print mechanical properties to output window
void TransIsotropic::PrintMechanicalProperties(void)
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

#ifdef MPM_CODE
// print transport properties to output window
void TransIsotropic::PrintTransportProperties(void)
{
	// Diffusion constants
	if(DiffusionTask::active)
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
	{	PrintProperty("ka",kcondA,"W/(m-K)");
		PrintProperty("kt",kcondT,"W/(m-K)");
		PrintProperty("Cp",1000.*heatCapacity,"J/(kg-K)");
		cout << endl;
	}
}
#endif

// Read material properties
char *TransIsotropic::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"EA")==0)
    {	read[EA_PROP]=1;
        return((char *)&EA);
    }
    
    else if(strcmp(xName,"ET")==0)
    {	read[ET_PROP]=1;
        return((char *)&ET);
    }
    
    else if(strcmp(xName,"GA")==0)
    {	read[GA_PROP]=1;
        return((char *)&GA);
    }
    
    else if(strcmp(xName,"GT")==0)
    {	read[GT_PROP]=1;
        return((char *)&GT);
    }
    
    else if(strcmp(xName,"nuA")==0)
    {	read[NUA_PROP]=1;
        return((char *)&nuA);
    }
    
    else if(strcmp(xName,"nuT")==0)
    {	read[NUT_PROP]=1;
        return((char *)&nuT);
    }
    
    else if(strcmp(xName,"alphaA")==0)
    {	read[AA_PROP]=1;
        return((char *)&aA);
    }
    
    else if(strcmp(xName,"alphaT")==0)
    {	read[AT_PROP]=1;
        return((char *)&aT);
    }

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
        return((char *)&kcondA);
    
    else if(strcmp(xName,"kCondT")==0)
        return((char *)&kcondT);
#endif
    
    return MaterialBase::InputMat(xName,input);
}

// calculate properties used in analyses
const char *TransIsotropic::VerifyProperties(int np)
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

    // set properties
	const char *err;
    if(MaterialTag()==TRANSISO1)
    {	err=SetAnalysisProps(np,1.e6*ET,1.e6*ET,1.e6*EA,nuT,ET*nuA/EA,ET*nuA/EA,
                1.e6*GT,1.e6*GA,1.e6*GA,1.e-6*aT,1.e-6*aT,1.e-6*aA,betaT*concSaturation,betaT*concSaturation,betaA*concSaturation);
    }
    else
    {	err=SetAnalysisProps(np,1.e6*ET,1.e6*EA,1.e6*ET,ET*nuA/EA,nuT,nuA,
                1.e6*GA,1.e6*GA,1.e6*GT,1.e-6*aT,1.e-6*aA,1.e-6*aT,betaT*concSaturation,betaA*concSaturation,betaT*concSaturation);
    }
	if(err!=NULL) return err;
	
	// superclass call
	return MaterialBase::VerifyProperties(np);
}

#pragma mark TransIsotropic::Methods

// fill in stiffness matrix if necessary
// makeSpecific divides by density, but only works in MPM code
// Used by TranIsoptropic 1 and 2 and by Orthotropic
// angle is in radians
void TransIsotropic::LoadMechProps(int makeSpecific,double angle,int np)
{
    if(MaterialTag()==TRANSISO1)
	{	// isotropic in the plane thus force angle to be zero
    	if(hasMatProps) return;
        angle=0.;
    }
    else
    {	if(hasMatProps && DbleEqual(angle,lastMatAngle)) return;
    }
    lastMatAngle=angle;
    hasMatProps=TRUE;

    // If angle not zero do rotation
    if(!DbleEqual(angle,0.))
	{	// analysis axes are ccw from material axes (note the -sin(ang))
        double cs=cos(angle);
        double c2=cs*cs;
        double c3=c2*cs;
        double c4=c2*c2;
        double sn=-sin(angle);
        double s2=sn*sn;
        double s3=s2*sn;
        double s4=s2*s2;
        double c2s2=c2*s2;
        double cs3=cs*s3;
        double c3s=c3*sn;
        double c2ms2=c2-s2;
		double cssn=cs*sn;

        // Rotated stiffness matrix ccw to return to analysis axes
        mdm[1][1]=c4*C11+2.*c2s2*(C12+2.*C66)+s4*C22;
        mdm[1][2]=c2s2*(C11+C22-4*C66)+(c4+s4)*C12;
        mdm[1][3]=c3s*(C11-C12-2*C66)+cs3*(C12-C22+2.*C66);
        mdm[2][2]=s4*C11+2.*c2s2*(C12+2.*C66)+c4*C22;
        mdm[2][3]=cs3*(C11-C12-2.*C66)+c3s*(C12-C22+2.*C66);
        mdm[3][3]=c2s2*(C11+C22-2.*C12)+c2ms2*c2ms2*C66;

#ifdef MPM_CODE
		mdm[4][1]=c2*C13+s2*C23;
		mdm[4][2]=s2*C13+c2*C23;
		mdm[4][3]=cssn*(C13-C23);
		mdm[4][4]=C33;
		
		me0[5]=c2*prop1+s2*prop2;
		me0[6]=s2*prop1+c2*prop2;
		me0[7]=2.*cssn*(prop1-prop2);
		
		// concentration strains
		mc0[1]=c2*CME1+s2*CME2;
		mc0[2]=s2*CME1+c2*CME2;
		mc0[3]=2.*cssn*(CME1-CME2);
		mc0[4]=CME3;
#endif

		// initial strains - all thermal and strain per temperature change
		me0[1]=c2*CTE1+s2*CTE2;
		me0[2]=s2*CTE1+c2*CTE2;
		me0[3]=2.*cssn*(CTE1-CTE2);
		me0[4]=CTE3;
		
    }
    
    else
    {	// Stiffness matrix
        mdm[1][1]=C11;
        mdm[1][2]=C12;
        mdm[1][3]=0.;
        mdm[2][2]=C22;
        mdm[2][3]=0.;
        mdm[3][3]=C66;
		
#ifdef MPM_CODE
		mdm[4][1]=C13;
		mdm[4][2]=C23;
		mdm[4][3]=0.;
		mdm[4][4]=C33;
		
		me0[5]=prop1;
		me0[6]=prop2;
		me0[7]=0.;
		
		// concentration strains 
		mc0[1]=CME1;
		mc0[2]=CME2;
		mc0[3]=0.;
		mc0[4]=CME3;
#endif

		// initial strains - all thermal and strain per temperature change
		me0[1]=CTE1;
		me0[2]=CTE2;
		me0[3]=0.;
		me0[4]=CTE3;
    }

#ifdef MPM_CODE
    // for MPM (units N/m^2 cm^3/g)
    if(makeSpecific)
    {	double rrho=1./rho;
    	mdm[1][1]*=rrho;
        mdm[1][2]*=rrho;
        mdm[1][3]*=rrho;
        mdm[2][2]*=rrho;
        mdm[2][3]*=rrho;
        mdm[3][3]*=rrho;
		if(np==PLANE_STRAIN_MPM)
    	{	mdm[4][1]*=rrho;
			mdm[4][2]*=rrho;
			mdm[4][3]*=rrho;
			mdm[4][4]*=rrho;
		}
    }
#endif

    // Fill bottom half of stiffness matrix
    mdm[2][1]=mdm[1][2];
    mdm[3][1]=mdm[1][3];
    mdm[3][2]=mdm[2][3];
	
#ifdef FEA_CODE
    /* For Plane Strain analysis, save term for strain energy (4 = E3*a3^2*T^2))
          or energy +=0.5*mdm[4][4]*volume*(delta T)^2 */
    if(np==PLANE_STRAIN)
        mdm[4][4]=prop3*CTE3;

    /* For axisymmetric, put theta direction properties in 4th
			column and row */
    else if(np==AXI_SYM)
    {	mdm[4][1]=mdm[1][4]=C13;
        mdm[4][2]=mdm[2][4]=C23;
        mdm[4][3]=mdm[3][4]=0.;
        mdm[4][4]=C33;
    }
#endif
}

#ifdef MPM_CODE
// fill in transport tensor if necessary
// Used by TranIsoptropic 1 and 2 and by Orthotropic
void TransIsotropic::LoadTransportProps(MPMBase *mptr)
{
	double angle=mptr->GetRotation();
	
    if(MaterialTag()==TRANSISO1)
	{	// isotropic in the plane thus force angle to be zero
    	if(hasTransProps) return;
        angle=0.;
    }
    else
    {	if(hasTransProps && DbleEqual(angle,lastTransAngle)) return;
    }
    lastTransAngle=angle;
    hasTransProps=TRUE;

    // If angle not zero do rotation
    if(!DbleEqual(angle,0.))
	{	// analysis axes are ccw from material axes (note the -sin(ang))
        double cs=cos(angle);
        double c2=cs*cs;
        double sn=-sin(angle);
        double s2=sn*sn;
        double cssn=cs*sn;

		// diffusion and conductivity tensors
		diffusionTensor.xx=diffA*s2 + diffT*c2;
		diffusionTensor.yy=diffA*c2 + diffT*s2;
		diffusionTensor.xy=(diffT-diffA)*cssn;
		kCondTensor.xx=kcondA*s2 + kcondT*c2;
		kCondTensor.yy=kcondA*c2 + kcondT*s2;
		kCondTensor.xy=(kcondT-kcondA)*cssn;
    }
    
    else
	{	if(MaterialTag()==TRANSISO1)
		{	diffusionTensor.xx=diffT;
			diffusionTensor.yy=diffT;
			kCondTensor.xx=kcondT;
			kCondTensor.yy=kcondT;
		}
		else
		{	diffusionTensor.xx=diffT;
			diffusionTensor.yy=diffA;
			kCondTensor.xx=kcondT;
			kCondTensor.yy=kcondA;
		}
		diffusionTensor.xy=0.;
		kCondTensor.xy=0.;
    }
}
#endif


#pragma mark TransIsotropic::Accessors

// Return the material tag
int TransIsotropic::MaterialTag(void) { return tiType; }

// return material type
const char *TransIsotropic::MaterialType(void)
{	if(MaterialTag()==TRANSISO1)
        return "Tranversely isotropic (A axis normal to x-y plane)";
    else
        return "Tranversely isotropic (A axis in x-y plane)";
}

#ifdef MPM_CODE
/*	calculate maximum wave speed in mm/sec (moduli in MPa, rho in g.cm^3
        TRANSISO1
            wave speeds are GT/rho and (KT+GT)/rho - return larger one
        TRANSISO2
            wave speeds are GT/rho, GA/rho, (KT+GT)/rho, and (EA + 4KT nuA^2)/rho
            return largest (assumes shear ones are not largest)
*/
double TransIsotropic::WaveSpeed(void)
{
    if(MaterialTag()==TRANSISO1)
        return sqrt(1.e9*(KT+GT)/rho);
    else
        return sqrt(1.e9*fmax(KT+GT,EA+4.*KT*nuA*nuA)/rho);
}

// maximum diffusion coefficient in cm^2/sec
double TransIsotropic::MaximumDiffusion(void) { return max(diffA,diffT)/100.; }

// maximum diffusivity in cm^2/sec
double TransIsotropic::MaximumDiffusivity(void) { return max(kcondA,kcondT)/(rho*heatCapacity*100.); }

// remove when 3D law written
bool TransIsotropic::ThreeDMaterial(void) { return false; }
#endif


