/********************************************************************************
    Elastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Elastic.hpp"
#ifdef MPM_CODE
	#include "Global_Quantities/ThermalRamp.hpp"
#endif

#pragma mark Elastic::Constructors and Destructors

// Constructors
Elastic::Elastic() {}

Elastic::Elastic(char *matName) : MaterialBase(matName)
{
#ifdef MPM_CODE
	useLargeRotation=0;
#endif
}

#pragma mark Elastic::Initialization

// Read material properties
char *Elastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
#ifdef MPM_CODE
	if(strcmp(xName,"largeRotation")==0)
	{	input=INT_NUM;
		return((char *)&useLargeRotation);
	}
#endif
	
    return MaterialBase::InputMaterialProperty(xName,input,gScaling);
}

#pragma mark Elastic::Methods

// Fill ElasticProperties variable with current particle state
void Elastic::FillUnrotatedElasticProperties(ElasticProperties *p,int np)
{
#ifdef MPM_CODE
	bool makeSpecific = true;
	
	if(np==THREED_MPM)
	{	double rrho=1./rho;
		p->C[0][0]=C11*rrho;
		p->C[0][1]=C12*rrho;
		p->C[0][2]=C13*rrho;
		p->C[0][3]=p->C[0][4]=p->C[0][5]=0.;
		
		p->C[1][0]=p->C[0][1];
		p->C[1][1]=C22*rrho;
		p->C[1][2]=C23*rrho;
		p->C[1][3]=p->C[1][4]=p->C[1][5]=0.;
		
		p->C[2][0]=p->C[0][2];
		p->C[2][1]=p->C[1][2];
		p->C[2][2]=C33*rrho;
		p->C[2][3]=p->C[2][4]=p->C[2][5]=0.;
		
		p->C[3][0]=p->C[3][1]=p->C[3][2]=0.;
		p->C[3][3]=C44*rrho;
		p->C[3][4]=p->C[3][5]=0.;
		
		p->C[4][0]=p->C[4][1]=p->C[4][2]=p->C[4][3]=0.;
		p->C[4][4]=C55*rrho;
		p->C[4][5]=0.;
		
		p->C[5][0]=p->C[5][1]=p->C[5][2]=p->C[5][3]=p->C[5][4]=0.;
		p->C[5][5]=C66*rrho;
		
		// need p->alpha[] and p->beta[] for thermal and moisture expansion
		p->alpha[0]=CTE1;
		p->alpha[1]=CTE2;
		p->alpha[2]=CTE3;
		p->alpha[3]=p->alpha[4]=p->alpha[5]=0.;
		
		p->beta[0]=CME1;
		p->beta[1]=CME2;
		p->beta[2]=CME3;
		p->beta[3]=p->beta[4]=p->beta[5]=0.;
		
		return;
	}
#endif
    
    // Stiffness matrix
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
#endif
	
    // initial strains - all thermal and strain per temperature change
    p->alpha[1]=CTE1;
    p->alpha[2]=CTE2;
    p->alpha[3]=0.;
    p->alpha[4]=CTE3;
	
#ifdef MPM_CODE
	p->alpha[5]=prop1;
	p->alpha[6]=prop2;
	p->alpha[7]=0.;
	
    // initial strains for moisture
    p->beta[1]=CME1;
    p->beta[2]=CME2;
    p->beta[3]=0.;
    p->beta[4]=CME3;
	
    // for MPM (N/m^2 mm^3/g = (g-mm^2/sec^2)/g when props in MPa and rho in g/mm^3)
    if(makeSpecific)
    {	double rrho=1./rho;
    	p->C[1][1]*=rrho;
        p->C[1][2]*=rrho;
        p->C[2][2]*=rrho;
        p->C[3][3]*=rrho;
		if(np==PLANE_STRAIN_MPM || np==AXISYMMETRIC_MPM)
		{	p->C[4][1]*=rrho;
			p->C[4][2]*=rrho;
			p->C[4][3]*=rrho;
			p->C[4][4]*=rrho;
		}
    }
#endif
	
    // Fill bottom half of stiffness matrix
    p->C[2][1]=p->C[1][2];
    p->C[3][1]=0.;
    p->C[3][2]=0.;
	
#ifdef FEA_CODE
    /* For Plane Strain analysis, save term for strain energy (4 = E3*a3^2*T^2))
	 or energy +=0.5*p->C[4][4]*volume*(delta T)^2 */
    if(np==PLANE_STRAIN)
        p->C[4][4]=prop3*CTE3;
	
    /* For axisymmetric, put theta direction properties in 4th
	 column and row */
    else if(np==AXI_SYM)
    {	p->C[4][1]=C13;
        p->C[4][2]=C23;
        p->C[4][3]=0.;
        p->C[4][4]=C33;
		
        p->C[1][4]=p->C[4][1];
        p->C[2][4]=p->C[4][2];
        p->C[3][4]=p->C[4][3];
    }
#endif
}

/*  Store properties constants used to fill elastic properties later
    Also check for invalid properties
        Moduli now in Pa
        CTE in temp^-1
		CME in strain per concentration potential
    Pays attention to plane stress vs plane strain 
    
    Plane stress and plain strain
        Cij give upper half of symmetric stiffness matrix (used to be prop[0] to prop[5])
			(i,j=3 is shear for planar, but theta stress to axisymmetric)
		CTEi gives diagonal element of thermal expansion tensor (used to be prop[9] to prop[11])
        prop1 to prop3 used to find s(zz) in plane strain or e(zz) in plane stress in FEA
		prop1 and prop2 hold Poisson's rations need for plain strain ss(zz) calculation
    
    WARNING: No verification that have valid properties such as Poisson ratios
*/
const char *Elastic::SetAnalysisProps(int np,double e1,double e2,double e3,double v12,
        double v13,double v23,double g12,double g13,double g23,double a1,
        double a2,double a3,double beta1,double beta2,double beta3)
{
    double v32,v31,v21,temp,xx;
    
    // other Poisson ratios
    v32=v23*e3/e2;
    v31=v13*e3/e1;
    v21=v12*e2/e1;
    
    /* Verify legal material properties, i.e. positive moduli
            and legal Poisson ratios */
    if(e1<0 || e2<0. || e3<0. || g13<0. || g12<0. || g23<0.)
		return "Modulus less than zero";
    temp=(1.-v12*v21-v31*v13-v23*v32)/2.;
    if(v12*v21>1)
		return "Illegal Poisson's ratio (v12*v21>1)";
    else if(v23*v32>1.)
		return "Illegal Poisson's ratio (v23*v32>1)";
    else if(v13*v31>1.)
		return "Illegal Poisson's ratio (v13*v31>1)";
    else if(v21*v32*v13>temp)
		return "Illegal Poisson's ratio (v21*v32*v13>(1-v12*v21-v31*v13-v23*v32)/2)";
#ifdef MPM_CODE
    else if(rho<=0.0)
		return "Density less than or equal to zero";
#endif

    // 2D plane stress properties
    if(np==PLANE_STRESS || np==PLANE_STRESS_MPM)
    {	xx=1.-v12*v21;
        C11=e1/xx;
        C12=e2*v12/xx;
        C22=e2/xx;
        C66=g12;
        CTE1=a1;
        CTE2=a2;
        CTE3=a3;

#ifdef MPM_CODE
        CME1=beta1;
        CME2=beta2;
        CME3=beta3;
#else
		// for calculatting eps(zz) in FEA
        prop1=v13/e1;
        prop2=v23/e2;
        prop3=-a3;
#endif
		// for calculatting eps(zz) in MPM
		xx=1.-v13*v31-v23*v32-v12*v21-2.*v13*v32*v21;
        C13=-(v13+v12*v23)/(1.-v21*v12);			// equal to -C13/C33
        C23=-(v23+v21*v13)/(1.-v21*v12);			// equal to -C23/C33
        C33=e3*(1.-v21*v12)/xx;
		
#ifdef MPM_CODE
		// Cadota in nJ/(g-K^2)
        double C113D = e1*(1.-v23*v32)/xx;
        double C123D = e2*(v12+v13*v32)/xx;
        double C223D = e2*(1.-v13*v31)/xx;
		Cadota = (C113D*a1*a1+C223D*a2*a2+C33*a3*a3 + 2.*(C123D*a2*a1-C33*(C13*a3*a1+C23*a3*a2)))/rho;
#endif
    }
	
    // 2D plane strain properties
    else if(np==PLANE_STRAIN || np==PLANE_STRAIN_MPM)
    {	xx=1.-v13*v31-v23*v32-v12*v21-2.*v12*v23*v31;
        C11=e1*(1.-v23*v32)/xx;
        C12=e2*(v12+v13*v32)/xx;
        C22=e2*(1.-v13*v31)/xx;
        C66=g12;
        CTE1=a1+v31*a3;
        CTE2=a2+v32*a3;
        CTE3=a3;

#ifdef MPM_CODE
		// for calculating sigma(zz) in MPM - excess cte in reduced cte
		prop1=v31;
		prop2=v32;
		
        CME1=beta1+v31*beta3;
        CME2=beta2+v32*beta3;
        CME3=beta3;
#else
		// for calculating sigma(zz) in FEA
        prop1=-e3*v13/e1;
        prop2=-e3*v23/e2;
        prop3=e3*a3;
#endif
		
		// for calculating sigma(zz) in MPM
        C13=e3*(v13+v12*v23)/xx;
        C23=e3*(v23+v21*v13)/xx;
        C33=e3*(1.-v21*v12)/xx;
		
#ifdef MPM_CODE
		// Cadota in nJ/(g-K^2)
		Cadota = (C11*a1*a1+C22*a2*a2+C33*a3*a3+2.*(C12*a2*a1+C13*a3*a1+C23*a3*a2))/rho;
#endif
    }
    
    /* Full 3D stiffness matrix or axisymmetric matrix components
            For axisymmetric, E1 along r, E2 along z, and E3 along hoop */
    else if(np==AXI_SYM || np==THREED_MPM || np==AXISYMMETRIC_MPM)
    {	xx=1.-v13*v31-v23*v32-v12*v21-2.*v13*v32*v21;
        C11=e1*(1.-v23*v32)/xx;
        C12=e2*(v12+v13*v32)/xx;
        C13=e3*(v13+v12*v23)/xx;
        C22=e2*(1.-v13*v31)/xx;
        C23=e3*(v23+v21*v13)/xx;
        C33=e3*(1.-v21*v12)/xx;
		C66=g12;
		C44=g23;
		C55=g13;
        CTE1=a1;
        CTE2=a2;
        CTE3=a3;
#ifdef MPM_CODE
        CME1=beta1;
        CME2=beta2;
        CME3=beta3;
		// Cadota in nJ/(g-K^2)
		Cadota = (C11*a1*a1+C22*a2*a2+C33*a3*a3+2.*(C12*a2*a1+C13*a3*a1+C23*a3*a2)/rho);
#endif
    }
    
    return NULL;
}

#ifdef FEA_CODE
/* get transverse stress in plane strain analyses or get transverse strain
    in plane stress analysis (angle in radians)
*/
double Elastic::GetStressStrainZZ(double sxx,double syy,double sxy,double theTemp,double angle,int np)
{	
    double c,s,c2,s2;
	double deltaT=theTemp;
    
    if(DbleEqual(prop1,prop2))
		return -prop1*sxx-prop2*syy-prop3*deltaT;
    else
	{	// ccw around z axis (hence the -sin(ang))
        c=cos(angle);
        s=-sin(angle);
        c2=c*c;
        s2=s*s;
        return -(prop1*c2+prop2*s2)*sxx-(prop1*s2+prop2*c2)*syy
                -2.*(prop1-prop2)*c*s*sxy-prop3*deltaT;
    }
}
#endif


