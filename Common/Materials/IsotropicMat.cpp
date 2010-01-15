/********************************************************************************
    IsotropicMat.cpp
    NairnMPM
    
    Created by John Nairn on Thu Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/IsotropicMat.hpp"
#include "Read_XML/CommonReadHandler.hpp"

#pragma mark IsotropicMat::Constructors and Destructors

// Constructors
IsotropicMat::IsotropicMat() {}

// Constructors
IsotropicMat::IsotropicMat(char *matName) : Elastic(matName)
{
    int i;
    
    for(i=0;i<ISO_PROPS;i++)
        read[i]=0;
}

#pragma mark IsotropicMat::Initialization

// Read material properties
char *IsotropicMat::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"E")==0)
    {	read[E_PROP]=1;
        return((char *)&E);
    }
    
    else if(strcmp(xName,"G")==0)
    {	read[G_PROP]=1;
        return((char *)&G);
    }
    
    else if(strcmp(xName,"nu")==0)
    {	read[NU_PROP]=1;
        return((char *)&nu);
    }
    
    else if(strcmp(xName,"alpha")==0)
    {	read[A_PROP]=1;
        return((char *)&aI);
    }

#ifdef MPM_CODE
    else if(strcmp(xName,"beta")==0)
        return((char *)&betaI);
#endif

    return MaterialBase::InputMat(xName,input);
}

// print to output window
void IsotropicMat::PrintMechanicalProperties(void)
{	
	PrintProperty("E",E,"");
	PrintProperty("v",nu,"");
	PrintProperty("G",G,"");
	cout << endl;
	
	PrintProperty("a",aI,"");
    cout << endl;
}

// calculate properties used in analyses
const char *IsotropicMat::VerifyProperties(int np)
{
    // finish input and verify all there
    if(!read[G_PROP])
    {	G=E/(2.*(1.+nu));
        read[G_PROP]=1;
    }
    else if(!read[E_PROP])
    {	E=2.*G*(1.+nu);
        read[E_PROP]=1;
    }
    else if(!read[NU_PROP])
    {	nu=E/(2.*G)-1.;
        read[NU_PROP]=1;
    }
    else
		return "E, nu, and G all specified. Only two allowed";
		
    int i;
    for(i=0;i<ISO_PROPS;i++)
    {	if(!read[i])
			return "A required material property is missing";
    }
    
    // analysis properties
    const char *err=SetAnalysisProps(np,1.e6*E,1.e6*E,1.e6*E,nu,nu,nu,1.e6*G,1.e6*G,1.e6*G,
			1.e-6*aI,1.e-6*aI,1.e-6*aI,betaI*concSaturation,betaI*concSaturation,betaI*concSaturation);
	if(err!=NULL) return err;
	
	// superclass call
	return MaterialBase::VerifyProperties(np);
}

// fill in stiffness matrix
// makeSpecific divides by density, but only works in MPM code
// Here gets mdm[6][6], me0[4], mc0[4]
void IsotropicMat::InitialLoadMechProps(int makeSpecific,int np)
{
    hasMatProps=TRUE;
	
#ifdef MPM_CODE
	if(np==THREED_MPM)
	{	double rrho=1./rho;
		mdm[0][0]=C11*rrho;
		mdm[0][1]=C12*rrho;
		mdm[0][2]=C13*rrho;
		mdm[0][3]=mdm[0][4]=mdm[0][5]=0.;
		mdm[1][0]=mdm[0][1];
		mdm[1][1]=C22*rrho;
		mdm[1][2]=C23*rrho;
		mdm[1][3]=mdm[1][4]=mdm[1][5]=0.;
		mdm[2][0]=mdm[0][2];
		mdm[2][1]=mdm[1][2];
		mdm[2][2]=C33*rrho;
		mdm[2][3]=mdm[2][4]=mdm[2][5]=0.;
		mdm[3][0]=mdm[3][1]=mdm[3][2]=0.;
		mdm[3][3]=C44*rrho;
		mdm[3][4]=mdm[3][5]=0.;
		mdm[4][0]=mdm[4][1]=mdm[4][2]=mdm[4][3]=0.;
		mdm[4][4]=C55*rrho;
		mdm[4][5]=0.;
		mdm[5][0]=mdm[5][1]=mdm[5][2]=mdm[5][3]=mdm[5][4]=0.;
		mdm[5][5]=C66*rrho;
		
		// need me0[] and mc0[] for thermal and moisture expansion
		me0[0]=me0[1]=me0[2]=CTE1;
		me0[3]=me0[4]=me0[5]=0.;
		
		mc0[0]=mc0[1]=mc0[2]=CME1;
		mc0[3]=mc0[4]=mc0[5]=0.;
		
		return;
	}
#endif
    
    // Stiffness matrix
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
#endif

    // initial strains - all thermal and strain per temperature change
    me0[1]=CTE1;
    me0[2]=CTE2;
    me0[3]=0.;
    me0[4]=CTE3;

#ifdef MPM_CODE
	me0[5]=prop1;
	me0[6]=prop2;
	me0[7]=0.;
	
    // initial strains for moisture
    mc0[1]=CME1;
    mc0[2]=CME2;
    mc0[3]=0.;
    mc0[4]=CME3;
	
    // for MPM (units N/m^2 cm^3/g)
    if(makeSpecific)
    {	double rrho=1./rho;
    	mdm[1][1]*=rrho;
        mdm[1][2]*=rrho;
        mdm[2][2]*=rrho;
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
    mdm[3][1]=0.;
    mdm[3][2]=0.;

#ifdef FEA_CODE
    /* For Plane Strain analysis, save term for strain energy (4 = E3*a3^2*T^2))
          or energy +=0.5*mdm[4][4]*volume*(delta T)^2 */
    if(np==PLANE_STRAIN)
        mdm[4][4]=prop3*CTE3;

    /* For axisymmetric, put theta direction properties in 4th
			column and row */
    else if(np==AXI_SYM)
    {	mdm[4][1]=C13;
        mdm[4][2]=C23;
        mdm[4][3]=0.;
        mdm[4][4]=C33;
            
        mdm[1][4]=mdm[4][1];
        mdm[2][4]=mdm[4][2];
        mdm[3][4]=mdm[4][3];
    }
#endif
}

#ifdef MPM_CODE
// Call base class to set diffusion and moisture tensors, but set flag thay
// they are set first
void IsotropicMat::InitialLoadTransProps(void)
{
	hasTransProps=TRUE;
	MaterialBase::InitialLoadTransProps();
}
#endif

#pragma mark IsotropicMat::Methods

#ifdef MPM_CODE
/* d -- crack opening displacement near crack tip, d.y--opening, d.x--shear
     C -- crack propagating velocity
    J0 -- J integral components
*/
Vector IsotropicMat::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{ 
    double Cs2,Cd2,C2;
    double B1,B2,A1,A2,A3,A4,DC;
    double kf=0.,term1,term2;
    double dx,dy,J0x,J0y;
    Vector SIF;

    dx=d.x;
    dy=d.y;
    J0x=fabs(J0.x);                        // J0.x should be always positive
    J0y=J0.y;

    if(np==PLANE_STRAIN_MPM)
        kf=3.-4.*nu;
    else if(np==PLANE_STRESS_MPM)
        kf=(3.-nu)/(1.+nu);

    C2=C.x*C.x+C.y*C.y;                    // square of crack velocity
    if(!DbleEqual(sqrt(C2),0.0)) {         // dynamic crack
        Cs2=1.e+3*G/rho;                   // now in m^2/sec^2
        Cd2=Cs2*(kf+1.)/(kf-1.);
        B1=sqrt(1.-C2/Cd2);
        B2=sqrt(1.-C2/Cs2);
        DC=4*B1*B2-(1.+B2*B2)*(1.+B2*B2);
        A1=B1*(1.-B2*B2)/DC;
        A2=B2*(1.-B2*B2)/DC;
        A3=1./B2;
        term1=0.5*(4.*B1*B2+(1.+B2*B2)*(1.+B2*B2))*(2.+B1+B2)/sqrt((1.+B1)*(1.+B2));
        A4=(B1-B2)*(1.-B2*B2)*(term1-2.*(1.+B2*B2))/DC/DC;
    }
    else {                                // stationary crack
        B1=B2=1.0;
        A3=1.;
        A1=A2=A4=(kf+1.)/4.;
    }

    term2=dy*dy*B2+dx*dx*B1;
    if(DbleEqual(term2,0.0)) {            // COD=0
       SIF.x=0.0;
       SIF.y=0.0;
    }
    else {
       SIF.x=dy*sqrt(2*G*J0x*B2/A1/term2);
       SIF.y=dx*sqrt(2*G*J0x*B1/A2/term2);
    }
    
    return SIF;
}

#endif

#pragma mark IsotropicMat::Accessors

// Return the material tag
int IsotropicMat::MaterialTag(void) { return ISOTROPIC; }

// return material type
const char *IsotropicMat::MaterialType(void) { return "Isotropic"; }

#ifdef MPM_CODE
/*	calculate wave speed in mm/sec (because G in MPa and rho in g/cm^3)
	Uses sqrt((K +4G/3)/rho) which is dilational wave speed
	Identity also: K + 4G/3 = Lambda + 2G = 2G(1-nu)/(1-2 nu)
	Another form: E(1-nu)/((1+nu)(1-2*nu)rho)
*/
double IsotropicMat::WaveSpeed(void)
{
    return sqrt(2.e9*G*(1.-nu)/(rho*(1.-2.*nu)));
}

/*	calculate shear wave speed in mm/sec (because G in MPa and rho in g/cm^3)
*/
double IsotropicMat::ShearWaveSpeed(void) { return sqrt(1.e9*G/rho); }

#endif



