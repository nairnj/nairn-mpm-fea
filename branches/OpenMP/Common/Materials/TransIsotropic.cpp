/********************************************************************************
    TransIsotropic.cpp
    NairnMPM
    
    Created by John Nairn on Tue Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/TransIsotropic.hpp"
#ifdef MPM_CODE
	#include "Custom_Tasks/ConductionTask.hpp"
	#include "Custom_Tasks/DiffusionTask.hpp"
	#include "MPM_Classes/MPMBase.hpp"
	#include "Exceptions/CommonException.hpp"
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
    
    nuA=0.33;
    aA=40;
    aT=40;

#ifdef MPM_CODE
	diffA=0.;
	diffT=0.;
	kcondA=0.;
	kcondT=0.;
#else
    hasMatProps=FALSE;
#endif
	betaA=0.;
	betaT=0.;
}

#pragma mark TransIsotropic::Initialization

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

#ifdef MPM_CODE
// print transport properties to output window
void TransIsotropic::PrintTransportProperties(void) const
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
	{	PrintProperty("ka",rho*kcondA/1000.,"W/(m-K)");
		PrintProperty("kt",rho*kcondT/1000.,"W/(m-K)");
		PrintProperty("C",heatCapacity,"J/(kg-K)");
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
        return((char *)&kcondA);
    
    else if(strcmp(xName,"kCondT")==0)
        return((char *)&kcondT);
#endif
    
    return MaterialBase::InputMat(xName,input);
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

#ifdef MPM_CODE
    // make conductivity specific (N mm^3/(sec-K-g))
    kcondA *= (1000./rho);
    kcondT *= (1000./rho);
#endif
	
	// load elastic properties with constant values when istropic in x-y plane
	if(MaterialTag()==TRANSISO1 && np!=THREED_MPM)
		FillElasticProperties2D(&pr,np>BEGIN_MPM_TYPES,0.,np);
	
	// superclass call
	return MaterialBase::VerifyAndLoadProperties(np);
}

#pragma mark TransIsotropic::Methods

#ifdef MPM_CODE

// buffer size for mechanical properties
int TransIsotropic::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = 0;
    return sizeof(ElasticProperties);
}

// Isotropic material can use read-only initial properties
void *TransIsotropic::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	// if isotropic in 2D plane, use initial properties
	if(MaterialTag()==TRANSISO1 && np!=THREED_MPM)
		return (void *)&pr;
	
	// create new elastic properties
	ElasticProperties *p = (ElasticProperties *)matBuffer;
	if(np!=THREED_MPM)
		FillElasticProperties2D(p,TRUE,mptr->GetRotationZ(),np);
	else
		FillElasticProperties3D(mptr,p,np);
	return p;
}

// Fill ElasticProperties variable with current particle state
void TransIsotropic::FillElasticProperties3D(MPMBase *mptr,ElasticProperties *p,int np) const
{
	/* Rotation of the stiffness matrix requires Rz(-z).Ry(-y).Rx(-z).C.Rx^T(-x).Ry^T(-y).Rz^T(-z)
		Doing matrix math here would be 7*6*6 = 252 multiplications.
	 
		To improve performance, the transformation was expanded in Mathematica and each term
		of the matrix convert to an expression (using CForm) to paste. The 252 multiplications
		are reduced to 21 complex expression. Also trigonometric terms are evaluated once first.
	*/
	
	double z=mptr->GetRotationZ();
	double cz=cos(z);
	double sz=sin(z);
	double cz2=cz*cz;
	double sz2=sz*sz;
	double cz3=cz2*cz;
	double sz3=sz2*sz;
	double cz4=cz2*cz2;
	double sz4=sz2*sz2;
	double c2z=cos(2.*z);
	double s2z=sin(2.*z);
	double c2z2=c2z*c2z;
	double s2z2=s2z*s2z;
	
	double y=mptr->GetRotationY();
	double cy=cos(y);
	double sy=sin(y);
	double cy2=cy*cy;
	double sy2=sy*sy;
	double cy3=cy2*cy;
	double sy3=sy2*sy;
	double sy4=sy2*sy2;
	double cy4=cy2*cy2;
	double c2y=cos(2.*y);
	double s2y=sin(2.*y);
	double s4y=sin(4.*y);
	double c2y2=c2y*c2y;
	double s2y2=s2y*s2y;
	
	double x=mptr->GetRotationX();
	double cx=cos(x);
	double sx=sin(x);
	double cx2=cx*cx;
	double sx2=sx*sx;
	double cx4=cx2*cx2;
	double sx4=sx2*sx2;
	double c2x=cos(2.*x);
	double s2x=sin(2.*x);
	double s4x=sin(4.*x);
	double c2x2=c2x*c2x;
	double s2x2=s2x*s2x;
	
	p->C[0][0] = (2*cz2*s2x*s2z*sy*(2*(C12 - C13 - 2*C55 + 2*C66)*cy2 - (C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*sy2) +
				 s2z2*(4*C66*cx2*cy2 + 4*C55*cy2*sx2 + (4*c2x2*C44 + (C22 - 2*C23 + C33)*s2x2)*sy2) + 
				 4*cz4*(C11*cy4 + C55*cx2*s2y2 + C66*s2y2*sx2 + (C12 + C13 - C12*c2x + C13*c2x)*cy2*sy2 + (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy4) + 
				 2*((-2*C23*c2x + C22*(1 + c2x) + (-1 + c2x)*C33 - 4*c2x*C44)*s2x*s2z*sy + 
					4*cz2*(C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2))*sz2 + 4*(C22*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C33*sx4)*sz4)/4.;
	p->C[0][1] = cz2*(((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2z*sy)/4. + cz2*(C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2) +
					 (C22*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C33*sx4)*sz2) - s2z*(-(cz2*s2x*sy*
					(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2))/4. + 
					s2z*(cy2*(C66*cx2 + C55*sx2) + ((4*c2x2*C44 + (C22 - 2*C23 + C33)*s2x2)*sy2)/4.) + ((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*sy*sz2)/4.) + 
					sz2*(-(s2x*s2z*sy*(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2))/4. + 
						cz2*(C11*cy4 + C55*cx2*s2y2 + C66*s2y2*sx2 + (C12 + C13 - C12*c2x + C13*c2x)*cy2*sy2 + (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy4) + 
						 (C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2)*sz2);
	p->C[0][2] = s2z*(-((C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*cy2*s2x*sy)/4. + (C12 - C13)*cx*sx*sy3) +
				cz2*(-(C55*cx2*s2y2) - C66*s2y2*sx2 + cy2*(C11 + C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(cy4 + sy4) + C12*sx2*(cy4 + sy4)) + 
				(-(C44*cy2*s2x2) + (C22 + C33)*cx2*cy2*sx2 + C23*cy2*(cx4 + sx4) + (C12*cx2 + C13*sx2)*sy2)*sz2;
	p->C[0][3] = (cy*cz*(s2z*(-4*c2x2*C44 + 4*C66*cx2 - (C22 - 2*C23 + C33)*s2x2 + 4*C55*sx2)*sy +
				cz2*s2x*(-2*C12*cy2 + 2*C13*cy2 + (C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*sy2) - 
				(C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*sz2) - 2*sz*
				 (2*s2z*(c2y*(C55 - C66)*cx*cy*sx + ((-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2y*sy)/8. + (C12 - C13)*cx*cy*sx*sy2) + 
				  cz2*s2y*(C11*cy2 - 2*c2y*(C55*cx2 + C66*sx2) - (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(-cy2 + sy2) + C12*sx2*(-cy2 + sy2)) + 
				  s2y*(C12*cx2 + C44*s2x2 + (C13 - (C22 + C33)*cx2)*sx2 - C23*(cx4 + sx4))*sz2))/4.;
	p->C[0][4] = (cy*sz*(s2z*(-4*c2x2*C44 + 4*C66*cx2 - (C22 - 2*C23 + C33)*s2x2 + 4*C55*sx2)*sy +
				cz2*s2x*(-2*C12*cy2 + 2*C13*cy2 + (C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*sy2) - 
				(C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*sz2) + 2*cz*
				 (2*s2z*(c2y*(C55 - C66)*cx*cy*sx + ((-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2y*sy)/8. + (C12 - C13)*cx*cy*sx*sy2) + 
				  cz2*s2y*(C11*cy2 - 2*c2y*(C55*cx2 + C66*sx2) - (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(-cy2 + sy2) + C12*sx2*(-cy2 + sy2)) + 
				  s2y*(C12*cx2 + C44*s2x2 + (C13 - (C22 + C33)*cx2)*sx2 - C23*(cx4 + sx4))*sz2))/4.;
	p->C[0][5] = (s2z*(((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2z*sy)/4. + cz2*(C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2) +
					(C22*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C33*sx4)*sz2))/2. + c2z*
					(-(cz2*s2x*sy*(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2))/4. + 
					 s2z*(cy2*(C66*cx2 + C55*sx2) + ((4*c2x2*C44 + (C22 - 2*C23 + C33)*s2x2)*sy2)/4.) + ((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*sy*sz2)/4.) - 
					(s2z*(-(s2x*s2z*sy*(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2))/4. + 
					cz2*(C11*cy4 + C55*cx2*s2y2 + C66*s2y2*sx2 + (C12 + C13 - C12*c2x + C13*c2x)*cy2*sy2 + (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy4) + 
					(C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2)*sz2))/2.;
	
	p->C[1][1] = cx4*(C22*cz4 + 2*C23*cz2*sy2*sz2 + C33*sy4*sz4) + (4*cz4*(C44*s2x2 + C33*sx4) + 4*c2x2*C44*s2z2*sy2 + C22*s2x2*s2z2*sy2 - 2*C23*s2x2*s2z2*sy2 + C33*s2x2*s2z2*sy2 -
					4*(C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*cz3*s2x*sy*sz + 8*cz2*(C13*cy2*sx2 - C44*s2x2*sy2 + C23*sx4*sy2)*sz2 + 
					4*cz*s2x*sy*(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2)*sz3 + 4*C11*cy4*sz4 + 4*C66*s2y2*sx2*sz4 + 
					4*C44*s2x2*sy4*sz4 + 4*C22*sx4*sy4*sz4 + 4*cy2*(C55*s2z2*sx2 + (C12 + C13 - C12*c2x + C13*c2x)*sy2*sz4))/4. + 
					cx2*(C66*cy2*s2z2 + 2*cz2*(C12*cy2 + (C22 + C33)*sx2*sy2)*sz2 + C55*s2y2*sz4 + 2*C23*sx2*(cz4 + sy4*sz4));
	p->C[1][2] = cz2*(-(C44*cy2*s2x2) + (C22 + C33)*cx2*cy2*sx2 + C23*cy2*(cx4 + sx4) + (C12*cx2 + C13*sx2)*sy2) +
					((C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*cy2*s2x*s2z*sy + 4*(-C12 + C13)*cx*s2z*sx*sy3)/4. + 
					(-(C55*cx2*s2y2) - C66*s2y2*sx2 + cy2*(C11 + C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(cy4 + sy4) + C12*sx2*(cy4 + sy4))*sz2;
	p->C[1][3] = (cy*cz*(-((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*cz2*s2x) + s2z*(4*c2x2*C44 - 4*C66*cx2 + (C22 - 2*C23 + C33)*s2x2 - 4*C55*sx2)*sy +
				s2x*(-2*C12*cy2 + 2*C13*cy2 + (C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*sy2)*sz2) - 
				 2*sz*(cz2*s2y*(C12*cx2 + C44*s2x2 + (C13 - (C22 + C33)*cx2)*sx2 - C23*(cx4 + sx4)) - 
				2*s2z*(c2y*(C55 - C66)*cx*cy*sx + ((-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2y*sy)/8. + (C12 - C13)*cx*cy*sx*sy2) + 
				s2y*(C11*cy2 - 2*c2y*(C55*cx2 + C66*sx2) - (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(-cy2 + sy2) + C12*sx2*(-cy2 + sy2))*sz2))/4.;
	p->C[1][4] = (cy*sz*(-((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*cz2*s2x) + s2z*(4*c2x2*C44 - 4*C66*cx2 + (C22 - 2*C23 + C33)*s2x2 - 4*C55*sx2)*sy +
				s2x*(-2*C12*cy2 + 2*C13*cy2 + (C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*sy2)*sz2) + 
				 2*cz*(cz2*s2y*(C12*cx2 + C44*s2x2 + (C13 - (C22 + C33)*cx2)*sx2 - C23*(cx4 + sx4)) - 
				2*s2z*(c2y*(C55 - C66)*cx*cy*sx + ((-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2y*sy)/8. + (C12 - C13)*cx*cy*sx*sy2) + 
				s2y*(C11*cy2 - 2*c2y*(C55*cx2 + C66*sx2) - (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(-cy2 + sy2) + C12*sx2*(-cy2 + sy2))*sz2))/4.;
	p->C[1][5] = (c2z*((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*cz2*s2x*sy - s2z*(4*C66*cx2*cy2 + 4*C55*cy2*sx2 + (4*c2x2*C44 + (C22 - 2*C23 + C33)*s2x2)*sy2) -
				s2x*sy*(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2)*sz2) + 
				 2*s2z*(cz2*(C22*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C33*sx4) - ((C22 - C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*s2x*s2z*sy)/4. + 
				(C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2)*sz2) - 
				 2*s2z*((s2x*s2z*sy*(2*(-C12 + C13 + 2*C55 - 2*C66)*cy2 + (-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*sy2))/4. + 
				cz2*(C12*cx2*cy2 + C13*cy2*sx2 + (-(C44*s2x2) + (C22 + C33)*cx2*sx2 + C23*(cx4 + sx4))*sy2) + 
				(C11*cy4 + C55*cx2*s2y2 + C66*s2y2*sx2 + (C12 + C13 - C12*c2x + C13*c2x)*cy2*sy2 + (C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy4)*sz2))/4.;
	
	p->C[2][2] = C33*cx4*cy4 + C44*cy4*s2x2 + C55*cx2*s2y2 + 2*C23*cx2*cy4*sx2 + C66*s2y2*sx2 + C22*cy4*sx4 + 2*cy2*(C13*cx2 + C12*sx2)*sy2 + C11*sy4;
	p->C[2][3] = (cy*cz*s2x*((C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*cy2 + 2*(-C12 + C13 + 2*C55 - 2*C66)*sy2) +
				 2*s2y*(-2*c2y*(C55*cx2 + C66*sx2) + cy2*(C33*cx4 + C44*s2x2 - C12*sx2 + 2*C23*cx2*sx2 + C22*sx4) - C11*sy2 + C12*sx2*sy2 + C13*cx2*(-cy2 + sy2))*sz)/4.;
	p->C[2][4] = (2*cz*s2y*(2*c2y*(C55*cx2 + C66*sx2) - cy2*(C33*cx4 + C44*s2x2 - C12*sx2 + 2*C23*cx2*sx2 + C22*sx4) + C13*cx2*(cy2 - sy2) + C11*sy2 - C12*sx2*sy2) +
				 cy*s2x*((-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44))*cy2 + 2*(-C12 + C13 + 2*C55 - 2*C66)*sy2)*sz)/4.;
	p->C[2][5] = (c2z*s2x*sy*(-((-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*cy2) + 2*(C12 - C13)*sy2) +
				 2*s2z*(-(C44*cy2*s2x2) + (C22 + C33)*cx2*cy2*sx2 + C23*cy2*(cx4 + sx4) + (C12*cx2 + C13*sx2)*sy2) - 
				 2*s2z*(-(C55*cx2*s2y2) - C66*s2y2*sx2 + cy2*(C11 + C33*cx4 + C44*s2x2 + 2*C23*cx2*sx2 + C22*sx4)*sy2 + C13*cx2*(cy4 + sy4) + C12*sx2*(cy4 + sy4)))/4.;
	
	p->C[3][3] = (2*cz2*(4*c2x2*C44*cy2 + (C22 - 2*C23 + C33)*cy2*s2x2 + 4*(C66*cx2 + C55*sx2)*sy2) +
				 cz*(2*C12*(1 + c2y + 2*cy2)*s2x - 2*C13*(1 + c2y + 2*cy2)*s2x - C22*(1 + c2y - 2*(-1 + c2x)*cy2)*s2x + C22*cy2*s4x + 
				C33*((1 + c2y + 2*(1 + c2x)*cy2)*s2x + cy2*s4x) - 2*((C23 + 2*C44)*cy2*(2*c2x*s2x + s4x) + 2*c2y*(C55 - C66)*(s2x + 2*cx*sx)))*sy*sz + 
				 2*(4*c2y2*(C55*cx2 + C66*sx2) + s2y2*(C11 - 2*C13*cx2 + C33*cx4 + C44*s2x2 - 2*C12*sx2 + 2*C23*cx2*sx2 + C22*sx4))*sz2)/8.;
	p->C[3][4] = -(cz*(cz*((2*C12*(1 + c2y) - 2*C13*(1 + c2y) - C22*(1 + c2y) + C33 + c2y*(C33 - 4*C55 + 4*C66))*s2x + (C22 - 2*C23 + C33 - 4*C44)*cy2*s4x)*sy +
				2*(4*c2y2*(C55*cx2 + C66*sx2) + s2y2*(C11 - 2*C13*cx2 + C33*cx4 + C44*s2x2 - 2*C12*sx2 + 2*C23*cx2*sx2 + C22*sx4))*sz))/8. + 
				(sz*(4*c2x2*C44*cy2*cz + (C22 - 2*C23 + C33)*cy2*cz*s2x2 + 4*cz*(C66*cx2 + C55*sx2)*sy2 + 
				((2*C12 - 2*C13 + C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*cy2*s2x + 4*c2y*(-C55 + C66)*cx*sx)*sy*sz))/4.;
	p->C[3][5] = (-((-2*C23*c2x + C22*(1 + c2x) + (-1 + c2x)*C33 - 4*c2x*C44)*cy*cz*s2x*s2z) + 2*s2y*s2z*(-(C12*cx2) - C44*s2x2 + (-C13 + (C22 + C33)*cx2)*sx2 + C23*(cx4 + sx4))*sz +
				 2*c2z*cy*(-(cz*(4*c2x2*C44 - 4*C66*cx2 + (C22 - 2*C23 + C33)*s2x2 - 4*C55*sx2)*sy) - 
				(4*c2y*(C55 - C66)*cx*sx + (2*C12 - 2*C13 + C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*s2x*sy2)*sz) + 
				 s2z*(-((C55 + c2x*C55 + C66 - c2x*C66)*s4y*sz) + 2*cy3*((C12 - C13)*cz*s2x + (2*C11 + C12*(-1 + c2x) - C13*(1 + c2x))*sy*sz) - 
				(cy*sy2*(2*(-C22 + C33 + c2x*(C22 - 2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*cz*s2x + 
				8*(-(C13*cx2) + C33*cx4 + C44*s2x2 - C12*sx2 + 2*C23*cx2*sx2 + C22*sx4)*sy*sz))/2.))/8.;
	
	p->C[4][4] = (2*cz2*(4*c2y2*(C55*cx2 + C66*sx2) + s2y2*(C11 - 2*C13*cx2 + C33*cx4 + C44*s2x2 - 2*C12*sx2 + 2*C23*cx2*sx2 + C22*sx4)) -
				 cz*((C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*cy*s2x*s2y + 
				2*((2*C12 - 2*C13 + C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*cy2*s2x + 4*cx*(2*c2y*(-C55 + C66) + (C12 - C13)*cy2)*sx)*sy)*sz + 
				 2*(4*c2x2*C44*cy2 + (C22 - 2*C23 + C33)*cy2*s2x2 + 4*(C66*cx2 + C55*sx2)*sy2)*sz2)/8.;
	p->C[4][5] = (2*cz*s2y*s2z*(C12*cx2 + C44*s2x2 + (C13 - (C22 + C33)*cx2)*sx2 - C23*(cx4 + sx4)) +
				 2*cz*s2z*(C55*cx2*s4y + C66*s4y*sx2 + (-2*C11 + C12 + C13 - C12*c2x + C13*c2x)*cy3*sy + 
				2*cy*(-(C13*cx2) + C33*cx4 + C44*s2x2 - C12*sx2 + 2*C23*cx2*sx2 + C22*sx4)*sy3) - (-2*C23*c2x + C22*(1 + c2x) + (-1 + c2x)*C33 - 4*c2x*C44)*cy*s2x*s2z*sz + 
				 cy*s2x*s2z*(2*C12*cy2 - 2*C13*cy2 - (C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44) - 4*C55 + 4*C66)*sy2)*sz + 
				 2*c2z*cy*(4*c2y*(C55 - C66)*cx*cz*sx + (2*C12 - 2*C13 + C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*cz*s2x*sy2 - 
				(4*c2x2*C44 - 4*C66*cx2 + (C22 - 2*C23 + C33)*s2x2 - 4*C55*sx2)*sy*sz))/8.;
	
	p->C[5][5] = (2*c2z2*(4*C66*cx2*cy2 + 4*C55*cy2*sx2 + (4*c2x2*C44 + (C22 - 2*C23 + C33)*s2x2)*sy2) +
				 c2z*s2x*s2z*(sy*(-2*(C33 + 2*(C12 - C13 - 2*C55 + 2*C66)*cy2) + C33*sy2 - 2*C23*c2x*(2 + sy2) + c2x*(C33 - 4*C44)*(2 + sy2) + C22*(2 - sy2 + c2x*(2 + sy2))) + 
				(C22*(-1 + c2x) + C33 + c2x*(-2*C23 + C33 - 4*C44))*sy3) + 2*s2z2*
				 (-2*C12*cx2*cy2 + C11*cy4 + C44*s2x2 + C55*cx2*s2y2 + 2*C23*cx2*sx2 - 2*C13*cy2*sx2 + C66*s2y2*sx2 + C33*sx4 + 
				  sy2*((C12 + C13 - C12*c2x + C13*c2x)*cy2 - 2*C33*cx2*sx2 - C23*(2*cx4 + (-1 + c2y)*cx2*sx2 + 2*sx4) + C44*s2x2*(2 + sy2)) + C33*cx4*sy4 + 
				  C22*(cx4 - 2*cx2*sx2*sy2 + sx4*sy4)))/8.;

	// make specific and symmetric
	double rrho=1./rho;
	int i,j;
	for(i=0;i<6;i++)
	{	p->C[i][i]*=rrho;
		for(j=i+1;j<6;j++)
		{	p->C[i][j]*=rrho;
			p->C[j][i]=p->C[i][j];
		}
	}

	// need p->alpha[] and p->beta[] too for thermal and moisture expansion
	p->alpha[0] = CTE1*cy2*cz2 + CTE2*(cx*s2z*sx*sy + cz2*sx2*sy2 + cx2*sz2) + CTE3*(-(cx*s2z*sx*sy) + cx2*cz2*sy2 + sx2*sz2);
	p->alpha[1] = CTE1*cy2*sz2 + CTE3*(cz2*sx2 + cx*s2z*sx*sy + cx2*sy2*sz2) + CTE2*(cx2*cz2 - cx*s2z*sx*sy + sx2*sy2*sz2);
	p->alpha[2] = CTE3*cx2*cy2 + CTE2*cy2*sx2 + CTE1*sy2;
	p->alpha[3] = -(CTE1*s2y*sz)/2. + CTE3*(cx*cy*cz*sx + (cx2*s2y*sz)/2.) + CTE2*(-(cx*cy*cz*sx) + (s2y*sx2*sz)/2.);
	p->alpha[4] = (CTE1*cz*s2y)/2. + CTE2*(-(cz*s2y*sx2)/2. - cx*cy*sx*sz) + CTE3*(-(cx2*cz*s2y)/2. + cx*cy*sx*sz);
	p->alpha[5] = -(CTE1*cy2*s2z)/2. + CTE3*((s2z*sx2)/2. - c2z*cx*sx*sy - (cx2*s2z*sy2)/2.) + CTE2*((cx2*s2z)/2. + c2z*cx*sx*sy - (s2z*sx2*sy2)/2.);
	
	p->alpha[0] = CME1*cy2*cz2 + CME2*(cx*s2z*sx*sy + cz2*sx2*sy2 + cx2*sz2) + CME3*(-(cx*s2z*sx*sy) + cx2*cz2*sy2 + sx2*sz2);
	p->alpha[1] = CME1*cy2*sz2 + CME3*(cz2*sx2 + cx*s2z*sx*sy + cx2*sy2*sz2) + CME2*(cx2*cz2 - cx*s2z*sx*sy + sx2*sy2*sz2);
	p->alpha[2] = CME3*cx2*cy2 + CME2*cy2*sx2 + CME1*sy2;
	p->alpha[3] = -(CME1*s2y*sz)/2. + CME3*(cx*cy*cz*sx + (cx2*s2y*sz)/2.) + CME2*(-(cx*cy*cz*sx) + (s2y*sx2*sz)/2.);
	p->alpha[4] = (CME1*cz*s2y)/2. + CME2*(-(cz*s2y*sx2)/2. - cx*cy*sx*sz) + CME3*(-(cx2*cz*s2y)/2. + cx*cy*sx*sz);
	p->alpha[5] =-(CME1*cy2*s2z)/2. + CME3*((s2z*sx2)/2. - c2z*cx*sx*sy - (cx2*s2z*sy2)/2.) + CME2*((cx2*s2z)/2. + c2z*cx*sx*sy - (s2z*sx2*sy2)/2.);

}

#else

// Fill the object in the material class
void TransIsotropic::LoadMechanicalPropertiesFEA(int makeSpecific,double angle,int np)
{	// if filled with same angle then no need to fill again
	if(hasMatProps && DbleEqual(angle,lastMatAngle)) return;
	lastMatAngle=angle;
	hasMatProps=TRUE;
	FillElasticProperties2D(&pr,FALSE,angle,np);
}

#endif

// fill in stiffness matrix if necessary
// makeSpecific divides by density, but only used in MPM code
// Used by TranIsoptropic 2 and by Orthotropic
// angle is in radians
void TransIsotropic::FillElasticProperties2D(ElasticProperties *p,int makeSpecific,double angle,int np) const
{
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
    // for MPM (units N/m^2 cm^3/g)
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
// Note: no angle, because can not depend on material angle
// Here fills in isotropic properties, materials with different anisotropic properties should override
void TransIsotropic::FillTransportProperties(TransportProperties *t)
{
	if(MaterialTag()==TRANSISO1)
	{	t->diffusionTensor.xx = diffT;
		t->diffusionTensor.yy = diffT;
		t->kCondTensor.xx = kcondT;
		t->kCondTensor.yy = kcondT;
	}
	else
	{	t->diffusionTensor.xx=diffT;
		t->diffusionTensor.yy=diffA;
		t->kCondTensor.xx = kcondT;
		t->kCondTensor.yy = kcondA;
	}
	t->diffusionTensor.xx = GetDiffZ();
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
		if(MaterialTag()==TRANSISO1)
		{	*t = tr;
			return;
		}

		// analysis axes are ccw from material axes (note the -sin(ang))
		double angle=mptr->GetRotationZ();
		double cs=cos(angle);
		double c2=cs*cs;
		double sn=-sin(angle);
		double s2=sn*sn;
		double cssn=cs*sn;

		// diffusion and conductivity tensors
		t->diffusionTensor.xx=diffA*s2 + diffT*c2;
		t->diffusionTensor.yy=diffA*c2 + diffT*s2;
		t->diffusionTensor.xy=(diffT-diffA)*cssn;
		t->kCondTensor.xx = kcondA*s2 + kcondT*c2;
		t->kCondTensor.yy = kcondA*c2 + kcondT*s2;
		t->kCondTensor.xy = (kcondT-kcondA)*cssn;
	}
	
	else
	{	/* Rotation of the transport tensor requires Rz(-z).Ry(-y).Rx(-z).T.Rx^T(-x).Ry^T(-y).Rz^T(-z)
			To improve performance, the transformation was expanded in Mathematica and each term
			of the matrix convert to an expression (using CForm) to paste.
			Also trigonometric terms are evaluated once first.
		*/
		
		double z=mptr->GetRotationZ();
		double cz=cos(z);
		double sz=sin(z);
		double cz2=cz*cz;
		double sz2=sz*sz;
		
		double y=mptr->GetRotationY();
		double cy=cos(y);
		double sy=sin(y);
		double cy2=cy*cy;
		double sy2=sy*sy;
		
		double x=mptr->GetRotationX();
		double cx=cos(x);
		double sx=sin(x);
		double sx2=sx*sx;
		double cx2=cx*cx;
		
		double diffz = GetDiffZ();
		
		t->diffusionTensor.xx=cy2*cz2*diffT + diffA*(-(cz2*sx2*sy2) + cx2*sz2) + diffz*(cx2*cz2*sy2 - sx2*sz2);
		t->diffusionTensor.xy=-(cy2*cz*diffT*sz) + diffz*(cx*sx*sy - cz*sx2*sz - cx2*cz*sy2*sz) + diffA*(cx*sx*sy + cx2*cz*sz + cz*sx2*sy2*sz);
		t->diffusionTensor.xz=cy*cz*diffT*sy + cy*diffz*(-(cx2*cz*sy) + cx*sx*sz) + cy*diffA*(cz*sx2*sy + cx*sx*sz);
		
		t->diffusionTensor.yy=cy2*diffT*sz2 + diffz*(-(cz2*sx2) + cx2*sy2*sz2) + diffA*(cx2*cz2 - sx2*sy2*sz2);
		t->diffusionTensor.yz=-(cy*diffT*sy*sz) + cy*diffz*(cx*cz*sx + cx2*sy*sz) + cy*diffA*(cx*cz*sx - sx2*sy*sz);
		
		t->diffusionTensor.zz=cx2*cy2*diffz - cy2*diffA*sx2 + diffT*sy2;
		
		double kz=GetKcondZ();
		
		t->kCondTensor.xx=cy2*cz2*kcondT + kcondA*(-(cz2*sx2*sy2) + cx2*sz2) + kz*(cx2*cz2*sy2 - sx2*sz2);
		t->kCondTensor.xy=-(cy2*cz*kcondT*sz) + kz*(cx*sx*sy - cz*sx2*sz - cx2*cz*sy2*sz) + kcondA*(cx*sx*sy + cx2*cz*sz + cz*sx2*sy2*sz);
		t->kCondTensor.xz=cy*cz*kcondT*sy + cy*kz*(-(cx2*cz*sy) + cx*sx*sz) + cy*kcondA*(cz*sx2*sy + cx*sx*sz);
		
		t->kCondTensor.yy=cy2*kcondT*sz2 + kz*(-(cz2*sx2) + cx2*sy2*sz2) + kcondA*(cx2*cz2 - sx2*sy2*sz2);
		t->kCondTensor.yz=-(cy*kcondT*sy*sz) + cy*kz*(cx*cz*sx + cx2*sy*sz) + cy*kcondA*(cx*cz*sx - sx2*sy*sz);
		
		t->kCondTensor.zz=cx2*cy2*kz - cy2*kcondA*sx2 + kcondT*sy2;
	}
}
#endif


#pragma mark TransIsotropic::Accessors

// Return the material tag
int TransIsotropic::MaterialTag(void) const { return tiType; }

// return material type
const char *TransIsotropic::MaterialType(void) const
{	if(MaterialTag()==TRANSISO1)
        return "Tranversely isotropic (A axis normal to x-y plane)";
    else
        return "Tranversely isotropic (A axis in x-y plane)";
}

#ifdef MPM_CODE

/* Calculate maximum wave speed in mm/sec (moduli in MPa, rho in g/cm^3)
	TRANSISO1
		wave speeds are GT/rho and (KT+GT)/rho - return larger one
	TRANSISO2
		wave speeds are GT/rho, GA/rho, (KT+GT)/rho, and (EA + 4KT nuA^2)/rho
		return largest (assumes shear ones are not largest)
*/
double TransIsotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{
    if(MaterialTag()==TRANSISO1 && !threeD)
        return sqrt(1.e9*(KT+GT)/rho);
    else
        return sqrt(1.e9*fmax(GA,fmax(KT+GT,EA+4.*KT*nuA*nuA))/rho);
}

// maximum diffusion coefficient in cm^2/sec (diff in mm^2/sec)
double TransIsotropic::MaximumDiffusion(void) const { return max(diffA,diffT)/100.; }

// maximum diffusivity in cm^2/sec
// specific k is mJ mm^2/(sec-K-g) and Cp is mJ/(g-K) so k/Cp = mm^2 /sec * 1e-2 = cm^2/sec
double TransIsotropic::MaximumDiffusivity(void) const { return 0.01*max(kcondA,kcondT)/heatCapacity; }

// diffusion and conductivity in the z direction
double TransIsotropic::GetDiffZ(void) const { return MaterialTag()==TRANSISO1 ? diffA : diffT; }
double TransIsotropic::GetKcondZ(void) const { return MaterialTag()==TRANSISO1 ? kcondA : kcondT; }

#endif


