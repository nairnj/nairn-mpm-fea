/********************************************************************************
    MaterialBaseMPM.cpp - more MaterialBase for MPM code
    NairnMPM
    
    Created by John Nairn on Mar 17 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include <vector>

// class statics for MPM - zero based material IDs when in multimaterial mode
vector<int> MaterialBase::fieldMatIDs;

#pragma mark MaterialBase::Initialization

// Read material properties common to all MPM materials
char *MaterialBase::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"rho")==0)
    {	input=DOUBLE_NUM;
        return((char *)&rho);
    }
    
    else if(strcmp(xName,"KIc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&KIc);
    }
    
	// crit 3 only
    else if(strcmp(xName,"p")==0)
    {	input=DOUBLE_NUM;
        return((char *)&pCrit3);
    }
    
	// crit 2 only
    else if(strcmp(xName,"maxLength")==0)
    {	input=DOUBLE_NUM;
        return((char *)&maxLength);
    }
    
    else if(strcmp(xName,"KIIc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&KIIc);
    }

	// gain in crit 3
    else if(strcmp(xName,"gain")==0)
    {	input=DOUBLE_NUM;
        return((char *)&gain);
    }

	// initTime in crit 2
    else if(strcmp(xName,"initTime")==0)
    {	input=DOUBLE_NUM;
        return((char *)&initTime);
    }

    else if(strcmp(xName,"KIexp")==0)
    {   input=DOUBLE_NUM;
        return((char *)&KIexp);
    }

    else if(strcmp(xName,"KIIexp")==0)
    {   input=DOUBLE_NUM;
        return((char *)&KIIexp);
    }
    
    else if(strcmp(xName,"JIc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&JIc);
    }
    
    else if(strcmp(xName,"JIIc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&JIIc);
    }
	
    else if(strcmp(xName,"nmix")==0)
    {	input=DOUBLE_NUM;
        return((char *)&nmix);
    }
	
	// speed in crit 2 and 3
    else if(strcmp(xName,"speed")==0)
    {	input=DOUBLE_NUM;
        return((char *)&initSpeed);
    }
	
	// growth direction in crit 2
    else if(strcmp(xName,"xGrow")==0)
    {	input=DOUBLE_NUM;
		constantDirection=TRUE;
        return((char *)&growDir.x);
    }
	
	// growth direction in crit 2
    else if(strcmp(xName,"yGrow")==0)
    {	input=DOUBLE_NUM;
		constantDirection=TRUE;
        return((char *)&growDir.y);
    }
	
	else if(strcmp(xName,"gamma")==0)
    {	input=DOUBLE_NUM;
        return((char *)&gamma);
    }
    
	// criterion 6 only
    else if(strcmp(xName,"delIc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&delIc);
    }
    
    else if(strcmp(xName,"delIIc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&delIIc);
    }
    
	else if(strcmp(xName,"Cp")==0)
    {	input=DOUBLE_NUM;
        return((char *)&heatCapacity);
    }

	else if(strcmp(xName,"Cv")==0)
    {	input=DOUBLE_NUM;
        return((char *)&heatCapacityVol);
    }

	else if(strcmp(xName,"csat")==0)
    {	input=DOUBLE_NUM;
        return((char *)&concSaturation);
    }

	else if(strcmp(xName,"D")==0)
    {	input=DOUBLE_NUM;
        return((char *)&diffusionCon);
	}
    
	else if(strcmp(xName,"kCond")==0)
    {	input=DOUBLE_NUM;
        return((char *)&kCond);
	}
    
    return((char *)NULL);
}

// print any properties common to all MPM material types
void MaterialBase::PrintCommonProperties(void)
{
	if(Rigid() || isTractionLaw()) return;
	
	// print density
	PrintProperty("rho",rho,"");
	cout << endl;
	
	// print growth criterion and relevant material properties for crack growth
	cout << "Crack Growth Criterion: ";
	PrintCriterion(criterion[0],matPropagateDirection[0]);
	
	// traction mat
	if(criterion[0]!=NO_PROPAGATION && tractionMat[0]>0)
		cout << "   New crack surface has traction material " << tractionMat[0] << endl;
	
	if(criterion[0]!=NO_PROPAGATION && criterion[1]!=NO_PROPAGATION)
	{	cout << "Alternate Crack Growth Criterion: ";
		PrintCriterion(criterion[1],matPropagateDirection[1]);
		
		// traction mat
		if(tractionMat[1]>0)
			cout << "   New crack surface has traction material " << tractionMat[1] << endl;
	}
	
	// convert MPa sqrt(m) to MPa sqrt(mm)
	KIc*=31.62277660168379;
	KIIc*=31.62277660168379;
	
	// convert other crack growth properties
	JIc*=1.e-3;			// convert J/m^2 to N/mm
	initTime*=1e-3;		// convert to sec
	if(criterion[0]==TOTALENERGYBALANCE)
		initSpeed*=0.01;	// convert % of WaveSpeed() to fraction of WaveSpeed()
	else
		initSpeed*=1.e3;	// convert m/sec to mm/sec
	gamma*=1.e-3;		// convert J/m^2 to N/mm
	// pCrit3 - dimensionless
	// gain - no units change
	
	// optional color
	if(red>=0.)
	{	char mline[200];
		sprintf(mline,"color= %g, %g, %g",red,green,blue);
		cout << mline << endl;
	}
}

// print fraction criterion
void MaterialBase::PrintCriterion(int thisCriterion,int thisDirection)
{
	char mline[200];
	
	switch(thisCriterion)
	{   case NO_PROPAGATION:
			cout << "No propagation" << endl;
			break;
			
		case MAXHOOPSTRESS:
			cout << "Maximum hoop stess" << PreferredDirection(thisDirection) << endl;
			PrintProperty("KIc",KIc,"MPa-sqrt(m)");
			PrintProperty("KIIc",KIIc,"MPa-sqrt(m)");
			cout << endl;
			break;
			
		case CRITICALERR:
			cout << "Critical Energy Release Rate" << PreferredDirection(thisDirection) << endl;
			PrintProperty("Jc",JIc,"J/m^2");
			cout << endl;
			break;
			
		case STEADYSTATEGROWTH:
			cout << "Constant crack speed" << PreferredDirection(thisDirection) << endl;
			if(initTime<0.)
			{	PrintProperty("Jc",JIc,"J/m^2");
				PrintProperty("initSpeed",initSpeed,"m/s");
			}
			else
			{	PrintProperty("ti",initTime,"ms");
				PrintProperty("initSpeed",initSpeed,"m/s");
			}
			cout << endl;
			if(maxLength>0.)
			{	PrintProperty("max length",maxLength,"mm");
				cout << endl;
			}
			if(constantDirection && thisDirection==DEFAULT_DIRECTION)
			{	double norm=sqrt(growDir.x*growDir.x + growDir.y*growDir.y);
				growDir.x/=norm;
				growDir.y/=norm;
				sprintf(mline,"direction = (%9.5f,%9.5f)",growDir.x,growDir.y);
				cout << mline << endl;
			}
			break;
			
		case TOTALENERGYBALANCE:
			// if gamma not provided or too large, set to JIc/2
			if(gamma<0. || JIc<2.*gamma) gamma=JIc/2.;
			
			cout << "Total energy balance" << PreferredDirection(thisDirection) << endl;
			sprintf(mline,"Jc =%12.3f J/m^2  vel=%12.3f%c wave speed",JIc,initSpeed,'%');
			cout << mline << endl;
			sprintf(mline,"gam=%12.3f J/m^2  p  =%12.3f        gain=%12.3g",gamma,pCrit3,gain);
			cout << mline << endl;
			break;
			
		case STRAINENERGYDENSITY:
			cout << "Minmum strain energy density" << PreferredDirection(thisDirection) << endl;  
			sprintf(mline,"KIc=%12.3f MPa-sqrt(m)  KIIc=%12.3f MPa-sqrt(m)",KIc,KIIc);
			cout << mline << endl;
			break;
			
		case EMPIRICALCRITERION:
			// If KIIc not set, set now to a default value
			if(KIIc<0.) KIIc=0.817*KIc;
			cout << "Empirical criterion" << PreferredDirection(thisDirection) << endl;  
			sprintf(mline,"KIc=%12.3f MPa-sqrt(m)  KIIc=%12.3f MPa-sqrt(m) KIexp=%12.3f KIIexp=%12.3f",
					KIc,KIIc,KIexp,KIIexp);
			cout << mline << endl;
			break;
			
		case MAXCTODCRITERION:
			cout << "Maximum CTOD" << PreferredDirection(thisDirection) << endl;
			sprintf(mline,"delIc=%12.6f mm  delIIc=%12.6f mm",delIc,delIIc);
			cout << mline << endl;
			break;
			
		default:
			cout << "Unknown" << endl;
	}
}	

// print transport properties to output window - default is isotropic properties
// aniostropic materials must override it
void MaterialBase::PrintTransportProperties(void)
{
	if(Rigid()) return;
	
	// Diffusion constants
	if(DiffusionTask::active)
	{	PrintProperty("D",diffusionCon,"mm^2/sec");
		PrintProperty("csat",concSaturation,"");
		cout << endl;
		PrintProperty("b",betaI,"1/wt fr");
		cout << endl;
	}
	// Conductivity constants
	if(ConductionTask::active)
	{	PrintProperty("k",kCond,"W/(m-K)");
		PrintProperty("Cp",1000.*heatCapacity,"J/(kg-K)");
		cout << endl;
	}
}

// return string to define crack propagation direction (a static class method)
const char *MaterialBase::PreferredDirection(int style)
{
	switch(style)
	{	case SELF_SIMILAR:
			return " in self-similar direction";
		case NORMAL_TO_COD:
			return " in direction normal to COD";
		case HOOP_FROM_COD:
			return " in hoop direction estimated from COD";
		case INITIAL_DIRECTION:
			return " in initial crack direction";
		case DEFAULT_DIRECTION:
		default:
			break;
	}
	return " in default criterion direction";
}

/* calculate properties used in analyses
	If superclass overrides this method, must call this too
*/
const char *MaterialBase::VerifyProperties(int np)
{
	// check which were set
	if(heatCapacityVol<0.) heatCapacityVol=fmax(heatCapacity,0.);
	if(heatCapacity<0.) heatCapacity=heatCapacityVol;
	
	// Convert to J/(g-K) such that (J/(sec-m-K)) / (rho (g/cm^3) Cp) has units mm^2/sec
	heatCapacity/=1000.;
	heatCapacityVol/=1000.;

	// in case only need to load some things once, load those mechanical properties now
	InitialLoadMechProps((int)(np>BEGIN_MPM_TYPES),np);
	InitialLoadTransProps();
	
	// check for unsupported alternate propagation criterion
	if(criterion[1]==TOTALENERGYBALANCE)
		return "The alternate propagation criterion cannot be energy balance method.";
	if(criterion[1]==STEADYSTATEGROWTH)
		return "The alternate propagation criterion cannot be steady state crack growth.";

	return NULL;
}

// Called before analysis, material can fill in things that never change during the analysis
// Note: no angle, because can not depend on material angle
// Here fills in isotropic properties, materials with different anisotropic properties should override
void MaterialBase::InitialLoadTransProps(void)
{
	// diffusion tensor (xx, yy, xy order)
	diffusionTensor.xx=diffusionCon;
	diffusionTensor.yy=diffusionCon;
	diffusionTensor.zz=diffusionCon;
	diffusionTensor.xy=0.;
	diffusionTensor.xz=0.;
	diffusionTensor.yz=0.;
	
	// conductivity tensor (xx, yy, xy order)
	kCondTensor.xx=kCond;
	kCondTensor.yy=kCond;
	kCondTensor.zz=kCond;
	kCondTensor.xy=0.;
	kCondTensor.xz=0.;
	kCondTensor.yz=0.;
}

// create and return pointer to material-specific data on a particle
//	using this material. Called once at start of analysis for each
//	particle
char *MaterialBase::MaterialData(void) { return NULL; }

// preliminary calculations (throw CommonException on problem)
void MaterialBase::PreliminaryMatCalcs(void)
{	int i;
	
	for(i=0;i<=1;i++)
	{	if(i==1 && criterion[i]==NO_PROPAGATION) break;		// skip if no alternate criterion
		if(tractionMat[i]>0)
		{	if(tractionMat[i]>nmat)
				throw CommonException("Material with undefined traction law material for propagation","MaterialBase::PreliminaryMatCalcs");
			if(!theMaterials[tractionMat[i]-1]->isTractionLaw())
				throw CommonException("Material with propagation material that is not a traction law","MaterialBase::PreliminaryMatCalcs");
		}
	}
}

// when set, return total number of materials if this is a new one, or 1 if not in multimaterial mode
int MaterialBase::SetField(int fieldNum,bool multiMaterials,int matid)
{	if(!multiMaterials)
	{	field=0;
		fieldNum=1;
	}
	else
	{	if(field<0)
		{	field=fieldNum;
			fieldNum++;
			fieldMatIDs.push_back(matid);
		}
	}
	return fieldNum;
}

// -1 if material not in use, otherwise zero-based field number
int MaterialBase::GetField(void) { return field; }

// maximum diffusion coefficient in cm^2/sec (anisotropic must override)
double MaterialBase::MaximumDiffusion(void) { return diffusionCon/100.; }

// maximum diffusivity in cm^2/sec  (anisotropic must override)
double MaterialBase::MaximumDiffusivity(void) { return kCond/(rho*heatCapacity*100.); }

// material-to-material contact
void MaterialBase::SetFriction(double friction,int matID,double Dn,double Dnc,double Dt)
{	
	if(lastFriction==NULL)
	{	lastFriction=(ContactDetails *)malloc(sizeof(ContactDetails));
		lastFriction->nextFriction=NULL;
	}
	else
	{	ContactDetails *newFriction=(ContactDetails *)malloc(sizeof(ContactDetails));
		newFriction->nextFriction=(char *)lastFriction;
		lastFriction=newFriction;
	}
	lastFriction->friction=friction;
	lastFriction->matID=matID;
	lastFriction->Dn=Dn;
	lastFriction->Dnc=Dnc;
	lastFriction->Dt=Dt;
}

// Look for contact to a given material and return contact details if fond
ContactDetails *MaterialBase::GetContactToMaterial(int otherMatID)
{	ContactDetails *currentFriction=lastFriction;
	while(currentFriction!=NULL)
	{	if(currentFriction->matID==otherMatID) return currentFriction;
		currentFriction=(ContactDetails *)currentFriction->nextFriction;
	}
	return NULL;
}	

// Print contact law settings for cracks and finalize variables
void MaterialBase::ContactOutput(int thisMatID)
{
	char hline[200];
	
	ContactDetails *currentFriction=lastFriction;
	
	if(currentFriction!=NULL)
		cout << "Custom contact between material " << thisMatID << " and" << endl;
	
	while(currentFriction!=NULL)
	{	// Custom Contact
		if(currentFriction->friction<-10.)
		{   currentFriction->law=NOCONTACT;
			sprintf(hline,"contact nodes revert to center of mass velocity field");
		}
		else if(currentFriction->friction<-.5)
		{   currentFriction->law=STICK;
			sprintf(hline,"stick conditions");
		}
		else if(DbleEqual(currentFriction->friction,(double)0.))
		{   currentFriction->law=FRICTIONLESS;
			sprintf(hline,"frictionless sliding");
		}
		else if(currentFriction->friction>10.)
		{   currentFriction->law=IMPERFECT_INTERFACE;
			if(currentFriction->Dnc<-100.) currentFriction->Dnc=currentFriction->Dn;
			sprintf(hline,"imperfect interface\n     Dn = %g MPa/mm, Dnc = %g MPa/mm, Dt = %g MPa/mm",
					currentFriction->Dn,currentFriction->Dnc,currentFriction->Dt);
		}
		else
		{   currentFriction->law=FRICTIONAL;
			sprintf(hline,"frictional with coefficient of friction: %.6f",currentFriction->friction);
		}
		
		cout << "     material " << currentFriction->matID << ": " << hline << endl;
		
		currentFriction=(ContactDetails *)currentFriction->nextFriction;
	}
}

#pragma mark MaterialBase::Methods

// if cannot be used in current analysis type throw an exception
void MaterialBase::MPMConstLaw(int np) {}

// MPM call to allow material to change properties depending on particle state
// The base method assumes angle is only variable and loads possible
//     rotated meechanical properties (which does nothing unless overridden)
void MaterialBase::LoadMechanicalProps(MPMBase *mptr,int np)
{	LoadMechProps(TRUE,mptr->GetRotationZ(),np);
}

// get transport property tensors (if change with particle state)
void MaterialBase::LoadTransportProps(MPMBase *mptr,int np) { return; }

// implemented in case heat capacity changes with particle state (Cp and Cv)
double MaterialBase::GetHeatCapacity(MPMBase *mptr) { return heatCapacity; }
double MaterialBase::GetHeatCapacityVol(MPMBase *mptr) { return heatCapacityVol; }

// Correct stress update for rotations using hypoelasticity approach
// rotD is -wxy of rotation tensor (or minus twice tensorial rotation, i.e. dvxy-dvyx = -(dvyx-dvxy))
// dsxx, dsyy, and dtxy are unrotated updates to sxx, syy, and txy
// This methods increments particle angle, rotational strain, and in-plane stresses
// This update by a midpoint rule update is documented in Mathematica notes (HypoUpdate.nb)
void MaterialBase::Hypo2DCalculations(MPMBase *mptr,double rotD,double dsxx,double dsyy,double dtxy)
{
	mptr->IncrementRotationStrain(-rotD);
	
	Tensor *sp=mptr->GetStressTensor();
    dsxx+=sp->xy*rotD;
	dsyy-=sp->xy*rotD;
    dtxy-=0.5*(sp->xx-sp->yy)*rotD;
					
    rotD*=0.5;					// comment out for midpoint rule to have rotD=D (endpoint rule) or D/2 (midpoint rule)
	double rotD2=rotD*rotD;		// D^2 (endpoint rule) or D^2/4 (midpoint rule)
    double det=1.+rotD2;		// 1+D^2 (endpoint rule) or 1+D^2/4 (midpoint rule)
    sp->xx+=((1.+0.5*rotD2)*dsxx + 0.5*rotD2*dsyy + rotD*dtxy)/det;
    sp->yy+=(0.5*rotD2*dsxx + (1.+0.5*rotD2)*dsyy - rotD*dtxy)/det;
	sp->xy+=(0.5*rotD*(dsyy-dsxx) + dtxy)/det;
}

// Correct stress update for rotations using hypoelasticity approach
// dwxy, dwxz, and dwyz are the engineering rotational strains
// dsig[6] on initial stress updates in order (xx,yy,zz,yz,xz,xy)
// This methods increments particle angles (will), rotational strains, and stresses
// This linear approximation to a midpoint rule update is documented in Mathematica notes (HypoUpdate.nb)
void MaterialBase::Hypo3DCalculations(MPMBase *mptr,double dwxy,double dwxz,double dwyz,double *dsig)
{
	// rotational strains
	mptr->IncrementRotationStrain(dwxy,dwxz,dwyz);
	
	// get stress tensor
	Tensor *sp=mptr->GetStressTensor();
	/*
	// nonhypo test
	sp->xx+=dsig[XX];
	sp->yy+=dsig[YY];
	sp->zz+=dsig[ZZ];
	sp->yz+=dsig[YZ];
	sp->xz+=dsig[XZ];
	sp->xy+=dsig[XY];
	return;
	*/
	
	// angles (for anisotropic materials)
	
	// stresss - first those involving current stress
	Tensor st;
	st.xx= -dwxy*sp->xy - dwxz*sp->xz;
	st.yy=  dwxy*sp->xy               - dwyz*sp->yz;
	st.zz=                dwxz*sp->xz + dwyz*sp->yz;
	st.yz= 0.5*( +dwxy*sp->xz          + dwxz*sp->xy          + dwyz*(sp->yy-sp->zz)  );
	st.xz= 0.5*( -dwxy*sp->yz          + dwxz*(sp->xx-sp->zz) + dwyz*sp->xy  );
	st.xy= 0.5*( +dwxy*(sp->xx-sp->yy) - dwxz*sp->yz          - dwyz*sp->xz );
	
	// now add all
	sp->xx+=dsig[XX] + st.xx - 0.5*dsig[XY]*dwxy - 0.5*dsig[XZ]*dwxz;
	sp->yy+=dsig[YY] + st.yy + 0.5*dsig[XY]*dwxy - 0.5*dsig[YZ]*dwyz;
	sp->zz+=dsig[ZZ] + st.zz + 0.5*dsig[XZ]*dwxz + 0.5*dsig[YZ]*dwyz;
	sp->yz+=dsig[YZ] + st.yz + 0.25*( +dsig[XZ]*dwxy + dsig[XY]*dwxz + (dsig[YY]-dsig[ZZ])*dwyz );
	sp->xz+=dsig[XZ] + st.xz + 0.25*( -dsig[YZ]*dwxy + (dsig[XX]-dsig[ZZ])*dwxz + dsig[XY]*dwyz  );
	sp->xy+=dsig[XY] + st.xy + 0.25*( +(dsig[XX]-dsig[YY])*dwxy - dsig[YZ]*dwxz - dsig[XZ]*dwyz  );
}

#pragma mark MaterialBase::Fracture Calculations

// set propagate and alternate propagate
MaterialBase *MaterialBase::SetFinalPropagate(void)
{	int i;
	for(i=0;i<=1;i++)
	{	if(criterion[i]==UNSPECIFIED)	criterion[i]=fmobj->propagate[i];
		if(matPropagateDirection[i]==UNSPECIFIED) matPropagateDirection[i]=fmobj->propagateDirection[i];
		if(tractionMat[i]<=0) tractionMat[i]=fmobj->propagateMat[i];
	}
	return (MaterialBase *)GetNextObject();
}

// stress intensity factor - zero if material does not know better
Vector MaterialBase::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{ 
   Vector SIF;
   SIF.x=0.0;
   SIF.y=0.0;
   return SIF;
}

// Determine what calculations are needed for the propagation criterion
// in this material - must match needs in ShouldPropagate() routine
int MaterialBase::CriterionNeeds(int critIndex)
{
	switch(criterion[critIndex])
	{	case MAXHOOPSTRESS:
		case STRAINENERGYDENSITY:
		case EMPIRICALCRITERION:
			return NEED_JANDK;
		
		case STEADYSTATEGROWTH:
			if(initTime<0.) return NEED_J;
			return FALSE;
		
		case TOTALENERGYBALANCE:
		case CRITICALERR:
			return NEED_J;
			
		case MAXCTODCRITERION:
			return FALSE;
		
		case NO_PROPAGATION:
			return FALSE;
		
		default:
			return NEED_JANDK;		// just to be sure
	}
}

// Crack Propagation Criterion
//	Input is crack tip segment (which will have needed parameters)
//		and crackDir is direction at crack tip (unit vector)
//	Output is
//		NOGROWTH - nothing more to do
//		GROWNOW - propagate and see if speed needs control
//	If GROWNOW, change crackDir to unit
//		vector in crack growth direction
//  If need J or K, must say so in CriterionNeeds() routine
int MaterialBase::ShouldPropagate(CrackSegment *crkTip,Vector &crackDir,CrackHeader *theCrack,int np,int critIndex)
{	
    double KI,KII,fCriterion,cosTheta0,sinTheta0;
    double deltaPotential,deltaPlastic,deltaLength,p;
    double balance,adjustSpeed,cosTheta2;
	Vector hoopDir;
    //double deltaHPlastic,deltaTime,dUirrda,dUirrdtCona,dUirrdaCont,avgSpeed;

    // retrieve fracture parameters
    switch(criterion[critIndex])
	{	// Criterion 1
        case MAXHOOPSTRESS:
            // Maximum hoop stress (or maximum principal stress) criterion only applies to
            // isotropic elastic materials
            KI=crkTip->sif.x;
            KII=crkTip->sif.y;
			HoopDirection(KI,KII,&hoopDir);
            
            // failure criterion
            cosTheta2=sqrt((1.+hoopDir.x)/2.);
            fCriterion=KI*pow(cosTheta2,3)-1.5*KII*cosTheta2*hoopDir.y-KIc;
        
            // growth, return direction (unless overridden)
            if(fCriterion>0.)
			{	if(!SelectDirection(crkTip,crackDir,theCrack,critIndex))
					RotateDirection(&crackDir,hoopDir.x,hoopDir.y);
                return GROWNOW;
            }
            break;
        
		// Criterion 7
		case CRITICALERR:
            // growth, direction by direction option
            if(crkTip->Jint.x>=JIc)
			{	SelectDirection(crkTip,crackDir,theCrack,critIndex);
				return GROWNOW;
            }
			break;
			
		// Criterion 2
        case STEADYSTATEGROWTH:
            switch(crkTip->steadyState)
            {	case STATIONARY:
                    // If J now > JIc start propagating at constant speed
					//   or if time was specified, go on that time
					if(initTime<0.)
                    {	if(crkTip->Jint.x<JIc) return NOGROWTH;
					}
					else
					{	if(mtime<initTime) return NOGROWTH;
					}
                    crkTip->steadyState=PROPAGATING;
                    crkTip->speed=initSpeed;		// constant speed
                    break;
                
                case PROPAGATING:
                    // continue propagating unless reached maximum length
                    if(maxLength>0)
                    {	if(theCrack->Length()>=maxLength)
                        {   crkTip->steadyState=ARRESTING;
                            return NOGROWTH;
                        }
                    }
					if(fmobj->dflag[0]==4)
					{	// stop in cutting simulation when reach edge of the sample
						if(crkTip->x<1.)
                        {   crkTip->steadyState=ARRESTING;
                            return NOGROWTH;
                        }
					}
                    break;
                
                case ARRESTING:
                    // ARRESTING applies to first interval after propagation stops
                    crkTip->steadyState=ARRESTED;
                case ARRESTED:
                    return NOGROWTH;
                
                default:
                    break;
            }
            
            // if not overriden self similar or a constant direction
			if(!SelectDirection(crkTip,crackDir,theCrack,critIndex))
			{	if(constantDirection) crackDir=growDir;		// if direction specified in material
			}
            return GROWNOW;
            break;
        
		// Criterion 3
        case TOTALENERGYBALANCE:
            switch(crkTip->steadyState)
            {	case STATIONARY:
                    // decide if it should now start
                    if(crkTip->Jint.x>=JIc)
                    {	// next check will just by to save plastic energy
						// (Note that JIc may differ from 2*gamma if desired
						// but must have JIc >= 2*gamma)
                    	crkTip->steadyState=BEGINPROPAGATING;
                        crkTip->speed=initSpeed*WaveSpeed(FALSE);		// fraction of wave speed in 2D
                        cout << "# Initiate t:" << 1000.*mtime <<
                                " s:" << crkTip->speed << endl;
						SelectDirection(crkTip,crackDir,theCrack,critIndex);
                        return GROWNOW;
                    }
                    break;
                
                case BEGINPROPAGATING:
                    cout << "# First Propagation t:" << 1000.*mtime << endl;
                    crkTip->steadyState=PROPAGATING;
					SelectDirection(crkTip,crackDir,theCrack,critIndex);
                    return GROWNOW;
                    break;
                
                case PROPAGATING:
                case SLOWLYPROPAGATING:
                    /* decide if crack speed should change
                        a. Find J-p*(dUplast/da)-Jc (in N/mm)
                        b. Change speed by gain*balance (within limits)
                    */
                    deltaPotential=crkTip->potential[0]-crkTip->potential[2];
                    deltaPlastic=crkTip->plastic[0]-crkTip->plastic[2];
                    deltaLength=crkTip->clength[0]-crkTip->clength[2];
                    p=pCrit3;
                    
                    // The energy change derviative is -balance
                    balance=-(deltaPotential+(1.+p)*deltaPlastic)/deltaLength - 2.*gamma;
                    //balance=crkTip->Jint.x - p*deltaPlastic/deltaLength - 2.*gamma;
                    adjustSpeed=gain*balance;
                    
                    // balance<0 means energy is increasing, this crack should not be growing
                    //    so here we slow it down
                    if(balance<0.)
                    {	if(crkTip->steadyState==SLOWLYPROPAGATING)
                        {   cout << "# No Growth t:" << 1000.*mtime <<
                            " J:" << crkTip->Jint.x <<
                            " J-DU:" << -(deltaPotential+deltaPlastic)/deltaLength <<
                            " DUirr:" <<  p*deltaPlastic/deltaLength <<
                            " -DU:" <<  balance  << endl;
                            return NOGROWTH;
                        }
                        adjustSpeed=-fmin(-adjustSpeed,0.5*crkTip->speed);
                    }
                    
                    // balance>0 means energy is decreasing. Here we speed up the
                    //    crack to try and use it up, but not above wave speed
                    else
                    {	adjustSpeed=fmin(adjustSpeed,WaveSpeed(FALSE)-crkTip->speed);
                        crkTip->steadyState=PROPAGATING;
                    }
                    
                    // adjust speed, but never below minimum
                    crkTip->speed+=adjustSpeed;
                    if(crkTip->speed<=0.5*initSpeed*WaveSpeed(FALSE))
                    {	crkTip->speed=0.5*initSpeed*WaveSpeed(FALSE);
                        crkTip->steadyState=SLOWLYPROPAGATING;
                    }
					
					// archive energy flow (in N/mm)
					crkTip->release=-(deltaPotential+deltaPlastic)/deltaLength;
					crkTip->absorb=p*deltaPlastic/deltaLength + 2*gamma;
					crkTip->crackIncrements++;
                    
                    cout << "# Growing t:" << 1000.*mtime << 
                            " J:" << crkTip->Jint.x <<
                            " J-DU:" << -(deltaPotential+deltaPlastic)/deltaLength <<
                            " DUirr:" << p*deltaPlastic/deltaLength <<
                            " -DU:" << balance <<
                            " Ds:" << adjustSpeed <<
                            " s:" << crkTip->speed << endl;
					SelectDirection(crkTip,crackDir,theCrack,critIndex);
                    return GROWNOW;
                
                default:
                    break;
            }
            break;
		
		// Criterion 4
        case STRAINENERGYDENSITY:
            // Minimum strain energy density criterion only applies to 
            // isotropic elastic materials
            KI=crkTip->sif.x;
            KII=crkTip->sif.y;

            // Poisson's ratio and the constant k
            double v,k;
            v=C12/C11;
            k=(np==PLANE_STRESS_MPM)? (3-v)/(1+v) : (3-4*v);

            // Crack propagation direction
            // pure mode I
            if(DbleEqual(KII,0.0))
            {   cosTheta0 = 1.0;
                sinTheta0 = 0.0;
            }

            // pure mode II or mixed mode
            // crack propagation direction by solving the criterion equation numerically.
            else
            {   double theta0=CrackPropagationAngleFromStrainEnergyDensityCriterion(k,KI,KII); 
                cosTheta0=cos(theta0);
                sinTheta0=sin(theta0);
            }

            // failure criterion of strain energy density
            double a11,a12,a22;
            a11=(1+cosTheta0)*(k-cosTheta0);
            a12=sinTheta0*(2*cosTheta0-k+1);
            a22=(k+1)*(1-cosTheta0)+(1+cosTheta0)*(3*cosTheta0-1);
            fCriterion=sqrt((a11*KI*KI+2*a12*KI*KII+a22*KII*KII)/2/(k-1))-KIc;

            // growth, return direction
            if(fCriterion>0.)
			{	if(!SelectDirection(crkTip,crackDir,theCrack,critIndex))
					RotateDirection(&crackDir,cosTheta0,sinTheta0);
                return GROWNOW;
            }
            break;

		// Criterion 5
        case EMPIRICALCRITERION:
            // Empirical criterion (KI/KIc)^p+(KII/KIIc)^q=1 only applies to
            // isotropic elastic materials.
            KI=crkTip->sif.x;
            KII=crkTip->sif.y;

            // Crack propagation direction
            // pure mode I
            if(DbleEqual(KII,0.0))
            {   cosTheta0 = 1.0;
                sinTheta0 = 0.0;
            }               
                            
            // For pure mode II or mixed mode, detremine crack propagation direction
            // by maximum principal stress criterion since the empirical criterion can only 
            // describe fracture locus. 
            else        
            {   double R=KI/KII;
                int sign = (KII>=0.)? -1 : 1;
                double theta0=2*atan((R+sign*sqrt(R*R+8))/4.);
                cosTheta0=cos(theta0);
                sinTheta0=sin(theta0);
            }       

            // failure criterion of strain energy density
            fCriterion=(pow(KI/KIc,KIexp)+pow(KII/KIIc,KIIexp))*KIc-KIc;
                    
            // growth, return direction
            if(fCriterion>0.)
			{	if(!SelectDirection(crkTip,crackDir,theCrack,critIndex))
					RotateDirection(&crackDir,cosTheta0,sinTheta0);
                return GROWNOW;         
            }
            break;
		
		// criterion 6 - fail if either CTOD exists, default by self simlar
		//	some critical value
		case MAXCTODCRITERION:
			double codn,cods;
			Vector cod;
			theCrack->GetCOD(crkTip,cod,true);
			cods=fabs(cod.x);
			codn=fabs(cod.y);
			if(cods>delIIc && delIIc>0.)
			{	SelectDirection(crkTip,crackDir,theCrack,critIndex);
				return GROWNOW;
			}
			if(codn>delIc && delIc>0.)
			{	SelectDirection(crkTip,crackDir,theCrack,critIndex);
				return GROWNOW;
			}
			break;

        default:
            break;
    }
    
    return NOGROWTH;
}

/* Pick crack growth direction according to selected option and install it into
	crackDir. If selected return TRUE or return FALSE for criterion to
	pick its own direction
*/

bool MaterialBase::SelectDirection(CrackSegment *crkTip,Vector &crackDir,CrackHeader *theCrack,int critIndex)
{
	Vector cod;
	double codamp;
	
	switch(matPropagateDirection[critIndex])
	{	case DEFAULT_DIRECTION:
			return FALSE;
		
		case NORMAL_TO_COD:
			theCrack->GetCOD(crkTip,cod,false);
			codamp=sqrt(cod.x*cod.x+cod.y*cod.y);
			crackDir.x=cod.y/codamp;
			crackDir.y=-cod.x/codamp;
			break;
			
		case HOOP_FROM_COD:
			theCrack->GetCOD(crkTip,cod,true);
			Vector hoopDir;
			HoopDirection(cod.y,cod.x,&hoopDir);
			RotateDirection(&crackDir,hoopDir.x,hoopDir.y);
			break;
		
		case INITIAL_DIRECTION:
			theCrack->GetInitialDirection(crkTip,crackDir);
			break;
		
		case SELF_SIMILAR:
		default:
			// no change but return TRUE
			break;
	}
	
	return TRUE;
}

/* find hoop direction unit vector relative to crack direction. Vector in hoop direction will by 
		(crackDir.x*hoopDir.x - crackDir.y*hoopDir.y, crackDir.x*hoopDir.y + crackDir.y*hoopDir.x)
	If theta is the ccw rotation of hoop direction relative to crack tip direction, then
		hooopDir = (cos(theta), sin(theta))
*/
void MaterialBase::HoopDirection(double KI,double KII,Vector *hoopDir)
{
	// pure mode II (70.5 degrees)
	if(DbleEqual(KI,0.0))
	{	hoopDir->x=1./3.;
		hoopDir->y=sqrt(8./9.);
		if(KII>=0.) hoopDir->y=-hoopDir->y;
	}
	
	// pure mode I or mixed mode
	else
	{	double KI2=KI*KI;
		double KII2=KII*KII;
		
		// these really obey hoopDir->x^2 + hoopDir->y^2 = 1
		hoopDir->x=(3*KII2 + KI*sqrt(KI2+8.*KII2))/(KI2+9*KII2);
		hoopDir->y=fabs(KII*(3*hoopDir->x-1.)/KI);
		if(KII>=0.) hoopDir->y=-hoopDir->y;
	}
}

// Rotate vector ccw by angle theta give cos(theta) and sin(theta)
void MaterialBase::RotateDirection(Vector *crackDir,double cosTheta,double sinTheta)
{
	double dx=crackDir->x;
	double dy=crackDir->y;
	crackDir->x=dx*cosTheta - dy*sinTheta;
	crackDir->y=dx*sinTheta + dy*cosTheta;
}

// Get crack propagation angle by solving the equation 
// of strain energy density criterion numerically
double 
MaterialBase::CrackPropagationAngleFromStrainEnergyDensityCriterion(double k,
                                                       double KI, double KII)
{
  double theta=0.0;
  double errF=1.e-6,errV=1.e-2;
  double a,b,c,fa,fb,fc;

  double A=-PI_CONSTANT, B=PI_CONSTANT;   // The region of the roots 
  int n=36;             // Divide [A,B] into n intervals
  double h=(B-A)/n;     // Subinterval length

  double theta0=atan(KI/KII);
  vector<double> root;  // Store the solutions of the equation
  // Solve the equation numerically
  for(int i=0; i<n; i++) { // Loop over the whole interval [A,B]
    a=A+i*h;
    b=A+(i+1)*h;
    fa=(k-1)*sin(a-2*theta0)-2*sin(2*(a-theta0))-sin(2*a);
    fb=(k-1)*sin(b-2*theta0)-2*sin(2*(b-theta0))-sin(2*b);

    // Find the root in [a,b)
    if(fabs(fa)<errF) { // Where f(a)=0
      root.push_back(a);
    }
    else if(fa*fb<0.) { // There is a root in (a,b)
      double cp=2*B;    // Set the value beyond [A,B]      
      for(int j=0; j<32768; j++) { // 32768=2^15 (a big int) 
        c=b-(a-b)*fb/(fa-fb);
        fc=(k-1)*sin(c-2*theta0)-2*sin(2*(c-theta0))-sin(2*c);
        if(fabs(fc)<errF || fabs(c-cp)<errV) { // c is the root
          root.push_back(c);
          break;
        }
        else { // Record the cross point with axis x
          cp=c;
        }

        // Narrow the region of the root
        if(fc*fa<0.) { // The root is in (a,c)
          fb=fc; b=c;
        }
        else if(fc*fb<0.) { // The root is in (c,b)
          fa=fc; a=c;
        }
      } // End of loop over j
    } // End of if(fa*fb<0.)
  } // End of loop over i

  // Select the direction from the solutions 
  // along which there exists the minimum strain energy density   
  int count=0;
  double S0=0.0;
  for(int i=0;i<(int)root.size();i++) {
    double r=root[i];

    // The signs of propagation angle and KII must be opposite 
    if(KII*r>0.) continue;

    // Calculate the second derivative of the strain energy density
    double sr=sin(r),cr=cos(r),sr2=sin(2*r),cr2=cos(2*r);
    double dsdr2=KI*KI*((1-k)*cr+2*cr2)-2*KI*KII*(4*sr2+(1-k)*sr)+
                 KII*KII*((k-1)*cr-6*cr2);
    if(dsdr2>0.) {
      // Determine propagation angle by comparison of strain energy density. 
      // Along the angle there exists the minimum strain energy density. 
      double S=(1+cr)*(k-cr)*KI*KI+2*sr*(2*cr-k+1)*KI*KII+
               ((k+1)*(1-cr)+(1+cr)*(3*cr-1))*KII*KII;
      if(count==0 || (count>0 && S<S0)) {
        theta=r;
        S0=S;
        count++;
      }
    }
  } // Enf of loop over i
  root.clear();

  return theta;
}

// If needed, adjust next time to check on propagation to control
//   crack speed
bool MaterialBase::ControlCrackSpeed(CrackSegment *crkTip,double &waitTime)
{
    // check if crack speed should be controlled
    switch(fmobj->propagate[0])
    {   case STEADYSTATEGROWTH:
            if(crkTip->steadyState==PROPAGATING)
            {	waitTime=crkTip->theGrowth/crkTip->speed;
                return TRUE;
            }
            break;
        
        case TOTALENERGYBALANCE:
            if(crkTip->steadyState==PROPAGATING)
            {	waitTime=crkTip->theGrowth/crkTip->speed;
                return TRUE;
            }
            break;
            
        default:
            break;
    }
    
    return FALSE;
}

#pragma mark MaterialBase::Accessors

// Provide something
double MaterialBase::ShearWaveSpeed(bool threeD) { return WaveSpeed(threeD)/sqrt(3.); }

// archive material data for this material type when requested.
double MaterialBase::GetHistory(int num,char *historyPtr) { return (double)0.; }

// return TRUE if rigid particle (for contact or for BC)
short MaterialBase::Rigid(void) { return FALSE; }

// return TRUE is rigid BC particle (not rigid for contact)
short MaterialBase::RigidBC(void) { return FALSE; }

// return TRUE is rigid BC particle (not rigid for contact)
short MaterialBase::RigidContact(void) { return FALSE; }

// check if traciton law material
bool MaterialBase::isTractionLaw(void) { return FALSE; }

// return pointer to k conduction tensor
Tensor *MaterialBase::GetkCondTensor(void) { return &kCondTensor; }

// return pointer to diffusion tensor
Tensor *MaterialBase::GetDiffusionTensor(void) { return &diffusionTensor; }

// get material rho (g/mm^3) for material velocity field (can only call in multimaterial mode)
double MaterialBase::GetMVFRho(int matfld) { return theMaterials[fieldMatIDs[matfld]]->rho*0.001; }

// see if material for a material velocity field is rigid (only rigid contact materials can be in a velocity field)
short MaterialBase::GetMVFIsRigid(int matfld)
{	return matfld<(int)fieldMatIDs.size() ? theMaterials[fieldMatIDs[matfld]]->Rigid() : FALSE ;
}

// convert field number (zero based) to material ID for that field (zero based)
int MaterialBase::GetFieldMatID(int matfld) { return fieldMatIDs[matfld]; }

// Get current relative volume change - only used to convert speific results to actual values when archiving
// Material with explicit treatment of large deformation might need it (e.g., Hyperelastic)
double MaterialBase::GetCurrentRelativeVolume(MPMBase *mptr,bool threeD) { return 1.; }





