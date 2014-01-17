/********************************************************************************
    MaterialBaseMPM.cpp - more MaterialBase for MPM code
    nairn-mpm-fea
    
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
#include "Cracks/CrackSurfaceContact.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Materials/LinearHardening.hpp"
#include "Materials/NonlinearHardening.hpp"
#include "Materials/JohnsonCook.hpp"
#include "Materials/SCGLHardening.hpp"
#include "Materials/SLMaterial.hpp"
#include "Materials/Nonlinear2Hardening.hpp"
#include "Materials/DDBHardening.hpp"
#include <vector>

// global
bool MaterialBase::isolatedSystemAndParticles = FALSE;

// class statics for MPM - zero based material IDs when in multimaterial mode
vector<int> MaterialBase::fieldMatIDs;
vector<int> MaterialBase::activeMatIDs;
int MaterialBase::incrementalDefGradTerms = 2;			// terms in exponential of deformation gradient
int MaterialBase::maxPropertyBufferSize = 0.;           // maximum buffer size needed among active materials to get copy of mechanical properties
int MaterialBase::maxAltBufferSize = 0.;                // maximum optional buffer size needed for more properties (e.g., hardenling law)


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
	
	else if(strcmp(xName,"constantTip")==0)
	{	input=INT_NUM;
		return((char *)&constantTip);
	}
    
	else if(strcmp(xName,"Cp")==0 || strcmp(xName,"Cv")==0)
    {	input=DOUBLE_NUM;
        return((char *)&heatCapacity);
    }

	else if(strcmp(xName,"csat")==0)
    {	input=DOUBLE_NUM;
        return((char *)&concSaturation);
    }

	else if(strcmp(xName,"D")==0)
    {	input=DOUBLE_NUM;
        return((char *)&diffusionCon);
	}
    
    else if(strcmp(xName,"beta")==0)
    {	input=DOUBLE_NUM;
        return((char *)&betaI);
    }
	
	else if(strcmp(xName,"kCond")==0)
    {	input=DOUBLE_NUM;
        return((char *)&kCond);
	}
    
	// check properties only for some materials
	if(SupportsArtificialViscosity())
	{	if(strcmp(xName,"ArtificialVisc")==0)
		{	artificialViscosity=TRUE;
			input=NOT_NUM;
			return((char *)&artificialViscosity);
		}
		
		else if(strcmp(xName,"avA1")==0)
		{	input=DOUBLE_NUM;
			return((char *)&avA1);
		}
		
		else if(strcmp(xName,"avA2")==0)
		{	input=DOUBLE_NUM;
			return((char *)&avA2);
		}
	}
    
    return((char *)NULL);
}

// Material that allow hardening laws must accept
// the final call and use if supported
// Newly created laws should be added here
void MaterialBase::SetHardeningLaw(char *lawName)
{
    HardeningLawBase *pLaw = NULL;
    int lawID = 0;
    
    // check options
    if(strcmp(lawName,"Linear")==0 || strcmp(lawName,"1")==0)
    {   pLaw = new LinearHardening(this);
        lawID = 1;
    }
    
    else if(strcmp(lawName,"Nonlinear")==0 || strcmp(lawName,"2")==0)
    {   pLaw = new NonlinearHardening(this);
        lawID = 2;
    }
    
    else if(strcmp(lawName,"Nonlinear2")==0 || strcmp(lawName,"6")==0)
    {   pLaw = new Nonlinear2Hardening(this);
        lawID = 6;
    }
    
    else if(strcmp(lawName,"JohnsonCook")==0 || strcmp(lawName,"3")==0)
    {   pLaw = new JohnsonCook(this);
        lawID = 3;
    }
    
    else if(strcmp(lawName,"SCGL")==0 || strcmp(lawName,"4")==0)
    {   pLaw = new SCGLHardening(this);
        lawID = 4;
    }

    else if(strcmp(lawName,"SL")==0 || strcmp(lawName,"5")==0)
    {   pLaw = new SLMaterial(this);
        lawID = 5;
    }
    
    else if(strcmp(lawName,"DDB-PPM")==0 || strcmp(lawName,"7")==0)
    {   pLaw = new DDBHardening(this);
      lawID = 7;
    }
    
    // was it found
    if(pLaw==NULL)
        ThrowSAXException("The hardening law '%s' is not valid",lawName);
    
    // pass on to the material
    if(!AcceptHardeningLaw(pLaw,lawID))
        ThrowSAXException("The hardening law '%s' is not allows by material",lawName);
}

// A Material that allows hardening laws should accept this one or it can
// veto the choice by returning FALSE
bool MaterialBase::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID) { return FALSE; }

/*	Verify andcCalculate properties used in analyses. If error, return string with an error message.
 This is called once at start of the calculation just before the material properties
 are printined to the output file and before any calculations. It is called for
 every material defined in the input file, even if it is not used by any
 particle.
 If superclass overrides this method, must call this too
 */
const char *MaterialBase::VerifyAndLoadProperties(int np)
{
    // make conductivity specific (N mm^3/(sec-K-g) = mJ mm^2/(sec-K-g))
    kCond *= (1000./rho);
	
	// in case only need to load some things once, load those mechanical properties now
	FillTransportProperties(&tr);
	
	// convert MPa sqrt(m) to MPa sqrt(mm)
	KIc*=31.62277660168379;
	KIIc*=31.62277660168379;
	
	// convert other crack growth properties
	JIc*=1.e-3;				// convert J/m^2 to N/mm
	initTime*=1e-3;			// convert to sec
	if(criterion[0]==TOTALENERGYBALANCE)
		initSpeed*=0.01;	// convert % of WaveSpeed() to fraction of WaveSpeed()
	else
		initSpeed*=1.e3;	// convert m/sec to mm/sec
	gamma*=1.e-3;			// convert J/m^2 to N/mm
	if(gamma<0. || JIc<2.*gamma) gamma=JIc/2.;
	// pCrit3 - dimensionless
	// gain - no units change
	
	double norm=sqrt(growDir.x*growDir.x + growDir.y*growDir.y);
	growDir.x/=norm;
	growDir.y/=norm;
	
	// may not work with two criteria
	if(matPropagateDirection[0]==EMPIRICALCRITERION)
	{	// If KIIc not set, set now to a default value
		if(KIIc<0.) KIIc=0.817*KIc;
	}
	
	return NULL;
}

// print any properties common to all MPM material types
void MaterialBase::PrintCommonProperties(void) const
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
	   
    // artificial visconsity
    if(artificialViscosity)
	{	PrintProperty("Artificial viscosity on",FALSE);
		PrintProperty("AV-A1",avA1,"");
		PrintProperty("AV-A2",avA2,"");
        cout << endl;
        if(ConductionTask::AVHeating)
            PrintProperty("    AV heating on",FALSE);
        else
            PrintProperty("    AV heating off",FALSE);
	}
	
	// optional color
	if(red>=0.)
	{	char mline[200];
		sprintf(mline,"color= %g, %g, %g, %g",red,green,blue,alpha);
		cout << mline << endl;
	}
}

// print fraction criterion
void MaterialBase::PrintCriterion(int thisCriterion,int thisDirection) const
{
	char mline[200];
	
	switch(thisCriterion)
	{   case NO_PROPAGATION:
			cout << "No propagation" << endl;
			break;
			
		case MAXHOOPSTRESS:
			cout << "Maximum hoop stess" << PreferredDirection(thisDirection) << endl;
			PrintProperty("KIc",KIc/31.62277660168379,"MPa-sqrt(m)");
			PrintProperty("KIIc",KIIc/31.62277660168379,"MPa-sqrt(m)");
			cout << endl;
			break;
			
		case CRITICALERR:
			cout << "Critical Energy Release Rate" << PreferredDirection(thisDirection) << endl;
			PrintProperty("Jc",1000.*JIc,"J/m^2");
			cout << endl;
			break;
			
		case STEADYSTATEGROWTH:
			cout << "Constant crack speed" << PreferredDirection(thisDirection) << endl;
			if(initTime<0.)
			{	PrintProperty("Jc",JIc,"J/m^2");
				PrintProperty("initSpeed",1.e-3*initSpeed,"m/s");
			}
			else
			{	PrintProperty("ti",1000.*initTime,"ms");
				PrintProperty("initSpeed",1.e-3*initSpeed,"m/s");
			}
			cout << endl;
			if(maxLength>0.)
			{	PrintProperty("max length",maxLength,"mm");
				cout << endl;
			}
			if(constantDirection && thisDirection==DEFAULT_DIRECTION)
			{	sprintf(mline,"direction = (%9.5f,%9.5f)",growDir.x,growDir.y);
				cout << mline << endl;
			}
			break;
			
		case TOTALENERGYBALANCE:
			cout << "Total energy balance" << PreferredDirection(thisDirection) << endl;
			sprintf(mline,"Jc =%12.3f J/m^2  vel=%12.3f%c wave speed",JIc,100.*initSpeed,'%');
			cout << mline << endl;
			sprintf(mline,"gam=%12.3f J/m^2  p  =%12.3f        gain=%12.3g",1000.*gamma,pCrit3,gain);
			cout << mline << endl;
			break;
			
		case STRAINENERGYDENSITY:
			cout << "Minmum strain energy density" << PreferredDirection(thisDirection) << endl;  
			sprintf(mline,"KIc=%12.3f MPa-sqrt(m)  KIIc=%12.3f MPa-sqrt(m)",KIc,KIIc);
			cout << mline << endl;
			break;
			
		case EMPIRICALCRITERION:
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
void MaterialBase::PrintTransportProperties(void) const
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
	{	PrintProperty("k",rho*kCond/1000.,"W/(m-K)");
		PrintProperty("Cv",heatCapacity,"J/(kg-K)");        // aka mJ/(g-K)
		PrintProperty("Cp",heatCapacity+GetCpMinusCv(NULL),"J/(kg-K)");        // aka mJ/(g-K)
		cout << endl;
	}
	else if(ConductionTask::adiabatic)
	{	PrintProperty("Cv",heatCapacity,"J/(kg-K)");        // aka mJ/(g-K)
        // Cp only used in conduction so not printed here when conduction is off
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

// Called before analysis, material can fill in things that never change during the analysis
// Note: no angle, because cannot depend on material angle
// Here fills in isotropic properties, materials with different anisotropic properties should override
void MaterialBase::FillTransportProperties(TransportProperties *t)
{
	// diffusion tensor (xx, yy, xy order)
	t->diffusionTensor.xx=diffusionCon;
	t->diffusionTensor.yy=diffusionCon;
	t->diffusionTensor.zz=diffusionCon;
	t->diffusionTensor.xy=0.;
	t->diffusionTensor.xz=0.;
	t->diffusionTensor.yz=0.;
	
	// conductivity tensor (xx, yy, xy order) normalized to rho
	t->kCondTensor.xx = kCond;
	t->kCondTensor.yy = kCond;
	t->kCondTensor.zz = kCond;
	t->kCondTensor.xy = 0.;
	t->kCondTensor.xz = 0.;
	t->kCondTensor.yz = 0.;
}

/* This is called after PreliminaryCalcs() and just before first MPM time step and it
		is only called if the material is actually in use by one or more particles
	If material cannot be used in current analysis type throw an exception
	Subclass that overrides must pass on to super class
 */
void MaterialBase::ValidateForUse(int np) const
{	int i;
	
	for(i=0;i<=1;i++)
	{	if(i==1 && criterion[i]==NO_PROPAGATION) break;		// skip if no alternate criterion
		if(tractionMat[i]>0)
		{	if(tractionMat[i]>nmat)
			{	throw CommonException("Material with undefined traction law material for propagation",
									  "MaterialBase::ValidateForUse");
			}
			if(!theMaterials[tractionMat[i]-1]->isTractionLaw())
			{	throw CommonException("Material with propagation material that is not a traction law",
									  "MaterialBase::ValidateForUse");
			}
		}
	}
	
	// check for unsupported alternate propagation criterion
	if(criterion[1]==TOTALENERGYBALANCE)
	{	throw CommonException("The alternate propagation criterion cannot be energy balance method.",
							  "MaterialBase::ValidateForUse");
	}
	if(criterion[1]==STEADYSTATEGROWTH)
	{	throw CommonException("The alternate propagation criterion cannot be steady state crack growth.",
							  "MaterialBase::ValidateForUse");
	}
	
	if(ConductionTask::active || ConductionTask::adiabatic)
	{	if(heatCapacity<=0. && !Rigid())
		{	throw CommonException("Thermal conduction and/or mechanical energy cannot be done using materials that have zero heat capacity.",
								  "MaterialBase::ValidateForUse");
		}
	}
}

// create and return pointer to material-specific data on a particle
//	using this material. Called once at start of analysis for each
//	particle.
// If a material subclass override this method, it is responsible for
//  creating space for super class history varibles too. There is
//  no mechanism now for calling super classes to prepend or
//  append the data.
char *MaterialBase::InitHistoryData(void) { return NULL; }

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
// Such a class must pass on the super class after its own initializations
void MaterialBase::SetInitialParticleState(MPMBase *mptr,int np) const
{
	if(isolatedSystemAndParticles)
    {   // need to add initial heat energy, because special cases in this mode
        // will ignore the Cv dT term
		double Cv = 1000.*GetHeatCapacity(mptr);
		mptr->AddHeatEnergy(Cv*(mptr->pTemperature-thermal.reference));
        // initial entropy kept to zero, could add something here if ever needed
	}
}

// when set, return total number of materials if this is a new one, or 1 if not in multimaterial mode
// not thread safe due to push_back()
int MaterialBase::SetField(int fieldNum,bool multiMaterials,int matid,int &activeNum)
{	if(!multiMaterials)
	{	if(field<0)
		{	field=0;
			activeField=activeNum;
			fieldNum=1;
			activeNum++;
			activeMatIDs.push_back(matid);
            int altBuffer,matBuffer = SizeOfMechanicalProperties(altBuffer);
            if(matBuffer > maxPropertyBufferSize) maxPropertyBufferSize = matBuffer;
            if(altBuffer > maxAltBufferSize) maxAltBufferSize = altBuffer;
		}
	}
	else
	{	if(field<0)
		{	field=fieldNum;
			fieldNum++;
			fieldMatIDs.push_back(matid);
			if(activeNum>=0)
			{	activeField=activeNum;
				activeNum++;
				activeMatIDs.push_back(matid);
                int altBuffer,matBuffer = SizeOfMechanicalProperties(altBuffer);
                if(matBuffer > maxPropertyBufferSize) maxPropertyBufferSize = matBuffer;
                if(altBuffer > maxAltBufferSize) maxAltBufferSize = altBuffer;
			}
		}
	}
	return fieldNum;
}

// -1 if material not in use, otherwise zero-based field number
int MaterialBase::GetField(void) const { return field; }

// -1 if material not in use, otherwise zero-based field number from 0 to numActiveMaterials
int MaterialBase::GetActiveField(void) const { return activeField; }

// maximum diffusion coefficient in cm^2/sec (anisotropic must override) (diff in mm^2/sec)
double MaterialBase::MaximumDiffusion(void) const { return diffusionCon/100.; }

// maximum diffusivity in cm^2/sec  (anisotropic must override)
// specific ks is mJ mm^2/(sec-K-g) and Cp is mJ/(g-K) so ks/Cp = mm^2 / sec * 1e-2 = cm^2/sec
double MaterialBase::MaximumDiffusivity(void) const { return 0.01*kCond/heatCapacity; }

// material-to-material contact
void MaterialBase::SetFriction(double friction,int matID,double Dn,double Dnc,double Dt)
{	
	if(lastFriction==NULL)
    {   // if this material did not have an friction settings, create one now
        // and tell the new object it is the only one (i.e., it's next is NULL)
		lastFriction=(ContactDetails *)malloc(sizeof(ContactDetails));
		lastFriction->nextFriction=NULL;
	}
	else
    {   // if this material already has a friction setting, create a new one,
        // set it to point to the prior one (in lastFriction), and then set
        // this material's lastFriction to this latest one
		ContactDetails *newFriction=(ContactDetails *)malloc(sizeof(ContactDetails));
		newFriction->nextFriction=(char *)lastFriction;
		lastFriction=newFriction;
	}
    // the law type is set later in MaterialBase::ContactOutput()
	lastFriction->friction=friction;
	lastFriction->matID=matID;
	lastFriction->Dn=Dn;
	lastFriction->Dnc=Dnc;
	lastFriction->Dt=Dt;
}

// Look for contact to a given material and return contact details
// by checking on friction settings for this material. If found return
// the ContactDetails, otherwise return NULL
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
	
    // Law determined by friction coefficient
    // <= -10 : Revert to single velocity field
    // -1 to -9 (actually -10 < friction < -.5 : stick conditions in contact, free in separation
    // 0 : frictionless
    // >=10 : imperfect interface
    // otherwise : friction with that coefficient of friction
	while(currentFriction!=NULL)
	{	// Custom Contact
		if(currentFriction->friction<=-10.)
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
			sprintf(hline,"imperfect interface\n         Dn = %g MPa/mm, Dnc = %g MPa/mm, Dt = %g MPa/mm",
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

// create buffer for material data with the requested number of
// doubles and set eash one to zero
double *MaterialBase::CreateAndZeroDoubles(int numDoubles) const
{
	// create buffer
	double *p=new double[numDoubles];
	
	// set all to zero
	int i;
	for(i=0;i<numDoubles;i++) p[i] = 0.0;
	
	// return pointer
	return p;
}


#pragma mark MaterialBase::Methods

// To handle elimination of old MPMConstLaw, this passes on to old one
// unless subclass overrides to use it directly
void MaterialBase::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 dv,double delTime,int np,void *properties,ResidualStrains *res) const
{   if(np==THREED_MPM)
    {   MPMConstLaw(mptr,dv(0,0),dv(1,1),dv(2,2),dv(0,1),dv(1,0),
                            dv(0,2),dv(2,0),dv(1,2),dv(2,1),delTime,np,properties,res);
    }
    else
    {   MPMConstLaw(mptr,dv(0,0),dv(1,1),dv(0,1),dv(1,0),dv(2,2),delTime,np,properties,res);
    }
}

// These methods are now deprecated. All new material should use the single matrix call (and check
//   np to see if 2D or 3D)
void MaterialBase::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
                                double dvzz,double delTime,int np,void *properties,ResidualStrains *res) const
{
}
void MaterialBase::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
                                double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np,void *properties,ResidualStrains *res) const
{
}

// buffer size for mechanical properties
int MaterialBase::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = 0;
    return 0;
}

// Get copy of properties that depend on material state
void *MaterialBase::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{	return NULL;
}

// Get transport property tensors (if change with particle state)
void MaterialBase::GetTransportProps(MPMBase *mptr,int np,TransportProperties *t) const { *t = tr; }

// Get Cv heat capacity
// Implemented in case heat capacity changes with particle state
// Units mJ/(g-K) = J/(kg-m)
double MaterialBase::GetHeatCapacity(MPMBase *mptr) const { return heatCapacity; }

// For Cp heat capacity
double MaterialBase::GetCpHeatCapacity(MPMBase *mptr) const { return GetHeatCapacity(mptr)+GetCpMinusCv(mptr); }

// A material can override to set Cp-Cv in mJ/(g-K) = J/(kg-m)
// From thermodyanamics Cp-Cv = (3K CTE ^2T/rho) where CTE is linear CTE
// (if mptr==NULL, can use stress free temperature instead)
double MaterialBase::GetCpMinusCv(MPMBase *mptr) const { return 0; }

// Increment heat energy using Cv(dT-dTq0) - dPhi, where Cv*dTq0 + dPhi = Cv dTad is total
//		dispated energy (it can by provided as either a temperature rise or an energy)
// dTq0 is temperature rise due to material mechanisms if the process was adiabatic
// dPhi is dissipated energy that is converted to temperature rise
void MaterialBase::IncrementHeatEnergy(MPMBase *mptr,double dT,double dTq0,double dPhi) const
{
	double Cv = 1000.*GetHeatCapacity(mptr);
	double dispEnergy = Cv*dTq0 + dPhi;                     // = Cv dTad
	
	// Isolated means no conduction and no thermal ramp (and in future if have other ways
	//		to change particle temperature, those are not active either)
	// In this mode, adiabatic has dq=0 and isothermal releases all as heat
    if(isolatedSystemAndParticles)
    {   // Here dText = 0
        // If adiabatic, dq = 0 (nothing to add)
        // If isothermal dq = -Cv dTad
		if(!ConductionTask::adiabatic)
        {   mptr->AddHeatEnergy(-dispEnergy);
            mptr->AddEntropy(-dispEnergy/mptr->pPreviousTemperature);
        }
    }
	else
	{	// For non isolated particle use dq = Cv(dT-dTad)
        // If adiabatic, the Cv dT term in next step will include Cv dTad from
        //    this step to give particle dq = 0. If isothermal, it will not and
        //    dq will be disspated heat
        double totalHeat = Cv*dT - dispEnergy;
		mptr->AddHeatEnergy(totalHeat);
        if(ConductionTask::adiabatic)
        {   double dTad = dispEnergy/Cv;
            mptr->AddEntropy(Cv*dT/mptr->pPreviousTemperature - dispEnergy/(mptr->pPreviousTemperature+dTad));
        }
        else
        {   mptr->AddEntropy(totalHeat/mptr->pPreviousTemperature);
        }
	}
    
	// The dispated energy is added here, but only if adiabatic, in which case it will be
    // converted to particle temperature rise later in the calculations. This temperature
    // change works with conduction on or off
    if(ConductionTask::adiabatic)
        mptr->AddDispEnergy(dispEnergy);
}

// Correct stress update for rotations using hypoelasticity approach
// rotD is -wxy of rotation tensor (or minus twice tensorial rotation, i.e. dvxy-dvyx = -(dvyx-dvxy))
// dsxx, dsyy, and dtxy are unrotated updates to sxx, syy, and txy
// This methods increments particle angle, rotational strain, and in-plane stresses
// This update by a midpoint rule update is documented in Mathematica notes (HypoUpdate.nb)
void MaterialBase::Hypo2DCalculations(MPMBase *mptr,double rotD,double dsxx,double dsyy,double dtxy) const
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
void MaterialBase::Hypo3DCalculations(MPMBase *mptr,double dwxy,double dwxz,double dwyz,double *dsig) const
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

// Convert J to K assuming an isotropic material
// d -- crack opening displacement near crack tip, d.y--opening, d.x--shear
// C -- crack propagating velocity
// J0 -- J integral components in J0.x and J0.y
// np -- PLANE_STRESS_MPM or PLANE_STRAIN_MPM (asysymtry not certain, current reverts to plane strain)
// nuLS and GLS -- low strain Poisson's ratio and shear modulus (in MPa)
Vector MaterialBase::IsotropicJToK(Vector d,Vector C,Vector J0,int np,double nuLS,double GLS)
{
    double Cs2,Cd2,C2;
    double B1,B2,A1,A2,A3,A4,DC;
    double term1,term2;
    Vector SIF;
	
    double dx = d.x;
    double dy = d.y;
    double J0x = fabs(J0.x);                        // J0.x should be always positive
    double J0y = J0.y;
	
	double kf=0.;
    if(np==PLANE_STRESS_MPM)
        kf=(3.-nuLS)/(1.+nuLS);
	else
        kf=3.-4.*nuLS;
	
    C2 = C.x*C.x+C.y*C.y;				// square of crack velocity
	// dynamic or stationary crack
    if(!DbleEqual(sqrt(C2),0.0)) 
	{	Cs2=1.e+3*GLS/rho;				// now in m^2/sec^2
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
    else
	{	B1=B2=1.0;
        A3=1.;
        A1=A2=A4=(kf+1.)/4.;
    }
	
    term2=dy*dy*B2+dx*dx*B1;
	// special case for zero COD
    if(DbleEqual(term2,0.0))
	{	SIF.x = 0.0;
		SIF.y = 0.0;
    }
    else
	{	SIF.x = dy*sqrt(2*GLS*J0x*B2/A1/term2);
		SIF.y = dx*sqrt(2*GLS*J0x*B1/A2/term2);
    }
    
    return SIF;
}

// Determine what calculations are needed for the propagation criterion
// in this material - must match needs in ShouldPropagate() routine
int MaterialBase::CriterionNeeds(int critIndex,bool &usesEnergyBalance)
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
            usesEnergyBalance = TRUE;
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
                        crkTip->speed=initSpeed*WaveSpeed(FALSE,NULL);		// fraction of wave speed in 2D
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
                    {	adjustSpeed=fmin(adjustSpeed,WaveSpeed(FALSE,NULL)-crkTip->speed);
                        crkTip->steadyState=PROPAGATING;
                    }
                    
                    // adjust speed, but never below minimum
                    crkTip->speed+=adjustSpeed;
                    if(crkTip->speed<=0.5*initSpeed*WaveSpeed(FALSE,NULL))
                    {	crkTip->speed=0.5*initSpeed*WaveSpeed(FALSE,NULL);
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

/* find hoop direction unit vector relative to crack direction. 
	If theta is the ccw rotation of hoop direction relative to crack tip direction, then
		hooopDir = (cos(theta), sin(theta))
	To get vector in hoop direction later use;
		(crackDir.x*hoopDir.x - crackDir.y*hoopDir.y, crackDir.x*hoopDir.y + crackDir.y*hoopDir.x)
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
// not thread safe due to push_back()
double MaterialBase::CrackPropagationAngleFromStrainEnergyDensityCriterion(double k,
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

/* Calculate maximum wave speed for material in mm/sec. WaveSpeed() is called
	once for each material point at beginning of calculation. If variable wave
	speed, be conservative and return the maximum possible save speed.
	It is also called by crack propagation, for which threeD is always FALSE
		and mptr is NULL
	It is also called by silent boundary conditions if ShearWaveSpeed() is not
		overridden. These boundary conditions only work (and not even sure of that)
		for isotropic materials.
	The method is abstract (in MaterialBase.hpp) so all sub classes must implement
*/
//double NewMaterial::WaveSpeed(bool threeD,MPMBase *mptr) { }

/* Calculate shear wave speed for material in mm/sec. This is only called for silent
	boundary conditions. This base class return WaveSpeed()/sqrt(3). A new
	material only needs to override this method if it will implement silent boundary
	conditions and if the base class method is not correct.
*/
double MaterialBase::ShearWaveSpeed(bool threeD,MPMBase *mptr) const { return WaveSpeed(threeD,mptr)/sqrt(3.); }

/* This method is only called by the custom task to adjust the wave speed depending on
    particle stress state. The default implementation calls the WaveSpeed() method used
    at the start of the time step. If wave speed depends on particle state in a known way,
    this method can be overridden to get the current wave speed. The custom task to
    adjust time step will only work for materials that override this method for
    particle-dependent results.
*/
double MaterialBase::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const { return WaveSpeed(threeD,mptr); }

// archive material data for this material type when requested.
double MaterialBase::GetHistory(int num,char *historyPtr) const { return (double)0.; }

// return TRUE if rigid particle (for contact or for BC)
bool MaterialBase::Rigid(void) const { return FALSE; }

// return TRUE is rigid BC particle (not rigid for contact)
short MaterialBase::RigidBC(void) const { return FALSE; }

// return TRUE is rigid BC particle (not rigid for contact)
short MaterialBase::RigidContact(void) const { return FALSE; }

// check if traciton law material
bool MaterialBase::isTractionLaw(void) const { return FALSE; }

// check if keeps crack tip
int MaterialBase::KeepsCrackTip(void) const { return constantTip; }

// see if material for a material velocity field is rigid (only rigid contact materials can be in a velocity field)
short MaterialBase::GetMVFIsRigid(int matfld)
{	return matfld<(int)fieldMatIDs.size() ? theMaterials[fieldMatIDs[matfld]]->Rigid() : FALSE ;
}

// convert field number (zero based) to material ID for that field (zero based)
int MaterialBase::GetFieldMatID(int matfld) { return fieldMatIDs[matfld]; }

// convert active material number (zero based) to material ID for that field (zero based)
int MaterialBase::GetActiveMatID(int matfld) { return activeMatIDs[matfld]; }

// Get current relative volume change - only used to convert speific results to actual values when archiving
// Materials with explicit treatment of large deformation will need it (e.g., Hyperelastic)
double MaterialBase::GetCurrentRelativeVolume(MPMBase *mptr) const { return 1.; }

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor MaterialBase::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    return stress;
}

// Calculate artficial damping where Dkk is the relative volume change rate = (V(k+1)-V(k))/(V(k+1)dt)
// and c is the current wave speed in m/sec
double MaterialBase::GetArtificalViscosity(double Dkk,double c) const
{
    double avred = 0.;
    if(Dkk<0 && artificialViscosity)
    {   double divuij = fabs(Dkk);                                  // sec^-1
        double dcell = mpmgrid.GetAverageCellSize();                // mm
        avred = dcell*divuij*(avA1*c + 1.e-3*avA2*dcell*divuij);    // Pa cm^3/g
    }
    return avred;
}

// If material partitions total strain into elastic and plastic strain saved in ep and eplast, it'
// should override this method and return TRUE. It is only used when material point is asked for
// its GetDeformationGradient(). When this is TRUE, gradient uses total strain, otherwise it
// uses only the terms in ep and wrot.
bool MaterialBase::PartitionsElasticAndPlasticStrain(void) { return FALSE; }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool MaterialBase::SupportsArtificialViscosity(void) const { return FALSE; }





