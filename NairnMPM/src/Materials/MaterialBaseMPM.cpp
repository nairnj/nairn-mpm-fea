/********************************************************************************
    MaterialBaseMPM.cpp - more MaterialBase for MPM code
    nairn-mpm-fea
    
    Created by John Nairn on Mar 17 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
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
#include "System/UnitsController.hpp"
#include "Materials/ContactLaw.hpp"
#include <vector>

// global
bool MaterialBase::isolatedSystemAndParticles = FALSE;

// class statics for MPM - zero based material IDs when in multimaterial mode
vector<int> MaterialBase::fieldMatIDs;					// material ID in velocity field [i]
vector<int> MaterialBase::activeMatIDs;					// list of active material velocity fields
int MaterialBase::incrementalDefGradTerms = 2;			// terms in exponential of deformation gradient
int MaterialBase::maxPropertyBufferSize = 0;            // maximum buffer size needed among active materials to get copy of mechanical properties
int MaterialBase::maxAltBufferSize = 0;                 // maximum optional buffer size needed for more properties (e.g., hardenling law)
bool MaterialBase::extrapolateRigidBCs = false;			// rigid BCs extrapolated (new) or projected (old)

#pragma mark MaterialBase::Initialization (required)

// Read material properties common to all MPM materials
char *MaterialBase::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"rho")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&rho,gScaling,0.001);
    }
    
    else if(strcmp(xName,"KIc")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&KIc,gScaling,31.62277660168379e6);
    }
    
    else if(strcmp(xName,"KIIc")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&KIIc,gScaling,31.62277660168379e6);
    }
	
	// crit 2 only
    else if(strcmp(xName,"maxLength")==0)
    {	input=DOUBLE_NUM;
        return((char *)&maxLength);
    }
    
	// initTime in crit 2
    else if(strcmp(xName,"initTime")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&initTime,gScaling,1.e-3);
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
		return UnitsController::ScaledPtr((char *)&JIc,gScaling,1000.);
    }
    
    else if(strcmp(xName,"JIIc")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&JIIc,gScaling,1000.);
    }
	
    else if(strcmp(xName,"nmix")==0)
    {	input=DOUBLE_NUM;
        return((char *)&nmix);
    }
	
	// speed in crit 2 (scale Legacy m/sec to mm/sec)
    else if(strcmp(xName,"speed")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&initSpeed,gScaling,1.e3);
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
    
	else if(strcmp(xName,"shareMatField")==0)
	{	input=INT_NUM;
		return((char *)&shareMatField);
	}
    
	else if(strcmp(xName,"Cp")==0 || strcmp(xName,"Cv")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&heatCapacity,gScaling,1.e6);
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
		return UnitsController::ScaledPtr((char *)&kCond,gScaling,1.e6);
	}
    
	// check properties only for some materials
	if(SupportsArtificialViscosity())
	{	if(strcmp(xName,"ArtificialVisc")==0)
		{	artificialViscosity = true;
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
    
    return (char *)NULL;
}

/*	Verify and calculate properties used in analyses. If error, return string with an error message.
 This is called once at start of the calculation just before the material properties
 are printined to the output file and before any calculations. It is called for
 every material defined in the input file, even if it is not used by any
 particle.
 If superclass overrides this method, must call this too
 */
const char *MaterialBase::VerifyAndLoadProperties(int np)
{
	// make conductivity specific (nJ mm^2/(sec-K-g))
	kCond /= rho;
	
	// in case only need to load some things once, load those mechanical properties now
	FillTransportProperties(&tr);
	
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

#pragma mark MaterialBase::Initialization (optional)

/* This is called after PreliminaryCalcs() and just before first MPM time step and it
 is only called if the material is actually in use by one or more particles
	If material cannot be used in current analysis type throw an exception
	Subclass that overrides must pass on to super class
	throws CommonException()
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
			if(theMaterials[tractionMat[i]-1]->MaterialStyle()!=TRACTION_MAT)
			{	throw CommonException("Material with propagation material that is not a traction law",
									  "MaterialBase::ValidateForUse");
			}
		}
	}
	
	// check for unsupported alternate propagation criterion
	if(criterion[1]==STEADYSTATEGROWTH)
	{	throw CommonException("The alternate propagation criterion cannot be steady state crack growth.",
							  "MaterialBase::ValidateForUse");
	}
	
	if(ConductionTask::active || ConductionTask::adiabatic)
	{	if(heatCapacity<=0. && !IsRigid())
		{	throw CommonException("Thermal conduction and/or mechanical energy cannot be done using materials that have zero heat capacity.",
								  "MaterialBase::ValidateForUse");
		}
	}
	
	// Check diffusion options, but don't ask if supports diffusion unless activated
	if(fmobj->HasDiffusion())
	{	if(!SupportsDiffusion())
		{	throw CommonException("Diffusion activated for a material that does not support it.",
								  "MaterialBase::ValidateForUse");
		}
	}
}

// Called before analysis, material can fill in things that never change during the analysis
// Note: no angle, because cannot depend on material angle
// Here fills in isotropic properties, materials with different anisotropic properties should override
void MaterialBase::FillTransportProperties(TransportProperties *t)
{
	// diffusion tensor
	t->diffusionTensor.xx=diffusionCon;
	t->diffusionTensor.yy=diffusionCon;
	t->diffusionTensor.zz=diffusionCon;
	t->diffusionTensor.xy=0.;
	t->diffusionTensor.xz=0.;
	t->diffusionTensor.yz=0.;
	
	// conductivity tensor normalized to rho
	t->kCondTensor.xx = kCond;
	t->kCondTensor.yy = kCond;
	t->kCondTensor.zz = kCond;
	t->kCondTensor.xy = 0.;
	t->kCondTensor.xz = 0.;
	t->kCondTensor.yz = 0.;
}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
// Such a class must pass on the super class after its own initializations
// The offset is to allow relocatable history data in child materials (e.g., phases)
void MaterialBase::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
	if(isolatedSystemAndParticles)
	{   // need to add initial heat energy, because special cases in this mode
		// will ignore the Cv dT term
		double Cv = GetHeatCapacity(mptr);            // in nJ/(g-K)
		mptr->AddHeatEnergy(Cv*(mptr->pTemperature-thermal.reference));
		// initial entropy kept to zero, could add something here if ever needed
	}
}

#pragma mark MaterialBase::Printing Properties

// print any properties common to all MPM material types
void MaterialBase::PrintCommonProperties(void) const
{
	if(MaterialStyle()==TRACTION_MAT) return;
	
	// rigid only allow cracks option
	if(IsRigid())
	{	if(AllowsCracks())
		cout << "Sees cracks" << endl;
	else
		cout << "Ignores cracks" << endl;
		return;
	}
	
	// print density
	PrintProperty("rho",rho*UnitsController::Scaling(1000.),"");
	cout << endl;
	
	if(AllowsCracks())
	{	// print growth criterion and relevant material properties for crack growth
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
	}
	else
		cout << "Ignores cracks" << endl;
	
	// artificial visconsity
	if(artificialViscosity)
	{	PrintProperty("Artificial viscosity on",FALSE);
		PrintProperty("AV-A1",avA1,"");
		PrintProperty("AV-A2",avA2,"");
		cout << endl;
		if(ConductionTask::AVHeating)
		{	PrintProperty("    AV heating on",FALSE);
			cout << endl;
		}
		else
		{	PrintProperty("    AV heating off",FALSE);
			cout << endl;
		}
	}
	
	// particle damping
	if(matUsePDamping)
	{	cout << "Particle damping: " << matPdamping << " /sec";
		if(matUsePICDamping)
		{	cout << ", PIC damping fraction: " << matFractionPIC;
		}
		cout << endl;
	}
	else if(matUsePICDamping)
	{	cout << "Particle PIC damping fraction: " << matFractionPIC;
		cout << endl;
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
			PrintProperty("KIc",KIc*UnitsController::Scaling(31.62277660168379e-9),UnitsController::Label(STRESSINTENSITY_UNITS));
			cout << endl;
			break;
			
		case CRITICALERR:
			cout << "Critical Energy Release Rate" << PreferredDirection(thisDirection) << endl;
			PrintProperty("Jc",JIc*UnitsController::Scaling(0.001),UnitsController::Label(ERR_UNITS));
			cout << endl;
			break;
			
		case STEADYSTATEGROWTH:
			cout << "Constant crack speed" << PreferredDirection(thisDirection) << endl;
			if(initTime<0.)
			{	PrintProperty("Jc",JIc*UnitsController::Scaling(0.001),UnitsController::Label(ERR_UNITS));
				PrintProperty("initSpeed",initSpeed*UnitsController::Scaling(0.001),"m/s");
			}
			else
			{	PrintProperty("ti",initTime*UnitsController::Scaling(1.e3),"ms");
				PrintProperty("initSpeed",initSpeed*UnitsController::Scaling(0.001),"m/s");
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
			
		case STRAINENERGYDENSITY:
			cout << "Minmum strain energy density" << PreferredDirection(thisDirection) << endl;
			PrintProperty("KIc",KIc*UnitsController::Scaling(31.62277660168379e-9),UnitsController::Label(STRESSINTENSITY_UNITS));
			cout << endl;
			break;
			
		case EMPIRICALCRITERION:
			cout << "Empirical criterion" << PreferredDirection(thisDirection) << endl;
			sprintf(mline,"KIc=%12.3f %s  KIIc=%12.3f %s KIexp=%12.3f KIIexp=%12.3f",
					KIc*UnitsController::Scaling(31.62277660168379e-9),UnitsController::Label(STRESSINTENSITY_UNITS),
					KIIc*UnitsController::Scaling(31.62277660168379e-9),UnitsController::Label(STRESSINTENSITY_UNITS),
					KIexp,KIIexp);
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

// print transport properties to output window - default is isotropic properties
// aniostropic materials must override it
// This is called just after VerifyAndLoadProperties
void MaterialBase::PrintTransportProperties(void) const
{
	if(IsRigid()) return;
	
	// Diffusion constants
	if(fmobj->HasDiffusion())
	{	PrintProperty("D",diffusionCon,UnitsController::Label(DIFFUSION_UNITS));
		PrintProperty("csat",concSaturation,"");
		cout << endl;
		PrintProperty("b",betaI,"1/wt fr");
		cout << endl;
	}
	
	// Conductivity constants
	if(ConductionTask::active)
	{	PrintProperty("k",rho*kCond*UnitsController::Scaling(1.e-6),UnitsController::Label(CONDUCTIVITY_UNITS));
		PrintProperty("Cv",heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		PrintProperty("Cp",(heatCapacity+GetCpMinusCv(NULL))*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		cout << endl;
	}
	else if(ConductionTask::adiabatic)
	{	PrintProperty("Cv",heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		// Cp only used in conduction so not printed here when conduction is off
		cout << endl;
	}
}

#pragma mark MaterialBase::Setting Internal Laws

// Material that allow hardening laws must accept
// the final call and use if supported
// Newly created laws should be added here
// throws std::bad_alloc, SAXException()
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
	{	ThrowSAXException("The hardening law '%s' is not valid",lawName);
		return;
	}
    
    // pass on to the material
    if(!AcceptHardeningLaw(pLaw,lawID))
	{	delete pLaw;
		ThrowSAXException("The hardening law '%s' is not allowed by material",lawName);
		return;
	}
}

// A Material that allows hardening laws should accept this one or it can
// veto the choice by returning FALSE
bool MaterialBase::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID) { return false; }

#pragma mark MaterialBase:velocities fields

// when set, return total number of materials if this is a new one, or 1 if not in multimaterial mode
// not thread safe due to push_back()
// throws CommonException()
int MaterialBase::SetField(int fieldNum,bool multiMaterials,int matid,int &activeNum)
{	if(!multiMaterials)
	{	if(!AllowsCracks())
		{	// Cannot ignore cracks in single material mode. Rather than an error, just
			// reset to see cracks here
			allowsCracks = true;
		}
		
		// for first particle using this material, add to active material IDs and check required
		// material buffer sizes
		if(field<0)
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
		{	// if sharing a field, look up shared material.
			if(shareMatField>0)
			{	if(shareMatField>nmat)
				{	throw CommonException("Material class trying to share velocity field with an undefined material type",
										  "MaterialBase::SetField");
				}
				
				// must match for rigid, allowing cracks, and membrane features
				MaterialBase *matRef = theMaterials[shareMatField-1];
				if(matRef->IsRigid() != IsRigid())
				{	throw CommonException("Material class trying to share velocity field with an incompatible material type",
										  "MaterialBase::SetField");
				}
				if(matRef->AllowsCracks() != AllowsCracks())
				{	throw CommonException("Material class trying to share velocity field with an incompatible material type",
										  "MaterialBase::SetField");
				}
				if(matRef->MaterialStyle() != MaterialStyle())
				{	throw CommonException("Material class trying to share velocity field with an incompatible material type",
										  "MaterialBase::SetField");
				}
				
				// base material cannot share too
				if(matRef->GetShareMatField()>=0)
				{	throw CommonException("Material class trying to share velocity field with a material that share's its field",
										  "MaterialBase::SetField");
				}

				// set field to other material (and set other material if needed
				field = matRef->GetField();
				if(field<0)
				{	fieldNum = matRef->SetField(fieldNum,multiMaterials,shareMatField-1,activeNum);
					field = fieldNum-1;
				}
			}
			else
			{	field=fieldNum;
				fieldNum++;
				
				// fieldMatIDs[i] for i=0 to # materials is material index for that material velocity field
				// when materials share fields, it points to the based shared material
				fieldMatIDs.push_back(matid);
			}
			
			// for first particle using this material, add to active material IDs and check required
			// material buffer sizes
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

// material ID (zero based) for field to share its velocity field (-1 if not sharing or not in multimaterial mode)
int MaterialBase::GetShareMatField(void) const { return shareMatField-1; }

// -1 if material not in use, otherwise zero-based field number from 0 to numActiveMaterials
int MaterialBase::GetActiveField(void) const { return activeField; }

#pragma mark MaterialBase:Contact and damping

// material-to-material contact
void MaterialBase::SetFriction(int lawID,int matID)
{	
	if(lastFriction==NULL)
    {   // if this material did not have any friction settings, create one now
        // and tell the new object it is the only one (i.e., it's next is NULL)
		lastFriction = new ContactPair;
		lastFriction->nextFriction = NULL;
	}
	else
    {   // if this material already has a friction setting, create a new one,
        // set it to point to the prior one (in lastFriction), and then set
        // this material's lastFriction to this latest one
		ContactPair *newFriction = new ContactPair;
		newFriction->nextFriction = (char *)lastFriction;
		lastFriction = newFriction;
	}
    // the law type is set later in MaterialBase::ContactOutput()
	lastFriction->lawID = lawID;
	lastFriction->matID = matID;
}

// Custom material damping.
// Note that matPIC will be -1 unless this material changed it with PIC attribute on the
//		the Damping command. Thus particle can do FLIP my matPIC=0 resulting
//		in matFractionPIC=0 too.
void MaterialBase::SetDamping(double matDamping,double matPIC)
{	if(matDamping>-1e12)
	{	matPdamping = matDamping;
		matUsePDamping = true;
	}
	else
		matUsePDamping = false;
    if(matPIC>=0.)
	{   matFractionPIC = matPIC;
        matUsePICDamping = true;
    }
    else
        matUsePICDamping = false;
}

// Change damping if this material requests it
// On input, particleAlpha and gridAlpha should be the global settings
// if material has particle damping (matUsePDamping)
//	 Changes PIC fraction (including to zero) (matUsePICDamping), change to
//      particleAlpha   = matFractionPIC/dt + matPdamping(t)
//      gridAlpha       = -m*matFractionPIC/dt + damping(t)
//		localXPIC		= m*matFractionPIC/dt
//   Uses global PIC fraction, change to
//      particleAlpha   = globalPIC + matPdamping(t)
//		localXPIC		= m*globalPIC
// if material has no damping, but wants to change PIC fraction (including to zero) (matUsePICDamping)
//      particleAlpha   = matFractionPIC/dt + pdamping(t) = matFractionPIC/dt + (particleAlpha - globalPIC)
//      gridAlpha       = -m*matFractionPIC/dt + damping(t)
//		localXPIC		= m*matFractionPIC/dt
// if material has no damping return
//		localXPIC		= m*globalPIC
double MaterialBase::GetMaterialDamping(double &particleAlpha,double &gridAlpha,double nonPICGridAlpha,double globalPIC) const
{
	double localXPIC = 0.;
	
	// has either type
	if(matUsePDamping)
	{	if(matUsePICDamping)
		{	particleAlpha = matFractionPIC/timestep + matPdamping;
			gridAlpha = -matFractionPIC/timestep + nonPICGridAlpha;
		}
		else
		{   particleAlpha = globalPIC + matPdamping;
		}
	}
	
	// if has particle PIC only
	else if(matUsePICDamping)
	{	particleAlpha += matFractionPIC/timestep - globalPIC;
		gridAlpha = -matFractionPIC/timestep + nonPICGridAlpha;
	}
	
	return localXPIC;
}

// Look for contact to a given material and return contact law ID
// by checking on friction settings for this material.
int MaterialBase::GetContactToMaterial(int otherMatID)
{	ContactPair *currentFriction=lastFriction;
	while(currentFriction!=NULL)
	{	if(currentFriction->matID==otherMatID) return currentFriction->lawID;
		currentFriction=(ContactPair *)currentFriction->nextFriction;
	}
	return -1;
}	

// Print contact law settings for cracks and finalize variables
void MaterialBase::ContactOutput(int thisMatID)
{
	ContactPair *currentFriction=lastFriction;
	
	while(currentFriction!=NULL)
	{	cout << "Custom contact between material " << thisMatID << " and material "
				<< currentFriction->matID << ": material "
				<< (MaterialBase::GetContactLawNum(currentFriction->lawID)+1)
				<< endl;
		currentFriction=(ContactPair *)currentFriction->nextFriction;
	}
}

// Convert readID (1 based) to material number for contact law (0 based)
int MaterialBase::GetContactLawNum(int readID)
{
	// check if already a real ID
	if(readID<1000) return readID-1;
	
	// look for other ID in defined amterials
	int clmat=nmat-1;
	while(clmat>=0)
	{	ContactLaw *matLaw = (ContactLaw *)theMaterials[clmat];
		if(matLaw->autoID == readID) return clmat;
		clmat--;
	}
	return -1;
}

#pragma mark MaterialBase::History Data Methods

// create and return pointer to material-specific data on a particle
//	using this material. Called once at start of analysis for each
//	particle.
// If a material subclass overrides this method, it is responsible for
//  creating space for super class history variables too. There is
//  no mechanism now for calling super classes to prepend or
//  append the data.
// First is material points and second is for traction laws
char *MaterialBase::InitHistoryData(char *pchr,MPMBase *mptr) { return NULL; }
char *MaterialBase::InitHistoryData(char *pchr) { return NULL; }

// Number of history variables (used by some materials when creating array of doubles)
int MaterialBase::NumberOfHistoryDoubles(void) const { return 0; }

// archive material data for this material type when requested
// assumes materials using simple array of doubles
double MaterialBase::GetHistory(int num,char *historyPtr) const
{	double history=0.;
    if(num>0 && num<=NumberOfHistoryDoubles())
    {	double *soft = (double *)historyPtr;
        history = soft[num-1];
    }
    return history;
}

// return number of bytes needed for history data
// a negative number means this material does not support combined history
//	data that might be offset from particle history pointer
// Only used by Phase Transition Material
int MaterialBase::SizeOfHistoryData(void) const { return -1; }

// if pchr==NULL, create buffer for material data with the requested number of double
//		otherwise assume it already exists and is big enough
// set each double in the history data to zero
// throws std::bad_alloc
double *MaterialBase::CreateAndZeroDoubles(char *pchr,int numDoubles) const
{
	// create buffer (if needed)
	if(pchr==NULL)
	{	// allocate history data space
		int historySize = numDoubles*sizeof(double);
		pchr = new char[historySize];
	}
	
	// cast to double *
	double *p = (double *)pchr;
	
	// set all to zero
	int i;
	for(i=0;i<numDoubles;i++) p[i] = 0.0;
	
	// return pointer
	return p;
}

#pragma mark MaterialBase::Methods

// All maternal classes must override to handle their constitutive law
void MaterialBase::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 dv,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
}

// Using small-strain, large rotation method to find incremental strain
// Tasks
// 1. find dF, increment and save new deformation on the particle
// 2. Find incremental rotation tensor
// R0 must be provided for MATERIAL_CONFIGURATION, but otherwise can be NULL
// If Rnm1 or Rn are not NULL, that rotation matrix is returned (rotation to
//		initial, but not to material axes)
Matrix3 MaterialBase::LRGetStrainIncrement(int axes,MPMBase *mptr,Matrix3 du,Matrix3 *dR,Matrix3 *R0,
                                           Matrix3 *Rnm1out,Matrix3 *Rnout) const
{
    // current previous deformation gradient and stretch
    Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
    
    // get incremental deformation gradient and decompose it
    const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    
    // Update total deformation gradient
    Matrix3 pF = dF*pFnm1;
    mptr->SetDeformationGradientMatrix(pF);
    
    // decompose to get previous Rn and Rn-1
    Matrix3 Rnm1,Rn;
    Matrix3 *Rnm1Ptr = (Rnm1out==NULL) ? &Rnm1 : Rnm1out ;
    Matrix3 *RnPtr = (Rnout==NULL) ? &Rn : Rnout ;
    
    // two decompositions
    pFnm1.RightDecompose(Rnm1Ptr,NULL);
    pF.LeftDecompose(RnPtr,NULL);
    *dR = (*RnPtr)*Rnm1Ptr->Transpose();
    
    // get strain increments in
    // 1. current configuration (dF-dR)F(n-1) Rn^T
    // 2. initial configuration Rn^T (dF-dR)F(n-1)
	// 3. material configuration R0^T (initial) R0
    Matrix3 dFmdR = dF - *dR;
    Matrix3 de;
    switch(axes)
    {	case CURRENT_CONFIGURATION:
            de = dFmdR*(pFnm1*RnPtr->Transpose());
            break;
        case INITIAL_CONFIGURATION:
            de = RnPtr->Transpose()*(dFmdR*pFnm1);
            break;
		default:
            de = RnPtr->Transpose()*(dFmdR*pFnm1);
            de = de.RTMR(*R0);
            break;
    }
    
    return de;
}

// get incremental residual voluje change (needed for phase change materials
double MaterialBase::GetIncrementalResJ(MPMBase *mptr,ResidualStrains *res) const { return 0.; }

// buffer size for mechanical properties
int MaterialBase::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = 0;
    return 0;
}

// Get copy of properties that depend on material state
void *MaterialBase::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{	return NULL;
}

// Get transport property tensors (if change with particle state)
// only called when there are transport tasks
void MaterialBase::GetTransportProps(MPMBase *mptr,int np,TransportProperties *t) const
{	*t = tr;
}

// Get Cv heat capacity (with option for parent to add excess heat capacity)
// Implemented in case heat capacity changes with particle state
// Legacy Units nJ/(g-K)
double MaterialBase::GetHeatCapacity(MPMBase *mptr) const
{
    return heatCapacity;
}

// For Cp heat capacity in nJ/(g-K)
double MaterialBase::GetCpHeatCapacity(MPMBase *mptr) const { return GetHeatCapacity(mptr)+GetCpMinusCv(mptr); }

// A material can override to set Cp-Cv in nJ/(g-K)
// From thermodyanamics Cp-Cv = (3K CTE ^2T/rho) where CTE is linear CTE
// (if mptr==NULL, can use stress free temperature instead)
double MaterialBase::GetCpMinusCv(MPMBase *mptr) const { return 0; }

// For eta/Q in poroelasticity, otherwise it is 1
double MaterialBase::GetDiffusionCT(void) const { return diffusionCT; }

// Increment heat energy - call whenever dTq0 or dPhi changes
// dTq0 is temperature rise due to material mechanisms if the process was adiabatic and reversible
// dPhi is dissipated energy that is converted to temperature rise (it is irreverisble)
// If adiabatic: add dTq0+dPhi/CV to adiabatice temperature rise and dPhi/T to entropy
// If isothermal: heat energy decrease by Cv*dTqo+dPhi to mantain constant temperature
//			Entropy change -Cv dTq0/T which is reverisble exchange with exterior
void MaterialBase::IncrementHeatEnergy(MPMBase *mptr,double dTq0,double dPhi) const
{
	double Cv = GetHeatCapacity(mptr);					// in nJ/(g-K)
 
    // Adiabatic temperature change dT_{ad} = dTq0+dPhi/Cv
    if(ConductionTask::adiabatic)
    {   // temperature increased to buffered dT_{ad}
		double dTad = dTq0+dPhi/Cv;
        mptr->Add_dTad(dTad);
		
        // no heat here because adiabatic
        
        // entropy
		mptr->AddEntropy(dPhi,mptr->pPreviousTemperature);
    }
    else
    {   // no temperature change because isothermal
        
        // temperature rise switched to heat
        double baseHeat = -Cv*dTq0;
        mptr->AddHeatEnergy(baseHeat-dPhi);
        
        // for entropy dS = -Cv dTad/T + dPhi/T = - Cv dTq0/T
        mptr->AddEntropy(baseHeat,mptr->pPreviousTemperature);
    }
}

// Correct stress update for rotations using hypoelasticity approach
// dwrotxy of rotation tensor (or twice tensorial rotation, i.e. dvyx-dvxy)
// dsxx, dsyy, and dtxy are unrotated updates to sxx, syy, and txy
// This methods increments in-plane stresses only
void MaterialBase::Hypo2DCalculations(MPMBase *mptr,double dwrotxy,double dsxx,double dsyy,double dtxy) const
{
	Tensor *sp=mptr->GetStressTensor();
	double dnorm = dwrotxy*sp->xy;
	double dshear =  0.5*dwrotxy*(sp->xx-sp->yy);
	sp->xx += dsxx - dnorm;
	sp->yy += dsyy + dnorm;
	sp->xy += dtxy + dshear;
}

// Correct stress update for rotations using hypoelasticity approach
// dwxy, dwxz, and dwyz are the engineering rotational strains (lower diagonal of tensor)
// dsig[6] on initial stress updates in order (xx,yy,zz,yz,xz,xy)
// This methods increments stresses
void MaterialBase::Hypo3DCalculations(MPMBase *mptr,double dwxy,double dwxz,double dwyz,double *dsig) const
{
	// get stress tensor
	Tensor *sp=mptr->GetStressTensor();
	
	// stress increments involving current stress
	Tensor st;
	st.xx =      -dwxy*sp->xy - dwxz*sp->xz;
	st.yy =       dwxy*sp->xy               - dwyz*sp->yz;
	st.zz =                     dwxz*sp->xz + dwyz*sp->yz;
	st.yz = 0.5*(  dwxy*sp->xz              + dwxz*sp->xy          + dwyz*(sp->yy-sp->zz)  );
	st.xz = 0.5*( -dwxy*sp->yz              + dwxz*(sp->xx-sp->zz) + dwyz*sp->xy  );
	st.xy = 0.5*(  dwxy*(sp->xx-sp->yy)     - dwxz*sp->yz          - dwyz*sp->xz );
	
	// now add all
	sp->xx += dsig[XX] + st.xx;
	sp->yy += dsig[YY] + st.yy;
	sp->zz += dsig[ZZ] + st.zz;
	sp->yz += dsig[YZ] + st.yz;
	sp->xz += dsig[XZ] + st.xz;
	sp->xy += dsig[XY] + st.xy;
}

// Damage zero normal vector - damage material will override
Vector MaterialBase::GetDamageNormal(MPMBase *mptr,bool threeD) const
{	Vector dnorm;
	ZeroVector(&dnorm);
	return dnorm;
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
// d -- crack opening displacement near crack tip in mm, d.y--opening, d.x--shear
// C -- crack propagating velocity in mm/sec
// J0 -- J integral components in J0.x and J0.y in uN/mm
// np -- PLANE_STRESS_MPM or PLANE_STRAIN_MPM (axisymmetry not certain, current reverts to plane strain)
// nuLS and GLS -- low strain Poisson's ratio and shear modulus (in Pa = uN/mm^2)
Vector MaterialBase::IsotropicJToK(Vector d,Vector C,Vector J0,int np,double nuLS,double GLS)
{
    //double Cs2,Cd2,C2;
    //double B1,B2,A1,A2,A3,A4,DC;
    //double term1,term2;
    Vector SIF;
	
    double dx = d.x;
    double dy = d.y;
	
	// J0.x should be always positive. If not set K's to zero
	if(J0.x<=0.)
	{	SIF.x = 0.0;
		SIF.y = 0.0;
		return SIF;
    }
    double J0x = J0.x;                        
	
	double kf=0.;
    if(np==PLANE_STRESS_MPM)
        kf=(3.-nuLS)/(1.+nuLS);
	else
        kf=3.-4.*nuLS;
	
	// This code currently ignores crack velocity because it is not calculation
	// The effect in generally small
	// Future could evaulation velocity is use this code
	
	/*
    C2 = C.x*C.x+C.y*C.y;				// square of crack velocity
	// dynamic or stationary crack
    if(!DbleEqual(sqrt(C2),0.0))
	{	Cs2=GLS/rho;
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
	
    term2=dy*dy*B2+dx*dx*B1;			// mm^2
	
	// special case for zero COD
    if(DbleEqual(term2,0.0))
	{	SIF.x = 0.0;
		SIF.y = 0.0;
    }
    else
	{	// Units mm sqrt(uN/mm^2 uN/mm 1/mm^2) = uN/mm^2 sqrt(mm)
		SIF.x = dy*sqrt(2*GLS*J0x*B2/A1/term2);
		SIF.y = dx*sqrt(2*GLS*J0x*B1/A2/term2);
    }
    
	return SIF;
	*/

	
	// cod inmm^2
	double cod2 = dy*dy + dx*dx;
	double Aterm = cod2*(kf+1.)/4.;
	
	// special case for zero COD
	if(DbleEqual(cod2,0.0))
	{	SIF.x = 0.0;
		SIF.y = 0.0;
	}
	else
	{	// Units mm sqrt(uN/mm^2 uN/mm 1/mm^2) = uN/mm^2 sqrt(mm)
		SIF.x = dy*sqrt(2*GLS*J0x/Aterm);			// KI
		SIF.y = dx*sqrt(2*GLS*J0x/Aterm);			// KII
	}
	
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
			return false;
		
		case CRITICALERR:
			return NEED_J;
			
		case MAXCTODCRITERION:
			return false;
		
		case NO_PROPAGATION:
			return false;
		
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
    double KI,KII,fCriterion,cosTheta0,sinTheta0,cosTheta2;
	Vector hoopDir;

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
        
		// Criterion 3 no longer used
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
                return true;
            }
            break;
        
        default:
            break;
    }
    
    return false;
}

#pragma mark MaterialBase::Accessors (optional)

/* Calculate shear wave speed for material in mm/sec. This is only called for silent
	boundary conditions. This base class return CurrentWaveSpeed()/sqrt(3). A new
	material only needs to override this method if it will implement silent boundary
	conditions and if the base class method is not correct.
	offset allows relocatable history data (if needed)
*/
double MaterialBase::ShearWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{	return CurrentWaveSpeed(threeD,mptr,offset)/sqrt(3.);
}

/* This method is called by the custom task to adjust the wave speed depending on
    particle stress state. The default implementation calls the WaveSpeed() method used
    at the start of the time step. If wave speed depends on particle state in a known way,
    this method can be overridden to get the current wave speed. The custom task to
    adjust time step will only work for materials that override this method for
    particle-dependent results.
	It is also called by silent boundary conditions if ShearWaveSpeed() is not
		overridden. These boundary conditions only work (and not even sure of that)
		for isotropic materials.
	offset allows relocatable history data (if needed)
*/
double MaterialBase::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{	return WaveSpeed(threeD,mptr);
}

// maximum diffusion coefficient in mm^2/sec (anisotropic must override) (diff in mm^2/sec)
double MaterialBase::MaximumDiffusion(void) const { return diffusionCon; }

// maximum diffusivity in mm^2/sec  (anisotropic must override)
// specific ks is nJ mm^2/(sec-K-g) and Cp is nJ/(g-K) so ks/Cp = mm^2 / sec
double MaterialBase::MaximumDiffusivity(void) const { return kCond/heatCapacity; }

// Get current relative volume change - only used to convert speific results to actual values when archiving
// Materials with explicit treatment of large deformation will need it (e.g., Hyperelastic)
// offset version if for phase transition materials
double MaterialBase::GetCurrentRelativeVolume(MPMBase *mptr,int offset) const { return 1.; }

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor MaterialBase::GetStress(Tensor *sp,double pressure,MPMBase *) const
{   Tensor stress = *sp;
	return stress;
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor MaterialBase::GetStressPandDev(Tensor *sp,double pressure,MPMBase *mptr) const
{	Tensor stress = *sp;
	stress.xx -= pressure;
	stress.yy -= pressure;
	stress.zz -= pressure;
	return stress;
}

// store a new total stress on a particle's stress and pressure variables
void MaterialBase::SetStress(Tensor *spnew,MPMBase *mptr) const
{	Tensor *sp = mptr->GetStressTensor();
	*sp = *spnew;
}

// store a new total stress on a particle's stress and pressure variables
void MaterialBase::SetStressPandDev(Tensor *spnew,MPMBase *mptr) const
{	double newP = -(spnew->xx+spnew->yy+spnew->zz)/3.;
	Tensor *sp = mptr->GetStressTensor();
	*sp = *spnew;
	sp->xx += newP;
	sp->yy += newP;
	sp->zz += newP;
	mptr->SetPressure(newP);
}

// For generalized plane stress, increment thickness (zz) stress
void MaterialBase::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{	Tensor *sp = mptr->GetStressTensor();
	sp->zz += dszz;
}

// For generalized plane stress, increment thickness (zz) stress through deviatoric stress and pressure
void MaterialBase::IncrementThicknessStressPandDev(double dszz,MPMBase *mptr) const
{	double delP = -dszz/3.;
	Tensor *sp = mptr->GetStressTensor();
	sp->xx += delP;
	sp->yy += delP;
	sp->zz -= 2.*delP;
	double newP = mptr->GetPressure() + delP;
	mptr->SetPressure(newP);
}

// if a subclass material supports artificial viscosity, override this and return true
bool MaterialBase::SupportsArtificialViscosity(void) const { return false; }

// if a material does not support diffusion or poroelasticity, override and return false
// only called when diffusion is active
bool MaterialBase::SupportsDiffusion(void) const { return true; }

// return code indicated what is store in the "plastic strain" on the material, which can
// only be a symmetrix tensor
int MaterialBase::AltStrainContains(void) const { return NOTHING; }

// concentration saturation (mptr cannot be NULL)
double MaterialBase::GetMaterialConcentrationSaturation(MPMBase *mptr) const { return concSaturation; }

// density (NULL is allowed as long as material not using child materials)
// in MPM use a MPMBase *, in FEA always use NULL
double MaterialBase::GetRho(MPMBase *mptr) const { return rho; }

#pragma mark Material Base:Other Accessors

// return true if any rigid particle
bool MaterialBase::IsRigid(void) const { return false; }

// return true is rigid BC particle
bool MaterialBase::IsRigidBC(void) const { return false; }

// return true is rigid contact particle
bool MaterialBase::IsRigidContact(void) const { return false; }

// return true is rigid black particle
bool MaterialBase::IsRigidBlock(void) const { return false; }

// check if membrane material
int MaterialBase::MaterialStyle(void) const { return SOLID_MAT; }

// check if keeps crack tip
int MaterialBase::KeepsCrackTip(void) const { return constantTip; }

// Set flags for this material being rigid and for if it ignores cracks
// Not that ignoring cracks only works for multimaterial mode. Rigid contact
//    materials default to ignoring cracks while non rigid default to see them
int MaterialBase::GetMVFFlags(int matfld)
{	// check number errors
	if(matfld>=(int)fieldMatIDs.size()) return 0;
	
	MaterialBase *matID = theMaterials[fieldMatIDs[matfld]];
	int flags = 0;
	if(matID->IsRigid())
	{	flags = RIGID_FIELD_BIT;
		if(matID->IsRigidBlock()) flags += RIGID_BLOCK_BIT;
	}
	if(!matID->AllowsCracks()) flags += IGORE_CRACKS_BIT;
	return flags;
}

// convert field number (zero based) to material ID for that field (zero based)
// When materials share a field, it is material ID for the base field
int MaterialBase::GetFieldMatID(int matfld) { return fieldMatIDs[matfld]; }

// convert active material number (zero based) to material ID for that field (zero based)
// These are no rigid only
int MaterialBase::GetActiveMatID(int matfld) { return activeMatIDs[matfld]; }

// Calculate artficial damping where Dkk is the relative volume change rate = (V(k+1)-V(k))/(V(k+1)dt)
// and c is the current wave speed in current units
// WARNING - only call this method after verifyting Dkk<0 && artificialViscosity
double MaterialBase::GetArtificalViscosity(double Dkk,double c,MPMBase *mptr) const
{
    double divuij = fabs(Dkk);                                  // T^-1
	double dcell = mpmgrid.GetAverageCellSize(mptr);            // L
	return dcell*divuij*(avA1*c + avA2*dcell*divuij);			// L^2/T^2 = F L/M = Pressure/rho = Kirchoff Pressure/rho0
}

// True to see cracks or false to ignore cracks and use single velocity field
// Can only ignore cracks in multimaterial mode.
int MaterialBase::AllowsCracks(void) const { return allowsCracks; }

// Find the compression needed for bulk modulus to increase Kmax fold
double MaterialBase::GetMGEOSXmax(double gamma0,double S1,double S2,double S3,double &Kmax)
{
	// increase until passes kmax or become negative
	double xl = 0.,xr = 0.;
	double step = 0.05,denom,kr,kl=1.;
	while(true)
	{	xr += step;
		denom = 1./(1. - xr*(S1 + xr*(S2 + xr*S3)));
		
		// if less then zero, went too far
		if(denom<=0.)
		{	xr -= step;
			step /= 2.;
			continue;
		}
		
		// check ratio
		kr = (1.-0.5*gamma0*xr)*denom*denom;
		if(kr>Kmax) break;
		
		// did it pass a maximum
		if(kr<kl)
		{	Kmax = kl;
			return xl;
		}
		
		// update and continue
		kl = kr;
		xl = xr;
	}
	
	// bisect to improve
	double xmid;
	for(int i=0;i<15;i++)
	{	xmid = (xl+xr)/2.;
		denom = 1./(1. - xmid*(S1 + xmid*(S2 + xmid*S3)));
		kr = (1.-0.5*gamma0*xmid)*denom*denom;
		if(kr<Kmax)
			xl = xmid;
		else
			xr = xmid;
	}
	
	return xmid;
}

// lquids return shear-rate dependent viscosity. Only used in liquid contact.
double MaterialBase::GetViscosity(double shearRate) const { return -1.; }

// This method has several options (and only used in liquid contact):
//	1. If can solve x(1+k*eta(x*gmaxdot)) - 1 = 0 on interval 0 < x < 1, then return eta(g(dot))
//		Note: x = g(dot)/gmaxdot
//  2. If solution not possible, bracket the solution to y = x(1+k*eta(x*gmaxdot)) - 1 = 0
//		Such that y1(x1)<0 and y(x2)>0 and return -1
//	3. If can't help, set brackets to (0,-1) and (1,k*viscosity) and return -1.
double MaterialBase::BracketContactLawShearRate(double gmaxdot,double k,double &x1,double &y1,double &x2,double &y2) const
{	x1 = 0.;
	y1 = -1.;
	x2 = 1.;
	y2 = k*GetViscosity(gmaxdot);
	return -1;
}
