/****************************************************************************************************************************
    MaterialBase.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Material Classes
		+ indicates actual material (others cannot be materials)
		(MPM) or (FEA) indicate only that code, otherwise in both
	
	MaterialBase (base class)
	
	Subclasses
	
	  Elastic	->	IsotropicMat+	-> IsoPlasticity	-> VonMisesHardening+ (MPM)	
														-> MGSCGLMaterial+ (MPM) -> SLMaterial+ (MPM)
									-> BistableIsotropic+ (MPM)
													
				->	TransIsotropic+	->	Orthotropic+	-> AnisoPlasticity (MPM)	-> HillPlastic+ (MPM)
					
	  Viscoelastic+ (MPM)
					
	  RubberElastic		->	Mooney+ (MPM)
					
	  ImperfectInterface+ (FEA)
					
	  RigidMaterial+
	  
	  TractionLaw (MPM)	->	CohesiveZone+ (MPM)	-> LinearTraction+ (MPM)
						->	CubicTraction+ (MPM)
****************************************************************************************************************************/

#include "Materials/MaterialBase.hpp"

// Materials globals
MaterialBase **theMaterials=NULL;		// list of materials
int nmat=0;								// number of materials

#pragma mark MaterialBase::Constructors and Destructors

// Constructors
MaterialBase::MaterialBase() {}

// Constructors
MaterialBase::MaterialBase(char *matName)
{
	name=new char[strlen(matName)+1];
    strcpy(name,matName);
    hasMatProps=FALSE;
	concSaturation=1.;
	betaI=0.;
	red=-1.;
#ifdef MPM_CODE
    rho=1.;
    KIc=KIIc=JIc=JIIc=gamma=-1.;		// traction laws assumes -1
	delIc=delIIc=-1.;					// traction laws assumes -1
	pCrit3=1.;
	maxLength=-1.;
	gain=1.e5;
	initTime=-1.;
	initSpeed=1.;
	diffusionCon=0.;
    // KIexp and KIIexp are the parameters of the empirical fracture criterion:
    // (KI/KIc)^KIexp+(KII/KIIc)^KIIexp=1.
    // Set the default values to 2 for elliptical fracture locus.
    KIexp=KIIexp=2.;
    criterion=UNSPECIFIED;
	constantDirection=FALSE;
	growDir.x=0.;
	growDir.y=0.;
	matPropagateDirection=UNSPECIFIED;
	tractionMat=0;
	field=-1;
	kCond=0.;
	heatCapacity=-1.;
	heatCapacityVol=-1.;
	lastFriction=NULL;
#endif
}

// Destructor (and it is virtual)
MaterialBase::~MaterialBase()
{	delete [] name;
}

#pragma mark MaterialBase::Initialization

// print to output window
void MaterialBase::PrintMaterial(int num)
{
	// material name
    cout << "Material " << num << ": " << name << endl;
    cout << "     " << MaterialType() << " with:" << endl;
	
	// call class to list mechanical properties
	PrintMechanicalProperties();
	
	// properties in command base class (differs for MPM and FEA)
	PrintCommonProperties();
    
#ifdef MPM_CODE
	// Transport Properties
	PrintTransportProperties();
#endif
}

// print mechanical properties to output window
void MaterialBase::PrintMechanicalProperties(void)
{	cout << "Undocumented mechanical properties" << endl;
}

#define PROP_LABEL_LENGTH 4
#define PROP_COLUMN 20

// print property with units (option) and align to columns for material properties
void MaterialBase::PrintProperty(const char *propName,double value,const char *unitsStr)
{
	char prop[200];
	
	// name
	/*
	strcpy(prop,"");
	int length=PROP_LABEL_LENGTH-strlen(propName);
	if(length>0)
	{	while(length>1)
		{	strcat(prop," ");
			length--;
		}
	}
	strcat(prop,propName);
	if(length>0) strcat(prop," ");
	*/
	strcpy(prop,propName);
	int length=strlen(prop);
	while(length<PROP_LABEL_LENGTH)
	{	strcat(prop," ");
		length++;
	}
	
	// value
	char valueStr[50];
	sprintf(valueStr,"= %g",value);
	strcat(prop,valueStr);
	
	// units
	if(strlen(unitsStr)>0)
	{	strcat(prop," ");
		strcat(prop,unitsStr);
	}
	
	// pad to column
	length=strlen(prop);
	int padLength = (length<PROP_COLUMN) ? PROP_COLUMN : 2*PROP_COLUMN ;
	while(length<padLength)
	{	strcat(prop," ");
		length++;
	}
	
	cout << prop;
}

// print text aligned with material columns, left or right justified
void MaterialBase::PrintProperty(const char *text,bool rightJustify)
{
	char prop[200];
	
	// up to two columns long
	int length=strlen(text);
	int padLength = (length<PROP_COLUMN) ? PROP_COLUMN : 2*PROP_COLUMN ;
	strcpy(prop,"");
	
	// right justify
	if(rightJustify)
	{	int pad=padLength-length-1;
		while(pad>0)
		{	strcat(prop," ");
			pad--;
		}
		strcat(prop,text);
		strcat(prop," ");
	}
	
	// left justify
	else
	{	strcat(prop,text);
		while(length<padLength)
		{	strcat(prop," ");
			length++;
		}
	}
	
	cout << prop;
}

    
// Called before analysis, material can fill in things that never change during the analysis
// Note: no angle, because cannot depend on material angle
void MaterialBase::InitialLoadMechProps(int makeSpecific,int np) { return; }

// return material type
const char *MaterialBase::MaterialType(void) { return "Unknown Material Type"; }

#pragma mark MaterialBase::Methods

// rotate meechanical properties if needed
void MaterialBase::LoadMechProps(int makeSpecific,double angle,int np) { return; }

