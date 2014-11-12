/****************************************************************************************************************************
    MaterialBase.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Material Classes
		+ indicates actual material (others cannot be materials)
		(MPM) or (FEA) indicate only that code, otherwise in both
	
	MaterialBase (base class)
****************************************************************************************************************************/

#include "Materials/MaterialBase.hpp"

// Constants for output fields for material properties
#define PROP_LABEL_LENGTH 4
#define PROP_COLUMN 20

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
	concSaturation=1.;
	betaI=0.;
	red=-1.;
#ifdef MPM_CODE
    rho=1.;
    KIc=KIIc=JIc=JIIc=gamma=-1.;		// traction laws assumes -1
	delIc=delIIc=-1.;					// traction laws assumes -1
	nmix=1.;							// mixed mode exponent (if used)
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
    criterion[0]=criterion[1]=UNSPECIFIED;
	constantDirection=FALSE;
	constantTip=FALSE;
	growDir.x=0.;
	growDir.y=0.;
	matPropagateDirection[0]=matPropagateDirection[1]=UNSPECIFIED;
	tractionMat[0]=tractionMat[1]=0;
	field=-1;
	shareMatField=0;					// share field with another material
	activeField=-1;
	kCond=0.;
	heatCapacity=1.;					// keep one because needed by ideal gas
	lastFriction=NULL;
    artificialViscosity=FALSE;
    avA1 = 0.2;
    avA2 = 2.0;
#endif
}

// Destructor (and it is virtual)
MaterialBase::~MaterialBase()
{	delete [] name;
}

#pragma mark MaterialBase::Initialization

// print to output window
void MaterialBase::PrintMaterial(int num) const
{
	// material name
    cout << "Material " << num << ": " << name << endl;
    cout << "     " << MaterialType() << " with:" << endl;
	
	// call class to list mechanical properties
	PrintMechanicalProperties();
	
	// properties in common base class (differs for MPM and FEA)
	PrintCommonProperties();
    
#ifdef MPM_CODE
	// Transport Properties
	PrintTransportProperties();
#endif
}

// print mechanical properties to output window
void MaterialBase::PrintMechanicalProperties(void) const {}

// print property with units (option) and align to columns for material properties
// In field of width PROP_COLUMN
void MaterialBase::PrintProperty(const char *propName,double value,const char *unitsStr)
{
	char prop[200];
	
	// name
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
// In field of width PROP_COLUMN, will use two if needed
void MaterialBase::PrintProperty(const char *text,bool rightJustify)
{
	char prop[200];
	
	// up to two columns wide
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

    
// return material type
const char *MaterialBase::MaterialType(void) const { return "Unknown Material Type"; }

