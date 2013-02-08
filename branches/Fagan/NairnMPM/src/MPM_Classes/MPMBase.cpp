/********************************************************************************
    MPMBase.cpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"

// globals
MPMBase **mpm;		// list of material points
int nmpms=0;		// number of material points

bool rhoScaling = FALSE;	// modiftf default to no mass scaling. Vincent 
double rhoScale = 1;		// modiftf default to a value of 1, though not used. Vincent

// class statics
int MPMBase::currentParticleNum=0;

#pragma mark MPMBase::Constructors and Destructors

// Constructors
MPMBase::MPMBase()
{
}

MPMBase::MPMBase(int elem,int theMatl,double angin)
{
    int i;
    
    inElem=elem;
	mp=-1.;						// calculated in PreliminaryCalcs, unless set in input file
    matnum=theMatl;
	SetAnglez0InDegrees(angin);
	SetAngley0InDegrees(0.0);
	SetAnglex0InDegrees(0.0);
    
    for(i=1;i<MaxShapeNds;i++)
        vfld[i]=NO_CRACK;
        
    // zero stresses and strains
	ZeroTensor(&sp);
	ZeroTensor(&ep);
	ZeroTensor(&eplast);
	ZeroTensorAntisym(&wrot);
    
    // zero work and energies
    extWork=0.;
    plastEnergy=0.;
	dispEnergy=0.;
    strainEnergy=0.;
	
	// for J integral if needed
	velGrad=NULL;
	
	// material data is needed
	matData=NULL;
	
	// concentration (c units) and gradient (c units/mm)
	pConcentration=0.;
	pDiffusion=NULL;
    
	// history variables for Microstructure Model
	// **Perhaps place if (Microstructure Model TRUE) here?
	yieldC=yieldP=0.;
	archiverhoC=archiverhoW=archiveTDL=archiveDSize=fr=0.;
 	rhoWDot=rhoCDot=0.;								
	initfr=-3.;	
	// temperature (degrees) and gradient (degrees/mm)
	SetTemperature(0.,0.);
	pTemp=NULL;
    
    // CPDI data
    cpdi = NULL;
    faceArea = NULL;

    // PS - when point created, velocity and position and ext force should be set too
	ZeroVector(&vel);
	
	elementCrossings=0;
}

// allocation diffusion data if need in this calculations
void MPMBase::AllocateTemperature(void)
{	pTemp=new TemperatureField;
	ZeroVector(&pTemp->DT);
}

// allocation diffusion data if need in this calculations
void MPMBase::AllocateDiffusion(void)
{	pDiffusion=new DiffusionField;
	pDiffusion->flux=0.;
	ZeroVector(&pDiffusion->Dc);
}

// allocation diffusion data if need in this calculations
void MPMBase::AllocateJStructures(void)
{	velGrad=new Tensor;
}

// allocate structures when needed for particle domains
bool MPMBase::AllocateCPDIStructures(int gimpType,bool isThreeD)
{
    int cpdiSize=0;
    
    if(gimpType==LINEAR_CPDI)
        cpdiSize = isThreeD ? 8 : 4 ;
    else if(gimpType==QUADRATIC_CPDI)
        cpdiSize = 9;
    else
        return TRUE;
    
    // create memory for cpdiSize pointers
    cpdi = (CPDIDomain **)malloc(sizeof(LinkedObject *)*(cpdiSize));
    if(cpdi == NULL) return FALSE;
    
    // create each one
    int i;
    for(i=0;i<cpdiSize;i++)
    {   cpdi[i] = new CPDIDomain;
        cpdi[i]->wg.z = 0.;             // set zero once for 2D calculations
        cpdi[0]->ncpos.z = 0.;          // set zero once for 2D calculations
    }
    
    // save face areas (or lengths in 2D)
    if(firstTractionPt!=NULL) faceArea = new Vector;
    
    return TRUE;
        
}

// Destructor (and it is virtual)
MPMBase::~MPMBase() { }

#pragma mark MPMBase::Methods

// reverse the direction (should only be done for rigid material particles)
void MPMBase::ReverseParticle(void)
{	vel.x=-vel.x;
    vel.y=-vel.y;
	vel.z=-vel.z;
}

// stop the particle (should only be done for rigid material particles)
void MPMBase::StopParticle(void)
{	vel.x=0.;
    vel.y=0.;
	vel.z=0.;
}

#pragma mark MPMBase::Accessors

// set mass in PreliminaryCalcs, but only if input file did not set it first
void MPMBase::InitializeMass(double mass)
{	if(mp<0.) mp=mass;
}

// set or average velocity gradient (2D) (only needed for second J Term)
// for USF and SZS, sets gradient on its only pass
// for USAVG, stores in on first pass and then averages on second (J calc follows the second pass)
void MPMBase::SetVelocityGradient(double dvxx,double dvyy,double dvxy,double dvyx,int secondPass)
{
    if(JTerms!=2) return;
    if(!secondPass)
    {	velGrad->xx=dvxx;
        velGrad->yy=dvyy;
        velGrad->xy=dvxy;
        velGrad->zz=dvyx;		// yx stored in zz
    }
    else
    {	velGrad->xx=(velGrad->xx+dvxx)*0.5;
        velGrad->yy=(velGrad->yy+dvyy)*0.5;
        velGrad->xy=(velGrad->xy+dvxy)*0.5;
        velGrad->zz=(velGrad->zz+dvyx)*0.5;		// yx stored in zz
    }
}

// set concentration (only, and for all particles, during initialization)
void MPMBase::SetConcentration(double pConc,double pRefConc)
{	pConcentration=pConc;
	pPreviousConcentration=pRefConc;
}

// set temperature (only and for all particles, during initialization)
void MPMBase::SetTemperature(double pTempSet,double pRefTemp)
{	pTemperature=pTempSet;
	pPreviousTemperature=pRefTemp;
}

// material ID (convert to zero based)
int MPMBase::MatID(void) { return matnum-1; }			// zero based material array in data storage
int MPMBase::ArchiveMatID(void) { return matnum; }		// one based for archiving

// element ID (convert to zero based)
int MPMBase::ElemID(void) { return inElem-1; }					// zero based element array in data storage
void MPMBase::ChangeElemID(int newElem)
{	inElem=newElem+1;		// set using zero basis
	elementCrossings++;		// count crossings
}
int MPMBase::ArchiveElemID(void) { return inElem; }			// one based for archiving

// return current element crossings for archiving and reset to zero
int MPMBase::GetResetElementCrossings(void)
{	int was=elementCrossings;
	elementCrossings=0;
	return was;
}

// get deformation gradient terms
double MPMBase::GetDuDy(void) { return (ep.xy+eplast.xy-wrot.xy)/2.; }
double MPMBase::GetDvDx(void) { return (ep.xy+eplast.xy+wrot.xy)/2.; }

// anglez0 is initial z cw orientation angle (2D and 3D z,y,x scheme)
double MPMBase::GetRotationZ(void) { return anglez0-0.5*wrot.xy; }
double MPMBase::GetParticleRotationZ(void) { return -0.5*wrot.xy; }		// rotation in simulation, ignoring initial angle
double MPMBase::GetRotationZInDegrees(void) { return 180.*(anglez0-0.5*wrot.xy)/PI_CONSTANT; }
void MPMBase::SetAnglez0(double angle) { anglez0=angle; }
void MPMBase::SetAnglez0InDegrees(double angle) { anglez0=PI_CONSTANT*angle/180.; }
double MPMBase::GetAnglez0InDegrees(void) { return 180.*anglez0/PI_CONSTANT; }

// angley0 is initial y cw orientation angle (3D z,y,x scheme)
double MPMBase::GetRotationY(void) { return angley0+0.5*wrot.xz; }
double MPMBase::GetParticleRotationY(void) { return +0.5*wrot.xz; }		// rotation in simulation, ignoring initial angle
double MPMBase::GetRotationYInDegrees(void) { return 180.*(angley0+0.5*wrot.xz)/PI_CONSTANT; }
void MPMBase::SetAngley0(double angle) { angley0=angle; }
void MPMBase::SetAngley0InDegrees(double angle) { angley0=PI_CONSTANT*angle/180.; }
double MPMBase::GetAngley0InDegrees(void) { return 180.*angley0/PI_CONSTANT; }

// anglex0 is initial x cw orientation angle (3D z,y,x scheme)
double MPMBase::GetRotationX(void) { return anglex0-0.5*wrot.yz; }
double MPMBase::GetParticleRotationX(void) { return -0.5*wrot.yz; }		// rotation in simulation, ignoring initial angle
double MPMBase::GetRotationXInDegrees(void) { return 180.*(anglex0-0.5*wrot.yz)/PI_CONSTANT; }
void MPMBase::SetAnglex0(double angle) { anglex0=angle; }
void MPMBase::SetAnglex0InDegrees(double angle) { anglex0=PI_CONSTANT*angle/180.; }
double MPMBase::GetAnglex0InDegrees(void) { return 180.*anglex0/PI_CONSTANT; }

// 2D and 3D rotational strain tensor
void MPMBase::IncrementRotationStrain(double rotXY) { wrot.xy+=rotXY; }
void MPMBase::IncrementRotationStrain(double rotXY,double rotXZ,double rotYZ)
{	wrot.xy+=rotXY;
	wrot.xz+=rotXZ;
	wrot.yz+=rotYZ;
}
TensorAntisym *MPMBase::GetRotationStrainTensor(void) { return & wrot; }

// energies
double MPMBase::GetPlastEnergy(void) { return plastEnergy; }
void MPMBase::AddPlastEnergy(double energyInc) { plastEnergy+=energyInc; }
double MPMBase::GetDispEnergy(void) { return dispEnergy; }
void MPMBase::AddDispEnergy(double energyInc) { dispEnergy+=energyInc; }
void MPMBase::SetDispEnergy(double energy) { dispEnergy=energy; }
double MPMBase::GetStrainEnergy(void) { return strainEnergy; }
void MPMBase::SetStrainEnergy(double energyTot) { strainEnergy=energyTot; }
void MPMBase::AddStrainEnergy(double energyInc) { strainEnergy+=energyInc; }
double MPMBase::GetExtWork(void) { return extWork; }

// pointers to variables
Vector *MPMBase::GetPFext(void) { return &pFext; }
Vector *MPMBase::GetNcpos(void) { return &ncpos; }
CPDIDomain **MPMBase::GetCPDIInfo(void) { return cpdi; }
Vector *MPMBase::GetAcc(void) { return &acc; }
Tensor *MPMBase::GetVelGrad(void) { return velGrad; }
Tensor *MPMBase::GetStressTensor(void) { return &sp; }
Tensor *MPMBase::GetStrainTensor(void) { return &ep; }

// These methods return pointer to the same symmetric tensor. A material class must
// decide how it will use the eplast tensor variable and use only that one purpose.
// The only reason for two names is to make the code more readable
Tensor *MPMBase::GetPlasticStrainTensor(void) { return &eplast; }
Tensor *MPMBase::GetElasticLeftCauchyTensor(void) { return &eplast; }

// history data
char *MPMBase::GetHistoryPtr(void) { return matData; }
void MPMBase::SetHistoryPtr(char *createMatData)
{	if(matData!=NULL) free(matData);
	matData=createMatData;
}
double MPMBase::GetHistoryDble(void)
{	// assumes matData starts with a double and gets that first one
	double *h=(double *)matData;
	return *h;
}
void MPMBase::SetHistoryDble(double history)
{	// assumes matData starts with a double and sets that first one
	double *h=(double *)matData;
	*h=history;
}
double MPMBase::GetHistoryDble(int index)
{	// assumes matData is array of doubles and gets one an index (0 based)
	double *h=(double *)matData;
	return h[index];
}
void MPMBase::SetHistoryDble(int index,double history)
{	// assumes matData is array of doubles and sets one an index (0 based)
	double *h=(double *)matData;
	h[index]=history;
}

#pragma mark MPMBase Class Methods

/**********************************************************
    Update Strains on all particles
    Must impose grid velocity BCs and velocity
		alterations due to contact first
	secondPass will be TRUE only for USAVG method
**********************************************************/
void MPMBase::FullStrainUpdate(double strainTime,int secondPass,int np)
{
    NodalPoint::GetGridVelocitiesForStrainUpdate();			// velocities needed for strain update
	
	// loop over particles
	for(MPMBase::currentParticleNum=0;MPMBase::currentParticleNum<nmpms;MPMBase::currentParticleNum++)
		mpm[MPMBase::currentParticleNum]->UpdateStrain(strainTime,secondPass,np);
}


