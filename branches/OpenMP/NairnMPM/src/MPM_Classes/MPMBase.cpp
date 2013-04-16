/********************************************************************************
    MPMBase.cpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Materials/MaterialBase.hpp"

// globals
MPMBase **mpm;		// list of material points
int nmpms=0;		// number of material points
int nmpmsNR=0;		// number for last non-rigid material point
int nmpmsRC=0;		// number for last rigid contact material

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
    
	// space to hold velocity fields
	// these only change if there are cracks
	vfld = (char *)malloc(maxShapeNodes*sizeof(char));
    for(i=1;i<maxShapeNodes;i++)
        vfld[i]=NO_CRACK;
	
    // zero stresses and strains
	ZeroTensor(&sp);
    pressure = 0.;
	ZeroTensor(&ep);
	ZeroTensor(&eplast);
	ZeroTensorAntisym(&wrot);
    
    // zero work and energies
    extWork=0.;
    plastEnergy=0.;
	dispEnergy=0.;
    strainEnergy=0.;
    heatEnergy=0.;
	
	// for J integral if needed (on non-rigid only)
	velGrad=NULL;
	
	// material data is needed
	matData=NULL;
	
	// concentration (c units) and gradient (c units/mm)
	pConcentration=0.;
	pDiffusion=NULL;
    
	// temperature (degrees) and gradient (degrees/mm)
	SetTemperature(0.,0.);
	pTemp=NULL;
    
    // CPDI data
    cpdi = NULL;
    faceArea = NULL;

    // PS - when point created, velocity and position and ext force should be set too
	ZeroVector(&vel);
	
	// counts crossing and sign is whether or not left the grid
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
	ZeroVector(&pDiffusion->Dc);
}

// allocation velGrad tensor data if need in this calculations (non rigid only)
void MPMBase::AllocateJStructures(void)
{	velGrad=new Tensor;
}

// allocate structures when needed for particle domains
bool MPMBase::AllocateCPDIStructures(int gimpType,bool isThreeD)
{
    int cpdiSize=0;
    
    if(gimpType==LINEAR_CPDI || gimpType==LINEAR_CPDI_AS)
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
		
		// weights constant except for axisymmetric CPDI
		if(gimpType==LINEAR_CPDI)
			cpdi[i]->ws = isThreeD ? 0.125 : 0.25 ;
		else if(gimpType==QUADRATIC_CPDI)
		{	if(i<4)
				cpdi[i]->ws = 1./36.;
			else if(i<8)
				cpdi[i]->ws = 1./9.;
			else
				cpdi[i]->ws = 4./9.;
		}
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
    {	velGrad->xx = dvxx;
        velGrad->yy = dvyy;
        velGrad->xy = dvxy;
        velGrad->zz = dvyx;						// yx stored in zz
    }
    else
    {	velGrad->xx = (velGrad->xx+dvxx)*0.5;
        velGrad->yy = (velGrad->yy+dvyy)*0.5;
        velGrad->xy = (velGrad->xy+dvxy)*0.5;
        velGrad->zz = (velGrad->zz+dvyx)*0.5;		// yx stored in zz
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

// zero the temperature gradient (non-rigid particles only)
void MPMBase::AddTemperatureGradient(void)
{   pTemp->DT.x=0.;
    pTemp->DT.y=0.;
	pTemp->DT.z=0.;
}

// zero the concentration gradient
void MPMBase::AddConcentrationGradient(void)
{	pDiffusion->Dc.x=0.;
    pDiffusion->Dc.y=0.;
    pDiffusion->Dc.z=0.;
}

// material ID (convert to zero based)
int MPMBase::MatID(void) { return matnum-1; }			// zero based material array in data storage
int MPMBase::ArchiveMatID(void) { return matnum; }		// one based for archiving

// element ID (convert to zero based)
int MPMBase::ElemID(void) { return inElem-1; }					// zero based element array in data storage
void MPMBase::ChangeElemID(int newElem)
{	inElem=newElem+1;				// set using zero basis
	IncrementElementCrossings();		// count crossing
}
int MPMBase::ArchiveElemID(void) { return inElem; }			// one based for archiving

// return current element crossings for archiving and reset to zero
int MPMBase::GetElementCrossings(void) { return elementCrossings>=0 ? elementCrossings : -elementCrossings; }
void MPMBase::SetElementCrossings(int ec) { elementCrossings = ec; }
void MPMBase::IncrementElementCrossings(void)
{	elementCrossings = elementCrossings>=0 ? elementCrossings+1 : elementCrossings-1;
}
bool MPMBase::HasLeftTheGrid(void) { return elementCrossings<0; }
void MPMBase::SetHasLeftTheGrid(bool setting) { elementCrossings = setting ? -GetElementCrossings() : GetElementCrossings(); }

// get unscaled volume for use only in contact and imperfect interface calculations
// return result in mm^3
double MPMBase::GetUnscaledVolume(void)
{	double rho=theMaterials[MatID()]->rho*0.001;			// in g/mm^3
	return mp/rho;                                          // in mm^3
}

// get mass when finding mass gradient for contact calculations
// Axisymmetric particles override to return mass/rp to get uniform particle mass
double MPMBase::GetMassForGradient(void) { return mp; }

// return if mterial for this particle includes plastic strainin gradient or if
// entire deformation is in the elastic strain
bool MPMBase::PartitionsElasticAndPlasticStrain(void)
{   return theMaterials[MatID()]->PartitionsElasticAndPlasticStrain();
}

// get deformation gradient terms from full strains
// used when archiving and by 2D J Integral calculations
double MPMBase::GetDuDx(void)
{   if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
		return ep.xx+eplast.xx;
	else
		return ep.xx;
}
double MPMBase::GetDuDy(void)
{   if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
        return (ep.xy+eplast.xy-wrot.xy)/2.;
    else
        return (ep.xy-wrot.xy)/2.;
}
double MPMBase::GetDvDx(void)
{   if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
        return (ep.xy+eplast.xy+wrot.xy)/2.;
    else
        return (ep.xy+wrot.xy)/2.;
}
double MPMBase::GetDvDy(void)
{   if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
		return ep.yy+eplast.yy;
	else
		return ep.yy;
}
double MPMBase::GetDwDz(void)
{   if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
		return ep.zz+eplast.zz;
	else
		return ep.zz;
}

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

// a material should never call this direction. It is only called in IncrementHeatEnergy.
void MPMBase::AddDispEnergy(double energyInc) { dispEnergy+=energyInc; }
void MPMBase::SetDispEnergy(double energy) { dispEnergy=energy; }
double MPMBase::GetStrainEnergy(void) { return strainEnergy; }
void MPMBase::SetStrainEnergy(double energyTot) { strainEnergy=energyTot; }
void MPMBase::AddStrainEnergy(double energyInc) { strainEnergy+=energyInc; }
double MPMBase::GetExtWork(void) { return extWork; }
double MPMBase::GetHeatEnergy(void) { return heatEnergy; }
void MPMBase::SetHeatEnergy(double energyTot) { heatEnergy = energyTot; }
void MPMBase::AddHeatEnergy(double energyInc) { heatEnergy += energyInc; }
double MPMBase::GetInternalEnergy(void) { return strainEnergy + heatEnergy; }

// pointers to variables
Vector *MPMBase::GetPFext(void) { return &pFext; }
Vector *MPMBase::GetNcpos(void) { return &ncpos; }
CPDIDomain **MPMBase::GetCPDIInfo(void) { return cpdi; }
Vector *MPMBase::GetAcc(void) { return &acc; }
Tensor *MPMBase::GetVelGrad(void) { return velGrad; }
Tensor *MPMBase::GetStrainTensor(void) { return &ep; }

// material classes only should call GetStressTensor and can change the stress
// all others should use ReadStressTensor() because that gives materials a chance to
//  customize the way stresses are stored, such as to separate pressure and deviatoric stress
Tensor *MPMBase::GetStressTensor(void) { return &sp; }
Tensor MPMBase::ReadStressTensor(void) { return theMaterials[MatID()]->GetStress(&sp,pressure); }
void MPMBase::IncrementPressure(double dP) { pressure += dP; }
void MPMBase::SetPressure(double P) { pressure = P; }
double MPMBase::GetPressure(void) { return pressure; }

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

// Decribe for debugging use or output on some errors
void MPMBase::Describe(void)
{	cout << "# pt: pos=(" << pos.x << "," << pos.y << "," << pos.z << ") mass=" << mp << 
                " matl=" << matnum << " elem=" << inElem << endl;
    cout << "#     vel=(" << vel.x*1e-3 << "," << vel.y*1e-3 << "," << vel.z*1e-3 << ") m/sec" << endl;
    Matrix3 pF = GetDeformationGradientMatrix();
    cout << "#       F=" << pF << ", |F|=" << pF.determinant() << endl;
    double rho0=theMaterials[MatID()]->rho;
    double rho = rho0/theMaterials[MatID()]->GetCurrentRelativeVolume(this);
    cout << "#       P= " << pressure*rho*1e-6 << " MPa" << endl;
    cout << "# sigmaii=(" << sp.xx*rho*1e-6 << "," << sp.yy*rho*1e-6 << "," << sp.zz*rho*1e-6 << ") MPa" << endl;
    cout << "#   tauij=(" << sp.xy*rho*1e-6 << "," << sp.xz*rho*1e-6 << "," << sp.yz*rho*1e-6 << ") MPa" << endl;
}



