/********************************************************************************
    VTKArchive.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 12/5/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
	
	To add new quantity:
		1. Add enum in ArchiveData.hpp
		2. Add parameter in InputParam() and set its size. If will not 
			extraplate set size to minus actual size (scalars only currently)
		3. If will extrapolate, add case in NodalExtrapolation()
		4. Add case to write the data in ArchiveVTKFile(), but if extrapolated
			and vtk is NULL (memory error) skip writing the data.
	
	Possible Quantities to Add
		material angle(s) or rotation strain, kinetic energy
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/VTKArchive.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "Exceptions/CommonException.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"

// globals
VTKArchive *vtkArchiveTask=NULL;
int dummyArg;

#pragma mark CREATION OF THE TASK

// Constructors
VTKArchive::VTKArchive() : GridArchive()
{
	bufferSize=0;
	vtk=NULL;
    intIndex=0;
}

// Return name of this task
const char *VTKArchive::TaskName(void) { return "Archive grid results to VTK files"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
// not thread safe due to push_back()
// throws std::bad_alloc
char *VTKArchive::InputParam(char *pName,int &input,double &gScaling)
{
	int q=-1,thisBuffer=0,pindex=-1;
    
    // default return value
    char *retPtr = (char *)&dummyArg;
	
    // check for archiving quantity
    if(strcmp(pName,"mass")==0)
    {	q=VTK_MASS;
		// no buffer since no need to extrapolate
		thisBuffer=-1;
    }
    
    else if(strcmp(pName,"numpoints")==0)
    {	q=VTK_NUMBERPOINTS;
		thisBuffer=-1;
    }
	
    else if(strcmp(pName,"velocity")==0)
    {	q=VTK_VELOCITY;
		thisBuffer=3;		// extrapolate for ctr mass velocity
    }
	
    else if(strcmp(pName,"stress")==0)
    {	q=VTK_STRESS;
		thisBuffer=6;
    }
	
    else if(strcmp(pName,"strain")==0 || strcmp(pName,"elasticstrain")==0)
    {	q=VTK_STRAIN;
		thisBuffer=6;
    }
	
    else if(strcmp(pName,"defgrad")==0)
    {	q=VTK_DEFGRAD;
		thisBuffer=9;
    }
	
    else if(strcmp(pName,"totalstrain")==0)
    {	q=VTK_TOTALSTRAIN;
		thisBuffer=6;
    }
	
    else if(strcmp(pName,"displacement")==0)
    {	q=VTK_DISPLACEMENT;
		thisBuffer=3;
    }
	
    else if(strcmp(pName,"plasticstrain")==0)
    {	q=VTK_PLASTICSTRAIN;
		thisBuffer=6;
    }
	
    else if(strcmp(pName,"temperature")==0)
    {	q=VTK_TEMPERATURE;
		// no buffer since no need to extrapolate
		thisBuffer=-1;
    }
	
    else if(strcmp(pName,"contactforce")==0)
    {	q=VTK_RIGIDCONTACTFORCES;
		// no buffer since no need to extrapolate
        archiver->SetDoingArchiveContact(TRUE);
		thisBuffer=-3;
    }
	
    else if(strcmp(pName,"volumegradient")==0)
    {   q=VTK_VOLUMEGRADIENT;
        // no buffer since no need to extrapolate, but need to read argument as material number
		if(intIndex>=MAX_INTEGER_ARGUMENTS)
		{	cout << "Too many parameter arguments in VTKArchive options" << endl;
			return NULL;
		}
        retPtr = (char *)&intArgs[intIndex];
        pindex = intIndex;
        intIndex++;
        thisBuffer=-3;
    }
	
    else if(strcmp(pName,"reactionforce")==0)
    {	q=VTK_BCFORCES;
		// no buffer since no need to extrapolate
        // currently not implement (or documented) and gets all zeros
		thisBuffer=-3;
    }
	
    else if(strcmp(pName,"pressure")==0)
    {	q=VTK_PRESSURE;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"deltav")==0)
    {	q=VTK_RELDELTAV;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"equivstrain")==0)
    {	q=VTK_EQUIVSTRAIN;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"equivstress")==0)
    {	q=VTK_EQUIVSTRESS;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"concentration")==0 || strcmp(pName,"porepressure")==0)
    {	q=VTK_CONCENTRATION;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"strainenergy")==0 || strcmp(pName,"workenergy")==0)
    {	q=VTK_WORKENERGY;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"plasticenergy")==0)
    {	q=VTK_PLASTICENERGY;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"heatenergy")==0)
    {	q=VTK_HEATENERGY;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"material")==0)
    {	q=VTK_MATERIAL;
		thisBuffer=1;
    }
    
	// if found one, add to arrays
	if(q>=0)
	{	quantity.push_back(q);
		quantitySize.push_back(thisBuffer);				// <0 is size for non-extrapolated quantities
        qparam.push_back(pindex);                       // index to argument array while reading (if>=0) (set to arguments when done)
		//if(thisBuffer<0) thisBuffer=-thisBuffer;
		if(thisBuffer>0) bufferSize+=thisBuffer;
		char *qname=new char[strlen(pName)+1];
		strcpy(qname,pName);
		quantityName.push_back(qname);
		input=INT_NUM;
        return retPtr;
	}
	
	// check remaining commands (e.g. archive times)
    return GridArchive::InputParam(pName,input,gScaling);
}

// called once at start of MPM analysis - initialize and print info
// throws CommonException()
CustomTask *VTKArchive::Initialize(void)
{
    cout << "Archive grid results to VTK files." << endl;
	
	// display quantities to be archived
	unsigned int q;
	cout << "   Archiving: " ;
	int len=14;
	for(q=0;q<quantity.size();q++)
	{	char *name=quantityName[q];
		if(len+strlen(name)>70)
		{	cout << "\n      ";
			len=3;
		}
		cout << name;
        if(quantity[q]==VTK_VOLUMEGRADIENT)
        {   if(qparam[q]<0)
            {   cout << endl;
                throw CommonException("VTKArchive volumegradient must be set to an available material","VTKArchive::Initialize");
            }
            int matnum = intArgs[qparam[q]];
            cout << " (material #" << matnum << ")";
            len+=14;		// this is for matnum<10, but exact length does not matter
            if(matnum<1 || matnum>nmat)
            {   cout << endl;
                throw CommonException("VTKArchive volumegradient must be for an available material","VTKArchive::Initialize");
            }
            qparam[q] = matnum;
        }
		if(q<quantity.size()-1) cout << ", ";
		len+=(int)strlen(name)+2;
	}
	cout << endl;
	
	// parent class prints the archive time
    return GridArchive::Initialize();
}

#pragma mark GRID ARCHVING TASKS

// It is called once (during PrepareForStep()) and doExport will be true or false is scheduled to archive.
//	  A subclass can cancel or force archiving by changind doExport.
// Return true or false if export will need extrapolations (but false if doExport is false)
bool VTKArchive::CheckExportForExtrapolations(void)
{	// task never picked something to export, cancel exports
	if(quantity.size()==0)
	{	doExport = false;
		return false;
	}
	
	// only need to extrapolate if about to export and has a buffer
	return doExport ? bufferSize>0 : false;
}

// Allocate buffers to store extrapolation results for each node
void VTKArchive::AllocateExtrapolationBuffers(void)
{
	// create list of pointers to buffer for each nodal point
	vtk = new (nothrow) double *[nnodes+1];
	
	// on error, print message and turn off rest of this task
	if(vtk==NULL)
	{	cout << "# memory error preparing data for vtk export" << endl;
		getExtrapolations = false;
		return;;
	}
	
	// create buffer for each node, on error print message and turn of extrapolations
	int i,j;
	for(i=1;i<=nnodes;i++)
	{	vtk[i] = new (nothrow) double[bufferSize];
		if(vtk[i]==NULL)
		{	// free previous buffers and main list
			for(j=1;j<i;j++) delete [] vtk[j];
			delete [] vtk;
			cout << "# memory error preparing data for vtk export" << endl;
			getExtrapolations = false;
			return;
		}
		
		// initialize all values to zero
		for(j=0;j<bufferSize;j++) vtk[i][j]=0.;
	}
}

// Extrapolate particle data for particle *mpnt to nodal point *ndmi
// vfld and matfld are the crack velocity and material velocity field for this particle-node paiur
// wt is extrapolation weight and equal to mp*Sip (or particle mass times the shape function)
// isRigid will be true or false if material for this particle is a rigid contact particle
CustomTask *VTKArchive::NodalExtrapolation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,double wt,short isRigid)
{
    // have to skip rigid because nodal masses ignore rigid materials
	if(!getExtrapolations || isRigid) return nextTask;
	
	unsigned int q;
	double *vtkquant=vtk[ndmi->num];
	double theWt=1.,rho,rho0,se;
	Tensor *ten=NULL,sp;
    
    // this loop for non-rigid particles
    for(q=0;q<quantity.size();q++)
    {	switch(quantity[q])
        {	case VTK_STRESS:
            case VTK_PRESSURE:
            case VTK_EQUIVSTRESS:
                rho0=theMaterials[mpnt->MatID()]->GetRho(mpnt);
                rho = rho0/theMaterials[mpnt->MatID()]->GetCurrentRelativeVolume(mpnt,0);
                theWt = wt*rho*UnitsController::Scaling(1.e-6);			// Legacy convert to MPa
                sp = mpnt->ReadStressTensor();
				switch(quantity[q])
				{	case VTK_PRESSURE:
						// pressure
						*vtkquant += -theWt*(sp.xx+sp.yy+sp.zz)/3.;
						vtkquant++;
						break;
					case VTK_EQUIVSTRESS:
						// equivalent or vonmises stress = sqrt(3 J2)
						se = pow(sp.xx-sp.yy,2.) + pow(sp.yy-sp.zz,2.) + pow(sp.xx-sp.zz,2.);
						se += 6.*sp.xy*sp.xy;
						if(fmobj->IsThreeD()) se += 6.*(sp.xz*sp.xz + sp.yz*sp.yz);
						*vtkquant += theWt*sqrt(0.5*se);
						vtkquant++;
						break;
					case VTK_STRESS:
						vtkquant[0]+=theWt*sp.xx;
						vtkquant[1]+=theWt*sp.yy;
						vtkquant[2]+=theWt*sp.zz;
						vtkquant[3]+=theWt*sp.xy;
						if(fmobj->IsThreeD())
						{	vtkquant[4]+=theWt*sp.xz;
							vtkquant[5]+=theWt*sp.yz;
						}
						vtkquant+=6;
						break;
				}
                break;
				
			// Elastic strain (absolute in Legacy)
			// New method small strain = Biot strain - archived plastic strain
			// New method hyperelastic = Biot strain from elastic B in plastic strain
			// Membranes = 0
            case VTK_STRAIN:
				if(theMaterials[mpnt->MatID()]->AltStrainContains()==ENG_BIOT_PLASTIC_STRAIN)
				{	// elastic strain = total strain minus plastic strain
					Matrix3 biot = mpnt->GetBiotStrain();
					ten = mpnt->GetAltStrainTensor();
					vtkquant[0] += wt*(biot(0,0)-ten->xx);
					vtkquant[1] += wt*(biot(1,1)-ten->yy);
					vtkquant[2] += wt*(biot(2,2)-ten->zz);
					vtkquant[3] += wt*(biot(0,1)-0.5*ten->xy);
					if(fmobj->IsThreeD())
					{	vtkquant[4] += wt*(biot(0,2)-0.5*ten->xz);
						vtkquant[5] += wt*(biot(1,2)-0.5*ten->yz);
					}
				}
				else if(theMaterials[mpnt->MatID()]->AltStrainContains()==LEFT_CAUCHY_ELASTIC_B_STRAIN)
				{	// get elastic strain from elastic B
					Matrix3 biot = mpnt->GetElasticBiotStrain();
					vtkquant[0] += wt*biot(0,0);
					vtkquant[1] += wt*biot(1,1);
					vtkquant[2] += wt*biot(2,2);
					vtkquant[3] += wt*biot(0,1);
					if(fmobj->IsThreeD())
					{	vtkquant[4] += wt*biot(0,2);
						vtkquant[5] += wt*biot(1,2);
					}
				}
				else
				{	// elastic strain = total strain for these materials
					Matrix3 biot = mpnt->GetBiotStrain();
					vtkquant[0] += wt*biot(0,0);
					vtkquant[1] += wt*biot(1,1);
					vtkquant[2] += wt*biot(2,2);
					vtkquant[3] += wt*biot(0,1);
					if(fmobj->IsThreeD())
					{	vtkquant[4] += wt*biot(0,2);
						vtkquant[5] += wt*biot(1,2);
					}
				}
                vtkquant+=6;
                break;
				
			// plastic strain (absolute in Legacy)
			// New method small strain = archived plastic strain
			// New method hyperelastic = Biot strain from F - Biot strain from elastic B in plastic strain
			// Membrane = 0
            case VTK_PLASTICSTRAIN:
				if(theMaterials[mpnt->MatID()]->AltStrainContains()==ENG_BIOT_PLASTIC_STRAIN)
				{	// plastic strain is in the tensor already
					ten = mpnt->GetAltStrainTensor();
					vtkquant[0] += wt*ten->xx;
					vtkquant[1] += wt*ten->yy;
					vtkquant[2] += wt*ten->zz;
					vtkquant[3] += 0.5*wt*ten->xy;
					if(fmobj->IsThreeD())
					{	vtkquant[4] += 0.5*wt*ten->xz;
						vtkquant[5] += 0.5*wt*ten->yz;
					}
				}
				else if(theMaterials[mpnt->MatID()]->AltStrainContains()==LEFT_CAUCHY_ELASTIC_B_STRAIN)
				{	// plastic strain = total strain minus elastic strain (from B)
					Matrix3 biotTot = mpnt->GetBiotStrain();
					Matrix3 biotElastic = mpnt->GetElasticBiotStrain();
					vtkquant[0] += wt*(biotTot(0,0)-biotElastic(0,0));
					vtkquant[1] += wt*(biotTot(1,1)-biotElastic(1,1));
					vtkquant[2] += wt*(biotTot(2,2)-biotElastic(2,2));
					vtkquant[3] += wt*(biotTot(0,1)-biotElastic(0,1));
					if(fmobj->IsThreeD())
					{	vtkquant[4] += wt*(biotTot(0,2)-biotElastic(0,2));
						vtkquant[5] += wt*(biotTot(1,2)-biotElastic(1,2));
					}
				}
				// others have zero plastic strain
                vtkquant+=6;
                break;
				
            case VTK_RELDELTAV:
                // Delta V/V0 - small or large strain
                se = mpnt->GetRelativeVolume();
                *vtkquant += wt*(se-1.);
                vtkquant++;
                break;
            
            case VTK_EQUIVSTRAIN:
			{	Matrix3 biot = mpnt->GetBiotStrain();
				double tre = (biot(0,0)+biot(1,1)+biot(2,2))/3.;
                double exx = biot(0,0) - tre;
                double eyy = biot(1,1) - tre;
                double ezz = biot(2,2) - tre;
                se = exx*exx + eyy*eyy * ezz*ezz + 2.0*biot(0,1)*biot(0,1);
                if(fmobj->IsThreeD()) se += 2.0*(biot(0,2)*biot(0,2) + biot(1,2)*biot(1,2));
                *vtkquant += wt*sqrt(2.*se/3.);
                vtkquant++;
                break;
            }
                
			// total strain (absolute in Legacy)
			// New method, small strain and hyperelastic = Biot strain from F
            case VTK_TOTALSTRAIN:
			{	Matrix3 biot = mpnt->GetBiotStrain();
				vtkquant[0] += wt*biot(0,0);
				vtkquant[1] += wt*biot(1,1);
				vtkquant[2] += wt*biot(2,2);
				vtkquant[3] += wt*biot(0,1);
				if(fmobj->IsThreeD())
				{	vtkquant[4] += wt*biot(0,2);
					vtkquant[5] += wt*biot(1,2);
				}
                vtkquant+=6;
                break;
			}
				
			case VTK_DEFGRAD:
			{	Matrix3 F = mpnt->GetDeformationGradientMatrix();
				vtkquant[0] += wt*F(0,0);
                vtkquant[1] += wt*F(0,1);
                vtkquant[2] += wt*F(0,2);
				vtkquant[3] += wt*F(1,0);
                vtkquant[4] += wt*F(1,1);
                vtkquant[5] += wt*F(1,2);
 				vtkquant[6] += wt*F(2,0);
                vtkquant[7] += wt*F(2,1);
                vtkquant[8] += wt*F(2,2);
                vtkquant+=9;
                break;
			}

            case VTK_DISPLACEMENT:
                vtkquant[0]+=wt*(mpnt->pos.x-mpnt->origpos.x);
                vtkquant[1]+=wt*(mpnt->pos.y-mpnt->origpos.y);
                vtkquant[2]+=wt*(mpnt->pos.z-mpnt->origpos.z);
                vtkquant+=3;
                break;
            
            case VTK_VELOCITY:
                vtkquant[0]+=wt*mpnt->vel.x;
                vtkquant[1]+=wt*mpnt->vel.y;
                vtkquant[2]+=wt*mpnt->vel.z;
                vtkquant+=3;
                break;
                
            case VTK_CONCENTRATION:
				if(fmobj->HasDiffusion())
                	theWt=wt*mpnt->GetConcSaturation();
				else
					theWt=wt*UnitsController::Scaling(1.e-6);
                *vtkquant+=theWt*mpnt->pConcentration;
                vtkquant++;
                break;
                
            case VTK_WORKENERGY:
                *vtkquant+=wt*mpnt->mp*mpnt->GetWorkEnergy()*UnitsController::Scaling(1.e-9);
                vtkquant++;
                break;
                
            case VTK_PLASTICENERGY:
                *vtkquant+=wt*mpnt->mp*mpnt->GetPlastEnergy()*UnitsController::Scaling(1.e-9);
                vtkquant++;
                break;
                
            case VTK_HEATENERGY:
                *vtkquant+=wt*mpnt->mp*mpnt->GetHeatEnergy()*UnitsController::Scaling(1.e-9);
                vtkquant++;
                break;
                
            case VTK_MATERIAL:
                *vtkquant+=wt*((double)mpnt->MatID()+1.);
                vtkquant++;
                break;
            
            default:
                // skip those not extrapolated
                break;
        }
    }
	
    return nextTask;
}

// When extrapolations are done, do any remaining calculations. Here divids by nodal mass
void VTKArchive::FinishExtrapolationCalculations(void)
{
	// divide all by nodal mass
	int i,j;
    for(i=1;i<=nnodes;i++)
	{	if(!nd[i]->NodeHasNonrigidParticles()) continue;
		double mnode=1./nd[i]->GetNodalMass(false);
		
		double *vtkquant=vtk[i];
		for(j=0;j<bufferSize;j++)
		{	*vtkquant*=mnode;
			vtkquant++;
		}
	}
}

// Archive VTK file now and free up buffers
void VTKArchive::ExportExtrapolationsToFiles(void)
{	archiver->ArchiveVTKFile(mtime+timestep,quantity,quantitySize,quantityName,qparam,vtk);

	// free buffer if used
	if(vtk!=NULL)
	{	int i;
		for(i=1;i<=nnodes;i++) delete [] vtk[i];
		delete [] vtk;
		vtk=NULL;
	}
}


