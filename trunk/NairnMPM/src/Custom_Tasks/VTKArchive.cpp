/********************************************************************************
    VTKArchive.cpp
    NairnMPM
    
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

#include "Custom_Tasks/VTKArchive.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "Exceptions/CommonException.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"

int dummyArg;

#pragma mark Constructors and Destructors

// Constructors
VTKArchive::VTKArchive()
{
	customArchiveTime=-1.;
	nextCustomArchiveTime=-1.;
	bufferSize=0;
	vtk=NULL;
}

// Return name of this task
const char *VTKArchive::TaskName(void) { return "Archive grid results to VTK files"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *VTKArchive::InputParam(char *pName,int &input)
{
	int q=-1,thisBuffer=0;
	
    // check for archiving quantity
    if(strcmp(pName,"mass")==0)
    {	q=VTK_MASS;
		// no buffer since no need to extrapolate
		thisBuffer=-1;
    }
    
    else if(strcmp(pName,"velocity")==0)
    {	q=VTK_VELOCITY;
		thisBuffer=3;		// extrapolate for cm velocity
    }
	
    else if(strcmp(pName,"stress")==0)
    {	q=VTK_STRESS;
		thisBuffer=6;
    }
	
    else if(strcmp(pName,"strain")==0)
    {	q=VTK_STRAIN;
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
	
    else if(strcmp(pName,"concentration")==0)
    {	q=VTK_CONCENTRATION;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"strainenergy")==0)
    {	q=VTK_STRAINENERGY;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"plasticenergy")==0)
    {	q=VTK_PLASTICENERGY;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"material")==0)
    {	q=VTK_MATERIAL;
		thisBuffer=1;
    }
	
    else if(strcmp(pName,"archiveTime")==0)
    {	input=DOUBLE_NUM;
        return (char *)&customArchiveTime;				// assumes in ms
    }
	
    else if(strcmp(pName,"firstArchiveTime")==0)
    {	input=DOUBLE_NUM;
        return (char *)&nextCustomArchiveTime;			// assumes in ms
    }
	
	// if found one, add to arrays
	if(q>=0)
	{	quantity.push_back(q);
		quantitySize.push_back(thisBuffer);				// <0 is size for non-extrapolated quantities
		if(thisBuffer<0) thisBuffer=-thisBuffer;
		bufferSize+=thisBuffer;
		char *qname=new char[strlen(pName)+1];
		strcpy(qname,pName);
		quantityName.push_back(qname);
		input=INT_NUM;
        return (char *)&dummyArg;
	}
	
	// check remaining commands
    return CustomTask::InputParam(pName,input);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *VTKArchive::Initialize(void)
{
    cout << "Archive grid results to VTK files." << endl;
	
	if(!mpmgrid.IsStructuredGrid())
		throw CommonException("VTKArchive task requires use of a generated grid","VTKArchive::Initialize");
	
	// time interval
	cout << "   Archive time: ";
	if(customArchiveTime>0)
	{	cout << customArchiveTime << " ms";
		customArchiveTime/=1000.;				// convert to sec
		if(nextCustomArchiveTime<0.)
		{	nextCustomArchiveTime=customArchiveTime;
			cout << endl;
		}
		else
		{	cout << ", starting at " << nextCustomArchiveTime << " ms" << endl;
			nextCustomArchiveTime/=1000.;				// convert to sec
		}
	}
	else
		cout << "same as particle archives" << endl;
	
	// quantities
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
		if(q<quantity.size()-1) cout << ", ";
		len+=strlen(name)+2;
	}
	cout << endl;
	
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
CustomTask *VTKArchive::PrepareForStep(bool &needExtraps)
{
	if(customArchiveTime>0.)
	{	if(mtime+timestep>=nextCustomArchiveTime)
		{	doVTKExport=TRUE;
			nextCustomArchiveTime+=customArchiveTime;
		}
		else
			doVTKExport=FALSE;
	}
	else
		doVTKExport=archiver->WillArchive();
	if(quantity.size()==0) doVTKExport=FALSE;
	getVTKExtraps = doVTKExport ? (bufferSize>0) : FALSE;
	if(getVTKExtraps) needExtraps=TRUE;
    return nextTask;
}

// Archive VTK file now
CustomTask *VTKArchive::StepCalculation(void)
{
	if(doVTKExport)
		archiver->ArchiveVTKFile(mtime+timestep,quantity,quantitySize,quantityName,vtk);
    return nextTask;
}

// Called when custom tasks are all done on a step
CustomTask *VTKArchive::FinishForStep(void)
{	// free buffer if used
	if(vtk!=NULL)
	{	int i;
		for(i=1;i<=nnodes;i++) free(vtk[i]);
		free(vtk);
		vtk=NULL;
	}
    return nextTask;
}

#pragma mark TASK EXTRAPOLATION METHODS

// initialize for crack extrapolations
CustomTask *VTKArchive::BeginExtrapolations(void)
{
	if(!getVTKExtraps) return nextTask;
	
	// create buffer for each nodalpoint
	vtk=(double **)malloc((nnodes+1)*sizeof(double *));
	if(vtk==NULL)
	{	cout << "# memory error preparing data for vtk export" << endl;
		getVTKExtraps=FALSE;
		return nextTask;
	}
	int i,j;
	for(i=1;i<=nnodes;i++)
	{	vtk[i]=(double *)malloc(bufferSize*sizeof(double));
		if(vtk[i]==NULL)
		{	for(j=1;j<i;j++) free(vtk[j]);
			free(vtk);
			cout << "# memory error preparing data for vtk export" << endl;
			getVTKExtraps=FALSE;
			return nextTask;
		}
		for(j=0;j<bufferSize;j++) vtk[i][j]=0.;
	}
	
    return nextTask;
}

// add particle data to a node
CustomTask *VTKArchive::NodalExtrapolation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,double wt)
{
	if(!getVTKExtraps) return nextTask;
	
	unsigned int q;
	double *vtkquant=vtk[ndmi->num];
	double theWt=1.;
	Tensor *ten=NULL;
	
	for(q=0;q<quantity.size();q++)
	{	switch(quantity[q])
		{	case VTK_STRESS:
				theWt=wt*theMaterials[mpnt->MatID()]->rho;
				ten=mpnt->GetStressTensor();
			case VTK_STRAIN:
				if(quantity[q]==VTK_STRAIN)
				{	theWt=wt;
					ten=mpnt->GetStrainTensor();
				}
			case VTK_PLASTICSTRAIN:
				if(quantity[q]==VTK_PLASTICSTRAIN)
				{	theWt=wt;
					ten=mpnt->GetPlasticStrainTensor();
				}
				vtkquant[0]+=theWt*ten->xx;
				vtkquant[1]+=theWt*ten->yy;
				vtkquant[2]+=theWt*ten->zz;
				vtkquant[3]+=theWt*ten->xy;
				if(fmobj->IsThreeD())
				{	vtkquant[4]+=theWt*ten->xz;
					vtkquant[5]+=theWt*ten->yz;
				}
				vtkquant+=6;
				break;
			
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
				theWt=wt*theMaterials[mpnt->MatID()]->concSaturation;
				*vtkquant+=theWt*mpnt->pConcentration;
				vtkquant++;
				break;
				
			case VTK_STRAINENERGY:
				*vtkquant+=wt*1.0e-6*mpnt->mp*mpnt->GetStrainEnergy();
				vtkquant++;
				break;
				
			case VTK_PLASTICENERGY:
				*vtkquant+=wt*1.0e-6*mpnt->mp*mpnt->GetPlastEnergy();
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

// initialize for crack extrapolations
CustomTask *VTKArchive::EndExtrapolations(void)
{
	if(!getVTKExtraps) return nextTask;
	
	// divide all by nodal mass
	int i,j;
    for(i=1;i<=nnodes;i++)
	{	if(nd[i]->NumberParticles()==0) continue;
		double mnode=1./nd[i]->mass;
		
		double *vtkquant=vtk[i];
		for(j=0;j<bufferSize;j++)
		{	*vtkquant*=mnode;
			vtkquant++;
		}
	}

    return nextTask;
}
        
