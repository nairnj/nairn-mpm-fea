/********************************************************************************
    ShapeController.cpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	This is base class for defining shapes to assigning boundary conditions
		in both MPM and FEA. It is never used as is, but always subclassed
		for a shape. The process is:
	
	1. Define new XML command for the shape with enough attributes to define
		the shape. In code for ne XML command, create new shape and set block to
		BCSHAPE. Read paramters using SetProperty() whenever possible
	
	2. Subclass needs only the following:
		a. Constructor to store enough parameters to define the shape (base
			class has xmin,xmax,ymin,ymax, and tolerance for use) or
			handle parameters in SetProperty().
		b. The function PtOnShape(Vector v) to decide if the point v is
			contained by the shape and should be assigned a BC.
		c. The optional FinishSetup() can be called after setting all
			parameters in case helpful to the object.
	
	3. Special case: the above works for any BC on nodes or particles and is
		found by enumerating over them all, depending on the type of BC. If the
		shape should not enumerate over all (e.g., a Path in FEA which enumerates
		only over nodes on that path, a single node, which enumerates over
		only that node), the class may instead overide
		a. nextNode() - return node number by index nodeNum (starting at 1)
		b. nextParticle() - return particle number by index particleNum (starting at 0)
		This type of subclass will not need PtOnShape() unless the overridden
		nextNode() or nextParticle() use it.
	
	4. For preventing use as particleBC, overide nextParticle() to return -1
	
	5. For preventing use as nodal BC, override nextNode() to return 0
	
	6. A special function GetContextInfo() returns a generic char * and NULL by
		default. Subclasses can override and return an object needed to set
		BCs.
********************************************************************************/

#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/ShapeController.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "MPM_Classes/MPMBase.hpp"
#endif

ShapeController *theShape=NULL;

/********************************************************************************
	ShapeController: Constructors and Destructor
********************************************************************************/

ShapeController::ShapeController(int block)
{
	sourceBlock=block;
	distScaling=1.;
	xmin=0.;
	xmax=-1.;
	ymin=0.;
	ymax=-1.;
	zmin=0.;
	zmax=-1.;
	nodeNum=1;
#ifdef MPM_CODE
	particleNum=0;
	numParticles=0;
#endif
}

ShapeController::ShapeController(int block,double x1,double x2,double y1,double y2,double tolerance)
{
	sourceBlock=block;
	distScaling=1.;
	xmin=x1;
	xmax=x2;
	ymin=y1;
	ymax=y2;
	nodeNum=1;
#ifdef MPM_CODE
	particleNum=0;
	numParticles=0;
#endif
}

ShapeController::~ShapeController() { }

/********************************************************************************
	ShapeController: methods
********************************************************************************/

// Deterime if on the shape (depending of the type of shape) 
bool ShapeController::PtOnShape(Vector v) { return FALSE; }

// the source block
int ShapeController::GetSourceBlock(void) { return sourceBlock; }
bool ShapeController::RequiredBlock(int block) { return block==sourceBlock; }

// reset nodeNum and return NULL (no errors except in other shapes)
void ShapeController::resetNodeEnumerator(void) { nodeNum=1; }

#ifdef MPM_CODE
// reset nodeNum and return NULL (no errors except in other shapes)
void ShapeController::resetParticleEnumerator(void) { particleNum=0; }
#endif

// reset nodeNum and return NULL (no errors except in other shapes)
const char *ShapeController::startNodeEnumerator(int command,int axis)
{	nodeNum=1;
	return NULL;
}

// return next node for this shape or 0 if no more
int ShapeController::nextNode(void)
{
	if(nodeNum>nnodes) return 0;
	int i;
	for(i=nodeNum;i<=nnodes;i++)
	{   if(PtOnShape(MakeVector(nd[i]->x,nd[i]->y,nd[i]->z)))
		{	nodeNum=i+1;
			return i;
		}
	}
	nodeNum=nnodes+1;
	return 0;
}

#ifdef MPM_CODE
// return next node for this shape or -1 if no more
int ShapeController::nextParticle(void)
{
	if(particleNum>=nmpms) return -1;
	int i;
	for(i=particleNum;i<nmpms;i++)
	{   if(PtOnShape(MakeVector(mpm[i]->pos.x,mpm[i]->pos.y,mpm[i]->pos.z)))
		{	particleNum=i+1;
			return i;
		}
	}
	particleNum=nmpms;
	return -1;
}

void ShapeController::setNetBC(bool setting)
{	// turn off is false
	if(!setting)
	{	numParticles=0;
		return;
	}
	
	// if already counted, then can exit
	if(numParticles>0) return;
	
	// count them
	int i;
	resetParticleEnumerator();
	while((i=nextParticle())>=0) numParticles++;
}

// return num of particles in current shape or 1 if on perParticle basis
double ShapeController::particleCount(void) { return numParticles>0 ? (double)numParticles : (double)1. ; }

#endif

// return no path
char *ShapeController::GetContextInfo(void) { return NULL; }

// set the scaling
void ShapeController::SetScaling(double scale) { distScaling=scale; }

// set a property
void ShapeController::SetProperty(char *aName,char *value,CommonReadHandler *reader)
{	if(strcmp(aName,"x1")==0 || strcmp(aName,"xmin")==0)
	{	xmin=reader->ReadX(value,distScaling);
	}
	else if(strcmp(aName,"y1")==0 || strcmp(aName,"ymin")==0)
	{	ymin=reader->ReadY(value,distScaling);
	}
	else if(strcmp(aName,"z1")==0 || strcmp(aName,"zmin")==0)
	{	zmin=reader->ReadZ(value,distScaling);
	}
	else if(strcmp(aName,"x2")==0 || strcmp(aName,"xmax")==0)
	{	xmax=reader->ReadX(value,distScaling);
	}
	else if(strcmp(aName,"y2")==0 || strcmp(aName,"ymax")==0)
	{	ymax=reader->ReadY(value,distScaling);
	}
	else if(strcmp(aName,"z2")==0 || strcmp(aName,"zmax")==0)
	{	zmax=reader->ReadZ(value,distScaling);
	}
}

// set a property
void ShapeController::SetProperty(char *aName,double value)
{
	if(strcmp(aName,"x1")==0 || strcmp(aName,"xmin")==0)
		xmin=value*distScaling;
	else if(strcmp(aName,"y1")==0 || strcmp(aName,"ymin")==0)
		ymin=value*distScaling;
	else if(strcmp(aName,"z1")==0 || strcmp(aName,"zmin")==0)
		zmin=value*distScaling;
	else if(strcmp(aName,"x2")==0 || strcmp(aName,"xmax")==0)
		xmax=value*distScaling;
	else if(strcmp(aName,"y2")==0 || strcmp(aName,"ymax")==0)
		ymax=value*distScaling;
	else if(strcmp(aName,"z2")==0 || strcmp(aName,"zmax")==0)
		zmax=value*distScaling;
}

// called after initialization is done
void ShapeController::FinishSetup(void) { }


