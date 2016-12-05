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
		b. The function ContainsPoint(Vector& v) to decide if the point v is
			contained by the shape and should be assigned a BC. The parent
			function (non-virtual) ShapeContainsPoint(Vector &v) decides
			if point in shape and not in any subordinate shapes
		c. The optional FinishSetup() can be called after setting all
			parameters in case helpful to the object.
	
	3. Special case: the above works for any BC on nodes or particles and is
		found by enumerating over them all, depending on the type of BC. If the
		shape should not enumerate over all (e.g., a Path in FEA which enumerates
		only over nodes on that path, a single node, which enumerates over
		only that node), the class may instead overide
		a. nextNode() - return node number by index nodeNum (starting at 1)
		b. nextParticle() - return particle number by index particleNum (starting at 0)
		This type of subclass will not need ContainsPoint() unless the overridden
		nextNode() or nextParticle() use it.
	
	4. For preventing use as particleBC, overide nextParticle() to return -1
	
	5. For preventing use as nodal BC, override nextNode() to return 0
	
	6. A special function GetContextInfo() returns a generic char * and NULL by
		default. Subclasses can override and return an object needed to set
		BCs.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/ShapeController.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Elements/ElementBase.hpp"
#ifdef MPM_CODE
	#include "MPM_Classes/MPMBase.hpp"
#endif

ShapeController *theShape=NULL;

#pragma mark ShapeController: constructors, destructors, and initializes

ShapeController::ShapeController(int block)
{
	sourceBlock=block;
	distScaling=1.;
	xmin=0.;
	xmax=0.;
	ymin=0.;
	ymax=0.;
	zmin=0.;
	zmax=0.;
	nodeNum=1;
    elemNum=0;
#ifdef MPM_CODE
	particleNum=0;
	numParticles=0;
#endif
	parentShape = NULL;
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
    elemNum=0;
#ifdef MPM_CODE
	particleNum=0;
	numParticles=0;
#endif
}

ShapeController::~ShapeController()
{
	// delete children shapes
	for(int i=0;i<children.size();i++)
	{	delete children[i];
	}
}

// set a property on reading for x, y, z, min and max
void ShapeController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
    if(strcmp(aName,"x1")==0 || strcmp(aName,"xmin")==0)
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

// set a property from value without read handler
void ShapeController::SetProperty(const char *aName,double value)
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

// to allow object to decode object-specific character data
// throw an exception if bad data
void ShapeController::SetProperty(char *bData,CommonReadHandler *reader) {}

// set the scaling
void ShapeController::SetScaling(double scale) { distScaling=scale; }

// set an object parameter (in subordinate command)
// called for attributes on XML objects subordinate to the shape command
void ShapeController::SetParameter(const char *aName,const char *value) { }

// called after finish attributes of subordinate command
// return FALSE if not set correctly, or TRUE is OK to continue
bool ShapeController::FinishParameter(void) { return true; }

// called after initialization is done, return TRUE if ready to use
// or FALSE if this object needs to wait for parameters
// This base class requires min and max (x, y and z) to differ and
//      reorders if needed. This it correct for rect, oval, box, sphere
//      and cylinder, but maybe not for others.
// throws SAXException()
bool ShapeController::FinishSetup(void)
{
    double temp;
    if(xmin>xmax)
	{	temp=xmax;
        xmax=xmin;
        xmin=temp;
    }
    if(DbleEqual(xmin,xmax))
        ThrowSAXException("%s: xmax cannot equal xmin in input parameters.",GetShapeName());
	
    if(ymin>ymax)
    {	temp=ymax;
        ymax=ymin;
        ymin=temp;
    }
    if(DbleEqual(ymin,ymax))
        ThrowSAXException("%s: ymax cannot equal ymin in input parameters.",GetShapeName());
	
    if(!Is2DShape())
    {	if(zmin>zmax)
        {	temp=zmax;
            zmax=zmin;
            zmin=temp;
        }
        if(DbleEqual(zmin,zmax))
            ThrowSAXException("%s: zmax cannot equal zmin in input parameters.",GetShapeName());
    }
	
	return true;
}

// some shapes might call this right be fore use. Return TRUE or FALSE
// if has all parameters. Normally only for shapes with subordinate commands.
bool ShapeController::HasAllParameters(void) { return TRUE; }

#pragma mark ShapeController: methods

// Determine if on the shape (depending of the type of shape)
bool ShapeController::ContainsPoint(Vector& v) { return FALSE; }

// Determine if on the shape (depending of the type of shape)
bool ShapeController::ShapeContainsPoint(Vector& v)
{	// check shape
	if(!ContainsPoint(v)) return false;
	
	// check subordinate cutouts, and false if any of them
	for(int i=0;i<children.size();i++)
	{	if(children[i]->ShapeContainsPoint(v))
			return false;
	}
	
	// finally true
	return true;
}

// reset nodeNum and return NULL (no errors except in other shapes)
void ShapeController::resetNodeEnumerator(void) { nodeNum=1; }

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
    {   Vector nv = MakeVector(nd[i]->x,nd[i]->y,nd[i]->z);
	    if(ShapeContainsPoint(nv))
		{	nodeNum=i+1;
			return i;
		}
	}
	nodeNum=nnodes+1;
	return 0;
}

// reset nodeNum and return NULL (no errors except in other shapes)
void ShapeController::resetElementEnumerator(void) { elemNum=0; }

// return next node for this shape or 0 if no more
int ShapeController::nextElement(void)
{
	if(elemNum>=nelems) return -1;
	int i;
    Vector ev;
    for(i=elemNum;i<nelems;i++)
    {   theElements[i]->GetXYZCentroid(&ev);
	    if(ShapeContainsPoint(ev))
		{	elemNum=i+1;
			return i;
		}
	}
	elemNum=nelems;
	return -1;
}

// Add a cutout shape
void ShapeController::AddCutoutShape(ShapeController *cutout)
{	// add child shape
	children.push_back(cutout);
}

#pragma mark ShapeController: MPM only methods

#ifdef MPM_CODE
// return next node for this shape or -1 if no more
int ShapeController::nextParticle(void)
{
	if(particleNum>=nmpms) return -1;
	int i;
	for(i=particleNum;i<nmpms;i++)
    {   Vector nv = MakeVector(mpm[i]->pos.x,mpm[i]->pos.y,mpm[i]->pos.z);
	    if(ShapeContainsPoint(nv))
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

// reset nodeNum and return NULL (no errors except in other shapes)
void ShapeController::resetParticleEnumerator(void) { particleNum=0; }

#endif

#pragma mark ShapeController: accessors

// type of object - used in some error messages
const char *ShapeController::GetShapeName(void) { return "Shape"; }

// type of object - used in some error messages
void ShapeController::DescribeShape(const char *prefix)
{	cout << prefix << "Shape: " << GetShapeName() << " (id: " << this << ")" << endl;
	cout << prefix << "   x range: " << xmin << " to " << xmax << endl;
	cout << prefix << "   y range: " << ymin << " to " << ymax << endl;
	if(!Is2DShape())
		cout << prefix << "   z range: " << zmin << " to " << zmax << endl;
	
	// cutouts
	if(children.size()>0)
	{	char cutPrefix[200];
		strcpy(cutPrefix,prefix);
		strcat(cutPrefix,"   ");
		for(int i=0;i<children.size();i++)
			children[i]->DescribeShape(cutPrefix);
	}
}

// override for 3D shapes and result false
bool ShapeController::Is2DShape(void) { return true; }

// the source block
int ShapeController::GetSourceBlock(void) { return sourceBlock; }

// check if sourceBlock is same as a required block (TRUE or FALSE)
bool ShapeController::RequiredBlock(int block) { return block==sourceBlock; }

// return no pointer
char *ShapeController::GetContextInfo(void) { return NULL; }

// parent shape
ShapeController *ShapeController::GetParentShape(void) const { return parentShape; }
void ShapeController::SetParentShape(ShapeController *obj) { parentShape=obj; }






