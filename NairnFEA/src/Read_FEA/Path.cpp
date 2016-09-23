/********************************************************************************
    PathsController.cpp
    NairnFEA
    
    Created by John Nairn on 6/23/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_FEA/Path.hpp"
#include "Read_FEA/KeypointsController.hpp"
#include "Read_FEA/Keypoint.hpp"
#include "Read_FEA/EdgeBCController.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"

// used in ConvertToRatio() and associated routines
double ratio2,argx0,argx1,argx2,argy0,argy1,argy2,comparea;

/********************************************************************************
	Path: Constructors and Destructor
********************************************************************************/

// Initialize, but assumes name already 32 characters or less and unique
// ... which is done by prior call to paths->ValidName()
// ... but internally created paths will have empty name
Path::Path(const char *pathName,int intvls,double rat)
{
	strcpy(pathID,pathName);
	intervals=intvls;
	ratio=rat;
	numKeypts=0;
	face=0;
	firstMainNode=-1;
	firstMidsideNode=-1;
	nodeIncrement=1;
	firstElem=-1;
	bcFlags[0]=bcFlags[1]=bcFlags[2]=FALSE;
	subPath[0]=subPath[1]=NULL;
}

/********************************************************************************
	Path: methods
********************************************************************************/

// add key point to the last path
int Path::AddKeypoint(char *keyName)
{
	// someday should allow more than 3 points and use spline paths
	if(numKeypts>=3) return FALSE;
	
	// find the keypoint object
	keys[numKeypts]=keyPts->FindKeypoint(keyName);
	if(keys[numKeypts]==NULL) return FALSE;
	
	// make sure keypoint not already in the path
	int i;
	for(i=0;i<numKeypts;i++)
	{	if(keys[i]==keys[numKeypts])
			return FALSE;
	}
	numKeypts++;
	
	return TRUE;
}

// assign edge boundary conditions to a path
// throws std::bad_alloc
void Path::AddEdgeBCsToPath(int dir,int nloads,double *stress)
{
	int i,j;
	
	// get three loads
	if(nloads==1)
		stress[1]=stress[2]=stress[0];
	else if(nloads==2)
	{	stress[2]=stress[1];
		stress[1]=(stress[0]+stress[2])/2.;
	}
	
	// Get line parametrization parameters
	double r,a,t=-1.,last;
	GetLineParameters(intervals,ratio,&r,&a);
	last=a;
	
	// element parameters
	int elementInc,switchPaths;
	int element=firstElem;
	if(element>=0)
	{	// use this path
		elementInc=elemIncrement;
		switchPaths=intervals+1;		// never switch paths
	}
	else
	{	// must look at sub paths
		element=subPath[0]->firstElem;
		elementInc=subPath[0]->elemIncrement;
		switchPaths=subPath[0]->intervals;		// when to switch to other subpath
	}
	
	// loop over intervals
	int faceLoc = (face<0) ? -face : face;		// negative face is when interface is present too
	double tn[3],values[3];
	for(i=1;i<=intervals;i++)
	{	tn[0]=t;
		t+=last;
		last*=r;
		tn[1]=(tn[0]+t)/2.;
		tn[2]=t;
		for(j=0;j<=2;j++)
		{	values[j]=0.5*(2*stress[1]+(stress[2]-stress[0])*tn[j]
				+(stress[2]+stress[0]-2*stress[1])*tn[j]*tn[j]);
		}
		
		// create new edge BC
		EdgeBC *newEdgeBC=new EdgeBC(element,faceLoc,dir);
		if(edgeBCCtrl->AddObject(newEdgeBC))
		{	if(firstMidsideNode>0)
				newEdgeBC->SetStress(values,3);
			else
			{	values[1]=values[2];
				newEdgeBC->SetStress(values,2);
			}
		}

		// should loop by switched to second subpath?
		if(i!=switchPaths)
			element+=elementInc;
		else
		{	element=subPath[1]->firstElem;
			elementInc=subPath[1]->elemIncrement;
		}
	}
}

// Reorient path prior to meshing an area
void Path::ReorientPath(void)
{
	int i,halfNum=numKeypts>>1,swapAddr;
	Keypoint *tempKey;
	
	// swap keypoints
	for(i=0;i<halfNum;i++)
	{	swapAddr=numKeypts-1-i;
		tempKey=keys[i];
		keys[i]=keys[swapAddr];
		keys[swapAddr]=tempKey;
	}
	
	// change ratio
	ratio=1./ratio;
	
	// node numbers (if there)
	if(firstMainNode>0)
		firstMainNode+=(intervals-2)*nodeIncrement;
	if(firstMidsideNode>0)
		firstMidsideNode+=(intervals-1)*nodeIncrement;
	if(firstMainNode>=0)
		nodeIncrement = (nodeIncrement==1) ? -1 : 1;
}

// return node along the path from 1 to intervals+1 (if no midside nodes)
//    or 1 to 2*intervals+1 if midside nodes
// return 0 if out of range or no nodes yet
int Path::nodeAtIndex(int nodeIndex)
{
	if(!IsMeshed()) return 0;
	
	// first key point
	if(nodeIndex<1) return 0;
	if(nodeIndex==1) return keys[0]->node;
	
	// find total number
	int numNodes = (firstMidsideNode>0) ? 2*intervals+1 : intervals+1;
	
	// last keypoint
	if(nodeIndex==numNodes) return keys[numKeypts-1]->node;
	if(nodeIndex>numNodes) return 0;
	
	// interior
	if(firstMidsideNode<0)
		return firstMainNode+(nodeIndex-2)*nodeIncrement;
	else
	{	if(IsEven(nodeIndex))
			return firstMidsideNode+((nodeIndex>>1)-1)*nodeIncrement;
		else
			return firstMainNode+(((nodeIndex-1)>>1)-1)*nodeIncrement;
	}
}

// check if allowed to set BC option on this path
const char *Path::VerifyBCOption(int command,int axis)
{
	// not a command
	if(command<0) return NULL;
	
	// if not meshed, can not set BC
	if(!IsMeshed())
		return "Can't set boundary conditions or rotate a path that is not in the mesh";
	
	// if path used twice, cannot use
	if(face==INTERNAL_PATH)
		return "Can't set boundary conditions or rotate an interior path";
	
	// no more checking if stress on edges
	if(command==LOAD_PATH_EDGES) return NULL;
	
	// can't set same DOF twice (dof=0 is for rotations)
	if(bcFlags[axis])
		return "Can't set same boundary condition direction twice or rotate twice on the same path";
	bcFlags[axis]=TRUE;
	
	if(command==ROTATE_PATH_NODES && (bcFlags[1] || bcFlags[2]))
		return "Can't rotate path that already has boundary condition set";
	
	return NULL;
}

// create node at the first point
void Path::FirstKeypointToNode(void)
{
	keys[0]->AssignToNode();
}

// Convert absolute distance input in ratio into a path ratio
void Path::ConvertToRatio(void)
{
	Vector lineNd[3];
	double a,equalInts,r,upper;
		
	// square of distance
	ratio2=ratio;
	ratio2*=ratio2;
	
	/* Get coordinates along line
		lineNd[0,1,2] are coordinate of three points along path */
	FirstKeypoint()->GetKeypointXY(&lineNd[0]);
	LastKeypoint()->GetKeypointXY(&lineNd[2]);
	if(numKeypts>2)
	{	keys[1]->GetKeypointXY(&lineNd[1]);
		
		/* Coefficients for quadratic parametrization of curve
				x = x[0] + argx0 + argx1*t + argx2*t^2
				y = y[0] + argy0 + argy1*t + argy2*t^2
			for t=-1 to 1
			
			Size of first element = sqrt((x-x[0])^ + (y-y[0])^2)
		*/
		argx0=lineNd[1].x-lineNd[0].x;
		argy0=lineNd[1].y-lineNd[0].y;
		argx1=(lineNd[2].x-lineNd[0].x)/2.;
		argy1=(lineNd[2].y-lineNd[0].y)/2.;
		argx2=(lineNd[2].x+lineNd[0].x)/2.-lineNd[1].x;
		argy2=(lineNd[2].y+lineNd[0].y)/2.-lineNd[1].y;
		
		// Find desired dimensionless a or absolute distance from -1 point
		a=1.+BinarySolve((double)(-1.),(double)1.,30,1);
	}
	else
	{	lineNd[1].x=(lineNd[0].x+lineNd[2].x)/2.;
		lineNd[1].y=(lineNd[0].y+lineNd[2].y)/2.;
		
		argx1=(lineNd[2].x-lineNd[0].x)/2.;
		argy1=(lineNd[2].y-lineNd[0].y)/2.;
		
		a=ratio/(sqrt(argx1*argx1+argy1*argy1));
	}
	
	/* Find r such that a + ar + ar^2 + ... a r^(n-1) = 2
		A. If a smaller than equal intervals, 1<r<upper (find this)
		B. If a bigger than equal intervals, 0<r<1
	*/
	equalInts=2./(double)intervals;
	comparea=2./a;
	if(a<equalInts)
	{	upper=2.;
		while(NormalSum(upper)<0) upper*=2.;
		r=BinarySolve((double)1.,(double)upper,30,2);
	}
	else if(a>equalInts)
		r=BinarySolve((double)0.,(double)1.,30,2);
	else
		r=1.;

	// find ratio of first to last element a/(a r^(n-1))
	ratio=1./(pow(r,(double)(intervals-1)));
}

// Binary Solution used for ConvertToRatio
double Path::BinarySolve(double smin,double smax,int times,short callID)
{
	double val,smid;
	int i;
	
	for(i=1;i<=times;i++)
	{	smid=(smin+smax)/2.;
		switch(callID)
		{	case 1:
				val=Distance(smid);
				break;
			case 2:
				val=NormalSum(smid);
				break;
			default:
				val=0.;
				break;
		}
		if(val<0.)
			smin=smid;
		else if(val>0.)
			smax=smid;
		else
			return smid;
	}
	return (smin+smax)/2.;
}

// Find distance for some bunary solutions
double Path::Distance(double a)
{
	double xdist=argx0+argx1*a+argx2*a*a;
	double ydist=argy0+argy1*a+argy2*a*a;
	return xdist*xdist+ydist*ydist-ratio2;
}

// Find sum of geometric series
double Path::NormalSum(double r)
{
	double sum;
	sum=(1.-pow(r,(double)intervals))/(1.-r);
	return sum-comparea;
}

/********************************************************************************
	Path: accessors
********************************************************************************/

// must be called before path used as it converts ratio<0 to actual ratio
int Path::ValidPath(void)
{
	if(numKeypts<2) return FALSE;
	if(ratio<0)
	{	ratio=-ratio;
		ConvertToRatio();
	}
	return TRUE;
}
Keypoint *Path::FirstKeypoint(void) { return numKeypts>0 ? keys[0] : NULL; }
Keypoint *Path::LastKeypoint(void) { return numKeypts>0 ? keys[numKeypts-1] : NULL; }
void Path::GetFirstAndMidXY(Vector *firstKey,Vector *midKey)
{
	keys[0]->GetKeypointXY(firstKey);
	if(numKeypts==2)
	{	Vector endKey;
		keys[1]->GetKeypointXY(&endKey);
		midKey->x=(firstKey->x+endKey.x)/2.;
		midKey->y=(firstKey->y+endKey.y)/2.;
	}
	else
		keys[1]->GetKeypointXY(midKey);
}
bool Path::IsMeshed(void) { return firstMainNode>=0; }
		
// add any non-NULL key to list, but return false if now too many
int Path::SetKeys(Keypoint *key1,Keypoint *key2,Keypoint *key3)
{
	if(key1!=NULL)
	{	if(numKeypts>=3) return FALSE;
		keys[numKeypts++]=key1;
	}
	if(key2!=NULL)
	{	if(numKeypts>=3) return FALSE;
		keys[numKeypts++]=key2;
	}
	if(key3!=NULL)
	{	if(numKeypts>=3) return FALSE;
		keys[numKeypts++]=key3;
	}
	return TRUE;
}
	
