/******************************************************************************** 
    Quad2D.cpp
    NairnFEA
    
    Created by John Nairn on Fri Oct 22 2002.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/Quad2D.hpp"
#include "Nodes/NodalPoint.hpp"

#pragma mark Quad2D: Constructors and Destructor

#ifdef MPM_CODE

// Main MPM constructor passes onto ElementBase constructor
Quad2D::Quad2D(int eNum,int *eNode) : ElementBase(eNum,eNode) {}

#else

// Main FEA constructor passes onto ElementBase constructor
Quad2D::Quad2D(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
			ElementBase(eNum,eNode,eMat,eAng)
{	thickness=eThick;
}

#endif

#pragma mark Quad2D: methods

/*	Use ray crossing algorithm to find out if a point (pt.x,pt.y) is in an element
	nd is pointer to 1-based array NodalPoints
*/
short Quad2D::PtInElement(Vector &pt)
{
	short i,ns=NumberSides(),crossings=0;
	double d,x1,y1,x2,y2;
	
	for(i=0;i<ns;i++)
	{	// start to mid side node
		x1=nd[nodes[i]]->x;
		y1=nd[nodes[i]]->y;
		x2=nd[nodes[i+ns]]->x;
		y2=nd[nodes[i+ns]]->y;
		d=(pt.y-y1)*(x2-x1) - (pt.x-x1)*(y2-y1);
		
		// get crossing unless both y's on same side of edge
		if((y1>=pt.y) != (y2>=pt.y))
		{	crossings+= (y2-y1>=0.) ? d>=0. : d<=0. ;
		}
		
		// if d is 0, check if point is on line (and thus in polygon)
		if(!d && fmin(x1,x2)<=pt.x && pt.x<=fmax(x1,x2) &&
				 fmin(y1,y2)<=pt.y && pt.y<=fmax(y1,y2))
		{	return(1);
		}
		
		// mid side node to end
		x1=x2;
		y1=y2;
		if(i+1<ns)
		{	x2=nd[nodes[i+1]]->x;
			y2=nd[nodes[i+1]]->y;
		}
		else
		{	x2=nd[nodes[0]]->x;
			y2=nd[nodes[0]]->y;
		}
		d=(pt.y-y1)*(x2-x1) - (pt.x-x1)*(y2-y1);
		
		// get crossing unless both y's on same side of egde
		if((y1>=pt.y) != (y2>=pt.y))
		{	crossings+= (y2-y1>=0.) ? d>=0. : d<=0. ;
		}
		
		// if d is 0, check if point is on line (and thus in polygon
		if(!d && fmin(x1,x2)<=pt.x && pt.x<=fmax(x1,x2) &&
				 fmin(y1,y2)<=pt.y && pt.y<=fmax(y1,y2))
		{	return(1);
		}
	}
	return crossings & 0x01;
}

#ifdef FEA_CODE

// Calculate edge loads
void Quad2D::CalcEdgeLoads(double *re,int iedge,int ndir,double *fload,int np)
{	// iedge is edge number (1 ... NumberSides())
	int ind1=iedge,ind2,ind3;
	
	// assumes corners at 1 ... NumberSides() and midside nodes at
	//    NumberSides+1 ... 2*NumberSides()
	ind2=ind1+NumberSides();
	if(ind1==NumberSides())
		ind3=1;
	else
		ind3=ind1+1;
	QuadEdgeLoad(ind1,ind2,ind3,ndir,fload,re,np);
}

// If this node has crack tip nodes, move neighboring nodes toward the crack tip
void Quad2D::MakeQuarterPointNodes(int crackTip,vector<int> &movedNodes)
{
	// check real corner nodes (in nodes[0] to nodes[NumberSides()]
	int i,ct=-1;;
    for(i=0;i<NumberSides();i++)
    {	if(nodes[i]==crackTip)
		{	ct=i;
			break;
		}
    }
	if(ct<0) return;
	
	// edge after ct (assumes midside nodes and ct+NumberSides())
	// note hav internal nodes after last midside node as 2*NumberSides()-1 (zero based)
	if(ct<NumberSides()-1)
		AdjustMidSideNode(ct,ct+NumberSides(),ct+1,movedNodes);
	else
		AdjustMidSideNode(NumberSides()-1,2*NumberSides()-1,0,movedNodes);
	
	// edge before ct
	if(ct>0)
		AdjustMidSideNode(ct,ct-1+NumberSides(),ct-1,movedNodes);
	else
		AdjustMidSideNode(0,2*NumberSides()-1,NumberSides()-1,movedNodes);
}

// move mid side node to quarter point, but only if was not moved before
// not thread safe due to push_back()
void Quad2D::AdjustMidSideNode(int tip,int mid,int end,vector<int> &movedNodes)
{
	// check if node was already moved
	int nd2=nodes[mid];
	unsigned i;
	for(i=0;i<movedNodes.size();i++)
	{	if(movedNodes[i]==nd2) return;
	}
	
	// mark this node now as moved
	movedNodes.push_back(nd2);
	
	// 1D extrapolation for these three nodes
	int nd1=nodes[tip],nd3=nodes[end];
	double dx,dmidx;
	
	dx=(nd[nd3]->x-nd[nd1]->x)/2.;
	dmidx=(nd[nd3]->x+nd[nd1]->x)/2. - nd[nd2]->x;
	nd[nd2]->x += (-0.5*dx + 0.25*dmidx);
	
	dx=(nd[nd3]->y-nd[nd1]->y)/2.;
	dmidx=(nd[nd3]->y+nd[nd1]->y)/2. - nd[nd2]->y;
	nd[nd2]->y += (-0.5*dx + 0.25*dmidx);
}

#endif

#pragma mark Quad2D: accessors

// thickness which may be in a subclass
double Quad2D::GetThickness(void) const { return thickness; }
void Quad2D::SetThickness(double thick) { thickness = thick; }

// face nodes
int Quad2D::FaceNodes(void) { return 3; }

/*	Calculate area of element (in mm^2 because nodes in mm)
	nodes is pointer to 0-based array NodalPoints
*/
double Quad2D::GetArea(void) const
{
	short i,ns=NumberSides(),nn=NumberNodes()-1;
	double area;
	
	area=0.;
	for(i=0;i<ns-1;i++)
	{	area+=nd[nodes[i]]->x*nd[nodes[i+ns]]->y
					-nd[nodes[i]]->y*nd[nodes[i+ns]]->x
					+nd[nodes[i+ns]]->x*nd[nodes[i+1]]->y
					-nd[nodes[i+ns]]->y*nd[nodes[i+1]]->x;
	}
	ns--;
	area+=nd[nodes[ns]]->x*nd[nodes[nn]]->y
				-nd[nodes[ns]]->y*nd[nodes[nn]]->x
				+nd[nodes[nn]]->x*nd[nodes[0]]->y
				-nd[nodes[nn]]->y*nd[nodes[0]]->x;
	return area/2.;
}

//	Calculate volume of element (in mm^3 because nodes in mm)
double Quad2D::GetVolume(void) const { return thickness*GetArea(); }

