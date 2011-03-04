/******************************************************************************** 
    Linear2D.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/Linear2D.hpp"
#include "Nodes/NodalPoint.hpp"

#pragma mark Linear2D::Constructors and Destructor

#ifdef MPM_CODE
// Main MPM constructor passes onto ElementBase constructor
Linear2D::Linear2D(int eNum,int *eNode) : ElementBase(eNum,eNode) {}
#else
// Main FEA constructor passes onto ElementBase constructor
Linear2D::Linear2D(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
                    ElementBase(eNum,eNode,eMat,eAng)
{ thickness=eThick; }
#endif

#pragma mark Linear2D::Methods

/*	Calculate area of element (in mm^2 because nodes in mm)
	node is pointer to 0-based array NodalPoints
	Assumes nodes[NumberNodes()] = nodes[0]
*/
double Linear2D::GetArea(void)
{
    int i,ns=NumberSides();
    double area;
    
    area=0.;
    for(i=0;i<ns;i++)
    {	area+=nd[nodes[i]]->x*nd[nodes[i+1]]->y
                        -nd[nodes[i]]->y*nd[nodes[i+1]]->x;
    }
    return area/2.;
}

//	Calculate volume of element (in mm^3 because nodes in mm)
double Linear2D::GetVolume(void) { return thickness*GetArea(); }

/*	Use ray crossing algorithm to find out if a point (pt.x,pt.y) is in an element
	nd is pointer to 1-based array NodalPoints
	Assumes nodes[NumberNodes()] = nodes[0]
*/
short Linear2D::PtInElement(Vector &pt)
{
    int i,ns=NumberSides(),crossings=0;
    double d,x1,y1,x2,y2;
    
    for(i=0;i<ns;i++)
    {	x1=nd[nodes[i]]->x;
        y1=nd[nodes[i]]->y;
        x2=nd[nodes[i+1]]->x;
        y2=nd[nodes[i+1]]->y;
        d=(pt.y-y1)*(x2-x1) - (pt.x-x1)*(y2-y1);
        
        // get crossing unless both y's on same side of egde
        if((y1>=pt.y) != (y2>=pt.y))
        {   crossings+= (y2-y1>=0.) ? d>=0. : d<=0. ;
        }
        
        // if d is 0, check if point is on line (and thus in polygon
        if(!d && fmin(x1,x2)<=pt.x && pt.x<=fmax(x1,x2) &&
                            fmin(y1,y2)<=pt.y && pt.y<=fmax(y1,y2))
        {   return 1;
        }
    }
    return crossings & 0x01;
}

#ifdef FEA_CODE
// Calculate edge load
void Linear2D::CalcEdgeLoads(double *re,int iedge,int ndir,double *fload,int np)
{
	int nd2;
	
	if(iedge==NumberSides())
		nd2=1;
	else
		nd2=iedge+1;
	LinearEdgeLoad(iedge,nd2,ndir,fload,re,np);
}
#endif

#pragma mark Linear2D::Acessors

// thickness which may be in a subclass
double Linear2D::GetThickness(void) { return thickness; }

// face nodes
int Linear2D::FaceNodes(void) { return 2; }

