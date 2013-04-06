/******************************************************************
	EightNodeIsoparamBrick.cpp
	NairnMPM

	Created by John Nairn on 7/20/06.
	Copyright 2006 RSAC Software. All rights reserved.
 ******************************************************************/

#include "Elements/EightNodeIsoparamBrick.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
#endif

// Local globals
static double xii[8]={-1.,1.,1.,-1.,-1.,1.,1.,-1.};
static double eti[8]={-1.,-1.,1.,1.,-1.,-1.,1.,1.};
static double zti[8]={-1.,-1.,-1.,-1.,1.,1.,1.,1.};

// globals for GIMP node locations
#ifdef MPM_CODE
static double gxii[64]={-1.,1.,1.,-1.,-1.,1.,1.,-1.,
						-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.,
						-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.,
						-1.,1.,1.,-1.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.,
						-1.,1.,1.,-1.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.};
static double geti[64]={-1.,-1.,1.,1.,-1.,-1.,1.,1.,
						-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,
						-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,
						-1.,-1.,1.,1.,-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,
						-1.,-1.,1.,1.,-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.};
static double gzti[64]={-1.,-1.,-1.,-1.,1.,1.,1.,1.,
						-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,
						1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
						-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,
						3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
#endif

#pragma mark EightNodeIsoparamBrick::Constructors and Destructor

// main constructor - pass up the chain
EightNodeIsoparamBrick::EightNodeIsoparamBrick(int eNum,int *eNode) : ElementBase3D(eNum,eNode)
{
	nodes[8]=eNode[0];
}

#pragma mark EightNodeIsoparamBrick::Methods

// Get shape functions, but does note get derivatives. All calls for derivative
// must use second method which assume orthogonal elements
void EightNodeIsoparamBrick::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe) const
{
    double temp1,temp2,temp3;
    int i;
    
    // shape function
    for(i=0;i<8;i++)
    {	temp1=(1.+xii[i]*xi->x);
        temp2=(1.+eti[i]*xi->y);
        temp3=(1.+zti[i]*xi->z);
        sfxn[i]=temp1*temp2*temp3/8.;
    }
}

// get just shape functions and optionally derivative wrt x and y
// General for shape funciton, but derivatives assumes a rectangular mesh
void EightNodeIsoparamBrick::ShapeFunction(Vector *xi,int getDeriv,double *sfxn,
											double *xDeriv,double *yDeriv,double *zDeriv) const
{
    double temp1,temp2,temp3;
    int i;
    
    // shape function
    for(i=0;i<8;i++)
    {	temp1=(1.+xii[i]*xi->x);
        temp2=(1.+eti[i]*xi->y);
        temp3=(1.+zti[i]*xi->z);
        sfxn[i]=0.125*temp1*temp2*temp3;
		if(getDeriv)
		{	xDeriv[i]=0.25*xii[i]*temp2*temp3/GetDeltaX();
			yDeriv[i]=0.25*eti[i]*temp1*temp3/GetDeltaY();
			zDeriv[i]=0.25*zti[i]*temp1*temp2/GetDeltaZ();
		}
	}
}

// see if point is in this element (assumes rectangular)
short EightNodeIsoparamBrick::PtInElement(Vector &pt)
{	if(pt.x<xmin || pt.x>xmax) return FALSE;
	if(pt.y<ymin || pt.y>ymax) return FALSE;
	if(pt.z<zmin || pt.z>zmax) return FALSE;
	return TRUE;
}

// see if this element is cube in cartesion coordinates and return TRUE or FALSE
// if cube, set dx, dy, dz to element dimensions
int EightNodeIsoparamBrick::Orthogonal(double *dx,double *dy,double *dz)
{
	int i;
	double xdel,ydel,zdel;
	
	*dx=0.;
	*dy=0.;
	*dz=0.;
	
	// front face
	for(i=0;i<3;i++)
    {	xdel=fabs(nd[nodes[i+1]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		zdel=fabs(nd[nodes[i+1]]->z-nd[nodes[i]]->z);
		if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
		*dz=fmax(zdel,*dz);
	}
    xdel=fabs(nd[nodes[0]]->x-nd[nodes[3]]->x);
	ydel=fabs(nd[nodes[0]]->y-nd[nodes[3]]->y);
	zdel=fabs(nd[nodes[0]]->z-nd[nodes[3]]->z);
	if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
	*dx=fmax(xdel,*dx);
	*dy=fmax(ydel,*dy);
	*dz=fmax(zdel,*dz);
	
	// back face
	for(i=4;i<7;i++)
    {	xdel=fabs(nd[nodes[i+1]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		zdel=fabs(nd[nodes[i+1]]->z-nd[nodes[i]]->z);
		if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
		*dz=fmax(zdel,*dz);
	}
    xdel=fabs(nd[nodes[4]]->x-nd[nodes[7]]->x);
	ydel=fabs(nd[nodes[4]]->y-nd[nodes[7]]->y);
	zdel=fabs(nd[nodes[4]]->z-nd[nodes[7]]->z);
	if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
	*dx=fmax(xdel,*dx);
	*dy=fmax(ydel,*dy);
	*dz=fmax(zdel,*dz);
	
	// side edges
	for(i=0;i<=3;i++)
    {	xdel=fabs(nd[nodes[i+4]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		zdel=fabs(nd[nodes[i+1]]->z-nd[nodes[i]]->z);
		if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
		*dz=fmax(zdel,*dz);
	}
	
	return TRUE;
}

// find dimensionless position, but assumes an orthogonal element
void EightNodeIsoparamBrick::GetXiPos(Vector *pos,Vector *xipos) const
{
	xipos->x=(2.*pos->x-xmin-xmax)/GetDeltaX();
	xipos->y=(2.*pos->y-ymin-ymax)/GetDeltaY();
	xipos->z=(2.*pos->z-zmin-zmax)/GetDeltaZ();
}

#pragma mark EightNodeIsoparamBrick::GIMP Methods

// Get GIMP nodes around an element #num, but only where shape functions are nonzero
// assumed to be properly numbered regular 3D array
// load nodes into nds[1]... and node ID (0-63) into ndIDs[0]...
void EightNodeIsoparamBrick::GetGimpNodes(int *numnds,int *nds,int *ndIDs,Vector *xipos) const
{
	// quadrant barriers assuming 8 particles
	double lp = mpmgrid.lp;
	double q1 = -1.+lp, q2 = 1.-lp;
	
	// nodes directly associated with the element
	nds[1]=nodes[0];
	nds[2]=nodes[1];
	nds[3]=nodes[2];
	nds[4]=nodes[3];
	nds[5]=nodes[4];
	nds[6]=nodes[5];
	nds[7]=nodes[6];
	nds[8]=nodes[7];
	ndIDs[0]=0;
	ndIDs[1]=1;
	ndIDs[2]=2;
	ndIDs[3]=3;
	ndIDs[4]=4;
	ndIDs[5]=5;
	ndIDs[6]=6;
	ndIDs[7]=7;

	// lower y quadrant
	if(xipos->y<q1)
	{	if(xipos->x<q1)
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;			// before element z direction
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.xplane;		// before element x direction
				nds[14]=nodes[3]-mpmgrid.xplane;
				nds[15]=nodes[4]-mpmgrid.xplane;
				nds[16]=nodes[7]-mpmgrid.xplane;
				nds[17]=nodes[0]-mpmgrid.yplane;		// before element y direction
				nds[18]=nodes[1]-mpmgrid.yplane;
				nds[19]=nodes[4]-mpmgrid.yplane;
				nds[20]=nodes[5]-mpmgrid.yplane;
				nds[21]=nds[13]-mpmgrid.zplane;
				nds[22]=nds[14]-mpmgrid.zplane;
				nds[23]=nds[17]-mpmgrid.zplane;
				nds[24]=nds[18]-mpmgrid.zplane;
				nds[25]=nds[17]-mpmgrid.xplane;
				nds[26]=nds[19]-mpmgrid.xplane;
				nds[27]=nds[25]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=19;
				ndIDs[13]=18;
				ndIDs[14]=31;
				ndIDs[15]=30;
				ndIDs[16]=9;
				ndIDs[17]=10;
				ndIDs[18]=21;
				ndIDs[19]=22;
				ndIDs[20]=47;
				ndIDs[21]=46;
				ndIDs[22]=37;
				ndIDs[23]=38;
				ndIDs[24]=8;
				ndIDs[25]=20;
				ndIDs[26]=36;
				*numnds=27;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[0]-mpmgrid.xplane;
				nds[10]=nodes[3]-mpmgrid.xplane;
				nds[11]=nodes[4]-mpmgrid.xplane;
				nds[12]=nodes[7]-mpmgrid.xplane;
				nds[13]=nodes[0]-mpmgrid.yplane;
				nds[14]=nodes[1]-mpmgrid.yplane;
				nds[15]=nodes[4]-mpmgrid.yplane;
				nds[16]=nodes[5]-mpmgrid.yplane;
				nds[17]=nds[13]-mpmgrid.xplane;
				nds[18]=nds[15]-mpmgrid.xplane;
				ndIDs[8]=19;
				ndIDs[9]=18;
				ndIDs[10]=31;
				ndIDs[11]=30;
				ndIDs[12]=9;
				ndIDs[13]=10;
				ndIDs[14]=21;
				ndIDs[15]=22;
				ndIDs[16]=8;
				ndIDs[17]=20;
				*numnds=18;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.xplane;
				nds[14]=nodes[3]-mpmgrid.xplane;
				nds[15]=nodes[4]-mpmgrid.xplane;
				nds[16]=nodes[7]-mpmgrid.xplane;
				nds[17]=nodes[0]-mpmgrid.yplane;
				nds[18]=nodes[1]-mpmgrid.yplane;
				nds[19]=nodes[4]-mpmgrid.yplane;
				nds[20]=nodes[5]-mpmgrid.yplane;
				nds[21]=nds[15]+mpmgrid.zplane;
				nds[22]=nds[16]+mpmgrid.zplane;
				nds[23]=nds[19]+mpmgrid.zplane;
				nds[24]=nds[20]+mpmgrid.zplane;
				nds[25]=nds[17]-mpmgrid.xplane;
				nds[26]=nds[19]-mpmgrid.xplane;
				nds[27]=nds[26]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=19;
				ndIDs[13]=18;
				ndIDs[14]=31;
				ndIDs[15]=30;
				ndIDs[16]=9;
				ndIDs[17]=10;
				ndIDs[18]=21;
				ndIDs[19]=22;
				ndIDs[20]=63;
				ndIDs[21]=62;
				ndIDs[22]=53;
				ndIDs[23]=54;
				ndIDs[24]=8;
				ndIDs[25]=20;
				ndIDs[26]=52;
				*numnds=27;
			}
		}
		else if(xipos->x<=q2)
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.yplane;
				nds[14]=nodes[1]-mpmgrid.yplane;
				nds[15]=nodes[4]-mpmgrid.yplane;
				nds[16]=nodes[5]-mpmgrid.yplane;
				nds[17]=nds[13]-mpmgrid.zplane;
				nds[18]=nds[14]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=9;
				ndIDs[13]=10;
				ndIDs[14]=21;
				ndIDs[15]=22;
				ndIDs[16]=37;
				ndIDs[17]=38;
				*numnds=18;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[0]-mpmgrid.yplane;
				nds[10]=nodes[1]-mpmgrid.yplane;
				nds[11]=nodes[4]-mpmgrid.yplane;
				nds[12]=nodes[5]-mpmgrid.yplane;
				ndIDs[8]=9;
				ndIDs[9]=10;
				ndIDs[10]=21;
				ndIDs[11]=22;
				*numnds=12;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.yplane;
				nds[14]=nodes[1]-mpmgrid.yplane;
				nds[15]=nodes[4]-mpmgrid.yplane;
				nds[16]=nodes[5]-mpmgrid.yplane;
				nds[17]=nds[15]+mpmgrid.zplane;
				nds[18]=nds[16]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=9;
				ndIDs[13]=10;
				ndIDs[14]=21;
				ndIDs[15]=22;
				ndIDs[16]=53;
				ndIDs[17]=54;
				*numnds=18;
			}
		}
		else
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[1]+mpmgrid.xplane;
				nds[14]=nodes[2]+mpmgrid.xplane;
				nds[15]=nodes[5]+mpmgrid.xplane;
				nds[16]=nodes[6]+mpmgrid.xplane;
				nds[17]=nodes[0]-mpmgrid.yplane;
				nds[18]=nodes[1]-mpmgrid.yplane;
				nds[19]=nodes[4]-mpmgrid.yplane;
				nds[20]=nodes[5]-mpmgrid.yplane;
				nds[21]=nds[13]-mpmgrid.zplane;
				nds[22]=nds[14]-mpmgrid.zplane;
				nds[23]=nds[17]-mpmgrid.zplane;
				nds[24]=nds[18]-mpmgrid.zplane;
				nds[25]=nds[18]+mpmgrid.xplane;
				nds[26]=nds[20]+mpmgrid.xplane;
				nds[27]=nds[25]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=12;
				ndIDs[13]=13;
				ndIDs[14]=24;
				ndIDs[15]=25;
				ndIDs[16]=9;
				ndIDs[17]=10;
				ndIDs[18]=21;
				ndIDs[19]=22;
				ndIDs[20]=40;
				ndIDs[21]=41;
				ndIDs[22]=37;
				ndIDs[23]=38;
				ndIDs[24]=11;
				ndIDs[25]=23;
				ndIDs[26]=39;
				*numnds=27;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[1]+mpmgrid.xplane;
				nds[10]=nodes[2]+mpmgrid.xplane;
				nds[11]=nodes[5]+mpmgrid.xplane;
				nds[12]=nodes[6]+mpmgrid.xplane;
				nds[13]=nodes[0]-mpmgrid.yplane;
				nds[14]=nodes[1]-mpmgrid.yplane;
				nds[15]=nodes[4]-mpmgrid.yplane;
				nds[16]=nodes[5]-mpmgrid.yplane;
				nds[17]=nds[14]+mpmgrid.xplane;
				nds[18]=nds[16]+mpmgrid.xplane;
				ndIDs[8]=12;
				ndIDs[9]=13;
				ndIDs[10]=24;
				ndIDs[11]=25;
				ndIDs[12]=9;
				ndIDs[13]=10;
				ndIDs[14]=21;
				ndIDs[15]=22;
				ndIDs[16]=11;
				ndIDs[17]=23;
				*numnds=18;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[1]+mpmgrid.xplane;
				nds[14]=nodes[2]+mpmgrid.xplane;
				nds[15]=nodes[5]+mpmgrid.xplane;
				nds[16]=nodes[6]+mpmgrid.xplane;
				nds[17]=nodes[0]-mpmgrid.yplane;
				nds[18]=nodes[1]-mpmgrid.yplane;
				nds[19]=nodes[4]-mpmgrid.yplane;
				nds[20]=nodes[5]-mpmgrid.yplane;
				nds[21]=nds[15]+mpmgrid.zplane;
				nds[22]=nds[16]+mpmgrid.zplane;
				nds[23]=nds[19]+mpmgrid.zplane;
				nds[24]=nds[20]+mpmgrid.zplane;
				nds[25]=nds[18]+mpmgrid.xplane;
				nds[26]=nds[20]+mpmgrid.xplane;
				nds[27]=nds[26]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=12;
				ndIDs[13]=13;
				ndIDs[14]=24;
				ndIDs[15]=25;
				ndIDs[16]=9;
				ndIDs[17]=10;
				ndIDs[18]=21;
				ndIDs[19]=22;
				ndIDs[20]=56;
				ndIDs[21]=57;
				ndIDs[22]=53;
				ndIDs[23]=54;
				ndIDs[24]=11;
				ndIDs[25]=23;
				ndIDs[26]=55;
				*numnds=27;
			}
		}
	}

	else if(xipos->y<=q2)
	{	if(xipos->x<q1)
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.xplane;
				nds[14]=nodes[3]-mpmgrid.xplane;
				nds[15]=nodes[4]-mpmgrid.xplane;
				nds[16]=nodes[7]-mpmgrid.xplane;
				nds[17]=nds[13]-mpmgrid.zplane;
				nds[18]=nds[14]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=19;
				ndIDs[13]=18;
				ndIDs[14]=31;
				ndIDs[15]=30;
				ndIDs[16]=47;
				ndIDs[17]=46;
				*numnds=18;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[0]-mpmgrid.xplane;
				nds[10]=nodes[3]-mpmgrid.xplane;
				nds[11]=nodes[4]-mpmgrid.xplane;
				nds[12]=nodes[7]-mpmgrid.xplane;
				ndIDs[8]=19;
				ndIDs[9]=18;
				ndIDs[10]=31;
				ndIDs[11]=30;
				*numnds=12;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.xplane;
				nds[14]=nodes[3]-mpmgrid.xplane;
				nds[15]=nodes[4]-mpmgrid.xplane;
				nds[16]=nodes[7]-mpmgrid.xplane;
				nds[17]=nds[15]+mpmgrid.zplane;
				nds[18]=nds[16]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=19;
				ndIDs[13]=18;
				ndIDs[14]=31;
				ndIDs[15]=30;
				ndIDs[16]=63;
				ndIDs[17]=62;
				*numnds=18;
			}
		}
		else if(xipos->x<=q2)
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				*numnds=12;
			}
			else if(xipos->z<=q2)
			{	*numnds=8;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				*numnds=12;
			}
		}
		else
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[1]+mpmgrid.xplane;
				nds[14]=nodes[2]+mpmgrid.xplane;
				nds[15]=nodes[5]+mpmgrid.xplane;
				nds[16]=nodes[6]+mpmgrid.xplane;
				nds[17]=nds[13]-mpmgrid.zplane;
				nds[18]=nds[14]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=12;
				ndIDs[13]=13;
				ndIDs[14]=24;
				ndIDs[15]=25;
				ndIDs[16]=40;
				ndIDs[17]=41;
				*numnds=18;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[1]+mpmgrid.xplane;
				nds[10]=nodes[2]+mpmgrid.xplane;
				nds[11]=nodes[5]+mpmgrid.xplane;
				nds[12]=nodes[6]+mpmgrid.xplane;
				ndIDs[8]=12;
				ndIDs[9]=13;
				ndIDs[10]=24;
				ndIDs[11]=25;
				*numnds=12;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[1]+mpmgrid.xplane;
				nds[14]=nodes[2]+mpmgrid.xplane;
				nds[15]=nodes[5]+mpmgrid.xplane;
				nds[16]=nodes[6]+mpmgrid.xplane;
				nds[17]=nds[15]+mpmgrid.zplane;
				nds[18]=nds[16]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=12;
				ndIDs[13]=13;
				ndIDs[14]=24;
				ndIDs[15]=25;
				ndIDs[16]=56;
				ndIDs[17]=57;
				*numnds=18;
			}
		}
	}

	else 
	{	if(xipos->x<q1)
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.xplane;
				nds[14]=nodes[3]-mpmgrid.xplane;
				nds[15]=nodes[4]-mpmgrid.xplane;
				nds[16]=nodes[7]-mpmgrid.xplane;
				nds[17]=nodes[3]+mpmgrid.yplane;
				nds[18]=nodes[2]+mpmgrid.yplane;
				nds[19]=nodes[7]+mpmgrid.yplane;
				nds[20]=nodes[6]+mpmgrid.yplane;
				nds[21]=nds[13]-mpmgrid.zplane;
				nds[22]=nds[14]-mpmgrid.zplane;
				nds[23]=nds[17]-mpmgrid.zplane;
				nds[24]=nds[18]-mpmgrid.zplane;
				nds[25]=nds[17]-mpmgrid.xplane;
				nds[26]=nds[19]-mpmgrid.xplane;
				nds[27]=nds[25]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=19;
				ndIDs[13]=18;
				ndIDs[14]=31;
				ndIDs[15]=30;
				ndIDs[16]=16;
				ndIDs[17]=15;
				ndIDs[18]=28;
				ndIDs[19]=27;
				ndIDs[20]=47;
				ndIDs[21]=46;
				ndIDs[22]=44;
				ndIDs[23]=43;
				ndIDs[24]=17;
				ndIDs[25]=29;
				ndIDs[26]=45;
				*numnds=27;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[0]-mpmgrid.xplane;
				nds[10]=nodes[3]-mpmgrid.xplane;
				nds[11]=nodes[4]-mpmgrid.xplane;
				nds[12]=nodes[7]-mpmgrid.xplane;
				nds[13]=nodes[3]+mpmgrid.yplane;
				nds[14]=nodes[2]+mpmgrid.yplane;
				nds[15]=nodes[7]+mpmgrid.yplane;
				nds[16]=nodes[6]+mpmgrid.yplane;
				nds[17]=nds[13]-mpmgrid.xplane;
				nds[18]=nds[15]-mpmgrid.xplane;
				ndIDs[8]=19;
				ndIDs[9]=18;
				ndIDs[10]=31;
				ndIDs[11]=30;
				ndIDs[12]=16;
				ndIDs[13]=15;
				ndIDs[14]=28;
				ndIDs[15]=27;
				ndIDs[16]=17;
				ndIDs[17]=29;
				*numnds=18;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[0]-mpmgrid.xplane;
				nds[14]=nodes[3]-mpmgrid.xplane;
				nds[15]=nodes[4]-mpmgrid.xplane;
				nds[16]=nodes[7]-mpmgrid.xplane;
				nds[17]=nodes[3]+mpmgrid.yplane;
				nds[18]=nodes[2]+mpmgrid.yplane;
				nds[19]=nodes[7]+mpmgrid.yplane;
				nds[20]=nodes[6]+mpmgrid.yplane;
				nds[21]=nds[15]+mpmgrid.zplane;
				nds[22]=nds[16]+mpmgrid.zplane;
				nds[23]=nds[19]+mpmgrid.zplane;
				nds[24]=nds[20]+mpmgrid.zplane;
				nds[25]=nds[17]-mpmgrid.xplane;
				nds[26]=nds[19]-mpmgrid.xplane;
				nds[27]=nds[26]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=19;
				ndIDs[13]=18;
				ndIDs[14]=31;
				ndIDs[15]=30;
				ndIDs[16]=16;
				ndIDs[17]=15;
				ndIDs[18]=28;
				ndIDs[19]=27;
				ndIDs[20]=63;
				ndIDs[21]=62;
				ndIDs[22]=60;
				ndIDs[23]=59;
				ndIDs[24]=17;
				ndIDs[25]=29;
				ndIDs[26]=61;
				*numnds=27;
			}
		}
		else if(xipos->x<=q2)
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[3]+mpmgrid.yplane;
				nds[14]=nodes[2]+mpmgrid.yplane;
				nds[15]=nodes[7]+mpmgrid.yplane;
				nds[16]=nodes[6]+mpmgrid.yplane;
				nds[17]=nds[13]-mpmgrid.zplane;
				nds[18]=nds[14]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=16;
				ndIDs[13]=15;
				ndIDs[14]=28;
				ndIDs[15]=27;
				ndIDs[16]=44;
				ndIDs[17]=43;
				*numnds=18;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[3]+mpmgrid.yplane;
				nds[10]=nodes[2]+mpmgrid.yplane;
				nds[11]=nodes[7]+mpmgrid.yplane;
				nds[12]=nodes[6]+mpmgrid.yplane;
				ndIDs[8]=16;
				ndIDs[9]=15;
				ndIDs[10]=28;
				ndIDs[11]=27;
				*numnds=12;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[3]+mpmgrid.yplane;
				nds[14]=nodes[2]+mpmgrid.yplane;
				nds[15]=nodes[7]+mpmgrid.yplane;
				nds[16]=nodes[6]+mpmgrid.yplane;
				nds[17]=nds[15]+mpmgrid.zplane;
				nds[18]=nds[16]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=16;
				ndIDs[13]=15;
				ndIDs[14]=28;
				ndIDs[15]=27;
				ndIDs[16]=60;
				ndIDs[17]=59;
				*numnds=18;
			}
		}
		else
		{	if(xipos->z<q1)
			{	nds[9]=nodes[0]-mpmgrid.zplane;
				nds[10]=nodes[1]-mpmgrid.zplane;
				nds[11]=nodes[2]-mpmgrid.zplane;
				nds[12]=nodes[3]-mpmgrid.zplane;
				nds[13]=nodes[1]+mpmgrid.xplane;
				nds[14]=nodes[2]+mpmgrid.xplane;
				nds[15]=nodes[5]+mpmgrid.xplane;
				nds[16]=nodes[6]+mpmgrid.xplane;
				nds[17]=nodes[3]+mpmgrid.yplane;
				nds[18]=nodes[2]+mpmgrid.yplane;
				nds[19]=nodes[7]+mpmgrid.yplane;
				nds[20]=nodes[6]+mpmgrid.yplane;
				nds[21]=nds[13]-mpmgrid.zplane;
				nds[22]=nds[14]-mpmgrid.zplane;
				nds[23]=nds[17]-mpmgrid.zplane;
				nds[24]=nds[18]-mpmgrid.zplane;
				nds[25]=nds[18]+mpmgrid.xplane;
				nds[26]=nds[20]+mpmgrid.xplane;
				nds[27]=nds[25]-mpmgrid.zplane;
				ndIDs[8]=32;
				ndIDs[9]=33;
				ndIDs[10]=34;
				ndIDs[11]=35;
				ndIDs[12]=12;
				ndIDs[13]=13;
				ndIDs[14]=24;
				ndIDs[15]=25;
				ndIDs[16]=16;
				ndIDs[17]=15;
				ndIDs[18]=28;
				ndIDs[19]=27;
				ndIDs[20]=40;
				ndIDs[21]=41;
				ndIDs[22]=44;
				ndIDs[23]=43;
				ndIDs[24]=14;
				ndIDs[25]=26;
				ndIDs[26]=42;
				*numnds=27;
			}
			else if(xipos->z<=q2)
			{	nds[9]=nodes[1]+mpmgrid.xplane;
				nds[10]=nodes[2]+mpmgrid.xplane;
				nds[11]=nodes[5]+mpmgrid.xplane;
				nds[12]=nodes[6]+mpmgrid.xplane;
				nds[13]=nodes[3]+mpmgrid.yplane;
				nds[14]=nodes[2]+mpmgrid.yplane;
				nds[15]=nodes[7]+mpmgrid.yplane;
				nds[16]=nodes[6]+mpmgrid.yplane;
				nds[17]=nds[14]+mpmgrid.xplane;
				nds[18]=nds[16]+mpmgrid.xplane;
				ndIDs[8]=12;
				ndIDs[9]=13;
				ndIDs[10]=24;
				ndIDs[11]=25;
				ndIDs[12]=16;
				ndIDs[13]=15;
				ndIDs[14]=28;
				ndIDs[15]=27;
				ndIDs[16]=14;
				ndIDs[17]=26;
				*numnds=18;
			}
			else
			{	nds[9]=nodes[4]+mpmgrid.zplane;
				nds[10]=nodes[5]+mpmgrid.zplane;
				nds[11]=nodes[6]+mpmgrid.zplane;
				nds[12]=nodes[7]+mpmgrid.zplane;
				nds[13]=nodes[1]+mpmgrid.xplane;
				nds[14]=nodes[2]+mpmgrid.xplane;
				nds[15]=nodes[5]+mpmgrid.xplane;
				nds[16]=nodes[6]+mpmgrid.xplane;
				nds[17]=nodes[3]+mpmgrid.yplane;
				nds[18]=nodes[2]+mpmgrid.yplane;
				nds[19]=nodes[7]+mpmgrid.yplane;
				nds[20]=nodes[6]+mpmgrid.yplane;
				nds[21]=nds[15]+mpmgrid.zplane;
				nds[22]=nds[16]+mpmgrid.zplane;
				nds[23]=nds[19]+mpmgrid.zplane;
				nds[24]=nds[20]+mpmgrid.zplane;
				nds[25]=nds[18]+mpmgrid.xplane;
				nds[26]=nds[20]+mpmgrid.xplane;
				nds[27]=nds[26]+mpmgrid.zplane;
				ndIDs[8]=48;
				ndIDs[9]=49;
				ndIDs[10]=50;
				ndIDs[11]=51;
				ndIDs[12]=12;
				ndIDs[13]=13;
				ndIDs[14]=24;
				ndIDs[15]=25;
				ndIDs[16]=16;
				ndIDs[17]=15;
				ndIDs[18]=28;
				ndIDs[19]=27;
				ndIDs[20]=56;
				ndIDs[21]=57;
				ndIDs[22]=60;
				ndIDs[23]=59;
				ndIDs[24]=14;
				ndIDs[25]=26;
				ndIDs[26]=58;
				*numnds=27;
			}
		}
	}
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular 3D array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-63)
void EightNodeIsoparamBrick::GimpShapeFunction(Vector *xi,int numnds,int *ndIDs,int getDeriv,double *sfxn,
						double *xDeriv,double *yDeriv,double *zDeriv) const
{
	int i;
	double xp,yp,zp,Svpx,Svpy,Svpz,dSvpx,dSvpy,dSvpz,xsign,ysign,zsign,argx=0.,argy=0.,argz=0.;
	
	// L is the cell spacing, 2*lp is the current particle size.
	// assuming the particle size is the same in x, y and z direction in element coordinate
	// the deformation of the particle is not considered yet.
	// assumes 8 particles per element
	double lp = mpmgrid.lp;
	double q1 = lp,q2 = 2.-lp, q3 = 2.+lp;

	for(i=0;i<numnds;i++)
	{	xp=fabs(xi->x-gxii[ndIDs[i]]);			// first quadrant (xp, yp)>=0
		yp=fabs(xi->y-geti[ndIDs[i]]);
		zp=fabs(xi->z-gzti[ndIDs[i]]);
		
		if(xp<=q1)
			Svpx = ((4.-lp)*lp-xp*xp)/(4.*lp);	// if lp=0.5: -(4.*xp*xp-7.)/8.;
		else if(xp<=q2)
			Svpx = (2.-xp)/2.;
		else if(xp<=q3)
		{	argx = (2.+lp-xp)/(4.*lp);			// if lp=0.5: (5.-2.*xp)/4
			Svpx = 2.*lp*argx*argx;				// if lp=0.5: (5.-2.*xp)^2/16
		}
		else
			Svpx=0.;
			
		if(yp<=q1)
			Svpy = ((4.-lp)*lp-yp*yp)/(4.*lp);	// if lp=0.5: -(4.*yp*yp-7.)/8.;
		else if(yp<=q2)
			Svpy = (2.-yp)/2.;
		else if(yp<=q3)
		{	argy = (2.+lp-yp)/(4.*lp);			// if lp=0.5: (5.-2.*yp)/4
			Svpy = 2.*lp*argy*argy;				// if lp=0.5: (5.-2.*yp)^2/16
		}
		else
			Svpy=0.;

		if(zp<=q1)
			Svpz = ((4.-lp)*lp-zp*zp)/(4.*lp);	// if lp=0.5: -(4.*zp*zp-7.)/8.;
		else if(zp<=q2)
			Svpz = (2.-zp)/2.;
		else if(zp<=q3)
		{	argz = (2.+lp-zp)/(4.*lp);			// if lp=0.5: (5.-2.*zp)/4
			Svpz = 2.*lp*argz*argz;				// if lp=0.5: (5.-2.*zp)^2/16
		}
		else
			Svpz=0.;
			
		sfxn[i] = Svpx*Svpy*Svpz;
				
		// find shape function at (xp,yp,zp) 		
		if(getDeriv)
		{	xsign = xi->x>gxii[ndIDs[i]] ? 1. : -1.;
			ysign = xi->y>geti[ndIDs[i]] ? 1. : -1.;
			zsign = xi->z>gzti[ndIDs[i]] ? 1. : -1.;

			if(xp<=q1)
				dSvpx = -xp/(2.*lp);			// if lp=0.5: -xp
			else if(xp<=q2)
				dSvpx = -0.5;
			else if(xp<=q3)
				dSvpx = -argx;
			else
				dSvpx = 0.;
				
			if(yp<=q1)
				dSvpy = -yp/(2.*lp);			// if lp=0.5: -yp
			else if(yp<=q2)
				dSvpy = -0.5;
			else if(yp<=q3)
				dSvpy = -argy;
			else
				dSvpy = 0.;

			if(zp<=q1)
				dSvpz = -zp/(2.*lp);			// if lp=0.5: -zp;
			else if(zp<=q2)
				dSvpz = -0.5;
			else if(zp<=q3)
				dSvpz = -argz;
			else
				dSvpz = 0.;

			xDeriv[i] = xsign*dSvpx*Svpy*Svpz*2.0/GetDeltaX();
			yDeriv[i] = ysign*Svpx*dSvpy*Svpz*2.0/GetDeltaY();
			zDeriv[i] = zsign*Svpx*Svpy*dSvpz*2.0/GetDeltaZ();
		}
	}


}

#pragma mark EightNodeIsoparamBrick::Accessors

// element name as an ID
short EightNodeIsoparamBrick::ElementName(void) { return(EIGHT_NODE_ISO_BRICK); }

// number of nodes in this element
int EightNodeIsoparamBrick::NumberNodes(void) const { return 8; }

// Get x-y area, z thickness, or volumne - all orthogonal brick
double EightNodeIsoparamBrick::GetArea(void) const { return (xmax-xmin)*(ymax-ymin); }
double EightNodeIsoparamBrick::GetVolume(void) const { return (xmax-xmin)*(ymax-ymin)*(zmax-zmin); }
double EightNodeIsoparamBrick::GetThickness(void) const { return zmax-zmin; }

