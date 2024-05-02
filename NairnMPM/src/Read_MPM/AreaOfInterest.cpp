/********************************************************************************
	AreaOfInterest.cpp
	NairnMPM

	Created by John Nairn on 2/19/2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif

#include "AreaOfInterest.hpp"
#include "System/MPMPrefix.hpp"

#pragma mark AreaOfInterest::Constructors and Destructors

// Constructors
AreaOfInterest::AreaOfInterest()
{
}

// Constructors
AreaOfInterest::AreaOfInterest(double ax1,double ax2,int anx,double ay1,double ay2,int any,double az1,double az2,int anz)
{
	x1 = ax1;
	x2 = ax2;
	nx = anx;
	xcell = (ax2-ax1)/(double)anx;
	
	y1 = ay1;
	y2 = ay2;
	ny = any;
	ycell = (ay2-ay1)/(double)any;
	
	z1 = az1;
	z2 = az2;
	nz = anz;
	zcell = (az2-az1)/(double)anz;

}

#pragma mark AreaOfInterest::Methods

// return true is axis iax (0,1,2) starts before that axis in aoi
bool AreaOfInterest::IsBefore(AreaOfInterest *aoi,int iax)
{
	switch(iax)
	{	case 0:
			return x1 < aoi->x1 ? true : false;
		case 1:
			return y1 < aoi->y1 ? true : false;
		case 2:
			return z1 < aoi->z1 ? true : false;
	}
	return false;
}

// Get cell counts on one axis
AreaOfInterest *AreaOfInterest::GetCellCounts(int ax,double &prevEdge,double lastEdge,double rmax,int *Ncells,int style)
{
	switch(ax)
	{	case 0:
			//cout << "x: " << x1 << " to " << x2 << " prev:" << prevEdge << " last:" << lastEdge;
			prevEdge = GetAxisCounts(0,x1,x2,xcell,nx,nxBefore,nxAfter,prevEdge,lastEdge,rmax,Ncells,style);
			//cout << " counts: " << nxBefore << "," << nx << "," << nxAfter << "," << *Ncells << endl;
			break;
		case 1:
			//cout << "y: " << y1 << " to " << y2 << " prev:" << prevEdge << " last:" << lastEdge;
			prevEdge = GetAxisCounts(1,y1,y2,ycell,ny,nyBefore,nyAfter,prevEdge,lastEdge,rmax,Ncells,style);
			//cout << " counts: " << nyBefore << "," << ny << "," << nyAfter << "," << *Ncells << endl;
			break;
		default:
			//cout << "z: " << z1 << " to " << z2 << " prev:" << prevEdge << " last:" << lastEdge;
			prevEdge = GetAxisCounts(2,z1,z2,zcell,nz,nzBefore,nzAfter,prevEdge,lastEdge,rmax,Ncells,style);
			//cout << " counts: " << nzBefore << "," << nz << "," << nzAfter << "," << *Ncells << endl;
			break;
	}
	return nextAOI[ax];
}

// Get counts on one axis
double AreaOfInterest::GetAxisCounts(int ax,double &a1,double &a2,double &acell,int &an,int &nbefore,int &nafter,
									 double prevEdge,double lastEdge,double rmax,int *Ncells,int style)
{
	// calculate only if not collapsed
	double delta;
	if(an>0)
	{	// cells before this area of interest
		delta = a1-prevEdge;
		
		// if less than cell, collapse is
		if(delta<acell)
		{	a1 = prevEdge;
			acell = (a2-a1)/(double)an;
			nbefore = 0;
		}
		else if(style==GEOMETRIC_STYLE)
		{	nbefore = int(log((rmax-1.)*(delta/acell)+rmax)/log(rmax) - 1.)+1;
		}
		else
		{	// linearly increasing
			double ada = 1./(rmax-1.);
			nbefore = int(sqrt(ada*(ada+1.+2.*delta/acell)+0.25)-ada-0.5)+1;
		}
	}
	
	// next area of interest
	if(nextAOI[ax] != NULL)
	{	double aoiEnd,aoiCell;
		double aoiStart = nextAOI[ax]->GetNextEdges(ax,aoiEnd,aoiCell);
		delta = 0.5*(aoiStart-a2);
		
		// if less than cell, combine them
		if(delta<acell)
		{	a2 = aoiEnd;
			double newCell = fmin(acell,aoiCell);
			an=(int)((a2-a1)/newCell+0.5);
			acell = (a2-a1)/(double)an;
			nextAOI[ax]->Collapse(ax,a1,a2);
			nafter = 0;
			*Ncells += nbefore + an;
			return a2;
		}
	}
	else
	{	delta = lastEdge - a2;
		if(delta<acell)
		{	a2 = lastEdge;
			acell = (a2-a1)/(double)an;
			nafter = 0;
			*Ncells += nbefore + an;
			return a2;
		}
	}
	
	// fill region after
	if(style==GEOMETRIC_STYLE)
		nafter = int(log((rmax-1.)*(delta/acell)+rmax)/log(rmax)-1.)+1;
	else
	{	// linearly increasing
		double ada = 1./(rmax-1.);
		nafter = int(sqrt(ada*(ada+1.+2.*delta/acell)+0.25)-ada-0.5)+1;
	}
	
	// increment count
	*Ncells += nbefore + an + nafter;
	
	// return new previous edge
	return a2 + delta;
}

// get edge and cell info along one ais
double AreaOfInterest::GetNextEdges(int ax,double &aoiEnd,double &aoiCell)
{
	switch(ax)
	{	case 0:
			aoiEnd = x2;
			aoiCell = xcell;
			return x1;
		case 1:
			aoiEnd = y2;
			aoiCell = ycell;
			return y1;
		default:
			aoiEnd = z2;
			aoiCell = zcell;
			return z1;
	}
}

// collapse aoi to combine with previous and none before or in the aoi
void AreaOfInterest::Collapse(int ax,double a1,double a2)
{
	switch(ax)
	{	case 0:
			x1=a1;
			x2=a2;
			nxBefore=0;
			nx=0;
		case 1:
			y1=a1;
			y2=a2;
			nyBefore=0;
			ny=0;
		default:
			z1=a1;
			z2=a2;
			nzBefore=0;
			nz=0;
	}
}

// Get nodes of before, in, and after this area of interest
AreaOfInterest *AreaOfInterest::MeshAreaOfInterest(int ax,double &prevEdge,double lastEdge,double rmax,int *num,double *pts,int style)
{
	switch(ax)
	{	case 0:
			prevEdge = MeshAxis(0,x1,x2,xcell,nx,nxBefore,nxAfter,prevEdge,lastEdge,rmax,num,pts,style);
			break;
		case 1:
			prevEdge = MeshAxis(1,y1,y2,ycell,ny,nyBefore,nyAfter,prevEdge,lastEdge,rmax,num,pts,style);
			break;
		default:
			prevEdge = MeshAxis(2,z1,z2,zcell,nz,nzBefore,nzAfter,prevEdge,lastEdge,rmax,num,pts,style);
			break;
	}
	return nextAOI[ax];
}

// Get nodes on one axis
double AreaOfInterest::MeshAxis(int ax,double a1,double a2,double acell,int an,int nbefore,int nafter,
									 double prevEdge,double lastEdge,double rmax,int *num,double *pts,int style)
{
	// calculate only if not collapsed
	double r,delta=0.,pt;
	
	// before the area
	if(nbefore>0)
	{	delta = a1-prevEdge;
		
		if(style==GEOMETRIC_STYLE)
		{	// find R
			if(nbefore==1)
			{	r = delta/acell;
			}
			else
			{	double r1=1;
				double r2=rmax;
				for(int i=0;i<10;i++)
				{	double rmid = 0.5*(r1+r2);
					double rtest = pow(rmid,(double)(nbefore+1)) - rmid - (rmid-1)*delta/acell;
					if(rtest<0.)
						r1 = rmid;
					else
						r2 = rmid;
				}
				r = 0.5*(r1+r2);
			}
			
			// start at prevEdge with cell size acell*pow(r,nbefore)
			double cellSize = acell*pow(r,(double)nbefore);
			pt = prevEdge;
			for(int i=0;i<nbefore;i++)
			{	pts[*num] = pt;
				*num = *num+1;
				pt += cellSize;
				cellSize /= r;
			}
		}
		else
		{	// find delta a
			double da = 2.*(delta-nbefore*acell)/((double)(nbefore*(nbefore+1)));
			
			// start at prevEdge with cell size acell+nbefore*da
			double cellSize = acell + nbefore*da;
			pt = prevEdge;
			for(int i=0;i<nbefore;i++)
			{	pts[*num] = pt;
				*num = *num+1;
				pt += cellSize;
				cellSize -= da;
			}
		}
	}
	
	// cells in area of interest
	if(an>0)
	{	pt = a1;
		for(int i=0;i<an;i++)
		{	pts[*num] = pt;
			*num = *num+1;
			pt += acell;
		}
	}
	
	// next area of interest and cells after this area of interest
	if(nafter>0)
	{	if(nextAOI[ax]!=NULL)
		{	double aoiEnd,aoiCell;
			double aoiStart = nextAOI[ax]->GetNextEdges(ax,aoiEnd,aoiCell);
			delta = 0.5*(aoiStart-a2);
		}
		else
			delta = lastEdge-a2;
		
		if(style==GEOMETRIC_STYLE)
		{	// find R
			if(nafter==1)
			{	r = delta/acell;
			}
			else
			{	double r1=1;
				double r2=rmax;
				for(int i=0;i<10;i++)
				{	double rmid = 0.5*(r1+r2);
					double rtest = pow(rmid,(double)(nafter+1)) - rmid - (rmid-1)*delta/acell;
					if(rtest<0.)
						r1 = rmid;
					else
						r2 = rmid;
				}
				r = 0.5*(r1+r2);
			}
			
			// start at prevEdge with cell size acell*pow(r,nbefore)
			double cellSize = acell*r;
			pt = a2;
			for(int i=0;i<nafter;i++)
			{	pts[*num] = pt;
				*num = *num+1;
				pt += cellSize;
				cellSize *= r;
			}
		}
		else
		{	// find delta a
			double da = 2.*(delta-nafter*acell)/((double)(nafter*(nafter+1)));
			
			// start at prevEdge with cell size acell+nbefore*da
			double cellSize = acell + da;
			pt = a2;
			for(int i=0;i<nafter;i++)
			{	pts[*num] = pt;
				*num = *num+1;
				pt += cellSize;
				cellSize += da;
			}
		}
	}
	
	// final point
	if(nextAOI[ax]==NULL)
	{	pts[*num] = lastEdge;
		*num = *num+1;
	}
	
	// return new previous edge
	return a2 + delta;
}

// linked objects (ax 0, 1, or 2, but not checked)
AreaOfInterest *AreaOfInterest::GetNextAOI(int ax) { return nextAOI[ax]; }
void AreaOfInterest::SetNextAOI(int ax,AreaOfInterest *aoi) { nextAOI[ax] = aoi; }


