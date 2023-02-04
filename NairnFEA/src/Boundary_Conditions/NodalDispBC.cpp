/********************************************************************************
    NodalDispBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/NodalDispBC.hpp"
#include "Nodes/NodalPoint.hpp"

// Nodal BC global
NodalDispBC *firstDispBC=NULL;

#pragma mark NodalDispBC: Constructors and Destructors

NodalDispBC::NodalDispBC(int num,int dof) : FEABoundaryCondition()
{
	DefaultSetup(num,dof,(double)0.);
}

NodalDispBC::NodalDispBC(int num,int dof,double displacement) : FEABoundaryCondition()
{
	DefaultSetup(num,dof,displacement);
}

void NodalDispBC::DefaultSetup(int num,int dof,double displacement)
{
    if(num>=0)
    {	nodeNum=num;
        direction=dof;
        axis=0;
    }
    else
    {	nodeNum=-num;
        direction=0;
        axis=dof;
    }
    angle=0.;
	bcValue=displacement;
}

#pragma mark NodalDispBC: Methods

// print it
NodalDispBC *NodalDispBC::PrintBC(ostream &os)
{
    char nline[200];
    size_t nlsize=200;
    
    if(direction>0)
        snprintf(nline,nlsize,"%5d  %2d  %15.7e",nodeNum,direction,bcValue);
    else
        snprintf(nline,nlsize,"%5d                          %2d  %15.7e",
                        nodeNum,axis,angle);
    os << nline << endl;
	
	return (NodalDispBC *)nextObject;
}

// remap if resequenced
NodalDispBC *NodalDispBC::MapNodes(int *revMap)
{
	nodeNum=revMap[nodeNum];
	return (NodalDispBC *)nextObject;
}

// fix nodal displacement or rotate coordinates and return nextBC
NodalDispBC *NodalDispBC::FixOrRotate(double **st,double *rm,int nsize,int nband,int nfree)
{
    int jend,j,i,ii,jj,rowj;
    double skew,cs,c2,s2,sn,cssn;
    double siii,sjjj,siiii,siijj,sjjjj,sijj,siij,rii,rjj;
    
    // Fixed displacement at nodeNum in direction of magnitude bcValue
    if(direction>0)
    {   // Subtract bcValue*col(j) from right side of system
        //     and set fixed displacement
        i=nfree*(nodeNum-1)+direction;		// fixed DOF
        if(bcValue!=0.)
        {   jend=fmin(nsize,i+nband-1);
            for(j=i+1;j<=jend;j++) rm[j]-=bcValue*st[i][j-i+1];
            jend=fmax(1,i-nband+1);
            for(j=i-1;j>=jend;j--) rm[j]-=bcValue*st[j][i-j+1];
        }
        rm[i]=bcValue;

        // Zero row and column i and put 1 on diagonal
        for(j=2;j<=nband;j++)
        {   jend=i-j+1;
            if(jend>=1) st[jend][j]=0.;		// zeros col i
            st[i][j]=0.;					// zeros row i
        }
        st[i][1]=1.;
    }

    /* If direction==0 then rotate nodeNum buy alpha degrees.
        This section just does the rotation (which is only allowed once per DOF)
        A later fixed displacement will set displacement at this DOF
    */
    else
    {	ii=nfree*(nodeNum-1)+1;
        jj=ii+1;
        skew=PI_CONSTANT*angle/180.;
        cs=cos(skew);
        c2=cs*cs;
        sn=sin(skew);
        s2=sn*sn;
        cssn=cs*sn;

        // Form TKTt in columns ii and jj (except near diagonals)
        jend=fmin(nband,jj);
        for(j=3;j<=jend;j++)
        {   rowj=jj-j+1;
            siii=st[rowj][j-1];
            sijj=st[rowj][j];
            st[rowj][j-1]=cs*siii-sn*sijj;
            st[rowj][j]=sn*siii+cs*sijj;
        }

        // Form TKTt in rows ii and jj (except near diagonals)
        jend=fmin(nband,nsize-ii+1);
        for(j=3;j<=jend;j++)
        {   siij=st[ii][j];
            sjjj=st[jj][j-1];
            st[ii][j]=cs*siij-sn*sjjj;
            st[jj][j-1]=sn*siij+cs*sjjj;
        }

        // Do diagonal terms
        siiii=st[ii][1];
        siijj=st[ii][2];
        sjjjj=st[jj][1];
        st[ii][1]=siiii*c2+sjjjj*s2-2*siijj*cssn;
        st[jj][1]=siiii*s2+sjjjj*c2+2*siijj*cssn;
        st[ii][2]=(siiii-sjjjj)*cssn+siijj*(c2-s2);

        // Transform right side vector
        rii=rm[ii];
        rjj=rm[jj];
        rm[ii]=rii*cs-rjj*sn;
        rm[jj]=rii*sn+rjj*cs;
    }
	
	return (NodalDispBC *)nextObject;
}

// unrotate any skewed nodes and return nextBC
NodalDispBC *NodalDispBC::Unrotate(double *rm,int nfree)
{
	int ii,jj;
	double alpha,cs,sn,rii,rjj;
	
	// unrotate if a skew BC
	if(direction<=0)
	{   // unrotate at nodeNum reaction
		ii=nfree*(nodeNum-1)+1;
		jj=ii+1;
		alpha=-PI_CONSTANT*angle/180.;
		cs=cos(alpha);
		sn=sin(alpha);
		rii=rm[ii];
		rjj=rm[jj];
		rm[ii]=rii*cs-rjj*sn;
		rm[jj]=rii*sn+rjj*cs;
	}
	
	return (NodalDispBC *)nextObject;
}

// Print reactivities at fixed node
void NodalDispBC::PrintReaction(void)
{
    char fline[200];
    size_t fsize=200;
    
    // skip if rotation
    if(direction<=0) return;
    
    snprintf(fline,fsize,"%5d     %15.7e     %15.7e",nodeNum,
            nd[nodeNum]->fs->force.x,nd[nodeNum]->fs->force.y);
    cout << fline << endl;
}

#pragma mark NodalDispBC: Accessors

// equal means applied to same node and dof, value does not matter
bool NodalDispBC::SameDofSetting(NodalDispBC *rhs)
{	return (nodeNum==rhs->nodeNum) && (direction==rhs->direction) && (axis==rhs->axis);
}
