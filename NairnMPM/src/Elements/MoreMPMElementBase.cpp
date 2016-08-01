/********************************************************************************
    MoreMPMElementBase.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Oct 22 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"

#define MAXITER 100

#pragma mark ElementBase: Constructors and Destructor MPM Only

/* Main MPM contructor when creating elements:
	input is node numbers (0-based array) but values
            are 1-based (length always MaxElNd, even if unused)
        subclass should set nodes[NumberNodes]=nodes[0]
*/
ElementBase::ElementBase(int eNum,int *eNode)
{	int i;
    num=eNum;
    for(i=0;i<MaxElNd;i++) nodes[i]=eNode[i];
	neighbors=NULL;
    filled=0;
}

ElementBase::~ElementBase()
{	if(neighbors!=NULL) delete [] neighbors;
}

#pragma mark ElementBase: Shape and Gradient calls for MPM


/*
 Find GIMP information or find CPDI info (one or the other, never both)
 Not called for rigid BC particles because do not use shape functions unless extrapolated
	and then use simple ones (non-GIMP)
 throws CommonException() if CPDI particle corner has left the grid
*/
void ElementBase::GetShapeFunctionData(MPMBase *mpmptr) const
{
    switch(useGimp)
    {
#ifdef LOAD_GIMP_INFO
		case POINT_GIMP:
			GetXiPos(&mpmptr->pos,mpmptr->GetNcpos());
			break;
			
        case UNIFORM_GIMP:
        case UNIFORM_GIMP_AS:
		{	// during initialization, find dimenless positiong, and nodes for non-zero shape function
			Vector *xipos = mpmptr->GetNcpos();
            GetXiPos(&mpmptr->pos,xipos);
			Vector lp;
            mpmptr->GetDimensionlessSize(lp);
			GIMPNodes *gimp = mpmptr->GetGIMPInfo();
            GetGimpNodes(gimp->nds,gimp->nds,gimp->ndIDs,xipos,lp);
			double fn[maxShapeNodes];
            GimpShapeFunction(xipos,gimp->nds[0],gimp->ndIDs,FALSE,&fn[1],NULL,NULL,NULL,lp);
			
			// compact the nodes and ndIDs
			int local=NumberNodes();
			int found=local;
			for(int i=local+1;i<=gimp->nds[0];i++)
			{	if(fn[i]<=0.0) continue;
				gimp->ndIDs[found] = gimp->ndIDs[i-1];
				found++;
				gimp->nds[found]=gimp->nds[i];
				fn[found]=fn[i];
			}
			gimp->nds[0]=found;
            break;
		}
			
#else
		case POINT_GIMP:
        case UNIFORM_GIMP:
        case UNIFORM_GIMP_AS:
			GetXiPos(&mpmptr->pos,mpmptr->GetNcpos());
			break;
			
#endif
        case LINEAR_CPDI:
		case LINEAR_CPDI_AS:
        case QUADRATIC_CPDI:
			mpmptr->GetCPDINodesAndWeights(useGimp);
            break;
    }
}

/* Find shape functions for crack particle using point GIMP method
	Finds nodes, dimensionless position, and then calls simple function
  NOTE: This is only called by crack update and do not need to save xipos.
*/
#define CRACK_POINT
void ElementBase::GetShapeFunctionsForCracks(int *numnds,double *fn,int *nds,Vector *pos) const
{
#ifdef CRACK_POINT
	Vector xipos;
	GetNodes(numnds,nds);
	GetXiPos(pos,&xipos);
	ShapeFunction(&xipos,FALSE,&fn[1],NULL,NULL,NULL);
#else
	This section no longer current. If switch to using it again, revise with method in following two routines
	/*
	Vector xipos,lp;
	unsigned char ndIDs[maxShapeNodes];
	
    switch(useGimp)
    {   case POINT_GIMP:
        	// Load element noodes, dimensionless position, and shape functinos
            GetNodes(numnds,nds);
            GetXiPos(pos,&xipos);
            ShapeFunction(&xipos,FALSE,&fn[1],NULL,NULL,NULL);
            break;
            
        case LINEAR_CPDI:
        case QUADRATIC_CPDI:
			// since no material point, CPDI uses GIMP method
        case UNIFORM_GIMP:
            GetXiPos(pos,&xipos);
            lp.x = mpmgrid.GetParticleSemiLength();		// no longer there
            lp.y = lp.x;
            lp.z = lp.x;
            GetGimpNodes(numnds,nds,ndIDs,&xipos,lp);
            GimpShapeFunction(&xipos,*numnds,ndIDs,FALSE,&fn[1],NULL,NULL,NULL,lp);
            GimpCompact(numnds,nds,fn,NULL,NULL,NULL);
            break;
            
		case LINEAR_CPDI_AS:
			// since no material point, CPDI uses GIMP method
        case UNIFORM_GIMP_AS:
            GetXiPos(pos,&xipos);
            lp.x = mpmgrid.GetParticleSemiLength();		// no longer there
            lp.y = lp.x;
            lp.z = lp.x;
            GetGimpNodes(numnds,nds,ndIDs,&xipos,lp);
            GimpShapeFunctionAS(&xipos,*numnds,ndIDs,FALSE,&fn[1],NULL,NULL,NULL,lp);
            GimpCompact(numnds,nds,fn,NULL,NULL,NULL);
            break;
     }
	 */
#endif
}

/*
 Return list of nodes (possible change array with address in ndsHandle)
	and the shape functions
 Load node numbers into nds[1]...
 Load shape functions into fn[1]...
 WARNING: This should neverbe called for Rigid BC particles
 NOTE: This is called at various places in the time step when shape functions are needed. It should
	recalculate the ones found at the begnning of the time step using precalculated xipos
    or CPDI info, which are found in initialization
 throws CommonException() if too many CPDI nodes
*/
void ElementBase::GetShapeFunctions(double *fn,int **ndsHandle,MPMBase *mpmptr) const
{
    Vector lp;
	int *nds = *ndsHandle;
    
    switch(useGimp)
    {   case POINT_GIMP:
			// load coordinates if not already done
            GetNodes(&nds[0],nds);
            ShapeFunction(mpmptr->GetNcpos(),false,&fn[1],NULL,NULL,NULL);
            break;

#ifdef LOAD_GIMP_INFO
        case UNIFORM_GIMP:
        {	// GIMP analysis
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
			GIMPNodes *gimp = mpmptr->GetGIMPInfo();
            GimpShapeFunction(xipos,gimp->nds[0],gimp->ndIDs,false,&fn[1],NULL,NULL,NULL,lp);
			*ndsHandle = gimp->nds;
            break;
        }
            
        case UNIFORM_GIMP_AS:
        {	// GIMP analysis
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
			GIMPNodes *gimp = mpmptr->GetGIMPInfo();
            GimpShapeFunctionAS(xipos,gimp->nds[0],gimp->ndIDs,false,&fn[1],NULL,NULL,NULL,lp);
			*ndsHandle = gimp->nds;
            break;
        }
            
#else
        case UNIFORM_GIMP:
        {	// GIMP analysis
            unsigned char ndIDs[maxShapeNodes];
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
            GetGimpNodes(&nds[0],nds,ndIDs,xipos,lp);
            GimpShapeFunction(xipos,nds[0],ndIDs,FALSE,&fn[1],NULL,NULL,NULL,lp);
            GimpCompact(&nds[0],nds,fn,NULL,NULL,NULL);
            break;
        }
            
        case UNIFORM_GIMP_AS:
        {	// GIMP analysis
            unsigned char ndIDs[maxShapeNodes];
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
			GetGimpNodes(&nds[0],nds,ndIDs,xipos,lp);
            GimpShapeFunctionAS(xipos,nds[0],ndIDs,FALSE,&fn[1],NULL,NULL,NULL,lp);
            GimpCompact(&nds[0],nds,fn,NULL,NULL,NULL);
			break;
        }
            
#endif
        case LINEAR_CPDI:
		case LINEAR_CPDI_AS:
        case QUADRATIC_CPDI:
            nds[0] = GetCPDIFunctions(nds,fn,NULL,NULL,NULL,mpmptr);
            break;
    }
}

/* 
 Return list of nodes (possibly change array with address in ndsHandle)
	and the shape functions and their gradients
 Load node numbers into nds[1]...
 Load shape functions into fn[1]...
 Load shape function derviatives into xDeriv[1]..., yDeriv[1]..., zDeriv[1]...
 WARNING: This should never be called for Rigid BC particles
 NOTE: This is called at various places in the time step when shape functions are needed. It should
	recalculate the ones found at the begnning of the time step using precalculated xipos
	or CPDI info, which are found in initialization
 throws CommonException() if too many CPDI nodes
*/
void ElementBase::GetShapeGradients(double *fn,int **ndsHandle,
                                    double *xDeriv,double *yDeriv,double *zDeriv,MPMBase *mpmptr) const
{
    Vector lp;
	int *nds = *ndsHandle;
    
    switch(useGimp)
    {   case POINT_GIMP:
			// load nodal numbers
            GetNodes(&nds[0],nds);
            
            // special case for regular mesh
            if(mpmgrid.GetCartesian()>0)
                ShapeFunction(mpmptr->GetNcpos(),true,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1]);
            else
            {	// Load element coordinates
                Vector ce[MaxElNd];
                double fnh[MaxElNd];
                GetCoordinates(ce,nds[0],nds);
                
                // find shape functions and derviatives
                ShapeFunction(mpmptr->GetNcpos(),BMATRIX,&fn[1],&xDeriv[1],&yDeriv[1],&ce[1],NULL,NULL,&fnh[1]);
            }
            break;
            
#ifdef LOAD_GIMP_INFO
        case UNIFORM_GIMP:
        {	// GIMP analysis
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
			GIMPNodes *gimp = mpmptr->GetGIMPInfo();
			GimpShapeFunction(xipos,gimp->nds[0],gimp->ndIDs,true,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1],lp);
			*ndsHandle = gimp->nds;
            break;
        }
            
        case UNIFORM_GIMP_AS:
        {	// GIMP analysis
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
			GIMPNodes *gimp = mpmptr->GetGIMPInfo();
            GimpShapeFunctionAS(xipos,gimp->nds[0],gimp->ndIDs,true,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1],lp);
			*ndsHandle = gimp->nds;
            break;
        }
            
#else
        case UNIFORM_GIMP:
        {	// uGIMP analysis
            unsigned char ndIDs[maxShapeNodes];
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
            GetGimpNodes(&nds[0],nds,ndIDs,xipos,lp);
            GimpShapeFunction(xipos,nds[0],ndIDs,TRUE,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1],lp);
            GimpCompact(&nds[0],nds,fn,xDeriv,yDeriv,zDeriv);
            break;
        }
            
        case UNIFORM_GIMP_AS:
        {	// uGIMP analysis
            unsigned char ndIDs[maxShapeNodes];
            Vector *xipos = mpmptr->GetNcpos();
            mpmptr->GetDimensionlessSize(lp);
            GetGimpNodes(&nds[0],nds,ndIDs,xipos,lp);
            GimpShapeFunctionAS(xipos,nds[0],ndIDs,TRUE,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1],lp);
            GimpCompact(&nds[0],nds,fn,xDeriv,yDeriv,zDeriv);
            break;
        }
            
#endif
        case LINEAR_CPDI:
		case LINEAR_CPDI_AS:
		case QUADRATIC_CPDI:
            nds[0] = GetCPDIFunctions(nds,fn,xDeriv,yDeriv,zDeriv,mpmptr);
            break;
    }
}

#pragma mark ElementBase: Shape and Gradient utility methods for MPM

/* Find dimensionless coordinates by numerical methods
	input: pos is position in the element
	output: xipos is dimensionless position
	only used in MPM and only here if non-rectangular elements
*/
void ElementBase::GetXiPos(Vector *pos,Vector *xipos) const
{
    double xt,yt,dxxi,dxeta,dyxi,dyeta;
    double deter,dxi,deta,dist;
    double gfn[MaxElNd],gdfnxi[MaxElNd],gdfnet[MaxElNd];
    int numnds=NumberNodes(),i,j;
	
	// initial guess
	GetCentroid(xipos);
	
	// nodal coordinates
	Vector eNode[MaxElNd];
	for(i=0;i<numnds;i++)
	{	eNode[i].x=nd[nodes[i]]->x;
		eNode[i].y=nd[nodes[i]]->y;
	}
    
    /* solve for xipos using Newton-Rapheson (see FEA Notes)
            using shape functions and their derivatives */
    for(j=1;j<=MAXITER;j++)
    {   ShapeFunction(xipos,TRUE,gfn,gdfnxi,gdfnet,NULL,NULL,NULL,NULL);
        xt=-pos->x;
        yt=-pos->y;
        dxxi=0.;
        dxeta=0.;
        dyxi=0.;
        dyeta=0.;
        for(i=0;i<numnds;i++)
        {   xt+=eNode[i].x*gfn[i];
            yt+=eNode[i].y*gfn[i];
            dxxi+=eNode[i].x*gdfnxi[i];
            dxeta+=eNode[i].x*gdfnet[i];
            dyxi+=eNode[i].y*gdfnxi[i];
            dyeta+=eNode[i].y*gdfnet[i];
        }
        deter=dxxi*dyeta-dxeta*dyxi;
        dxi=(-xt*dyeta+yt*dxeta)/deter;
        deta=(xt*dyxi-yt*dxxi)/deter;
        xipos->x+=dxi;
        xipos->y+=deta;
        dist=sqrt(dxi*dxi+deta*deta);
        if(dist<.001) break;
    }
}

// return dimensionless location for material points
void ElementBase::MPMPoints(short numPerElement,Vector *mpos) const
{
    int i,j,k;
    double fxn[MaxElNd],row,zrow;
    
    if(NumberSides()==4)
    {	switch(numPerElement)
        {   case 4:
                // ENI or FNI - 2D only
                mpos[0].x=-.5;
                mpos[0].y=-.5;
                mpos[0].z=0.;
                mpos[1].x=.5;
                mpos[1].y=-.5;
                mpos[1].z=0.;
                mpos[2].x=-.5;
                mpos[2].y=.5;
                mpos[2].z=0.;
                mpos[3].x=.5;
                mpos[3].y=.5;
                mpos[3].z=0.;
                break;
            case 1:
                // CM of square or brick
                mpos[0].x=0.;
                mpos[0].y=0.;
				mpos[0].z=0.;
				break;
			case 8:
				// 3D box
                mpos[0].x=-.5;
                mpos[0].y=-.5;
				mpos[0].z=-.5;
                mpos[1].x=.5;
                mpos[1].y=-.5;
				mpos[1].z=-.5;
                mpos[2].x=-.5;
                mpos[2].y=.5;
				mpos[2].z=-.5;
                mpos[3].x=.5;
                mpos[3].y=.5;
				mpos[3].z=-.5;
                mpos[4].x=-.5;
                mpos[4].y=-.5;
				mpos[4].z=.5;
                mpos[5].x=.5;
                mpos[5].y=-.5;
				mpos[5].z=.5;
                mpos[6].x=-.5;
                mpos[6].y=.5;
				mpos[6].z=.5;
                mpos[7].x=.5;
                mpos[7].y=.5;
				mpos[7].z=.5;
				break;
            case 9:
                // 2D
                k=0;
                row = -2./3.;
                for(j=0;j<3;j++)
                {   mpos[k].x=-2./3.;
                    mpos[k].y=row;
                    mpos[k].z=0.;
                    mpos[k+1].x=0.;
                    mpos[k+1].y=row;
                    mpos[k+1].z=0.;
                    mpos[k+2].x=2./3.;
                    mpos[k+2].y=row;
                    mpos[k+2].z=0.;
                    k += 3;
                    row += 2./3.;
                }
                break;
            case 16:
                // 2D
                k=0;
                row = -0.75;
                for(j=0;j<4;j++)
                {   mpos[k].x=-0.75;
                    mpos[k].y=row;
                    mpos[k].z=0.;
                    mpos[k+1].x=-.25;
                    mpos[k+1].y=row;
                    mpos[k+1].z=0.;
                    mpos[k+2].x=.25;
                    mpos[k+2].y=row;
                    mpos[k+2].z=0.;
                    mpos[k+3].x=.75;
                    mpos[k+3].y=row;
                    mpos[k+3].z=0.;
                    k += 4;
                    row += 0.5;
                }
                break;
            case 25:
                // 2D
                k=0;
                row = -0.8;
                for(j=0;j<5;j++)
                {   mpos[k].x=-0.8;
                    mpos[k].y=row;
                    mpos[k].z=0.;
                    mpos[k+1].x=-.4;
                    mpos[k+1].y=row;
                    mpos[k+1].z=0.;
                    mpos[k+2].x=0.;
                    mpos[k+2].y=row;
                    mpos[k+2].z=0.;
                    mpos[k+3].x=.4;
                    mpos[k+3].y=row;
                    mpos[k+3].z=0.;
                    mpos[k+4].x=.8;
                    mpos[k+4].y=row;
                    mpos[k+4].z=0.;
                    k += 5;
                    row += 0.4;
                }
                break;
            case 27:
                // 3D
                k=0;
                zrow = -2./3.;
                for(i=0;i<3;i++)
                {   row = -2./3.;
                    for(j=0;j<3;j++)
                    {   mpos[k].x=-2./3.;
                        mpos[k].y=row;
                        mpos[k].z=zrow;
                        mpos[k+1].x=0.;
                        mpos[k+1].y=row;
                        mpos[k+1].z=zrow;
                        mpos[k+2].x=2./3.;
                        mpos[k+2].y=row;
                        mpos[k+2].z=zrow;
                        k += 3;
                        row += 2./3.;
                    }
                    zrow += 2./3.;
                }
                break;
            default:
                throw CommonException("Invalid number of material points per element.","ElementBase::MPMPoints");
                break;
        }
    }
    
    // covert to x-y-z locations
    for(k=0;k<numPerElement;k++)
    {	ShapeFunction(&mpos[k],FALSE,fxn,NULL,NULL,NULL);
		ZeroVector(&mpos[k]);
        for(j=0;j<NumberNodes();j++)
        {   mpos[k].x+=nd[nodes[j]]->x*fxn[j];
            mpos[k].y+=nd[nodes[j]]->y*fxn[j];
            mpos[k].z+=nd[nodes[j]]->z*fxn[j];
        }
    }
}

// Get GIMP nodes around an element #num, but only where shape functions is non zero
// assumed to be properly numbered regular array
// load nodes into nds[1]... and node IDs (0-15) into ndIDs[0]...
// Elements that support GIMP must override
void ElementBase::GetGimpNodes(int *numnds,int *nds,unsigned char  *ndIDs,Vector *xipos,Vector &lp) const
{
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-15)
// Elements that support GIMP must override
void ElementBase::GimpShapeFunction(Vector *xi,int numnds,unsigned char *ndIDs,int getDeriv,double *sfxn,
										 double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-15)
// Only that elements that support axisymmetric GIMP must override
void ElementBase::GimpShapeFunctionAS(Vector *xi,int numnds,unsigned char *ndIDs,int getDeriv,double *sfxn,
									double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
}

// Ignore zero shape functions in GIMP, but only need to check remote nodes
void ElementBase::GimpCompact(int *numnds,int *nds,double *sfxn,double *xDeriv,double *yDeriv,double *zDeriv) const
{
	int local=NumberNodes();
	int i,found=local;
	for(i=local+1;i<=*numnds;i++)
	{	if(sfxn[i]>0.0)
		{	found++;
			nds[found]=nds[i];
			sfxn[found]=sfxn[i];
			if(xDeriv!=NULL)
			{	xDeriv[found]=xDeriv[i];
				yDeriv[found]=yDeriv[i];
				if(zDeriv!=NULL) zDeriv[found]=zDeriv[i];
			}
		}
	}
	*numnds=found;
}

// check if this GIMP element is on the edge of the grid
// assumes a generated 2D structured grid
bool ElementBase::OnTheEdge(void)
{	if(useGimp == POINT_GIMP) return FALSE;
	return mpmgrid.EdgeElement2D(num);
}

// Calculate CPDI nodes, optional shape functions and optional gradients
// numDnds - number of nodes in the particle domain
// cpdi[i]->ncpos - array of coordinates for corner nodes (numDnds of them)
// cpdi[i]->inElem - the element they are in
// cpdi[i]->ws - shape function weights (numNds of them)
// cpdi[i]->wg - gradient weights (numNds of them)
// throws CommonException() if too many CPDI nodes
int ElementBase::GetCPDIFunctions(int *nds,double *fn,double *xDeriv,double *yDeriv,double *zDeriv,MPMBase *mpmptr) const
{
	int i,j,numnds;
	CPDIDomain **cpdi = mpmptr->GetCPDIInfo();
    
	// Need 8X8 for linear CPDI in 3D, 8X4 for linear in 2D and 9X4 for quadratic in 2D (max is 64)
	int cnodes[64],ncnds=0;         // corner nodes and counter for those nodes
	double wsSi[64];				// hold ws * Si(xa)
	Vector wgSi[64];				// hold wg * Si(xa)
	
	// loop over the domain nodes
	for(i=0;i<numCPDINodes;i++)
	{	// get straight grid shape functions
		ElementBase *elem = theElements[cpdi[i]->inElem];
		elem->GetNodes(&numnds,nds);
		elem->ShapeFunction(&cpdi[i]->ncpos,false,&fn[1],NULL,NULL,NULL);
		
		// loop over shape grid shape functions and collect in arrays
		for(j=1;j<=numnds;j++)
		{   cnodes[ncnds] = nds[j];
			wsSi[ncnds] = cpdi[i]->ws*fn[j];
			if(xDeriv!=NULL) CopyScaleVector(&wgSi[ncnds], &cpdi[i]->wg, fn[j]);
			ncnds++;
		}
	}
	
    /*
    if(xDeriv!=NULL)
	{	cout << "Initial:" << endl;
		for(i=0;i<ncnds;i++)
		{   cout << "# node = " << cnodes[i] << ", ws*Si = " << wsSi[i] << ", wgx*Si = " << wgSi[i].x << ", wgy*Si = " << wgSi[i].y << endl;
		}
	}
    */
	
	// shell sort by node numbers in cnodes[] (always 16 for linear CPDI)
	int lognb2=(int)(log((double)ncnds)*1.442695022+1.0e-5);	// log base 2
	int k=ncnds,l,cmpNode;
	double cmpWsSi;
	Vector cmpWgSi;
	for(l=1;l<=lognb2;l++)
	{	k>>=1;		// divide by 2
		for(j=k;j<ncnds;j++)
		{	i=j-k;
			cmpNode = cnodes[j];
			cmpWsSi = wsSi[j];
			if(xDeriv!=NULL) cmpWgSi = wgSi[j];
			
			// Back up until find insertion point
			while(i>=0 && cnodes[i]>cmpNode)
			{	cnodes[i+k] = cnodes[i];
				wsSi[i+k] = wsSi[i];
				if(xDeriv!=NULL) wgSi[i+k] = wgSi[i];
				i-=k;
			}
			
			// Insert point
			cnodes[i+k]=cmpNode;
			wsSi[i+k]=cmpWsSi;
			if(xDeriv!=NULL) wgSi[i+k]=cmpWgSi;
		}
	}

    /*
 	if(xDeriv!=NULL)
    {   for(j=1;j<ncnds;j++)
        {   if(cnodes[j]<cnodes[j-1])
            {
#pragma omp critical (output)
                {   cout << "Not Sorted: " << endl;
                    for(i=0;i<ncnds;i++)
                    {   cout << "# node = " << cnodes[i] << ", ws*Si = " << wsSi[i] << ", wgx*Si = " << wgSi[i].x << ", wgy*Si = " << wgSi[i].y << endl;
                    }
                }
                break;
            }
        }
	}
    */
	
	// compact same node number
	int count = 0;
	nds[0] = -1;
	fn[0] = 1.;
	for(i=0;i<ncnds;i++)
	{   if(cnodes[i] == nds[count])
		{	// note that code never entered until all are initialized below
			fn[count] += wsSi[i];
			if(xDeriv!=NULL)
			{	xDeriv[count] += wgSi[i].x;
				yDeriv[count] += wgSi[i].y;
                zDeriv[count] += wgSi[i].z;
			}
		}
		else
		{	if(fn[count]>1.e-10)
			{	count++;		// keep previous one only if shape ends up nonzero (otherwise previous one is removed as too small)
				if(count>=maxShapeNodes)
                {	cout << "# Found " << count-1 << " nodes; need room for remaining nodes:" << endl;
                    for(j=i;j<ncnds;j++)
                    {   cout << "#   node = " << cnodes[j] << ", ws*Si = " << endl;
                    }
					throw CommonException("Too many CPDI nodes found; increase maxShapeNodes in source code by at least number of remaining nodes","ElementBase::GetCPDIFunctions");
				}
			}
			nds[count] = cnodes[i];
			fn[count] = wsSi[i];
			if(xDeriv!=NULL)
			{	xDeriv[count] = wgSi[i].x;
				yDeriv[count] = wgSi[i].y;
                zDeriv[count] = wgSi[i].z;
			}
		}
	}
	
	// remove last one if it was too small
	if(fn[count]<=1.e-10) count--;
	
	return count;
}

// by non-element methods that need access to grid shape functions only, and those methods are protected
void ElementBase::GridShapeFunctions(int *numnds,int *nds,Vector *xipos,double *fn) const
{
    GetNodes(numnds,nds);
    ShapeFunction(xipos,FALSE,&fn[1],NULL,NULL,NULL);
}

#pragma mark ElementBase: MPM Only Methods

/* Get list of nodes for this element and the number of nodes
	nodes will be in nds[1]...
*/
void ElementBase::GetNodes(int *numnds,int *nds) const
{
	int i;
	*numnds = NumberNodes();
	for(i=1;i<=*numnds;i++)
	   nds[i] = nodes[i-1];
}

/* Load nodal coordinates into eNodes[1]...
	input numnds and nds from previous call to GetNodes()
*/
void ElementBase::GetCoordinates(Vector *eNodes,int numnds,int *nds) const
{
	int i;
	for(i=1;i<=numnds;i++)
	{	eNodes[i].x=nd[nds[i]]->x;
		eNodes[i].y=nd[nds[i]]->y;
		eNodes[i].z=nd[nds[i]]->z;
	}
}

/*  Give dimensionless coordinates for central point
    of element to provide initial guess for
    numerical solution to find natural coordinates.
	Only needed for non-rectangular elements
*/
void ElementBase::GetCentroid(Vector *xipos) const
{
    xipos->x=0.;
    xipos->y=0.;
    xipos->z=0.;
}

/*  Find node of this element closest to specified point
    Only checks nodes of this element. In distorted mesh,
    could be closer node in neighboring element, even if
    this point is in this element
*/
int ElementBase::NearestNode(double x,double y,int *secondChoice)
{	int ind,nearest,nextNearest;
    double dist,distMin,distNextMin,dltx,dlty;
    int k,numnds=NumberNodes();
    
    // distance to first node
	nearest=nodes[0];
    dltx=nd[nearest]->x-x;
    dlty=nd[nearest]->y-y;
    distMin=dltx*dltx+dlty*dlty;
    
    // distance to second node
	nextNearest=nodes[1];
	dltx=nd[nextNearest]->x-x;
	dlty=nd[nextNearest]->y-y;
	distNextMin=dltx*dltx+dlty*dlty;
	if(distNextMin<distMin)
	{	dist=distMin;
		distMin=distNextMin;
		distNextMin=distMin;
		nearest=nodes[1];
		nextNearest=nodes[0];
	}
	
    // check the rest
    for(k=2;k<numnds;k++)
    {	ind=nodes[k];
    	dltx=nd[ind]->x-x;
        dlty=nd[ind]->y-y;
        dist=dltx*dltx+dlty*dlty;
        if(dist<distMin)
		{	distNextMin=distMin;
			nextNearest=nearest;
			distMin=dist;
            nearest=ind;
        }
		else if(dist<distNextMin)
		{	distNextMin=dist;
			nextNearest=ind;
		}
    }
    
    // return the node
	*secondChoice=nextNearest;
    return nearest;
}

/* find next node (in ccw direction) from specified node
    Assumes nodes[NumberNodes()]=nodes[0]
    return 0 if gridNode not in this element
*/
int ElementBase::NextNode(int gridNode)
{	int i;
    for(i=0;i<NumberNodes();i++)
    {	if(nodes[i]==gridNode)
            return nodes[i+1];
    }
    return 0;
}

// for structured grid, return 0-terminated list of neighbors
void ElementBase::GetListOfNeighbors(int *theList) { mpmgrid.ListOfNeighbors2D(num,theList); }

// If needed, create the neighbors array (only done for 2D with cracks)
void ElementBase::AllocateNeighborsArray(void)
{	neighbors = new int[NumberSides()];
	int i;
	for(i=0;i<NumberSides();i++)
		neighbors[i]=UNKNOWN_NEIGHBOR;
}

/* Find element (0 based) that touches the element face stating in gridNode
    First time called, search other elements for the edge, but then store result
    Assumes nodes[NumberNodes()]=nodes[0]
    return NO_NEIGHBOR (-1) if no neighbor (i.e., on edge of grid)
    return UNKNOWN_NEIGHBOR (-2) if gridNode not in this element
*/
int ElementBase::Neighbor(int gridNode)
{
    int i,edge=-1;
    
    // first find which edge (0 to NumberSides()-1) to search?
    for(i=0;i<NumberSides();i++)
    {	if(nodes[i]==gridNode)
        {   edge=i;
            break;
        }
    }
    
    // if edge not found and error because the gridNode is not even in this element
    if(edge<0) return UNKNOWN_NEIGHBOR;
    
    // search for neighbor if needed
    if(neighbors[edge]==UNKNOWN_NEIGHBOR)
    {	int nextNode=nodes[edge+1];
        neighbors[edge]=NO_NEIGHBOR;
        
        // search for element with nextNode,gridNode in ccw direction
        for(i=0;i<nelems;i++)
        {   if(theElements[i]->FindEdge(nextNode,gridNode)>=0)
            {	neighbors[edge]=i;
                break;
            }
        }
    }
    
    // return result
    return neighbors[edge];
}

/* find an edge (in ccw direction) with a pair of nodes
    Assumes nodes[NumberNodes()]=nodes[0]
    return edge number (0 based) or -1 if no such edge
*/
int ElementBase::FindEdge(int beginNode,int endNode)
{
    int i;
    for(i=0;i<NumberNodes();i++)
    {	if(nodes[i]==beginNode && nodes[i+1]==endNode)
            return i;
    }
    return -1;
}

// return TRUE or false if Orthogonal in x-y-z planes and if orthogonal then find width, height, and depth
// of the element. If 2D element, return depth=0. Only for MPM elements.
int ElementBase::Orthogonal(double *dx,double *dy,double *dz) { return FALSE; }

// Find Cartesion position from natural coordinates
void ElementBase::GetPosition(Vector *xipos,Vector *xyzpos)
{
	double fn[MaxElNd];
	ShapeFunction(xipos,FALSE,&fn[1],NULL,NULL,NULL);
	int i;
	ZeroVector(xyzpos);
	for(i=1;i<=NumberNodes();i++)
	{	xyzpos->x += fn[i]*nd[nodes[i-1]]->x;
		xyzpos->y += fn[i]*nd[nodes[i-1]]->y;
		xyzpos->z += fn[i]*nd[nodes[i-1]]->z;
	}
}

#pragma mark CLASS METHODS

// zero all velocity fields at start of time step
void ElementBase::AllocateNeighbors(void)
{	int i;
    for(i=0;i<nelems;i++)
        theElements[i]->AllocateNeighborsArray();
}

// on start up initialize number of CPDI nodes (if needed)
// done here in case need more initializations in the future
// throws CommonException() if CPDI type is not allowed
void ElementBase::InitializeCPDI(bool isThreeD)
{
    if(isThreeD)
    {   if(useGimp == LINEAR_CPDI)
        {   numCPDINodes = 8;
        }
        else if(useGimp == QUADRATIC_CPDI)
        {	throw CommonException("qCPDI is not yet implemented for 3D (use lCPDI instead).","ElementBase::InitializeCPDI");
        }
    }
    else
    {   if(useGimp == LINEAR_CPDI || useGimp==LINEAR_CPDI_AS)
        {   numCPDINodes = 4;
        }
        else if(useGimp == QUADRATIC_CPDI)
        {   numCPDINodes = 9;
        }
    }
}



