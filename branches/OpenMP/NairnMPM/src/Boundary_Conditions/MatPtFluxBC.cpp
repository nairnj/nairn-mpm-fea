/********************************************************************************
    MatPtFluxBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Mar 17 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_XML/mathexpr.hpp"
#include "Nodes/NodalPoint.hpp"

// global
MatPtFluxBC *firstFluxPt=NULL;

/*******************************************************************
	MatPtFluxBC: Constructors and Destructors
*******************************************************************/

// Constructors
MatPtFluxBC::MatPtFluxBC(int num,int dof,int sty,int edge) : MatPtLoadBC(num,dof,sty)
{
	face = edge;
}

/*******************************************************************
	MatPtFluxBC: Methods
*******************************************************************/

// print to output
BoundaryCondition *MatPtFluxBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,value,ftime);
    os << nline;
	PrintFunction(os);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
// (only called when diffusion is active)
MatPtFluxBC *MatPtFluxBC::AddMPFlux(double bctime)
{
    // condition value
	// flux BC in kg/(m^2-sec) - find J/rho in units of mm/sec
	MPMBase *mpmptr = mpm[ptNum-1];
    MaterialBase *matptr = theMaterials[mpmptr->MatID()];
	double csatrho = matptr->rho*matptr->concSaturation/mpmptr->GetRelativeVolume();
	double fluxMagX,fluxMagY=0.;
    int bcDir=X_DIRECTION;
	
	if(style==SILENT)
	{	// silent assumes isotropic material
		TransportProperties t;
		matptr->GetTransportProps(mpm[ptNum-1],fmobj->np,&t);
		Tensor *D = &(t.diffusionTensor);
        
        // get in mm/sec
		fluxMagX = D->xx*mpmptr->pDiffusion->Dc.x + D->xy*mpmptr->pDiffusion->Dc.y;
        fluxMagY = D->xy*mpmptr->pDiffusion->Dc.x + D->yy*mpmptr->pDiffusion->Dc.y;
        bcDir = N_DIRECTION;
	}
	else if(direction==EXTERNAL_FLUX)
	{	double mstime=1000.*bctime;
		fluxMagX = BCValue(mstime)/csatrho;
	}
	else
    {   // coupled surface flux and ftime is bath concentration
		varTime = mpmptr->pConcentration-ftime;
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		double currentValue = fabs(function->Val());
		if(varTime>0.) currentValue=-currentValue;
		fluxMagX = currentValue/csatrho;
	}
	
	// get corners and direction from material point
	int cElem[4],numDnds;
	Vector corners[4],tscaled;
	double ratio = mpmptr->GetTractionInfo(face,bcDir,cElem,corners,&tscaled,&numDnds);
	
    // compact CPDI nodes into list of nodes (nds) and final shape function term (fn)
    // May need up to 8 (in 3D) for each of the numDnds (2 in 2D or 4 in 3D)
    int nds[8*numDnds+1];
    double fn[8*numDnds+1];
    int numnds = CompactCornerNodes(numDnds,corners,cElem,ratio,nds,fn);
	
    // add force to each node
	int i;
    for(i=1;i<=numnds;i++)
    {   // skip empty nodes
        if(nd[nds[i]]->NumberNonrigidParticles())
		{	nd[nds[i]]->fdiff += (fluxMagX*tscaled.x + fluxMagY*tscaled.y)*fn[i];
        }
    }
	
    return (MatPtFluxBC *)GetNextObject();
}
