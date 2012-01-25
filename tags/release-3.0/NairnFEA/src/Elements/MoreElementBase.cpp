/********************************************************************************
    MoreElementBase.cpp
    NairnFEA
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnFEA_Class/NairnFEA.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/CommonException.hpp"

// better to not need this by moving the stuff to element method
#include "Elements/FourNodeIsoparam.hpp"
#include "Elements/EightNodeIsoparam.hpp"
#include "Elements/SixNodeTriangle.hpp"

// Quadrature points in various sets
// [0] is for 2X2 quadrilaterals (4 points)
// [1] is for 3 points in triangle
// [2] is for 3X3 quadrilaterals (9 points)
// define maximum gauss points + 1
#define MaxGaussPts 10
static double placeXi[3][9] = {
	{-0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.,0.,0.,0.,0.},
	{0.666666666666667,0.166666666666667,0.166666666666667,0.,0.,0.,0.,0.,0.},
	{-0.7745966692414834,-0.7745966692414834,-0.7745966692414834,0.,0.,0.,0.7745966692414834,0.7745966692414834,0.7745966692414834}
};
static double placeEta[3][9] = {
	{-0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,0.,0.,0.,0.,0.},
	{0.166666666666667,0.666666666666667,0.166666666666667,0.,0.,0.,0.,0.,0.},
	{-0.7745966692414834,0.,0.7745966692414834,-0.7745966692414834,0.,0.7745966692414834,-0.7745966692414834,0.,0.7745966692414834}
};
static double weight[3][9] = {
	{1.,1.,1.,1.,0.,0.,0.,0.,0.},
	{0.166666666666667,0.166666666666667,0.166666666666667,0.,0.,0.,0.,0.,0.},
	{0.3086419753086420,0.4938271604938272,0.3086419753086420,0.4938271604938272,0.7901234567901235,0.4938271604938272,
		0.3086419753086420,0.4938271604938272,0.3086419753086420}
};

// globals for common element storage while building stiffness matrix or
//    while evaluating results
Vector ce[MaxElNd];
double te[MaxElNd];
double re[MxFree*MaxElNd];
double se[MxFree*MaxElNd][MxFree*MaxElNd];

#pragma mark ElementBase: Constructors and Destructor FEA Only

/* Main FEA contructor when creating elements:
	input is node numbers (0-based array) but values
            are 1-based (length always MaxElNd, even if unused)
*/
ElementBase::ElementBase(int eNum,int *eNode,int eMat,double eAng)
{	int i;
    num=eNum;
    for(i=0;i<MaxElNd;i++) nodes[i]=eNode[i];
    filled=0;
    material=eMat;
	SetAngleInDegrees(eAng);
	numGauss=4;
	gaussSet=0;
}

ElementBase::~ElementBase()
{
}

#pragma mark More ElementBase Methods

// remap nodes if requenced
void ElementBase::MapNodes(int *revMap)
{
	int i;
	for(i=0;i<NumberNodes();i++)
		nodes[i]=revMap[nodes[i]];
}

// Element stiffness matrix (override if doesn't fit)
void ElementBase::Stiffness(int np) { IsoparametricStiffness(np); }

// Find forces and stresses in element (override if doesn't fit)
void ElementBase::ForceStress(double *rm,int np,int nfree)
{	IsoparametricForceStress(rm,np,nfree); }

/*	Calculate Stiffness Matrix
	Generalized here for any isoparametric element using
            Gaussian Quadrature integration with numgaus
            points stored in placeXi, placeEta and weight arrays.
*/
void ElementBase::IsoparametricStiffness(int np)
{
    int numnds=NumberNodes();
    int nx,ind1,ind2,i,j,irow,jcol,nst=2*numnds;
    double dv,bte[2*MaxElNd-1][5],temp;
    double xiDeriv[MaxElNd],etaDeriv[MaxElNd],asbe[MaxElNd],fn[MaxElNd];
    double detjac,asr,deltaT;
    MaterialBase *matl=theMaterials[material-1];
	Vector place;
	
    // Load nodal coordinates (ce[]), temperature (te[]), and
	//    material props (mdm[][] and me0[])
    GetProperties(np);
    
    // Zero upper se[][] and re[]
	ZeroUpperHalfStiffness();

    // Gaussian Quadrature integration over xi and neta
    for(nx=0;nx<numGauss;nx++)
    {	/* Call shape routine to calculate element B (in xiDeriv, etaDeriv, and asbe)
                matrix and the determinant of the Jacobian - both at pxi,pet */
		place.x=placeXi[gaussSet][nx];
		place.y=placeEta[gaussSet][nx];
        ShapeFunction(&place,BMATRIX,&fn[1],&xiDeriv[1],&etaDeriv[1],&ce[1],&detjac,&asr,&asbe[1]);

        /* Evaluate volume element */
        if(np!=AXI_SYM)
            dv=weight[gaussSet][nx]*GetThickness()*detjac/1000.;
        else
            dv=weight[gaussSet][nx]*asr*detjac;		// force per radian, hence no 2 pi

        /* Form matrix product BT E - exploit known sparcity of
                B matrix and only include multiplications by nonzero elements */
        ind1=-1;
		deltaT=0.;
        if(np!=AXI_SYM)
        {   for(i=1;i<=numnds;i++)
            {   ind1=ind1+2;
                ind2=ind1+1;
                for(j=1;j<=3;j++)
                {   bte[ind1][j]=xiDeriv[i]*matl->mdm[1][j]+etaDeriv[i]*matl->mdm[3][j];
                    bte[ind2][j]=etaDeriv[i]*matl->mdm[2][j]+xiDeriv[i]*matl->mdm[3][j];
                }
				deltaT+=te[i]*fn[i];
            }
        }
        else
        {   for(i=1;i<=numnds;i++)
            {   ind1=ind1+2;
                ind2=ind1+1;
                for(j=1;j<=4;j++)
                {   bte[ind1][j]=xiDeriv[i]*matl->mdm[1][j]+etaDeriv[i]*matl->mdm[3][j]
                                                    +asbe[i]*matl->mdm[4][j];
                    bte[ind2][j]=etaDeriv[i]*matl->mdm[2][j]+xiDeriv[i]*matl->mdm[3][j];
                }
				deltaT+=te[i]*fn[i];
            }
        }

        /* Form stiffness matrix by getting BT E B and add in initial
                strains into element load vector */
        for(irow=1;irow<=nst;irow++)
        {   for(j=1;j<=3;j++)
                re[irow]+=bte[irow][j]*matl->me0[j]*deltaT*dv;
            
            if(np!=AXI_SYM)
            {   for(jcol=irow;jcol<=nst;jcol++)
                {   if(IsEven(jcol))
                    {   ind1=jcol/2;
                        temp=bte[irow][2]*etaDeriv[ind1]
                                        +bte[irow][3]*xiDeriv[ind1];
                    }
                    else
                    {   ind1=(jcol+1)/2;
                        temp=bte[irow][1]*xiDeriv[ind1]
                                        +bte[irow][3]*etaDeriv[ind1];
                    }
                    se[irow][jcol]+=temp*dv;
                }
            }
            else
            {   re[irow]+=bte[irow][4]*matl->me0[4]*deltaT*dv;
                for(jcol=irow;jcol<=nst;jcol++)
                {   if(IsEven(jcol))
                    {   ind1=jcol/2;
                        temp=bte[irow][2]*etaDeriv[ind1]
                                        +bte[irow][3]*xiDeriv[ind1];
                    }
                    else
                    {   ind1=(jcol+1)/2;
                        temp=bte[irow][1]*xiDeriv[ind1]
                                        +bte[irow][3]*etaDeriv[ind1]
                                        +bte[irow][4]*asbe[ind1];
                    }
                    se[irow][jcol]+=temp*dv;
                }
            }
        }
    }   // End quadrature loop
    
    // Fill in lower half of stiffness matrix
	FillLowerHalfStiffness();
}

// Zero upper half of stiffness matrix (se[][]) reaction vector (re[])
void ElementBase::ZeroUpperHalfStiffness(void)
{
	int irow,jcol,nst=2*NumberNodes();
	
    for(irow=1;irow<=nst;irow++)
    {   re[irow]=0.;
        for(jcol=irow;jcol<=nst;jcol++)
            se[irow][jcol]=0.;
    }
}

// Fill lower half of stiffness matrix
void ElementBase::FillLowerHalfStiffness(void)
{
	int irow,jcol,nst=2*NumberNodes();
	
	for(irow=1;irow<=nst-1;irow++)
	{	for(jcol=irow+1;jcol<=nst;jcol++)
		{	se[jcol][irow]=se[irow][jcol];
		}
	}
}

/*	Calculate Forces and Stresses
	Generalized here for any isoparametric element using
            Gaussian Quadrature integration with numgaus
            points stored in placeXi, placeEta and weight arrays.
*/
void ElementBase::IsoparametricForceStress(double *rm,int np,int nfree)
{
    int numnds=NumberNodes(),ind,j,i,nst=2*numnds;
    int ngp,nx,ind1,ind2;
    int indg;
    double sgp[MaxGaussPts][5],etot[5];
    double xiDeriv[MaxElNd],etaDeriv[MaxElNd],asbe[MaxElNd],fn[MaxElNd];
    double detjac,asr,temp,dv,deltaT;
	
    double thck=GetThickness()/1000.;
    MaterialBase *matl=theMaterials[material-1];
	Vector place;
	
    // Load element coordinates (ce[]), noodal temperature (te[]), 
	//    and material props (mdm[][] and me0[])
    GetProperties(np);
    
    // Load nodal displacements into re[]
    ind=0;
    for(j=1;j<=numnds;j++)
    {	indg=nfree*(nodes[j-1]-1);
        for(i=1;i<=nfree;i++)
            re[++ind]=rm[indg+i];
    }

    // zero force at each degree of freedom (stored in se[i][7])
    for(i=1;i<=nst;i++)
        se[i][7]=0.;
        
    // element strainEnergy
    strainEnergy=0.;

    // Zero stress and force vectors to hold results at Gauss points
    for(i=1;i<=numGauss;i++)
    {	for(j=1;j<=4;j++)
            sgp[i][j]=0.;
    }

    // Gaussian Quadrature integration over xi and neta
    for(nx=0;nx<numGauss;nx++)
    {	ngp=nx+1;

        /* Call shape routine to calculate element B (in xiDeriv, etaDeriv, and asbe)
                    matrix and the determinant of the Jacobian - both at pxi,pet */
		place.x=placeXi[gaussSet][nx];
		place.y=placeEta[gaussSet][nx];
        ShapeFunction(&place,BMATRIX,&fn[1],&xiDeriv[1],&etaDeriv[1],&ce[1],&detjac,&asr,&asbe[1]);
			
        // Evaluate volume element
        if(np!=AXI_SYM)
            dv=weight[gaussSet][nx]*thck*detjac;
        else
            dv=weight[gaussSet][nx]*asr*detjac;

        /* Evaluate etot=(B d - e0). In forming
                B d, avoid multiplications by zero */
		deltaT=0.;
		for(i=1;i<=numnds;i++) deltaT+=te[i]*fn[i];
        etot[1]=-matl->me0[1]*deltaT;
        etot[2]=-matl->me0[2]*deltaT;
        etot[3]=-matl->me0[3]*deltaT;
        etot[4]=-matl->me0[4]*deltaT;
        ind1=-1;
        for(i=1;i<=numnds;i++)
        {   ind1=ind1+2;
            ind2=ind1+1;
            etot[1]+=xiDeriv[i]*re[ind1];
            etot[2]+=etaDeriv[i]*re[ind2];
            etot[3]+=etaDeriv[i]*re[ind1]+xiDeriv[i]*re[ind2];
            if(np==AXI_SYM)
                etot[4]+=asbe[i]*re[ind1];
        }

        // Multply by stiffness matrix: sig = dm * etot
        if(np!=AXI_SYM)
        {   for(i=1;i<=3;i++)
            {	for(j=1;j<=3;j++)
                    sgp[ngp][i]+=matl->mdm[i][j]*etot[j];
            }
        }
        else
        {   for(i=1;i<=4;i++)
            {	for(j=1;j<=4;j++)
                    sgp[ngp][i]+=matl->mdm[i][j]*etot[j];
            }
        }

        /* Add terms for getting integral(BT sig) which gives nodal forces
                2π added for axisymmetric forces at end of subroutine */
        ind1=-1;
        for(i=1;i<=numnds;i++)
        {   ind1=ind1+2;
            ind2=ind1+1;
            temp=xiDeriv[i]*sgp[ngp][1]+etaDeriv[i]*sgp[ngp][3];
            if(np==AXI_SYM) temp+=asbe[i]*sgp[ngp][4];
            se[ind1][7]+=temp*dv;
            temp=etaDeriv[i]*sgp[ngp][2]+xiDeriv[i]*sgp[ngp][3];
            se[ind2][7]+=temp*dv;
        }

        // Get initial/thermal strain contribution to strain energy
        temp=0.;
        for(i=1;i<=3;i++) temp+=sgp[ngp][i]*matl->me0[i]*deltaT;
        if(np==AXI_SYM) temp+=sgp[ngp][4]*matl->me0[4]*deltaT;
        strainEnergy-=0.5*temp*dv;

        /* When plane strain account for constrained 1D shrinkage effect
                        on strain energy */
        if(np==PLANE_STRAIN)
            strainEnergy+=0.5*matl->mdm[4][4]*dv*deltaT*deltaT;
            
    } // End of quadrature loop
	
    // Add 1/2 Fd to strain energy
    temp=0.;
    for(i=1;i<=nst;i++) temp+=re[i]*se[i][7];
    strainEnergy+=0.5*temp;

    // Extrapolate gaussian point stresses to nodal point stresses
	ExtrapolateGaussStressToNodes(sgp);

    /* For plane strain analysis, calculate sigz stress
            For axisymmetric, multiply force and energy by 2 pi
			FEA better have thermal.reference=0 */
    if(np==PLANE_STRAIN)
    {   for(i=1;i<=numnds;i++)
			se[i][4]=matl->GetStressStrainZZ(se[i][1],se[i][2],se[i][3],te[i],angle,np);
    }
    else if(np==AXI_SYM)
    {   for(i=1;i<=nst;i++)
            se[i][7]*=2.*PI_CONSTANT;
        strainEnergy*=2*PI_CONSTANT;
    }
}

// Take stress at gauss points and map them to element nodes
// Elements must override to support stress calculations
// sgp[i][j] is stress j (1 to 4) at Gauss point i
// se[i][j] is output stress j (1 to 4) at node i (1 to numnds) (externed variable)
void ElementBase::ExtrapolateGaussStressToNodes(double sgp[][5])
{
	int i,j,numnds = NumberNodes();
	for(i=1;i<=numnds;i++)
	{	for(j=1;j<=4;j++)
			se[i][j]=0.;
	}
}
	
// Determine if the element should be output
int ElementBase::WantElement(char thisFlag,const vector< int > &selected)
{
    // automatic result
    if(thisFlag=='N') return FALSE;
    if(thisFlag=='Y') return TRUE;
	if(selected.size()==0) return FALSE;

    // Check if any wanted node is present in element
	int i,look;
	unsigned int j;
    for(i=1;i<=NumberNodes();i++)
    {	look=nodes[i-1];
        for(j=0;j<selected.size();j++)
        {   if(look==selected[j]) return TRUE;
			if(look<selected[j]) break;
        }
    }
    return FALSE;
}

/* Calculate equivalent forces for edge load on linear edges only
	nd1 and nd2 are nodes at corners of the edge (1 based)
	ndir is 1 for normal, 2 for shear
	fload[0] and fload[1] are nodal stresses in MPa
	re is output resultants (1 based)
	(See theory nodes in CalcEdgeLoads())
*/
void ElementBase::LinearEdgeLoad(int nd1,int nd2,int ndir,double *fload,double *re,int np)
{
	int ind1,ind2;
	double r1,r2,delx,dely,arg1,arg2;
	
    // Load nodal coordinates (in m)
    for(ind2=1;ind2<=NumberNodes();ind2++)
    {   ind1=nodes[ind2-1];
        ce[ind2].x=nd[ind1]->x/1000.;
        ce[ind2].y=nd[ind1]->y/1000.;
    }
	
	// zero re[] vector
	int nst=2*NumberNodes();
	for(ind1=1;ind1<=nst;ind1++) re[ind1]=0.;
	
	// thickness (in m)
	double thck=GetThickness()/1000.;
	
	// Handle planar elements
	if(np!=AXI_SYM)
	{	delx=(ce[nd2].x-ce[nd1].x)*thck/2.;		// i.e., delx t/2
		dely=(ce[nd2].y-ce[nd1].y)*thck/2.;     // i.e., dely t/2
		arg1=(2.*fload[0]+fload[1])/3.e-6;
		arg2=(fload[0]+2.*fload[1])/3.e-6;
	}

	// Handle axisymmetric elements
	else
	{	r1=ce[nd1].x;
		r2=ce[nd2].x;
		delx=(ce[nd2].x-ce[nd1].x)/2.;		// i.e., delr/2
		dely=(ce[nd2].y-ce[nd1].y)/2.;		// i.e., delz/2
		arg1=r1*fload[0]/3.e-6;
		arg1+=(r1+r2)*(fload[0]+fload[1])/6.e-6;
		arg2=r2*fload[1]/3.e-6;
		arg2+=(r1+r2)*(fload[0]+fload[1])/6.e-6;
	}

	// DOF index into output array
	ind1=2*nd1-1;
	ind2=2*nd2-1;
	
	// Handle normal stress condition
	if(ndir==1)
	{	re[ind1]=dely*arg1;
		re[ind1+1]=-delx*arg1;
		re[ind2]=dely*arg2;
		re[ind2+1]=-delx*arg2;
	}

	// Handle shear stress condition
	else
	{	re[ind1]=delx*arg1;
		re[ind1+1]=dely*arg1;
		re[ind2]=delx*arg2;
		re[ind2+1]=dely*arg2;
	}
}

/* Calculate equivalent forces for edge load on quadratic edges only
	nd1, nd2, and nd3 are element node #s along the edges (1 based)
	ndir is 1 for normal, 2 for shear
	fload[0] to fload[2] are nodal stresses in MPa
	re is output resultants (1 based)
	(See theory nodes in CalcEdgeLoads())
*/
void ElementBase::QuadEdgeLoad(int nd1,int nd2,int nd3,int ndir,double *fload,double *re,int np)
{
	double delx,delxh,dely,delyh,delr,delrh,delz,delzh;
	double r1,r2,r3,stof[7][4];
	int i,j,ind[7];
	
    // Load nodal coordinates (in m) in 1-based in cd[]
    for(i=1;i<=NumberNodes();i++)
    {   j=nodes[i-1];
        ce[i].x=nd[j]->x/1000.;
        ce[i].y=nd[j]->y/1000.;
    }
	
	// zero re[] vector
	int nst=2*NumberNodes();
	for(i=1;i<=nst;i++) re[i]=0.;
	
	// thickness (in m)
	double thck=GetThickness()/1000.;
	
	// deltas
	if(ndir==1)
	{	delx=(ce[nd3].x-ce[nd1].x);
		delxh=(ce[nd3].x+ce[nd1].x-2.*ce[nd2].x);
		dely=(ce[nd3].y-ce[nd1].y);
		delyh=(ce[nd3].y+ce[nd1].y-2.*ce[nd2].y);
	}
	else
	{	dely=(ce[nd3].x-ce[nd1].x);
		delyh=(ce[nd3].x+ce[nd1].x-2.*ce[nd2].x);
		delx=-(ce[nd3].y-ce[nd1].y);
		delxh=-(ce[nd3].y+ce[nd1].y-2.*ce[nd2].y);
	}

	// Planar elements
	if(np!=AXI_SYM)
	{	delx=delx*thck/30.;
		dely=dely*thck/30.;
		delxh=delxh*thck/30.;
		delyh=delyh*thck/30.;
		stof[1][1]=4.*dely-6.*delyh;
		stof[1][2]=2.*dely-4*delyh;
		stof[1][3]=-dely;
		stof[2][1]=-4.*delx+6.*delxh;
		stof[2][2]=-2.*delx+4*delxh;
		stof[2][3]=+delx;
		stof[3][1]=2.*dely-4.*delyh;
		stof[3][2]=16*dely;
		stof[3][3]=2*dely+4*delyh;
		stof[4][1]=-2.*delx+4.*delxh;
		stof[4][2]=-16*delx;
		stof[4][3]=-2*delx-4*delxh;
		stof[5][1]=-dely;
		stof[5][2]=2.*dely+4.*delyh;
		stof[5][3]=4.*dely+6.*delyh;
		stof[6][1]=+delx;
		stof[6][2]=-2.*delx-4.*delxh;
		stof[6][3]=-4.*delx-6.*delxh;
	}

	// Handle axisymmetric elements
	else
	{	delr=delx/420.;
		delz=dely/420.;
		delrh=delxh/210.;
		delzh=delyh/210.;
		r1=ce[nd1].x;
		r2=ce[nd2].x;
		r3=ce[nd3].x;
		stof[1][1]=delz*(39.*r1+20.*r2-3.*r3)
					+delzh*(-33.*r1-12.*r2+3.*r3);
		stof[1][2]=delz*(20.*r1+16.*r2-8.*r3)
					+delzh*(-12.*r1-16.*r2);
		stof[1][3]=delz*(-3.*r1-8.*r2-3.*r3)
					+delzh*(3.*r1-3.*r3);
		stof[2][1]=-delr*(39.*r1+20.*r2-3.*r3)
					-delrh*(-33.*r1-12.*r2+3.*r3);
		stof[2][2]=-delr*(20.*r1+16.*r2-8.*r3)
					-delrh*(-12.*r1-16.*r2);
		stof[2][3]=-delr*(-3.*r1-8.*r2-3.*r3)
					-delrh*(3.*r1-3.*r3);
		stof[3][1]=delz*(20.*r1+16.*r2-8.*r3)
					+delzh*(-12.*r1-16.*r2);
		stof[3][2]=delz*(16.*r1+192.*r2+16.*r3)
					+delzh*(-16.*r1+16.*r3);
		stof[3][3]=delz*(-8.*r1+16.*r2+20.*r3)
					+delzh*(16.*r2+12.*r3);
		stof[4][1]=-delr*(20.*r1+16.*r2-8.*r3)
					-delrh*(-12.*r1-16.*r2);
		stof[4][2]=-delr*(16.*r1+192.*r2+16.*r3)
					-delrh*(-16.*r1+16.*r3);
		stof[4][3]=-delr*(-8.*r1+16.*r2+20.*r3)
					-delrh*(16.*r2+12.*r3);
		stof[5][1]=delz*(-3.*r1-8.*r2-3.*r3)
					+delzh*(3.*r1-3.*r3);
		stof[5][2]=delz*(-8.*r1+16.*r2+20.*r3)
					+delzh*(16.*r2+12.*r3);
		stof[5][3]=delz*(-3.*r1+20.*r2+39.*r3)
					+delzh*(-3.*r1+12.*r2+33.*r3);
		stof[6][1]=-delr*(-3.*r1-8.*r2-3.*r3)
					-delrh*(3.*r1-3.*r3);
		stof[6][2]=-delr*(-8.*r1+16.*r2+20.*r3)
					-delrh*(16.*r2+12.*r3);
		stof[6][3]=-delr*(-3.*r1+20.*r2+39.*r3)
					-delrh*(-3.*r1+12.*r2+33.*r3);
	}

	// Matrix multiply to get forces
	ind[1]=2*nd1-1;
	ind[2]=ind[1]+1;
	ind[3]=2*nd2-1;
	ind[4]=ind[3]+1;
	ind[5]=2*nd3-1;
	ind[6]=ind[5]+1;
	for(i=1;i<=6;i++)
	{	for(j=1;j<=3;j++)
		{	re[ind[i]]+=stof[i][j]*fload[j-1]*1.e6;
		}
	}
}

/*	Edge loads - override in elements that allow stress boundary conditions
	Theory: Let Ni(s) be shape functions to parameterize position and stresses
		along the edge with dimensionless coordinates from -1 to 1
		for i nodes on the edge, then
			(Fx)i = Integral(-1 to 1 ds,
		           Ni(s) t [y'(s) Sum(k,Nk(s)sigma(k)) + x'(s) Sum(k,Nk(s)tau(k))] )
			(Fy)i = Integral(-1 to 1 ds,
		           Ni(s) t [-x'(s) Sum(k,Nk(s)sigma(k)) + y'(s) Sum(k,Nk(s)tau(k))] )
		where sigma(k) and tau(k) are normal and tangential stresses at the nodes
		and the derivatives are from
			x'(s) = Sum(k, Nk'(s) x(k))  and  y'(s) = Sum(k, Nk'(s) y(k))
		For axisymmetric calculations get force per radian by replacing t by
			t -> (2 pi r(s))/(2 pi) = Sum(k, Nk(s) r(k))
*/
void ElementBase::CalcEdgeLoads(double *re,int iedge,int ndir,double *fload,int np)
{
	throw CommonException("Element does not allow stress boundary conditions.",
                                "ElementBase::ForcesOnEdges()");
}

// override if can have quarter point elements
void ElementBase::MakeQuarterPointNodes(int crackTip,vector<int> &movedNodes) {}

#pragma mark More ElementBase Accessors

/* Load element properties into globals for this element (for FEA only)
    1. Get nodal coordinates and temperature in ce[] and te[] arrays
    2. Make sure material mdm[][] and me0[] are defined
*/
void ElementBase::GetProperties(int np)
{
    int i,ind;
    
    // Load nodal coordinates (in m)
    for(i=1;i<=NumberNodes();i++)
    {   ind=nodes[i-1];
        ce[i].x=nd[ind]->x/1000.;
        ce[i].y=nd[ind]->y/1000.;
		te[i]=nd[ind]->gTemperature;
    }
    
    /* Get mechanical properties that depend on angle
            Currently special case for orthotropic */
    theMaterials[material-1]->LoadMechProps(FALSE,angle,np);
}

// does element have this node?
bool ElementBase::HasNode(int nodeNum)
{	// check real nodes
	int i;
    for(i=0;i<NumberNodes();i++)
    {	if(nodes[i]==nodeNum)
			return TRUE;
    }
	return FALSE;
}

// When deleting a node, decrement node numbers that are larger
void ElementBase::DecrementNodeNums(int nodeNum)
{	int i;
    for(i=0;i<NumberNodes();i++)
    {	if(nodes[i]>nodeNum)
			nodes[i]--;
    }
}

/* find maximum difference between nodes in the element
    This difference used to find FEA bandwidth
*/
void ElementBase::MaxMinNode(int *maxn,int *minn)
{
	// find new max and min (when 0) or build on previous one
	if(*maxn==0)
    {	*maxn=nodes[0];
		*minn=nodes[0];
	}
    
	// check real nodes
    int i;
    for(i=0;i<NumberNodes();i++)
    {	*maxn=fmax(nodes[i],*maxn);
        *minn=fmin(nodes[i],*minn);
    }
}

// Bulk FEA element (some elements may override with no)
bool ElementBase::BulkElement(void) { return TRUE; }

// element material angle
void ElementBase::SetAngleInDegrees(double eAng) { angle=PI_CONSTANT*eAng/180.; }
double ElementBase::GetAngleInDegrees(void) { return 180.*angle/PI_CONSTANT; }

#pragma mark CLASS METHODS

/* Move all mid side nodes near crack tip to the 1/4 location
	closer to the crack tip
*/
void ElementBase::MoveCrackTipNodes(int crackTip)
{	int i;
	vector<int> movedNodes;
    for(i=0;i<nelems;i++)
		theElements[i]->MakeQuarterPointNodes(crackTip,movedNodes);
}



