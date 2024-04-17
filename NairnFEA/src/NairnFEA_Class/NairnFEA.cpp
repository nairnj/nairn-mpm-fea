/*********************************************************************
    NairnFEA.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mar 12 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
*********************************************************************/

#include "stdafx.h"
#include "System/FEAArchiveData.hpp"
#include "System/UnitsController.hpp"
#include "NairnFEA_Class/NairnFEA.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/NodalLoad.hpp"
#include "Boundary_Conditions/NodalDispBC.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"
#include "Boundary_Conditions/Constraint.hpp"
#include "Exceptions/CommonException.hpp"
#include <time.h>

// global analysis object
NairnFEA *fmobj=NULL;

// prototypes to solve banded, symmetric matrix problem
int gelbnd(double **,int,int,double *,double *,int);

#pragma mark NairnFEA: Constructors and Destructor

// throws std::bad_alloc
NairnFEA::NairnFEA()
{
	version=9;					// NairnFEA version number
	subversion=0;				// subversion number
	buildnumber=0;				// build number
	xax='x';					// default axis names
	yax='y';
	zax='z';
	temperatureExpr=NULL;		// temperature expression
	stressFreeTemperature=0.;	// stress free temperature
	periodic.dof=0;				// periodic analysis
    
	// Default output flags
	int i;
	for(i=0;i<NUMBER_OUT;i++) outFlags[i]='Y';
	outFlags[REACT_OUT]='N';
    outFlags[NUMBER_OUT]=0;
	
	archiver=new FEAArchiveData();		// archiving object
}

#pragma mark NairnFEA: Run FEA Analysis

// Do the FEA analysis
void NairnFEA::CMAnalysis(bool abort)
{
	if(abort)
	{	cout << "\n***** NairnFEA RUN COMPLETED\n" << endl;
		return;
	}

	char nline[200];
    size_t nlsize=200;
    int result;
    int i;
    NodalDispBC *nextBC;
    double times[5];

#pragma mark --- TASK 0: INITIALIZE
    // start timer
    times[0]=CPUTime();
    
    // get problem size and other initialization
    nsize=nnodes*nfree + numConstraints;
    nband=GetBandWidth();
	
    if(np==AXI_SYM)
    {	xax='r';
        yax='z';
        zax='t';
    }
    
    // Stiffness matrix info
    PrintSection("TOTAL STIFFNESS MATRIX");
    snprintf(nline,nlsize,"Initial number of equations:%6d     Initial bandwidth:%6d",nsize,nband);
    cout << nline << endl << endl;

#pragma mark --- TASK 1: ALLOCATE R VECTOR
    // Allocate reaction vector and load with nodal loads and edge loads
    rm=(double *)malloc(sizeof(double)*(nsize+1));
    if(rm==NULL) throw CommonException("Memory error allocating reaction vector (rm)","NairnFEA::FEAAnalysis");
    for(i=1;i<=nsize;i++) rm[i]=0.;
    
    // add nodal loads to rm[] vector
    NodalLoad *nextLoad=firstLoadBC;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->Reaction(rm,np,nfree);
	
	// add stresses on element edges to rm[] vector
	ForcesOnEdges();

#pragma mark --- TASK 2: GET STIFFNESS MATRIX
    // allocate and fill global stiffness matrix, st[][], and reaction vector, rm[]
    times[1]=CPUTime();
    BuildStiffnessMatrix();

#pragma mark --- TASK 3: DISPLACEMENT BCs
    // Impose displacement boundary conditions and rotate
	//     nodes for skew boundary conditions
    nextBC=firstDispBC;
    while(nextBC!=NULL)
    	nextBC=nextBC->FixOrRotate(st,rm,nsize,nband,nfree);
    
#pragma mark --- TASK 4: INVERT STIFFNESS MATRIX
    // Solve linear system for nodal displacements
    times[2]=CPUTime();
    double *work=(double *)malloc(sizeof(double)*(nsize+1));
	if(work==NULL) throw CommonException("Memory error allocating work vector for linear solver",
                                        "NairnFEA::FEAAnalysis");
    result=gelbnd(st,nsize,nband,rm,work,0);
    if(result==1)
    {	throw CommonException("Linear solver error: matrix is singular. Check boundary conditions.\n  (Hint: turn on resequencing to check for mesh connectivity problem)",
                                "NairnFEA::FEAAnalysis");
    }
    else if(result==-1)
	{	cout << "Linear solver warning: solution process was close to singular. Results might be invalid." << endl;
    }
    free(work);
    free(st);
	free(stiffnessMemory);
    
#pragma mark --- TASK 5a: UNSKEW ROTATED NODES

    nextBC=firstDispBC;
    while(nextBC!=NULL)
    	nextBC=nextBC->Unrotate(rm,nfree);
    
#pragma mark --- TASK 6: OUTPUT RESULTS

	// time to here for performance evaluation
    double execTime=ElapsedTime();						// elpased time in secs
    times[3]=CPUTime();
    
    // Print Displacements
    DisplacementResults();
    
    // Calculate forces, stresses, and energy
    //	print element forces and stresses
    ForceStressEnergyResults();
    
    // Average nodal stresses
    AvgNodalStresses();
    
    // reactivities at fixed nodes
    ReactionResults();
    
    // strain energies
    EnergyResults();
    
    // execution times
    times[4]=CPUTime();
    PrintSection("EXECUTION TIMES AND MEMORY");
    cout << "1. Allocate Memory: " << (times[1]-times[0]) << " secs" << endl;		
    cout << "2. Build Stiffness Matrix: " << (times[2]-times[1]) << " secs" << endl;
    cout << "3. Solve Linear System: " << (times[3]-times[2]) << " secs" << endl;
    cout << "4. Write Results: " << (times[4]-times[3]) << " secs" << endl;
    cout << "5. Total Execution CPU Time: " << times[3] << " secs" << endl;
    cout << "6. Total Execution Elapsed Time: " << execTime << " secs" << endl;
	cout << "7. Scaling: " << times[3]/execTime << endl;
    
    //---------------------------------------------------
    // Trailer
    cout << "\n***** NAIRNFEA RUN COMPLETED\n";
}

/***********************************************************
    Calculate forces on edges
***********************************************************/

void NairnFEA::ForcesOnEdges(void)
{
	int i,mi0,ni0,ii,mi,ni,iel;
	int numnds;
	EdgeBC *nextEdge=firstEdgeBC;
	double re[2*MaxElNd+1];		// max free * max nodes in element + 1
	
	// loop of edge boundary conditions
	while(nextEdge!=NULL)
	{	// Have element calculate consistent loads in re[] vector
		nextEdge->GetConsistentLoads(re,np);
		
		/* Transfer consistent load vector re[] into global load vector rm[]
			m,n are adresses in global and element vectors */
		iel=nextEdge->ElementIndex();
		numnds=theElements[iel]->NumberNodes();
		for(i=1;i<=numnds;i++)
		{	mi0=nfree*(theElements[iel]->NodeIndex(i));
			ni0=nfree*(i-1);
			for(ii=1;ii<=nfree;ii++)
			{	mi=mi0+ii;
				ni=ni0+ii;
				rm[mi]+=re[ni];
			}
		}
		
		// next BC
		nextEdge=(EdgeBC *)nextEdge->GetNextObject();
	}
}

/***********************************************************************************
    Calculate stiffness matrix
	
	The stiffness matrix is symmetric, band diagonal with bandwidth nband
	If Kij is the (i,j) element of the stiffness matrix. It is stored in compact form
		that has only upper half of the matrix and only possible nonzero element
	    
	if j>=i
		Kij=st[i][j-i+1] for j-i+1 <= nband otherwise Kij=0
	if j<i
		Kij=st[j][i-j+1] for i-j+1 <= nband otherwise Kij=0
	
	Elements of upper half matrix stored as follows:
							Col i
							---------
	st[i-2][1]	st[i-2][2]	st[i-2][3]	st[i-2][4]	st[i-2][5]	st[i-2][6]	...
				st[i-1][1]	st[i-1][2]	st[i-1][3]	st[i-1][4]	st[i-1][5]	...
	Row i -------------->	st[i][1]	st[i][2]	st[i][3]	st[i][4]    ...
										st[i+1][1]	st[i+1][2]	st[i+1][3]	...
										
	For row i, first index=i and vary second from 1 to nband
                (note near bottom of matrix, non-zero elements stop at edge,
                    but memory beyond edge is there and zeroed)
	For col i, second index 1 to nband and first from i down to i-nband+1
				(note second index may go negative near top of matrix)
 
	throws CommonException()
***********************************************************************************/

void NairnFEA::BuildStiffnessMatrix(void)
{
    int i,j,iel,mi0,ni0,mi,ni,mj0,nj0,ind,mj,nj;
    int numnds,ii,jj;
    
    // allocate memory for stiffness matrix and zero it
	
	// st[] are pointers to rows of the stiffness matrix
    st = (double **)malloc(sizeof(double *)*(nsize+1));
    if(st==NULL) throw CommonException("Memory error creating stiffness matrix pointers (st)",
											"NairnFEA::BuildStiffnessMatrix");
											
	// allocate all in one contiguous block for better speed in algorithms
    stiffnessMemory = (double *)malloc(sizeof(double)*(nsize*nband));
    if(stiffnessMemory==NULL) throw CommonException("Memory error creating stiffness matrix memory block",
											"NairnFEA::BuildStiffnessMatrix");
											
	// allocate each row
	int baseAddr=0;
    for(i=1;i<=nsize;i++)
	{	st[i]=&stiffnessMemory[baseAddr];
		st[i]--;					// to make it 1 based
		baseAddr+=nband;
        for(j=1;j<=nband;j++) st[i][j]=0.;
    }
    
    // Loop over all elements
    for(iel=0;iel<nelems;iel++)
    {	// Get element stiffness matrix
        theElements[iel]->Stiffness(np);
        
        /*Transfer element stiffness matrix to global stiffness matrix
                m,n address in global and element matrices, respectively
                i,j refer to row and column
                0 refers to base address */
        numnds=theElements[iel]->NumberNodes();
        for(i=1;i<=numnds;i++)
        {   mi0=nfree*(theElements[iel]->NodeIndex(i));
            ni0=nfree*(i-1);
            for(j=1;j<=numnds;j++)
            {	if(theElements[iel]->NodeIndex(j)>=theElements[iel]->NodeIndex(i))
                {   mj0=nfree*(theElements[iel]->NodeIndex(j));
                    nj0=nfree*(j-1);
                    for(ii=1;ii<=nfree;ii++)
                    {	mi=mi0+ii;
                        ni=ni0+ii;
                        for(jj=1;jj<=nfree;jj++)
                        {   mj=mj0+jj;
                            ind=mj-mi+1;
                            if(ind>0)
                            {	nj=nj0+jj;
                                st[mi][ind]+=se[ni][nj];
                            }
                        }
                    }
                }
            }
        }
		
        /* Transfer element load vector into global load vector
                m,n are adresses in global and element vectors */
        for(i=1;i<=numnds;i++)
        {   mi0=nfree*(theElements[iel]->NodeIndex(i));
            ni0=nfree*(i-1);
            for(ii=1;ii<=nfree;ii++)
            {	mi=mi0+ii;
                ni=ni0+ii;
                rm[mi]+=re[ni];
            }
        }
    }
	
	// terms for contraints with Lagrange multiplier DOFs
	if(firstConstraint!=NULL)
	{	int nbase=nsize-numConstraints;
		Constraint *nextConstraint=firstConstraint;
		while(nextConstraint!=NULL)
        {	numnds=nextConstraint->NumberNodes();
			mj=nbase+nextConstraint->GetLambdaNum();		// mj > mi always
			for(i=1;i<=numnds;i++)
			{   mi=nextConstraint->NodalDof(i,nfree);
				st[mi][mj-mi+1]+=nextConstraint->GetCoeff(i);
			}
			rm[mj]+=nextConstraint->GetQ();
			nextConstraint=(Constraint *)nextConstraint->GetNextObject();
		}
	}
}

/***********************************************************
	Calculate bandwidth
***********************************************************/

int NairnFEA::GetBandWidth(void)
{
    int nband=1,i;
	int minn,maxn;
    
    for(i=0;i<nelems;i++)
	{	maxn=0;
		theElements[i]->MaxMinNode(&maxn,&minn);
		nband=fmax(nfree*(maxn-minn+1),nband);
	}
	
	// get bandwidth for Lagrange multiplier DOFs
	if(firstConstraint!=NULL)
	{	int nbase=nsize-numConstraints;
		Constraint *nextConstraint=firstConstraint;
		while(nextConstraint!=NULL)
		{	nband=fmax(nband,nextConstraint->GetBandWidth(nbase,nfree));
			nextConstraint=(Constraint *)nextConstraint->GetNextObject();
		}
	}
	
    return nband;
}

#pragma mark NairnFEA: Output Results

// Print Displacement Results
void NairnFEA::DisplacementResults(void)
{
    int i,ind=0,maxi,ii;
    char fline[200];
    size_t fsize=200;
    
    if(outFlags[DISPLACEMENT_OUT]=='N') return;
    
    // heading
	snprintf(fline,fsize,"NODAL DISPLACEMENTS (in %s)",UnitsController::Label(CULENGTH_UNITS));
    PrintSection(fline);
    if(np==AXI_SYM)
	cout << " Node        u               w               v" << endl;
    else
        cout << " Node        u               v               w" << endl;
    cout << "------------------------------------------------------" << endl;
    
    // each nodal displacement
    maxi = outFlags[DISPLACEMENT_OUT]=='Y' ? nnodes : (int)selectedNodes.size();
    for(ii=1;ii<=maxi;ii++)
    {	// find node and dof
    	if(outFlags[DISPLACEMENT_OUT]=='Y')
            i=ii;
        else
            i=selectedNodes[ii-1];
    	ind=nfree*i;
    
    	// 2D output
        if(nfree==2)
            snprintf(fline,fsize,"%5d %15.7e %15.7e",i,rm[ind-1],rm[ind]);
        
        // 3D output
        else if(nfree==3)
            snprintf(fline,fsize,"%5d %15.7e %15.7e %15.7e",i,rm[ind-2],rm[ind-1],rm[ind]);
        
        cout << fline << endl;
    }
    cout << endl;
}

// Calculate stresses and energies
// Print stress and forces
void NairnFEA::ForceStressEnergyResults(void)
{
    int i,j,iel,ind,kftemp=0,kstemp=0,numnds;
    int nodeNum;
    char gline[16],fline[200];
    size_t gsize=16,fsize=200;
	
    if(outFlags[FORCE_OUT]!='N' || outFlags[ELEMSTRESS_OUT]!='N')
	{	snprintf(fline,fsize,"NODAL FORCES (in %s) AND STRESSES (in %s) IN EACH ELEMENT",
						UnitsController::Label(FEAFORCE_UNITS),UnitsController::Label(PRESSURE_UNITS));
        PrintSection(fline);
	}
   
    /* The nodal stresses will store nodal point objects
            fs.stress.sig[], fs.force, fs.numElems
        Strain energy in element object strainEnergy
    */
    for(i=1;i<=nnodes;i++)
        nd[i]->InitForceField();

    // Loop over elements, calculating forces and averaging stresses
    for(iel=0;iel<nelems;iel++)
    {	theElements[iel]->ForceStress(rm,np,nfree);
        numnds=theElements[iel]->NumberNodes();
    
        // Print forces at nodes
        snprintf(gline,gsize,"%5d",iel+1);
        if(theElements[iel]->WantElement(outFlags[FORCE_OUT],selectedNodes))
        {   kftemp=1;
            cout << "--------------------------------------------------------------------------" << endl;
            cout << "Element  Node           F" << xax << "                  F" << yax
                    << "                  F" << zax << endl;
            cout << "--------------------------------------------------------------------------" << endl;
        }
        else
            kftemp=0;

        // print and sum forces
        for(j=1;j<=numnds;j++)
        {   ind=nfree*(j-1)+1;
            nodeNum=theElements[iel]->nodes[j-1];

            // print force
            if(kftemp==1)
            {	snprintf(fline,fsize,"%5s   %5d     %15.7e     %15.7e",gline,
                                nodeNum,-se[ind][7],-se[ind+1][7]);
                cout << fline << endl;
                strcpy(gline,"     ");
            }

            // Sum forces at nodes
            nd[nodeNum]->fs->force.x+=se[ind][7];
            nd[nodeNum]->fs->force.y+=se[ind+1][7];
        }

        // Print stresses at nodes in elements
        if(theElements[iel]->WantElement(outFlags[ELEMSTRESS_OUT],selectedNodes))
        {   kstemp=1;
            if(kftemp==1)
            {	cout << "                  --------------------------------------------------------" << endl;
                cout << "                       sig" << xax << "                sig" << yax
                        << "               sig" << xax << yax << endl;
                cout << "                  --------------------------------------------------------" << endl;
            }
            else
            {	cout << "--------------------------------------------------------------------------" << endl;
                cout << "Element  Node          sig" << xax << "                sig" << yax
                        << "               sig" << xax << yax << endl;
                cout << "--------------------------------------------------------------------------" << endl;
            }
        }
        else
            kstemp=0;

        for(j=1;j<=numnds;j++)
        {   nodeNum=theElements[iel]->nodes[j-1];
        
            if(kstemp==1)
            {   snprintf(fline,fsize,"%5s   %5d     %15.7e     %15.7e     %15.7e",
                            gline,nodeNum,se[j][1],se[j][2],se[j][3]);
                cout << fline << endl;
                strcpy(gline,"     ");
            }

            // Transfer stresses to global array - if bulk element
			if(theElements[iel]->BulkElement())
			{	nd[nodeNum]->fs->stress.xx+=se[j][1];
				nd[nodeNum]->fs->stress.yy+=se[j][2];
				nd[nodeNum]->fs->stress.xy+=se[j][3];
				nd[nodeNum]->fs->numElems++;
			}
        }

        // Print rest of stresses and add stresses into global array
        if(np==PLANE_STRAIN ||  np==AXI_SYM)
        {	if(kstemp==1)
            {   cout << "                  --------------------------------------------------------" << endl;
                cout << "                       sig" << zax << "               sig" << xax << zax
                        << "               sig" << yax << zax << endl;
                cout << "                  --------------------------------------------------------" << endl;
            }
            
            for(j=1;j<=numnds;j++)
            {   nodeNum=theElements[iel]->nodes[j-1];
                if(kstemp==1)
                {	snprintf(fline,fsize,"%5s   %5d     %15.7e     %15.7e     %15.7e",
                                    gline,nodeNum,se[j][4],(double)0.0,(double)0.0);
                    cout << fline << endl;
                }

                // Tranfer stresses to global array
				if(theElements[iel]->BulkElement())
					nd[nodeNum]->fs->stress.zz+=se[j][4];
            }
        }
    }
    
    // blank line
    if(outFlags[FORCE_OUT]!='N' || outFlags[ELEMSTRESS_OUT]!='N')
		cout << endl;
}

// Print average modal stresses
//  (Calculated in ForceStressEnergyResults())
void NairnFEA::AvgNodalStresses(void)
{
    int numshw,ii,i;
    
    if(outFlags[AVGSTRESS_OUT]=='N') return;
    
    // heading
	char fline[200];
    size_t fsize=200;
	snprintf(fline,fsize,"AVERAGE NODAL STRESSES (in %s)",UnitsController::Label(PRESSURE_UNITS));
	PrintSection(fline);
    cout << " Node       sig(" << xax << ")           sig(" << yax << ")           sig(" 
            << zax << ")          sig(" << xax << yax << ")" << endl;
    cout << "--------------------------------------------------------------------------" << endl;

    // print all or selected nodes
    numshw = (outFlags[AVGSTRESS_OUT]=='Y') ? nnodes : (int)selectedNodes.size();
    for(ii=1; ii<=numshw;ii++)
    {	i = (outFlags[AVGSTRESS_OUT]=='Y') ? ii : selectedNodes[ii-1];
        nd[i]->PrintAvgStress();
    }
    cout << endl;
}

// Print reactivities at fixed nodes
void NairnFEA::ReactionResults(void)
{
    NodalDispBC *nextBC=firstDispBC;
    
    if(outFlags[REACT_OUT]=='N') return;
    
	char fline[200];
    size_t fsize=200;
	snprintf(fline,fsize,"REACTIVITIES AT FIXED NODES (in %s)",UnitsController::Label(FEAFORCE_UNITS));
	PrintSection(fline);
    cout << " Node           F" << xax << "                  F" << yax
            << "                  F" << zax << endl;
    cout << "------------------------------------------------------------------" << endl;
    
    // each one
    while(nextBC!=NULL)
    {	nextBC->PrintReaction();
        nextBC=(NodalDispBC *)nextBC->GetNextObject();
    }
    cout << endl;
}

// Print strain energy results
void NairnFEA::EnergyResults(void)
{
    double temp;
    int incolm,i,ind;
    char fline[200];
    size_t fsize=200;
    
    if(outFlags[ENERGY_OUT]=='N') return;
    
    // heading
	snprintf(fline,fsize,"STRAIN ENERGIES IN ELEMENTS (in %s)",UnitsController::Label(FEAWORK_UNITS));
	PrintSection(fline);
    cout << " Elem      Strain Energy                 Elem      Strain Energy" << endl;
    cout << "------------------------------------------------------------------" << endl;
		
    temp=0.;
    if(IsEven(nelems))
        incolm=nelems/2;
    else
        incolm=(nelems+1)/2;
	double escale = UnitsController::Scaling(1.e-3);
    for(i=1;i<=incolm;i++)
    {	ind=i+incolm;
        if(ind<=nelems)
        {   snprintf(fline,fsize,"%5d     %15.7e               %5d     %15.7e",
                    i,escale*theElements[i-1]->strainEnergy,ind,escale*theElements[ind-1]->strainEnergy);
            temp+=theElements[ind-1]->strainEnergy;
        }
        else
            snprintf(fline,fsize,"%5d     %15.7e",i,escale*theElements[i-1]->strainEnergy);
        temp+=theElements[i-1]->strainEnergy;
        cout << fline << endl;
    }
    cout << "------------------------------------------------------------------" << endl;
    snprintf(fline,fsize,"Total     %15.7e",escale*temp);
    cout << fline << endl << endl;
}

#pragma mark NairnFEA: accessors

// verify analysis type
bool NairnFEA::ValidAnalysisType(void)
{	return np>=PLANE_STRAIN && np<END_FEA_TYPES;
}

#pragma NairnFEA: archiver pass alongs

// archiver pass throught while setting up analysis
void NairnFEA::ArchiveNodalPoints(int np) { archiver->ArchiveNodalPoints(np); }
void NairnFEA::ArchiveElements(int np) { archiver->ArchiveElements(np); }
void NairnFEA::SetInputDirPath(const char *xmlFIle,bool useWorkingDir)
{	archiver->SetInputDirPath(xmlFIle,useWorkingDir);
}
