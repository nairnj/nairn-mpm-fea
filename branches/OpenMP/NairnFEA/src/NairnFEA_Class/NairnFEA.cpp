/*********************************************************************
    NairnFEA.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mar 12 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
*********************************************************************/

#include "NairnFEA_Class/NairnFEA.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/NodalLoad.hpp"
#include "Boundary_Conditions/NodalDispBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"
#include "Boundary_Conditions/Constraint.hpp"
#include "System/FEAArchiveData.hpp"
#include <time.h>

// global analysis object
NairnFEA *fmobj=NULL;

// prototypes to solve banded, symmetric matrix problem
int gelbnd(double **,int,int,double *,double *,int);

/********************************************************************************
	NairnFEA: Constructors and Destructor
********************************************************************************/

NairnFEA::NairnFEA()
{
	version=3;					// NairnFEA version number
	subversion=0;				// subversion number
	buildnumber=0;				// build number
	xax='x';					// default axis names
	yax='y';
	zax='z';
	temperatureExpr=NULL;		// temperature expression
	stressFreeTemperature=0.;	// streess free temperature
	periodic.dof=0;				// periodic analysis
    
	// Default output flags
	int i;
	for(i=0;i<NUMBER_OUT;i++) outFlags[i]='Y';
	outFlags[REACT_OUT]='N';
    outFlags[NUMBER_OUT]=0;
	
	archiver=new FEAArchiveData();		// archiving object
}

/********************************************************************************
	NairnMPM: Start Analysis
********************************************************************************/

void NairnFEA::StartAnalysis(bool abort)
{
	// start analysis file
	StartResultsOutput();
	if(abort)
	{	cout << "\n***** NairnFEA RUN COMPLETED\n" << endl;
		return;
	}
    
	// Do FEA analysis
	FEAAnalysis();
}

/*********************************************************************
    NairnFEA: Main entry to read file and decode into objects
	On exception, throw CommonException() or char * exception
*********************************************************************/

void NairnFEA::FEAAnalysis()
{
    char nline[200];
    int result;
    int i;
    double *work;
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
    sprintf(nline,"Initial number of equations:%6d     Initial bandwidth:%6d",nsize,nband);
    cout << nline << endl << endl;

#pragma mark --- TASK 1: ALLOCATE R VECTOR
    // Allocate reaction vector and load with nodal loads and edge loads
    rm=new double[nsize+1];
    if(rm==NULL) throw CommonException("Out of memory creating reaction vector","NairnFEA::FEAAnalysis");
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
    work=new double[nsize+1];
    result=gelbnd(st,nsize,nband,rm,work,0);
    if(result==1)
    {	throw CommonException("Linear solver error: matrix is singular. Check boundary conditions.\n  (Hint: turn on resequencing to check for mesh connectivity problem)",
                                "NairnFEA::FEAAnalysis");
    }
    else if(result==-1)
	{	cout << "Linear solver warning: solution process was close to singular. Results might be invalid." << endl;
    }
    delete [] work;
    delete [] st;
	delete [] stiffnessMemory;
    
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
		Kij=st[j-i+1][i] for j-i+1 <= nband otherwise Kij=0
	if j<i
		Kij=st[i-j+1][j] for i-j+1 <= nband otherwise Kij=0
	
	Elements of upper half matrix stored as follows:
							Col i
							---------
	st[1][i-2]	st[2][i-2]	st[3][i-2]	st[4][i-2]	st[5][i-2]	st[6][i-2]	...
				st[1][i-1]	st[2][i-1]	st[3][i-1]	st[4][i-1]	st[5][i-1]	...
	Row i -------------->	st[1][i]	st[2][i]	st[3][i]	st[4][i]    ...
										st[1][i+1]	st[2][i+1]	st[3][i+1]	...
										
	For row i, second index=i and vary first from 1 to nband
	For col, first index 1 to nband and second from i to i-nband+1
				(note second index may go negative near edge of matrix)
***********************************************************************************/

void NairnFEA::BuildStiffnessMatrix(void)
{
    int i,j,iel,mi0,ni0,mi,ni,mj0,nj0,ind,mj,nj;
    int numnds,ii,jj;
    
    // allocate memory for stiffness matrix and zero it
	
	// st[] are points to diagonal columns of the stiffness matrix
    st=new double *[nband+1];
    if(st==NULL) throw CommonException("Out of memory creating stiffness matrix",
											"NairnFEA::BuildStiffnessMatrix");
											
	// allocate all in one contiguous block for better speed in algorithms
	stiffnessMemory=new double[nsize*nband];
    if(stiffnessMemory==NULL) throw CommonException("Out of memory creating stiffness matrix",
											"NairnFEA::BuildStiffnessMatrix");
											
	// allocate each column
	int baseAddr=0;
    for(i=1;i<=nband;i++)
	{	st[i]=&stiffnessMemory[baseAddr];
		st[i]--;					// to make it 1 based
		baseAddr+=nsize;
        for(j=1;j<=nsize;j++) st[i][j]=0.;
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
                                st[ind][mi]+=se[ni][nj];
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
				st[mj-mi+1][mi]+=nextConstraint->GetCoeff(i);
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

/***********************************************************
	Print Displacement Results
***********************************************************/

void NairnFEA::DisplacementResults(void)
{
    int i,ind=0,maxi,ii;
    char fline[200];
    
    if(outFlags[DISPLACEMENT_OUT]=='N') return;
    
    // heading
    PrintSection("NODAL DISPLACEMENTS (in mm)");
    if(np==AXI_SYM)
	cout << " Node        u               w               v" << endl;
    else
        cout << " Node        u               v               w" << endl;
    cout << "------------------------------------------------------" << endl;
    
    // each nodal displacement
    maxi = outFlags[DISPLACEMENT_OUT]=='Y' ? nnodes : selectedNodes.size();
    for(ii=1;ii<=maxi;ii++)
    {	// find node and dof
    	if(outFlags[DISPLACEMENT_OUT]=='Y')
            i=ii;
        else
            i=selectedNodes[ii-1];
    	ind=nfree*i;
    
    	// 2D output
        if(nfree==2)
            sprintf(fline,"%5d %15.7e %15.7e",i,1000.*rm[ind-1],1000.*rm[ind]);
        
        // 3D output
        else if(nfree==3)
        {   sprintf(fline,"%5d %15.7e %15.7e %15.7e",i,1000.*rm[ind-2],
                                            1000.*rm[ind-1],1000.*rm[ind]);
        }
        
        cout << fline << endl;
    }
    cout << endl;
}

/**********************************************************
	Calculate stresses and energies
	Print stress and forces
**********************************************************/

void NairnFEA::ForceStressEnergyResults(void)
{
    int i,j,iel,ind,kftemp=0,kstemp=0,numnds;
    int nodeNum;
    char gline[16],fline[200];
	
    if(outFlags[FORCE_OUT]!='N' || outFlags[ELEMSTRESS_OUT]!='N')
        PrintSection("NODAL FORCES (in N) AND STRESSES (in MPa) IN EACH ELEMENT");
   
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
        sprintf(gline,"%5d",iel+1);
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
            {	sprintf(fline,"%5s   %5d     %15.7e     %15.7e",gline,
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
            {   sprintf(fline,"%5s   %5d     %15.7e     %15.7e     %15.7e",
                            gline,nodeNum,1.e-6*se[j][1],1.e-6*se[j][2],1.e-6*se[j][3]);
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
                {	sprintf(fline,"%5s   %5d     %15.7e     %15.7e     %15.7e",
                                    gline,nodeNum,1.e-6*se[j][4],(double)0.0,(double)0.0);
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

/**********************************************************
    Print average modal stresses
        (Calculated in ForceStressEnergyResults())
**********************************************************/

void NairnFEA::AvgNodalStresses(void)
{
    int numshw,ii,i;
    
    if(outFlags[AVGSTRESS_OUT]=='N') return;
    
    // heading
    PrintSection("AVERAGE NODAL STRESSES (in MPa)");
    cout << " Node       sig(" << xax << ")           sig(" << yax << ")           sig(" 
            << zax << ")          sig(" << xax << yax << ")" << endl;
    cout << "--------------------------------------------------------------------------" << endl;

    // print all or selected nodes
    numshw = (outFlags[AVGSTRESS_OUT]=='Y') ? nnodes : selectedNodes.size();
    for(ii=1; ii<=numshw;ii++)
    {	i = (outFlags[AVGSTRESS_OUT]=='Y') ? ii : selectedNodes[ii-1];
        nd[i]->PrintAvgStress();
    }
    cout << endl;
}

/**********************************************************
	Print reactivities at fixed nodes
**********************************************************/

void NairnFEA::ReactionResults(void)
{
    NodalDispBC *nextBC=firstDispBC;
    
    if(outFlags[REACT_OUT]=='N') return;
    
    PrintSection("REACTIVITIES AT FIXED NODES (in N)");
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

/**********************************************************
	Print strain energy results
**********************************************************/

void NairnFEA::EnergyResults(void)
{
    double temp;
    int incolm,i,ind;
    char fline[200];
    
    if(outFlags[ENERGY_OUT]=='N') return;
    
    // heading
    PrintSection("STRAIN ENERGIES IN ELEMENTS (in J)");
    cout << " Elem      Strain Energy                 Elem      Strain Energy" << endl;
    cout << "------------------------------------------------------------------" << endl;
		
    temp=0.;
    if(IsEven(nelems))
        incolm=nelems/2;
    else
        incolm=(nelems+1)/2;
    for(i=1;i<=incolm;i++)
    {	ind=i+incolm;
        if(ind<=nelems)
        {   sprintf(fline,"%5d     %15.7e               %5d     %15.7e",
                    i,theElements[i-1]->strainEnergy,ind,theElements[ind-1]->strainEnergy);
            temp+=theElements[ind-1]->strainEnergy;
        }
        else
            sprintf(fline,"%5d     %15.7e",i,theElements[i-1]->strainEnergy);
        temp+=theElements[i-1]->strainEnergy;
        cout << fline << endl;
    }
    cout << "------------------------------------------------------------------" << endl;
    sprintf(fline,"Total     %15.7e",temp);
    cout << fline << endl << endl;
}

/********************************************************************************
	NairnFEA: accessors
********************************************************************************/

// return name, caller should delete
const char *NairnFEA::CodeName(void) const
{
	return "NairnFEA";
}

// Explain usage of this program
void NairnFEA::Usage()
{
    CoutCodeVersion();
    cout << "    Expects Xerces 3.1.1 or newer" << endl;
    cout << "\nUsage:\n"
            "    NairnFEA [-options] <InputFile>\n\n"
            "This program reads the <InputFile> and does an FEA analysis.\n"
            "The results summary is directed to standard output.\n\n"
            "Options:\n"
            "    -a          Abort after setting up problem but before\n"
			"                   FEA Analysis.\n"
            "    -H          Show this help.\n"
            "    -v          Validate input file if DTD is provided in !DOCTYPE\n"
            "                   (default is to skip validation)\n\n"
            "See http://oregonstate.edu/nairnj for full documentation.\n\n"
          <<  endl;
}

// verify analysis type
bool NairnFEA::ValidAnalysisType(void)
{	return np>=PLANE_STRAIN && np<END_FEA_TYPES;
}
