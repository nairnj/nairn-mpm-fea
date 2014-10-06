/********************************************************************************
    ArchiveData.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 06 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.    
********************************************************************************/

#include <fstream>
#include <errno.h>

#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/GlobalQuantity.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"

// archiver global
ArchiveData *archiver;

#pragma mark Constructors and Destructor and Initializers

ArchiveData::ArchiveData()
{
	archTime=1.;			// time interval between archives (sec)
	globalTime=-1.;			// time interval for archiving global results (sec)
	firstArchTime=-1.;		// first archive time (sec)
	nextArchTime=0.;		// next time to archive results (sec)
	nextGlobalTime=0.;		// next time to archive global results (sec)
	globalFile=NULL;		// path to global results file
	threeD=FALSE;			// three D calculations
    SetMPMOrder("mYYYYYNYYYNNNNNNNN");		// byte order + defaults + 16 items
    SetCrackOrder("mYNNN");					// byte order + defaults + 3 items
	timeStamp==NULL;
    
    // contact archive to coordinate with global contact archiving
	lastArchiveContactStep=0;               // last time contact force was archived
	doingArchiveContact=FALSE;
    ZeroVector(&lastContactForce);
}

// create archive folder
bool ArchiveData::MakeArchiveFolder(void)
{
	if(archiveRoot==NULL) return false;
	
    //---------------------------------------------------
    // Create directory for archived files (if needed)
	char syscmd[600];
	if(strlen(archiveParent)>0 || forceUnique)
	{	// find unique archiveParent if requested
		if(forceUnique)
		{	int folderID=1;
			while(folderID<1000)
			{	if(strlen(archiveParent)>0)
					sprintf(syscmd,"test -d '%s%s/%d'",outputDir,archiveParent,folderID);
				else
					sprintf(syscmd,"test -d '%s%d'",outputDir,folderID);
				int exists=system(syscmd);
				if(exists!=0) break;			// zero means it already exists
				folderID++;
			}
			
			// if not found, an error
			if(folderID>=1000) return false;
			
			// adjust archiveParent and archiveRoot
			int insertPos=strlen(archiveParent);
			if(insertPos>0)
			{	char fldrNum[10];
				sprintf(fldrNum,"/%d",folderID);			// max length is 4, and space was saved for it
				int endPos=strlen(archiveRoot);
				int numLength=strlen(fldrNum);
				int i;
				for(i=endPos+numLength;i>=insertPos;i--)
				{	if(i>=insertPos+numLength)
						archiveRoot[i]=archiveRoot[i-numLength];
					else
						archiveRoot[i]=fldrNum[i-insertPos];
				}
				strcat(archiveParent,fldrNum);
			}
			else
			{	sprintf(syscmd,"%d/%s",folderID,archiveRoot);			// space was saved for insertion
				strcpy(archiveRoot,syscmd);
				sprintf(syscmd,"%s/%d",archiveParent,folderID);
				strcpy(archiveParent,syscmd);
			}
		}
		
		// now make the folder
    	strcpy(syscmd,"mkdir -p '");
		strcat(syscmd,outputDir);
		strcat(syscmd,archiveParent);
		strcat(syscmd,"'");
		system(syscmd);
	}
	
	// copy input commands
	strcpy(syscmd,"cp '");
	strcat(syscmd,inputDir);							// input folder
	strcat(syscmd,&inputDir[strlen(inputDir)+1]);		// input name
	strcat(syscmd,"' '");
	strcat(syscmd,outputDir);
	strcat(syscmd,archiveRoot);
	strcat(syscmd,".fmcmd'");
    system(syscmd);	
	
	// test by creating dummy file because return value above may be system dependent
	// If logging progress, keep the temporary file for future use
    FILE *fp;
#ifdef LOG_PROGRESS
	logFile=new char[strlen(outputDir)+strlen(archiveRoot)+6];
	logStartTime=fmobj->CPUTime();
#else
	char *logFile=new char[strlen(outputDir)+strlen(archiveRoot)+6];
#endif
	sprintf(logFile,"%s%s.log",outputDir,archiveRoot);
    if((fp=fopen(logFile,"w"))==NULL) return false;
	fclose(fp);
#ifndef LOG_PROGRESS
	// delete file unless want to keep it for logging
	remove(logFile);
	delete [] logFile;
#endif
	
	return true;
}

// Create global file, get archive size, and print to output file archiving table
void ArchiveData::BeginArchives(bool isThreeD)
{
	// set flag if 3D calculations
	threeD=isThreeD;
	
	// global archiving
	CreateGlobalFile();
	
    // Archive file list headind
    CalcArchiveSize();
	SetArchiveHeader();
    PrintSection("ARCHIVED ANALYSIS RESULTS");
    cout << "Root file name: " << archiveRoot << "." << endl
        << "Archive format: " << mpmOrder << endl
        << "Crack archive format: " << crackOrder << endl << endl
        << "  Step    Time (msec)    Filename" << endl
        << "----------------------------------------------"
        << endl;
}
	
/**********************************************************
    get archive record size using current save orders
	
	Although these methods use sizeof and will therefore
	work on any computer, for the archive files to be
	same on all computers then must have following sizes:
		int: 4 bytes		short: 2 bytes
		double: 8 bytes
*/
void ArchiveData::CalcArchiveSize(void)
{
    int i;
    
    // pad if needed, and truncate if not recognized
    if(strlen(mpmOrder)<2) strcpy(mpmOrder,"mY");
    if(strlen(mpmOrder)<ARCH_MAXMPMITEMS)
    {	for(i=strlen(mpmOrder);i<ARCH_MAXMPMITEMS;i++)
            mpmOrder[i]='N';
    }
    mpmOrder[ARCH_MAXMPMITEMS]=0;
	mpmOrder[ARCH_Defaults]='Y';			// now required to be 'Y'
	mpmOrder[ARCH_OldOrigPosition]='N';		// now requried to be 'N'
	mpmOrder[ARCH_ver2Empty]='N';			// now required to be 'N'
	
    if(strlen(crackOrder)<2) strcpy(crackOrder,"mY");
    if(strlen(crackOrder)<ARCH_MAXCRACKITEMS)
    {	for(i=strlen(crackOrder);i<ARCH_MAXCRACKITEMS;i++)
            crackOrder[i]='N';
    }
    crackOrder[ARCH_MAXCRACKITEMS]=0;
	crackOrder[ARCH_Defaults]='Y';			// now required to be 'Y'
    
    /* byte order marker
		Intel chips use little endian in which int 1 will have 0x01 in first byte.
		G5, etc, use big endian in which int 1 will have 0x01 in last byte.
		When done, byte order will be set to 'm' if output file is big endian
			or to 'i' if output file is little endian.
	*/
    int test=1;
    char *testPtr=(char *)&test;			// point to first byte
    if(!fmobj->GetReverseBytes())			// if no -r flag, point to last byte
        testPtr+=sizeof(int)-1;
    if(*testPtr==1)
	{	// either -r flag and running on little endian or no -r and running on big endian
    	mpmOrder[ARCH_ByteOrder]='m';
        crackOrder[ARCH_ByteOrder]='m';
    }
    else
	{	// either -r flag and running on big endian or no -r and running on little endian
    	mpmOrder[ARCH_ByteOrder]='i';
        crackOrder[ARCH_ByteOrder]='i';
    }
    
    // initialize
    mpmRecSize=0;
    crackRecSize=0;
	
	// sizes
	int vectorSize,tensorSize;
	if(threeD)
	{	vectorSize=3*sizeof(double);
		tensorSize=6*sizeof(double);
	}
	else
	{	vectorSize=2*sizeof(double);
		tensorSize=4*sizeof(double);
	}
    
    // check what will be there for material points
	
	/* ARCH_Defaults are
		2D: elemID (int), mass (double), matId (short) angle (double), thickness (double),
							pos (Vector), origPos (Vector) (64)
		3D: thickness replaced by two angles and Vectors are longer (88)
	*/
	mpmRecSize+=sizeof(int)+3*sizeof(double)+2*vectorSize+sizeof(short)+2;
	if(threeD) mpmRecSize+=sizeof(double);
	
    if(mpmOrder[ARCH_Velocity]=='Y')
        mpmRecSize+=vectorSize;
    if(mpmOrder[ARCH_Stress]=='Y')
        mpmRecSize+=tensorSize;
    if(mpmOrder[ARCH_Strain]=='Y')
        mpmRecSize+=tensorSize;
    if(mpmOrder[ARCH_PlasticStrain]=='Y')
        mpmRecSize+=tensorSize;
    if(mpmOrder[ARCH_WorkEnergy]=='Y')
        mpmRecSize+=sizeof(double);
    if(mpmOrder[ARCH_DeltaTemp]=='Y')
        mpmRecSize+=sizeof(double);
    if(mpmOrder[ARCH_PlasticEnergy]=='Y')
        mpmRecSize+=sizeof(double);
    if(mpmOrder[ARCH_ShearComponents]=='Y')
        mpmRecSize+=2*sizeof(double);
    if(mpmOrder[ARCH_StrainEnergy]=='Y')
        mpmRecSize+=sizeof(double);
    if(mpmOrder[ARCH_History]=='Y')
        mpmRecSize+=sizeof(double);
	else if(mpmOrder[ARCH_History]!='N')
	{	if(mpmOrder[ARCH_History]&0x01) mpmRecSize+=sizeof(double);
		if(mpmOrder[ARCH_History]&0x02) mpmRecSize+=sizeof(double);
		if(mpmOrder[ARCH_History]&0x04) mpmRecSize+=sizeof(double);
		if(mpmOrder[ARCH_History]&0x08) mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_Concentration]=='Y')
        mpmRecSize+=sizeof(double)+vectorSize;
    if(mpmOrder[ARCH_HeatEnergy]=='Y')
        mpmRecSize+=sizeof(double);
    if(mpmOrder[ARCH_ElementCrossings]=='Y')
        mpmRecSize+=sizeof(int);
    if(mpmOrder[ARCH_RotStrain]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
            
    // check what will be there for crack segments
	
	/* ARCH_Defaults are
		planeInElem (int), empty (double), newCrack (short+2), pos (Vector), origPos (Vector),
			aboveInElem (int), above (Vector) belowInElem (int), below (Vector)
	*/
    crackRecSize+=sizeof(int)+sizeof(double)+sizeof(short)+2;
	crackRecSize+=2*sizeof(int)+8*sizeof(double);
    if(crackOrder[ARCH_JIntegral]=='Y')
        crackRecSize+=2*sizeof(double);
    if(crackOrder[ARCH_StressIntensity]=='Y')
        crackRecSize+=2*sizeof(double);
    if(crackOrder[ARCH_BalanceResults]=='Y')
        crackRecSize+=sizeof(int)+2*sizeof(double);
    
    // record is max of these two sizes
    recSize=mpmRecSize>crackRecSize ? mpmRecSize : crackRecSize;
}

/**********************************************************
    Create archive file header once at start of the
	calculations
		ver3: new defaults and oldest one support by current tools
				header was only 4 bytes with ID string
		ver4: add 60 more bytes to make 64 byte header
				0-3: version ID
				4: number of mpm archive items
				5 on: mpmOrder string  (>=18)
				next: number of crack archive items
				next+1  on: crackOrder (>=5)
				next: '2' or '3' for dimensionality
		ver5:	next: '0' or '1' for structured grid
				next 4: archive time in ms (float, using Endian of archive items setting)
				rest: 0
		ver6: changed default properties to have 3 angles in 3D results
	
	Current length 35 (if mpmOrder is 18 and crackOrder is 5)
*/
void ArchiveData::SetArchiveHeader(void)
{
	unsigned i;
	for(i=0;i<HEADER_LENGTH;i++) archHeader[i]=0;
	
	// version ID
	strcpy(archHeader,"ver6");
	
	// mpmOrder
	archHeader[strlen(archHeader)]=strlen(mpmOrder);
	strcat(archHeader,mpmOrder);
	
	// crackOrder
	archHeader[strlen(archHeader)]=strlen(crackOrder);
	strcat(archHeader,crackOrder);
	
	// 2 or 3 dimensions
	archHeader[strlen(archHeader)] = (threeD) ? '3' : '2' ;
	
	// structured
	archHeader[strlen(archHeader)] = mpmgrid.IsStructuredGrid() ? '1' : '0' ;
	
	// save offset to 4 byte time for time float
	timeStamp=(float *)&archHeader[strlen(archHeader)];
	for(i=0;i<sizeof(float);i++) archHeader[strlen(archHeader)]='-';
}

// Create global archive file
void ArchiveData::CreateGlobalFile(void)
{
    FILE *fp;
	char fline[1000];
	GlobalQuantity *nextGlobal;
	
	// skip if none, but use archTime if wants global archive
	if(globalTime<0.)
	{	if(firstGlobal==NULL) return;
		globalTime = archTime;
	}
	
	// get relative path name to the file
	globalFile=new char[strlen(outputDir)+strlen(archiveRoot)+8];
	sprintf(globalFile,"%s%s.global",outputDir,archiveRoot);
	
    // create and open the file
    if((fp=fopen(globalFile,"w"))==NULL) goto abort;
	
	// write color and count archives
	strcpy(fline,"#setColor");
	nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
	    nextGlobal=nextGlobal->AppendColor(fline);
	strcat(fline,"\n");
	if(fwrite(fline,strlen(fline),1,fp)!=1) goto abort;
	
	// write name
	strcpy(fline,"#setName");
	nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
		nextGlobal=nextGlobal->AppendName(fline);
	strcat(fline,"\n");
	if(fwrite(fline,strlen(fline),1,fp)!=1) goto abort;
	
	// close the file
    if(fclose(fp)!=0) goto abort;
	
	// section in output file
    PrintSection("ARCHIVED GLOBAL RESULTS");
    cout << "Global data file: " << archiveRoot << ".global" << endl;
	cout << endl;
	
	return;

abort:
	FileError("Global archive file creation failed",globalFile,"ArchiveData::CreateGlobalFile");
}

#pragma mark ARCHIVING METHODS

// Nodal velocity conditions
void ArchiveData::ArchiveVelocityBCs(BoundaryCondition *firstBC)
{
	BoundaryCondition *nextBC=firstBC;
	
	if(archiveMesh && fmobj->IsThreeD())
	{	char fname[500];
		ofstream outfile;
		sprintf(fname,"%s%s_VelBCs.txt",outputDir,archiveRoot);
		outfile.open(fname);
		if(outfile)
    	{	sprintf(fname,"%s_VelBCs.txt",archiveRoot);
			cout << "File: " << fname << endl << endl;
			while(nextBC!=NULL)
				nextBC=nextBC->PrintBC(outfile);
			outfile.close();
			return;
		}
	}
	
	// list in output file (by request or if file error)
    cout << " Node    DOF ID  Vel (mm/sec)   Arg (ms/ms^-1)  Angle1  Angle2  Function\n"
   	     << "--------------------------------------------------------------------------\n";
    nextBC=(BoundaryCondition *)firstBC;
    while(nextBC!=NULL)
		nextBC=nextBC->PrintBC(cout);
    cout << endl;
}

// Archive the results if it is time
void ArchiveData::ArchiveResults(double atime)
{
	double rho,rho0;
    double sxx,syy,sxy;
    char fname[500],fline[500];
    int i,p;
    CrackHeader *nextCrack;
	
	// test global archiving
	GlobalArchive(atime);
    
    // see if desired
    if(atime<nextArchTime || timeStamp==NULL) return;
    nextArchTime+=archTime;
	
	// exit if using delayed archiving
	if(atime>0.9*timestep && atime<firstArchTime) return;
    
    // get relative path name to the file
    sprintf(fname,"%s%s.%d",outputDir,archiveRoot,fmobj->mstep);
    
    // output step number, time, and file name to results file
    for(i=strlen(fname);i>=0;i--)
    {	if(fname[i]=='/') break;
    }
    sprintf(fline,"%7d %15.7e  %s",fmobj->mstep,1000.*atime,&fname[i+1]);
    cout << fline << endl;

    // open the file
	ofstream afile;
	try
	{	afile.open(fname, ios::out | ios::binary);
		if(!afile.is_open())
			FileError("Cannot open an archive file",fname,"ArchiveData::ArchiveResults");
		
		// write header created in SetArchiveHeader
		*timeStamp=1000.*atime;
		afile.write(archHeader,HEADER_LENGTH);
		if(afile.bad())
			FileError("File error writing archive file header",fname,"ArchiveData::ArchiveResults");
	}
	catch(CommonException err)
	{   // give up on hopefully temporary file problem
		cout << "# " << err.Message() << endl;
		cout << "# Will try to continue" << endl;
		if(afile.is_open()) afile.close();
		return;
	}
		
	// allocate space for one material point
	long blen=recSize;
	char *aptr=(char *)malloc(blen);
	if(aptr==NULL)
		throw CommonException("Out of memory allocating buffer for archive file","ArchiveData::ArchiveResults");
    
    // all material points
    for(p=0;p<nmpms;p++)
	{	// buffer is for one particle
    	char *app=aptr;
		
		// must have these defaults
	   
		// element ID
        *(int *)app=mpm[p]->ArchiveElemID();
        app+=sizeof(int);
        
		// mass (g)
        *(double *)app=mpm[p]->mp;
        app+=sizeof(double);
        
		// material ID
        *(short *)app=mpm[p]->ArchiveMatID();
        app+=sizeof(short);
		// fill in two zeros for byte alignment
		*app=0;
		app+=1;
		*app=0;
		app+=1;
        
		// 3D has three angles, 2D has one angle and thickness
		if(threeD)
		{	// 3 material rotation angles in degrees
			*(double *)app=mpm[p]->GetRotationZInDegrees();
			app+=sizeof(double);
			*(double *)app=mpm[p]->GetRotationYInDegrees();
			app+=sizeof(double);
			*(double *)app=mpm[p]->GetRotationXInDegrees();
			app+=sizeof(double);
		}
		else
		{	// material rotation angle in degrees
			*(double *)app=mpm[p]->GetRotationZInDegrees();
			app+=sizeof(double);
			
			// thickness (2D) in mm
			*(double *)app=mpm[p]->thickness();
			app+=sizeof(double);
		}
        
		// (x,y,z) position (mm)
        *(double *)app=mpm[p]->pos.x;
        app+=sizeof(double);
        
        *(double *)app=mpm[p]->pos.y;
        app+=sizeof(double);
		
		if(threeD)
		{	*(double *)app=mpm[p]->pos.z;
			app+=sizeof(double);
		}

		// original (x,y,z) position (mm)
        *(double *)app=mpm[p]->origpos.x;
        app+=sizeof(double);
                
        *(double *)app=mpm[p]->origpos.y;
        app+=sizeof(double);

		if(threeD)
		{	*(double *)app=mpm[p]->origpos.z;
			app+=sizeof(double);
		}

        // velocity (mm/sec)
        if(mpmOrder[ARCH_Velocity]=='Y')
        {   *(double *)app=mpm[p]->vel.x;
            app+=sizeof(double);
                
            *(double *)app=mpm[p]->vel.y;
            app+=sizeof(double);
			
			if(threeD)
			{	*(double *)app=mpm[p]->vel.z;
				app+=sizeof(double);
			}
        }

        // stress - internally it is N/m^2 cm^3/g, output in N/m^2 or Pa
        //          internal SI units are kPa/(kg/m^3)
		// For large deformation, need to convert Kirchoff Stress/rho0 to Cauchy stress
        int matid = mpm[p]->MatID();
        rho0=theMaterials[matid]->rho;
        rho = rho0/theMaterials[matid]->GetCurrentRelativeVolume(mpm[p]);
        Tensor sp = mpm[p]->ReadStressTensor();
        sxx=rho*sp.xx;
        syy=rho*sp.yy;
        sxy=rho*sp.xy;
        if(mpmOrder[ARCH_Stress]=='Y')
        {   *(double *)app=sxx;
            app+=sizeof(double);
                
            *(double *)app=syy;
            app+=sizeof(double);
            
			*(double *)app=rho*sp.zz;
            app+=sizeof(double);
                
            *(double *)app=sxy;
            app+=sizeof(double);
			
			if(threeD)
            {	*(double *)app=rho*sp.xz;
				app+=sizeof(double);
				
            	*(double *)app=rho*sp.yz;
				app+=sizeof(double);
			}
        }

        // elastic strain (absolute)
        if(mpmOrder[ARCH_Strain]=='Y')
		{	Tensor *ep=mpm[p]->GetStrainTensor();
			*(double *)app=ep->xx;
            app+=sizeof(double);
                
            *(double *)app=ep->yy;
            app+=sizeof(double);
            
			*(double *)app=ep->zz;
            app+=sizeof(double);
                
            *(double *)app=ep->xy;
            app+=sizeof(double);
			
			if(threeD)
            {	*(double *)app=ep->xz;
				app+=sizeof(double);
				
            	*(double *)app=ep->yz;
				app+=sizeof(double);
			}
        }
        
        // plastic strain (absolute)
        if(mpmOrder[ARCH_PlasticStrain]=='Y')
		{	Tensor *eplast=mpm[p]->GetPlasticStrainTensor();
            *(double *)app=eplast->xx;
            app+=sizeof(double);
                
            *(double *)app=eplast->yy;
            app+=sizeof(double);
                
            *(double *)app=eplast->zz;
            app+=sizeof(double);
                
            *(double *)app=eplast->xy;
            app+=sizeof(double);
			
			if(threeD)
            {	*(double *)app=eplast->xz;
				app+=sizeof(double);
				
            	*(double *)app=eplast->yz;
				app+=sizeof(double);
			}
        }
        
        // external work (cumulative) in J
        if(mpmOrder[ARCH_WorkEnergy]=='Y')
		{	*(double *)app=1.0e-6*mpm[p]->mp*mpm[p]->GetWorkEnergy();
            app+=sizeof(double);
        }
                
        // temperature
        if(mpmOrder[ARCH_DeltaTemp]=='Y')
        {   *(double *)app=mpm[p]->pTemperature;
            app+=sizeof(double);
        }
        
        // total plastic energy (Volume*energy) in J
        // energies in material point based on energy per unit mass
        // here need mass * U/(rho0 V0)
        if(mpmOrder[ARCH_PlasticEnergy]=='Y')
        {   *(double *)app=1.0e-6*mpm[p]->mp*mpm[p]->GetPlastEnergy();
            app+=sizeof(double);
        }
                
        // shear components (absolute)
        if(mpmOrder[ARCH_ShearComponents]=='Y')
        {   *(double *)app=mpm[p]->GetDuDy();
            app+=sizeof(double);
                
            *(double *)app=mpm[p]->GetDvDx();
            app+=sizeof(double);
        }

        // total strain energy (Volume*energy) in J
        // energies in material point based on energy per unit mass
        // here need mass * U/(rho0 V0)
        // internal units are same as stress: N/m^2 cm^3/g = microJ/g = mJ/kg
		// note that rho*energy has units J/m^3 = N/m^2 (if rho in g/cm^3)
        if(mpmOrder[ARCH_StrainEnergy]=='Y')
        {   *(double *)app=1.0e-6*mpm[p]->mp*mpm[p]->GetStrainEnergy();
            app+=sizeof(double);
        }
        
        // material history data on particle (whatever units the material chooses)
        if(mpmOrder[ARCH_History]=='Y')
        {   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(1,mpm[p]->GetHistoryPtr());
            app+=sizeof(double);
        }
		else if(mpmOrder[ARCH_History]!='N')
		{	if(mpmOrder[ARCH_History]&0x01)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(1,mpm[p]->GetHistoryPtr());
				app+=sizeof(double);
			}
			if(mpmOrder[ARCH_History]&0x02)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(2,mpm[p]->GetHistoryPtr());
				app+=sizeof(double);
			}
			if(mpmOrder[ARCH_History]&0x04)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(3,mpm[p]->GetHistoryPtr());
				app+=sizeof(double);
			}
			if(mpmOrder[ARCH_History]&0x08)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(4,mpm[p]->GetHistoryPtr());
				app+=sizeof(double);
			}
		}
		
		// concentration and gradients convert to wt fraction units using csat for this material
        if(mpmOrder[ARCH_Concentration]=='Y')
		{	double csat=theMaterials[mpm[p]->MatID()]->concSaturation;
		
           *(double *)app=mpm[p]->pConcentration*csat;
            app+=sizeof(double);
			
			if(mpm[p]->pDiffusion!=NULL)
			{	*(double *)app=mpm[p]->pDiffusion->Dc.x*csat;
				app+=sizeof(double);
				
				*(double *)app=mpm[p]->pDiffusion->Dc.y*csat;
				app+=sizeof(double);
				
				if(threeD)
				{	*(double *)app=mpm[p]->pDiffusion->Dc.z*csat;
					app+=sizeof(double);
				}
			}
			else
 			{	*(double *)app=0.;
				app+=sizeof(double);
				
				*(double *)app=0.;
				app+=sizeof(double);
				
				if(threeD)
				{	*(double *)app=0.;
					app+=sizeof(double);
				}
			}
       }
		
        // total heat energy (Volume*energy) in J
        // energies in material point based on energy per unit mass
        // here need mass * U/(rho0 V0)
        // internal units are same as stress: N/m^2 cm^3/g = microJ/g = mJ/kg
		// note that rho*energy has units J/m^3 = N/m^2 (if rho in g/cm^3)
        if(mpmOrder[ARCH_HeatEnergy]=='Y')
        {   *(double *)app=1.0e-6*mpm[p]->mp*mpm[p]->GetHeatEnergy();
            app+=sizeof(double);
        }
		
		// element crossings since last archive - now cumulative
        if(mpmOrder[ARCH_ElementCrossings]=='Y')
        {	*(int *)app=mpm[p]->GetElementCrossings();
			app+=sizeof(int);
		}
		
		// here=initial angle z (degrees) while angle(above)=here-0.5*180*wxy/PI (degrees)
		//		Thus 0.5*180*wxy/PI = here-angle(above) or wxy = (PI/90)*(here-angle(above))
		// here=initial angle y (degrees) while angle(above)=here+0.5*180*wrot.xz/PI_CONSTANT (degrees)
		//		Thus 0.5*180*wxz/PI = -(here-angle(above)) or wxz = -(PI/90)*(here-angle(above))
		// here=initial angle x (degrees) while angle(above)=here-0.5*180.*wrot.yz/PI_CONSTANT; (degrees)
		//		Thus 0.5*180*wyz/PI = (here-angle(above)) or wyz = (PI/90)*(here-angle(above))
        if(mpmOrder[ARCH_RotStrain]=='Y')
		{	if(threeD)
			{	*(double *)app=mpm[p]->GetAnglez0InDegrees();
				app+=sizeof(double);
				*(double *)app=mpm[p]->GetAngley0InDegrees();
				app+=sizeof(double);
				*(double *)app=mpm[p]->GetAnglex0InDegrees();
				app+=sizeof(double);
			}
			else
        	{	*(double *)app=mpm[p]->GetAnglez0InDegrees();
				app+=sizeof(double);
			}
		}
		
        // padding
        if(mpmRecSize<recSize)
            app+=recSize-mpmRecSize;
        
        // reversing bytes?
        if(fmobj->GetReverseBytes())
        {   app-=recSize;
        
            // defaults
            app+=Reverse(app,sizeof(int));			// element
            app+=Reverse(app,sizeof(double));		// mass
            app+=Reverse(app,sizeof(short))+2;		// material ID
            app+=Reverse(app,sizeof(double));		// angle
            app+=Reverse(app,sizeof(double));		// thickness or dihedral
            app+=Reverse(app,sizeof(double));		// position
            app+=Reverse(app,sizeof(double));
			if(threeD) app+=Reverse(app,sizeof(double));
            app+=Reverse(app,sizeof(double));		// orig position
            app+=Reverse(app,sizeof(double));
			if(threeD) app+=Reverse(app,sizeof(double));
    
            // velocity (mm/sec)
            if(mpmOrder[ARCH_Velocity]=='Y')
            {	app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
				if(threeD) app+=Reverse(app,sizeof(double));
            }
    
            // stress 2D (in N/m^2)
            if(mpmOrder[ARCH_Stress]=='Y')
            {	app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
				if(threeD)
				{	app+=Reverse(app,sizeof(double));
					app+=Reverse(app,sizeof(double));
				}
            }

            // strain 2D (absolute)
            if(mpmOrder[ARCH_Strain]=='Y')
            {	app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
				if(threeD)
				{	app+=Reverse(app,sizeof(double));
					app+=Reverse(app,sizeof(double));
				}
            }
            
            // plastic strain (absolute)
            if(mpmOrder[ARCH_PlasticStrain]=='Y')
            {	app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
				if(threeD)
				{	app+=Reverse(app,sizeof(double));
					app+=Reverse(app,sizeof(double));
				}
            }
            
            // work energy (cumulative) in J
            if(mpmOrder[ARCH_WorkEnergy]=='Y')
                app+=Reverse(app,sizeof(double));
                    
            // temperature
            if(mpmOrder[ARCH_DeltaTemp]=='Y')
                app+=Reverse(app,sizeof(double));
            
            // total plastic energy
            if(mpmOrder[ARCH_PlasticEnergy]=='Y')
                app+=Reverse(app,sizeof(double));
                    
            // shear components
            if(mpmOrder[ARCH_ShearComponents]=='Y')
            {	app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
            }
    
            // total strain energy
            if(mpmOrder[ARCH_StrainEnergy]=='Y')
                app+=Reverse(app,sizeof(double));
                    
            // material history data on particle
            if(mpmOrder[ARCH_History]=='Y')
                app+=Reverse(app,sizeof(double));
			else if(mpmOrder[ARCH_History]!='N')
			{	if(mpmOrder[ARCH_History]&0x01) app+=Reverse(app,sizeof(double));
				if(mpmOrder[ARCH_History]&0x02) app+=Reverse(app,sizeof(double));
				if(mpmOrder[ARCH_History]&0x04) app+=Reverse(app,sizeof(double));
				if(mpmOrder[ARCH_History]&0x08) app+=Reverse(app,sizeof(double));
			}
                    
            // concentration and gradients
            if(mpmOrder[ARCH_Concentration]=='Y')
            {	app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
                app+=Reverse(app,sizeof(double));
				if(threeD)
					app+=Reverse(app,sizeof(double));
            }
			
            // total strain energy
            if(mpmOrder[ARCH_HeatEnergy]=='Y')
                app+=Reverse(app,sizeof(double));

            // element crossings
            if(mpmOrder[ARCH_ElementCrossings]=='Y')
                app+=Reverse(app,sizeof(int));
			
			// rotational strain
			if(mpmOrder[ARCH_RotStrain]=='Y')
             {	app+=Reverse(app,sizeof(double));
 				if(threeD)
				{	app+=Reverse(app,sizeof(double));
 					app+=Reverse(app,sizeof(double));
				}
            }
			
			// padding
            if(mpmRecSize<recSize)
                app+=recSize-mpmRecSize;
        }
		
		// write this particle (ofstream should buffer for us)
		try
		{	afile.write(aptr,blen);
			if(afile.bad())
				FileError("File error writing material point data",fname,"ArchiveData::ArchiveResults");
		}
		catch(CommonException err)
		{   // give up on hopefully temporary file problem
			cout << "# " << err.Message() << endl;
			cout << "# Will try to continue" << endl;
			afile.close();
			free(aptr);
			return;
		}
    }
    
	// clear material point record buffer
	free(aptr);
    
    // add the cracks
    nextCrack=firstCrack;
    while(nextCrack!=NULL)
    {	nextCrack->Archive(afile);
        nextCrack=(CrackHeader *)nextCrack->GetNextObject();
    }
    
    // close the file
	try
	{	afile.close();
		if(afile.bad())
			FileError("File error closing an archive file",fname,"ArchiveData::ArchiveResults");
	}
	catch(CommonException err)
	{   // give up on hopefully temporary file problem
		cout << "# " << err.Message() << endl;
		cout << "# Will try to continue" << endl;
	}
}

// Archive global results if it is time
void ArchiveData::GlobalArchive(double atime)
{
    // see if desired
	if(globalTime<0.) return;
    if(atime<nextGlobalTime) return;
    nextGlobalTime+=globalTime;
	
	// clear previous ones
	lastArchived.clear();
	
	// each global quantity
	GlobalQuantity *nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
	    nextGlobal=nextGlobal->AppendQuantity(lastArchived);
    
	// time in ms
	char fline[1000],numStr[100];
	sprintf(fline,"%g",1000.*atime);
	int i;
	for(i=0;i<lastArchived.size();i++)
	{	sprintf(numStr,"\t%e",lastArchived[i]);
		strcat(fline,numStr);
	}
	
	// append to global results file
	ofstream global;
	try
	{	global.open(globalFile,ios::out | ios::app);
		if(!global.is_open())
			FileError("File error opening global results",globalFile,"ArchiveData::GlobalArchive");
		global << fline << endl;
		if(global.bad())
			FileError("File error writing global results",globalFile,"ArchiveData::GlobalArchive");
		global.close();
		if(global.bad())
			FileError("File error closing global results",globalFile,"ArchiveData::GlobalArchive");
	}
	
    catch(CommonException err)
	{   // divert to standard output and try to continue
		cout << "# " << err.Message() << endl;
		cout << "# Data: " << fline << endl;
		if(global.is_open()) global.close();
	}
}

// Archive the results if it is time
void ArchiveData::ArchiveVTKFile(double atime,vector< int > quantity,vector< int > quantitySize,
											vector< char * > quantityName,vector< int > qparam,double **vtk)
{
    char fname[300],fline[300];
	
    // get relative path name to the file
    sprintf(fname,"%s%s_%d.vtk",outputDir,archiveRoot,fmobj->mstep);
    
    // open the file
	ofstream afile;
	afile.open(fname, ios::out);
	if(!afile.is_open())
        FileError("Cannot open a vtk archive file",fname,"ArchiveData::ArchiveVTKFile");
	
    // required header line
	afile << "# vtk DataFile Version 4.2" << endl;
	
	// title
	sprintf(fline,"step:%d time:%15.7e ms",fmobj->mstep,1000.*atime);
    afile << fline << endl;
	
	// header
	afile << "ASCII" << endl;
	afile << "DATASET STRUCTURED_POINTS" << endl;
	
	int ptx,pty,ptz;
	mpmgrid.GetGridPoints(&ptx,&pty,&ptz);
	afile << "DIMENSIONS " << ptx << " " << pty;
	if(fmobj->IsThreeD())
		afile << " " << ptz << endl;
	else
		afile << " 1" << endl;
	
	afile << "ORIGIN " << mpmgrid.xmin << " " << mpmgrid.ymin;
	if(fmobj->IsThreeD())
		afile << " " << mpmgrid.zmin << endl;
	else
		afile << " 0" << endl;
	
	afile << "SPACING "  << mpmgrid.gridx << " " << mpmgrid.gridy;
	if(fmobj->IsThreeD())
		afile << " " << mpmgrid.gridz << endl;
	else
		afile << " " << mpmgrid.gridx << endl;
	
	afile << "POINT_DATA " << nnodes << endl;
	
	// export selected data
	int i,offset=0;
	unsigned int q;
	double *vtkquant=NULL;
	double scale=1.;
    
    // contact force special case
    int archiveStepInterval=1;
    if(GetDoingArchiveContact())
    {   archiveStepInterval=fmobj->mstep-lastArchiveContactStep;
        lastArchiveContactStep=fmobj->mstep;
        ZeroVector(&lastContactForce);
    }
	
	for(q=0;q<quantity.size();q++)
	{	// header for next quantity
		switch(quantitySize[q])
		{	case 1:
				if(vtk==NULL) break;
			case -1:
				afile << "SCALARS ";
				afile << quantityName[q];
				afile << " double 1" << endl;
				afile << "LOOKUP_TABLE default" << endl;
				break;
			
			case 3:
				if(vtk==NULL) break;
			case -3:
				afile << "VECTORS ";
				afile << quantityName[q];
                if(quantity[q]==VTK_VOLUMEGRADIENT)
                    afile << qparam[q];
				afile << " double" << endl;
				break;
			
			case 6:
				if(vtk==NULL) break;
				afile << "TENSORS ";
				afile << quantityName[q];
				afile << " double" << endl;
				break;
				
			default:
				break;
		}
	
		// the quantity for each node
		for(i=1;i<=nnodes;i++)
		{	if(vtk!=NULL) vtkquant=vtk[i];
			switch(quantity[q])
			{	case VTK_MASS:
					// mass in g
					afile << nd[i]->GetNodalMass() << endl;
					break;
				
				case VTK_NUMBERPOINTS:
                    // number of points (no rigid will show up)
                    afile << nd[i]->NumberParticles() << endl;
                    break;
					
				case VTK_TEMPERATURE:
					afile << nd[i]->gTemperature << endl;
					break;
				
				case VTK_RIGIDCONTACTFORCES:
                {   Vector fcontact=nd[i]->GetTotalContactForce(TRUE);
					ScaleVector(&fcontact,-1./(double)archiveStepInterval);		// force of rigid particles on the object
					// average over steps since last archive
					afile << fcontact.x << " " << fcontact.y << " " << fcontact.z << endl;
					AddVector(&lastContactForce,&fcontact);
					break;
				}
				
                case VTK_VOLUMEGRADIENT:
                {   Vector grad;
                    nd[i]->GetMatVolumeGradient(qparam[q],&grad);
					afile << grad.x << " " << grad.y << " " << grad.z << endl;
                    break;
                }
					
				case VTK_BCFORCES:
					if(nd[i]->fixedDirection&XYZ_SKEWED_DIRECTION)
					{	//Vector fbc = nd[i]->GetCMatFtot();
						//afile << fbc.x << " " << fbc.y << " " << fbc.z << endl;
					}
					else
						afile << "0. 0. 0." << endl;
					break;
				
				case VTK_CONCENTRATION:
				case VTK_WORKENERGY:
				case VTK_PLASTICENERGY:
				case VTK_MATERIAL:
                case VTK_HEATENERGY:
                case VTK_PRESSURE:
                case VTK_EQUIVSTRESS:
                case VTK_RELDELTAV:
                case VTK_EQUIVSTRAIN:
					if(vtk==NULL) break;
					afile << vtkquant[offset] << endl;
					break;
				
				case VTK_VELOCITY:				// extraplated to always get cm velocity
				case VTK_DISPLACEMENT:
					// Displacement in mm
					if(vtk==NULL) break;
					afile << vtkquant[offset] << " " << vtkquant[offset+1] << " " << vtkquant[offset+2] << endl;
					break;
				
				case VTK_PLASTICSTRAIN:
					scale=1.;
				case VTK_STRESS:
					if(quantity[q]==VTK_STRESS) scale=1.e-6;
				case VTK_STRAIN:
				case VTK_TOTALSTRAIN:
					if(quantity[q]==VTK_STRAIN || quantity[q]==VTK_TOTALSTRAIN) scale=1.;
					// stress in MPa, Strains absolute
					if(vtk==NULL) break;
					afile << scale*vtkquant[offset] << " " << scale*vtkquant[offset+3] << " " << scale*vtkquant[offset+4] << endl;
					afile << scale*vtkquant[offset+3] << " " << scale*vtkquant[offset+1] << " " << scale*vtkquant[offset+5] << endl;
					afile << scale*vtkquant[offset+4] << " " << scale*vtkquant[offset+5] << " " << scale*vtkquant[offset+2] << endl;
					break;
				
				default:
					break;
			}
		}
		
		// offset for next quantity (if in the buffer)
		if(quantitySize[q]>0) offset+=quantitySize[q];
	}
    
    // close the file
	afile.close();
	if(afile.bad())
        FileError("File error closing a vtk archive file",fname,"ArchiveData::ArchiveVTKFile");
}

// Archive the results if it is time
void ArchiveData::ArchiveHistoryFile(double atime,vector< int > quantity)
{
    char fname[300],fline[600],subline[100];
	
    // get relative path name to the file
    sprintf(fname,"%s%s_History_%d.txt",outputDir,archiveRoot,fmobj->mstep);
    
    // open the file
	ofstream afile;
	afile.open(fname, ios::out);
	if(!afile.is_open())
        FileError("Cannot open a particle history archive file",fname,"ArchiveData::ArchiveHistoryFile");
	
    // header line
	afile << "Particle History Data File" << endl;
	
	// title
	sprintf(fline,"step:%d time:%15.7e ms",fmobj->mstep,1000.*atime);
    afile << fline << endl;
	
	strcpy(fline,"#\tx\ty");
    if(threeD) strcat(fline,"\tz");
	unsigned int q;
	for(q=0;q<quantity.size();q++)
	{	sprintf(subline,"\t%d",quantity[q]);
		strcat(fline,subline);
	}
	afile << fline << endl;
    
    // each particle
    int p;
    for(p=0;p<nmpms;p++)
    {   // number and position
        afile << p+1 << "\t" << mpm[p]->pos.x << "\t" << mpm[p]->pos.y ;
        if(threeD) afile << "\t" << mpm[p]->pos.z ;
        
        // history data
        MaterialBase *matptr = theMaterials[mpm[p]->MatID()];
        char *hptr = mpm[p]->GetHistoryPtr();
        for(q=0;q<quantity.size();q++)
        {   afile << "\t" << matptr->GetHistory(quantity[q],hptr);
        }
        afile << endl;
    }
    
    // close the file
	afile.close();
	if(afile.bad())
        FileError("File error closing a particle history archive file",fname,"ArchiveData::ArchiveHistoryFile");
}

// force archive now, but stay on archiving schedule after that
void ArchiveData::ForceArchiving(void)
{	nextArchTime-=archTime;
	if(globalTime>0.) nextGlobalTime-=globalTime;
}

// report a file error to some file
void ArchiveData::FileError(const char *msg,const char *filename,const char *method)
{	char errNo[50];
	sprintf(errNo,"%d",errno);
	char *errMsg=new char[strlen(msg)+strlen(filename)+strlen(errNo)+15];
	strcpy(errMsg,msg);
	strcat(errMsg," (file: ");
	strcat(errMsg,filename);
	strcat(errMsg,", #");
	strcat(errMsg,errNo);
	strcat(errMsg,").");
	throw CommonException(errMsg,method);
}

#pragma mark ArchiveData: Log File Methods

#ifdef LOG_PROGRESS

// empty the file
void ArchiveData::ClearLogFile(void)
{
	// append to global results file
	ofstream logstream;
	try
	{	logstream.open(logFile,ios::trunc);
		if(!logstream.is_open())
			FileError("File error opening log file",globalFile,"ArchiveData::ClearLogFile");
		logstream.close();
		if(logstream.bad())
			FileError("File error closing log file",globalFile,"ArchiveData::ClearLogFile");
		logStartTime=fmobj->CPUTime();
	}
	
    catch(CommonException err)
	{   // divert to standard output and try to continue
		cout << "# " << err.Message() << endl;
		if(logstream.is_open()) logstream.close();
	}
}

// empty the file
void ArchiveData::WriteLogFile(const char *logLine,double *num)
{
	// append to global results file
	ofstream logstream;
	try
	{	logstream.open(logFile,ios::out | ios::app);
		if(!logstream.is_open())
			FileError("File error opening log file",logFile,"ArchiveData::WriteLogFile");
		double logWriteTime=fmobj->CPUTime();
		logstream << logLine;
		double recentTime=logWriteTime-logStartTime;
		if(recentTime>1.e-4)
		{	logstream << " (time: " << logWriteTime-logStartTime;
			if(num!=NULL) logstream << ", number: " << *num;
			logstream << ")";
		}
		else
		{	if(num!=NULL)
				logstream << "(number: " << *num << ")";
		}
		logstream << endl;
		logStartTime=logWriteTime;
		if(logstream.bad())
			FileError("File error writing line to log file",logFile,"ArchiveData::WriteLogFile");
		logstream.close();
		if(logstream.bad())
			FileError("File error closing log file",logFile,"ArchiveData::WriteLogFile");
	}
	
    catch(CommonException err)
	{   // divert to standard output and try to continue
		cout << "# " << err.Message() << endl;
		cout << "# Log Line: " << logLine;
		if(num!=NULL) cout << "(number: " << *num << ")";
		logstream << endl;
		if(logstream.is_open()) logstream.close();
	}
}

#endif

#pragma mark ArchiveData: Accessors

// Return if need to find J and/or K or FALSE if not time to archive or not needed
// if rightNow is FALSE, return if ever needed to get J and K
int ArchiveData::WillArchiveJK(bool rightNow)
{
	if(mtime+timestep<nextArchTime && rightNow) return FALSE;
	if(crackOrder[ARCH_StressIntensity]=='Y')
		return NEED_J+NEED_K;
	else if(crackOrder[ARCH_JIntegral]=='Y')
		return NEED_J;
	else
		return FALSE;
}

// return TRUE if archiving will happen on the next time step
int ArchiveData::WillArchive(void)
{ 	if(mtime+timestep<nextArchTime) return FALSE;
	return TRUE;
}

// Record size after it is calculated
int ArchiveData::GetRecordSize(void) { return recSize; }
 
// set archive orders
void ArchiveData::SetMPMOrder(const char *xData) { strcpy(mpmOrder,xData); }
void ArchiveData::SetCrackOrder(const char *xData) { strcpy(crackOrder,xData); }
bool ArchiveData::PointArchive(int orderBit) { return mpmOrder[orderBit]=='Y'; }
bool ArchiveData::CrackArchive(int orderBit) { return crackOrder[orderBit]=='Y'; }
void ArchiveData::SetDoingArchiveContact(bool setting) { doingArchiveContact=setting; }
bool ArchiveData::GetDoingArchiveContact(void) { return doingArchiveContact; }

// get interval since last VTK archive. Reset if not doing VTK archiving, otherwise reset in ArchiveVTKFile
int ArchiveData::GetArchiveContactStepInterval(void)
{	int archiveStepInterval=fmobj->mstep-lastArchiveContactStep;
	if(!doingArchiveContact) lastArchiveContactStep=fmobj->mstep;
	return archiveStepInterval;
}

// store recent contact force in case needed for global archiving
Vector ArchiveData::GetLastContactForce(void) { return lastContactForce; }
void ArchiveData::SetLastContactForce(Vector fcontact) { lastContactForce=fcontact; }

// check if passed last archve
bool ArchiveData::PassedLastArchived(int qIndex,double criticalValue)
{
	// if no archived yet say false or invalid
	if(qIndex<0 || qIndex>=lastArchived.size()) return false;
	
	//cout << qIndex << "," << criticalValue << "," << lastArchived[qIndex] << "," << (lastArchived[qIndex]>=criticalValue) << endl;
	
	if(criticalValue>=0.)
		return lastArchived[qIndex] >= criticalValue;
	else
		return lastArchived[qIndex] <= criticalValue;
	
}



