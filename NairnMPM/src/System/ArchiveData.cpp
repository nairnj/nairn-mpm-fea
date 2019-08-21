/********************************************************************************
    ArchiveData.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 06 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "stdafx.h"
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
#include "System/UnitsController.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

// archiver global
ArchiveData *archiver;

#pragma mark Constructors and Destructor and Initializers

ArchiveData::ArchiveData()
{
	// clear blocks of archive times
	archTimes.clear();		
	firstArchTimes.clear();
	maxProps.clear();
	
	globalTime=-1.;			// time interval for archiving global results (sec)
	nextGlobalTime=0.;		// next time to archive global results (sec)
	
	globalFile=NULL;		// path to global results file
	decohesionFile=NULL;	// path to decohesion file
	decohesionModes[0]=0;	// observed decohesion modes
	threeD=FALSE;			// three D calculations
	
	// default archive has byte order and defaults only
	char defaultOrder[ARCH_MAXMPMITEMS+1];
	strcpy(defaultOrder,"mY");
	for(int i=2;i<ARCH_MAXMPMITEMS;i++)
		strcat(defaultOrder,"N");
	SetMPMOrder(defaultOrder);
	
	// default archive has byte order and defaults
	strcpy(defaultOrder,"mY");
	for(int i=2;i<ARCH_MAXCRACKITEMS;i++)
		strcat(defaultOrder,"N");
    SetCrackOrder(defaultOrder);
	
	timeStamp = NULL;			// pointer to header
	propgationCounter = 0;					// counts crack propagation
    
    // contact archive to coordinate with global contact archiving
	lastArchiveContactStep = 0;               // last time contact force was archived
	doingArchiveContact = false;
 	contactForce = NULL;
 }

// create archive folder - true if works or false if fails
// throws std::bad_alloc
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
			{
#ifdef WINDOWS_EXE
				if(strlen(archiveParent)>0)
					sprintf(syscmd,"if not exist \"%s%s\\%d\" exit 1",outputDir,archiveParent,folderID);
				else
					sprintf(syscmd,"if not exist \"%s%d\" exit 1",outputDir,folderID);
				int exists=system(syscmd);
#else
				if(strlen(archiveParent)>0)
					sprintf(syscmd,"test -d '%s%s/%d'",outputDir,archiveParent,folderID);
				else
					sprintf(syscmd,"test -d '%s%d'",outputDir,folderID);
				int exists=system(syscmd);
#endif
				if(exists!=0) break;			// zero means it already exists
				folderID++;
			}
			
			// if not found, an error
			if(folderID>=1000) return false;
			
			// adjust archiveParent and archiveRoot if changed by unique subfolder
			int insertPos=(int)strlen(archiveParent);
			if(insertPos>0)
			{	char fldrNum[10];
				sprintf(fldrNum,"/%d",folderID);			// max length is 4, and space was saved for it
				int endPos=(int)strlen(archiveRoot);
				int numLength=(int)strlen(fldrNum);
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
#ifdef WINDOWS_EXE
			// update changed in for DOS files
			MakeDOSPath(archiveParent);
			strcpy(archiveDosRoot, archiveRoot);
			MakeDOSPath(archiveDosRoot);
#endif
		}
		
		// now make the folder
#ifdef WINDOWS_EXE
		strcpy(syscmd,"if not exist \"");
		char *dosPath = new char[strlen(outputDir) + strlen(archiveParent) + 1];
		strcpy(dosPath, outputDir);
		strcat(dosPath, archiveParent);
		strcat(syscmd,MakeDOSPath(dosPath));
		strcat(syscmd,"\"");
		strcat(syscmd, " MD \"");
		strcat(syscmd,dosPath);
		strcat(syscmd,"\"");
		system(syscmd);
		delete[] dosPath;
#else
    	strcpy(syscmd,"mkdir -p '");
		strcat(syscmd,outputDir);
		strcat(syscmd,archiveParent);
		strcat(syscmd,"'");
		system(syscmd);
#endif
	}
	
	// copy input commands
#ifdef WINDOWS_EXE
	strcpy(syscmd,"COPY /Y \"");
	strcat(syscmd, inputDir);		// empty or ends in backslash
	strcat(syscmd, &inputDir[strlen(inputDir) + 1]);		// name only
	strcat(syscmd, "\" \"");
	strcat(syscmd, outputDir);
	strcat(syscmd, archiveDosRoot);
	strcat(syscmd, ".fmcmd\" > NUL");
	system(syscmd);
#else
	strcpy(syscmd,"cp '");
	strcat(syscmd,inputDir);							// input folder
	strcat(syscmd,&inputDir[strlen(inputDir)+1]);		// input name
	strcat(syscmd,"' '");
	strcat(syscmd,outputDir);
	strcat(syscmd,archiveRoot);
	strcat(syscmd,".fmcmd'");
    system(syscmd);
#endif
	
	// test by creating dummy file because return value above may be system dependent
	// If logging progress, keep the temporary file for future use
    FILE *fp;
#ifdef LOG_PROGRESS
	logFile=new char[strlen(outputDir)+strlen(archiveRoot)+6];
	logStartTime=fmobj->CPUTime();
#else
	char *logFile=new char[strlen(outputDir)+strlen(archiveRoot)+6];
#endif
	GetFilePath(logFile,"%s%s.log");
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
// true if works or false if fails
// throws std::bad_alloc
bool ArchiveData::BeginArchives(bool isThreeD,int maxMats)
{
	// set up archive times
	int blocks = (int)archTimes.size();
	if(blocks==0) return false;

	// archve the first one always
	archBlock=0;
	nextArchTime=0.;
	
	// fill in blank initial phase
	if(blocks>1 && archTimes[0]<0)
	{	archTimes[0] = firstArchTimes[1];
		firstArchTimes[0] = 0.;
		
		// but skip if second block starts at zero
		if(firstArchTimes[1]<=0.) archBlock=1;
	}
	
	// check that blocks are ordered in time
	for(int i=archBlock;i<blocks-1;i++)
	{	if(firstArchTimes[i]>=firstArchTimes[i+1]) return false;
	}
	
	// set flag if 3D calculations
	threeD=isThreeD;
	
	// prepare for rigid material contact if multimaterial mode
	if(maxMats>0)
		contactForce = new Vector[maxMats];
	
	// global archiving
	CreateGlobalFile();
	
    // Archive file list headind
    CalcArchiveSize();
	SetArchiveHeader();
    PrintSection("ARCHIVED ANALYSIS RESULTS");
    cout << "Root file name: " << archiveRoot << "." << endl
        << "Archive format: " << mpmOrder << endl
        << "Crack archive format: " << crackOrder << endl << endl
		<< "  Step     Time (" << UnitsController::Label(ALTTIME_UNITS) << ")     Filename" << endl
        << "----------------------------------------------"
        << endl;
	
	return true;
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
    {	for(i=(int)strlen(mpmOrder);i<ARCH_MAXMPMITEMS;i++)
            mpmOrder[i]='N';
    }
    mpmOrder[ARCH_MAXMPMITEMS]=0;
	mpmOrder[ARCH_Defaults]='Y';			// now required to be 'Y'
	mpmOrder[ARCH_OldOrigPosition]='N';		// now requried to be 'N'
	mpmOrder[ARCH_ver2Empty]='N';			// now required to be 'N'
	
	// now forces ARCH_RotStrain if has strain and disable if does not have strain
	mpmOrder[ARCH_RotStrain] = mpmOrder[ARCH_Strain];
	
    if(strlen(crackOrder)<2) strcpy(crackOrder,"mY");
    if(strlen(crackOrder)<ARCH_MAXCRACKITEMS)
    {	for(i=(int)strlen(crackOrder);i<ARCH_MAXCRACKITEMS;i++)
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
    if(mpmOrder[ARCH_DamageNormal]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=2*sizeof(double);
	}
    if(mpmOrder[ARCH_SpinMomentum]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_SpinVelocity]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
    
    // check what will be there for crack segments
	
	/* ARCH_Defaults are
		plane element (int), empty (double), newCrack (short+2), pos (Vector), origPos (Vector),
			above element (int), above (Vector) below element (int), below (Vector)
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
	archHeader[strlen(archHeader)]=(char)strlen(mpmOrder);
	strcat(archHeader,mpmOrder);
	
	// crackOrder
	archHeader[strlen(archHeader)]=(char)strlen(crackOrder);
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
// throws std::bad_alloc
void ArchiveData::CreateGlobalFile(void)
{
    FILE *fp;
	char fline[1000];
	GlobalQuantity *nextGlobal;
	
	// no file if no quantities created
	if(firstGlobal==NULL) return;
	
	// see what files are needed
	bool hasDecohesion = false;
	bool hasGlobalValues = false;
	bool hasTracers = 0;
	nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
	{	switch(nextGlobal->GetQuantity())
		{	case DECOHESION:
				hasDecohesion = true;
				break;
			case UNKNOWN_QUANTITY:
				break;
			default:
				if(nextGlobal->IsTracerParticle())
					hasTracers = true;
				hasGlobalValues = true;
				break;
		}
		nextGlobal = nextGlobal->GetNextGlobal();
	}
	if(!hasDecohesion && !hasGlobalValues) return;

	// create global file and add headers
	if(hasGlobalValues)
	{	// find tracer particles
		if(hasTracers)
		{	for(int p=0;p<nmpms;p++)
			{	nextGlobal=firstGlobal;
				while(nextGlobal!=NULL)
					nextGlobal = nextGlobal->FindTracerParticle(p,&mpm[p]->origpos);
			}
			nextGlobal=firstGlobal;
			while(nextGlobal!=NULL)
				nextGlobal = nextGlobal->SetTracerParticle();
		}

		// get relative path name to the file
		globalFile = new char[strlen(outputDir)+strlen(archiveRoot)+8];
		GetFilePath(globalFile,"%s%s.global");
	
		// create and open the file
		if((fp=fopen(globalFile,"w"))==NULL)
		{	FileError("Global archive file creation failed",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}
	
		// write color and count archives
		strcpy(fline,"#setColor");
		nextGlobal=firstGlobal;
		while(nextGlobal!=NULL)
			nextGlobal = nextGlobal->AppendColor(fline);
		strcat(fline,"\n");
		if(fwrite(fline,strlen(fline),1,fp)!=1)
		{	FileError("Global archive file failed to add colors",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}

		// write name
		strcpy(fline,"#setName");
		nextGlobal=firstGlobal;
		while(nextGlobal!=NULL)
			nextGlobal = nextGlobal->AppendName(fline);
		strcat(fline,"\n");
		if(fwrite(fline,strlen(fline),1,fp)!=1)
		{	FileError("Global archive file failed to add quantity names",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}

		// close the file
		if(fclose(fp)!=0)
		{	FileError("Global archive file failed to close",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}
		
	}
	
	// Create Decohsion file
	if(hasDecohesion)
	{	// get relative path name to the file
		decohesionFile = new char[strlen(outputDir)+strlen(archiveRoot)+8];
		GetFilePath(decohesionFile,"%s%s.decohn");
		
		// create and open the file
		if((fp=fopen(decohesionFile,"w"))==NULL)
		{	FileError("Decohesion file creation failed",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}
		
		// write heading
		strcpy(fline,"t\tmat\tID\tXp\tYp\tZp\tAng1\tAng2\tAng3\tGI\tGII1\tGII2\tGtot\n");
		if(fwrite(fline,strlen(fline),1,fp)!=1)
		{	FileError("Decohesion file failed to add header",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}
		
		// close the file
		if(fclose(fp)!=0)
		{	FileError("Decohesion file failed to close",globalFile,"ArchiveData::CreateGlobalFile");
			return;
		}
	}

	// section in output file
    PrintSection("ARCHIVED GLOBAL RESULTS");
	if(hasGlobalValues)
    	cout << "Global data file: " << archiveRoot << ".global" << endl;
	if(hasDecohesion)
		cout << "Decohesion data file: " << archiveRoot << ".decohn" << endl;
	cout << endl;
}

// Create file in archive folder (outputDir)/(rootName)(fileName)
// throws std::bad_alloc
char *ArchiveData::CreateFileInArchiveFolder(char *fileName)
{
	FILE *fp;
	
	// get relative path name to the file
	char *newFile = new char[strlen(outputDir)+strlen(archiveRoot)+strlen(fileName)+1];
	GetFilePath(newFile,"%s%s");
	strcat(newFile,fileName);
	
	// create and open the file
	if((fp=fopen(newFile,"w"))==NULL) goto abort;
	
	// close the file
	if(fclose(fp)!=0) goto abort;
	
	return newFile;
	
abort:
	FileError("Create of file in acrhive folder failed",newFile,"ArchiveData::CreateFileInArchiveFolder");
	return NULL;
}


#pragma mark ARCHIVING METHODS

// Nodal velocity conditions
void ArchiveData::ArchiveVelocityBCs(BoundaryCondition *firstBC)
{
	BoundaryCondition *nextBC=firstBC;
	
	if(archiveMesh)
	{	char fname[500];
		ofstream outfile;
		GetFilePath(fname,"%s%s_VelBCs.txt");
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
    cout << " Node    DOF ID  Vel (" << UnitsController::Label(CUVELOCITY_UNITS) << ")"
	<< "   Arg (" << UnitsController::Label(BCARG_UNITS) << ")  Angle1  Angle2  Function(:side:grad:pos)\n"
   	     << "---------------------------------------------------------------------------------------------\n";
    nextBC=(BoundaryCondition *)firstBC;
    while(nextBC!=NULL)
		nextBC=nextBC->PrintBC(cout);
    cout << endl;
}

// Particle sizes - file has actual sizes in current dimensions (e.g. mm)
// Full particle size (i.e. "diameter")
void ArchiveData::ArchivePointDimensions(void)
{	char fname[500];
	ofstream outfile;
	GetFilePath(fname,"%s%s_PtDims.txt");
	outfile.open(fname);
	
	if(outfile)
	{	outfile << "Particle dimensions" << endl;
		outfile << "  Pt#      dX      dY      dZ" << endl;
		outfile << "--------------------------------" << endl;
		
		int p;
		Vector lp;
		char nline[200];
		
		if(mpmgrid.IsStructuredEqualElementsGrid())
		{	Vector grid = mpmgrid.GetCellSize();
			for(p=0;p<nmpms;p++)
			{
				mpm[p]->GetDimensionlessSize(lp);
				sprintf(nline,"%7d %g %g %g",p+1,lp.x*grid.x,lp.y*grid.y,lp.z*grid.z);
				outfile << nline << endl;
			}
		}
		outfile.close();
	}
}

// Archive the results if it is time
// throws CommonException()
void ArchiveData::ArchiveResults(double atime)
{
	double rho,rho0;
    double sxx,syy,sxy;
    char fname[500],fline[500];
    int i,p;
    CrackHeader *nextCrack;
	
	// test global archiving based on specified time
	if(firstGlobal!=NULL && globalTime>=0.)
    {	if(atime>nextGlobalTime)
	 	{	GlobalArchive(atime);
			nextGlobalTime+=globalTime;
		}
	}
    
    // If not ready, not archive, unless forces by propagation counts
	if(atime<nextArchTime)
	{	// not ready to archive, unless propagations have happened
		if(maxProps[archBlock]<=0 || propgationCounter<maxProps[archBlock]) return;
		propgationCounter=0;
		nextArchTime = atime;
	}
	
	// increment for next archive time
	nextArchTime += archTimes[archBlock];
	if(archBlock+1<firstArchTimes.size())
	{	if(nextArchTime > firstArchTimes[archBlock+1])
		{	archBlock++;
			nextArchTime = atime + archTimes[archBlock];
		}
	}
	
	// global archive too, if using archive time
	if(firstGlobal!=NULL && globalTime<0.)
		GlobalArchive(atime);
	
    // get relative path name to the file
	GetFilePathNum(fname,"%s%s.%d",fmobj->mstep);
    
    // output step number, time, and file name to results file
    for(i=(int)strlen(fname);i>=0;i--)
    {	if(fname[i]=='/' || fname[i]=='\\') break;
    }
    sprintf(fline,"%7d %15.7e  %s",fmobj->mstep,atime*UnitsController::Scaling(1.e3),&fname[i+1]);
    cout << fline << endl;

    // open the file
	ofstream afile;
	try
	{	afile.open(fname, ios::out | ios::binary);
		if(!afile.is_open())
			FileError("Cannot open an archive file",fname,"ArchiveData::ArchiveResults");
		
		// write header created in SetArchiveHeader
		*timeStamp=(float)(atime*UnitsController::Scaling(1.e3));
		afile.write(archHeader,HEADER_LENGTH);
		if(afile.bad())
			FileError("File error writing archive file header",fname,"ArchiveData::ArchiveResults");
	}
	catch(CommonException& err)
	{   // give up on hopefully temporary file problem
		cout << "# File error - check disk for amount of free space" << endl;
		cout << "# " << err.Message() << endl;
		cout << "# Will skip this file and try to continue" << endl;
		if(afile.is_open()) afile.close();
		return;
	}
	catch(...)
	{	cout << "Unknown exception ArchiveResults() and failed to archive current results" << endl;
		return;
	}
	
	// allocate space for one material point
	long blen=recSize;
	char *aptr = new (std::nothrow) char[blen];
	if(aptr==NULL)
		throw CommonException("Out of memory allocating buffer for archive file","ArchiveData::ArchiveResults");
    
    // all material points
    for(p=0;p<nmpms;p++)
	{	// buffer is for one particle
    	char *app=aptr;
		
		// must have these defaults
	   
		// ------- element ID
        *(int *)app=mpm[p]->ArchiveElemID();
        app+=sizeof(int);
        
		// ------- mass (Legacy units g)
        *(double *)app=mpm[p]->mp;
        app+=sizeof(double);
        
		// ------- material ID
        *(short *)app=mpm[p]->ArchiveMatID();
        app+=sizeof(short);
		// fill in two zeros for byte alignment
		*app=0;
		app+=1;
		*app=0;
		app+=1;
        
		// ------- 3D has three angles, 2D has one angle and thickness
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
        
		// ------- (x,y,z) position (Legacy units mm)
        *(double *)app=mpm[p]->pos.x;
        app+=sizeof(double);
        
        *(double *)app=mpm[p]->pos.y;
        app+=sizeof(double);
		
		if(threeD)
		{	*(double *)app=mpm[p]->pos.z;
			app+=sizeof(double);
		}

		// ------- original (x,y,z) position (Legacy units mm)
        *(double *)app=mpm[p]->origpos.x;
        app+=sizeof(double);
                
        *(double *)app=mpm[p]->origpos.y;
        app+=sizeof(double);

		if(threeD)
		{	*(double *)app=mpm[p]->origpos.z;
			app+=sizeof(double);
		}

        // ------- velocity (Legacy units mm/sec)
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

        // ------- stress
		// Tracked stress is (Kirchoff Stress)/rho0 = (Cauchy stress)/rho
		// Convert to actual stress (Legacy units Pa)
        int matid = mpm[p]->MatID();
        rho0=theMaterials[matid]->GetRho(mpm[p]);
        rho = rho0/theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
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

        // ------- elastic strain (absolute)
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
        
        // ------- plastic strain (absolute)
        if(mpmOrder[ARCH_PlasticStrain]=='Y')
		{	Tensor *eplast=mpm[p]->GetAltStrainTensor();
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
        
        // ------- external work (cumulative) (Legacy units J)
        if(mpmOrder[ARCH_WorkEnergy]=='Y')
		{	*(double *)app = UnitsController::Scaling(1.e-9)*mpm[p]->mp*mpm[p]->GetWorkEnergy();
            app+=sizeof(double);
        }
                
        // ------- temperature (K)
        if(mpmOrder[ARCH_DeltaTemp]=='Y')
        {   *(double *)app=mpm[p]->pTemperature;
			//*(double *)app=mpm[p]->pPreviousTemperature;
            app+=sizeof(double);
        }
        
        // ------- total plastic energy (Volume*energy) (Legacy units J)
        // energies in material point based on energy per unit mass
         if(mpmOrder[ARCH_PlasticEnergy]=='Y')
        {   *(double *)app = UnitsController::Scaling(1.e-9)*mpm[p]->mp*mpm[p]->GetPlastEnergy();
            app+=sizeof(double);
        }
                
        // ------- shear components (absolute)
        if(mpmOrder[ARCH_ShearComponents]=='Y')
		{	Matrix3 gradU = mpm[p]->GetDisplacementGradientMatrix();
            *(double *)app=gradU(0,1);
            app+=sizeof(double);
                
            *(double *)app=gradU(1,0);
            app+=sizeof(double);
        }

        // ------- total energy (Volume*energy) (Legacy units J)
        // energies in material point based on energy per unit mass
        if(mpmOrder[ARCH_StrainEnergy]=='Y')
        {   *(double *)app = UnitsController::Scaling(1.e-9)*mpm[p]->mp*mpm[p]->GetStrainEnergy();
            app+=sizeof(double);
        }
        
        // ------- material history data on particle (whatever units the material chooses)
        if(mpmOrder[ARCH_History]=='Y')
        {   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(1,mpm[p]->GetHistoryPtr(0));
            app+=sizeof(double);
        }
		else if(mpmOrder[ARCH_History]!='N')
		{	if(mpmOrder[ARCH_History]&0x01)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(1,mpm[p]->GetHistoryPtr(0));
				app+=sizeof(double);
			}
			if(mpmOrder[ARCH_History]&0x02)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(2,mpm[p]->GetHistoryPtr(0));
				app+=sizeof(double);
			}
			if(mpmOrder[ARCH_History]&0x04)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(3,mpm[p]->GetHistoryPtr(0));
				app+=sizeof(double);
			}
			if(mpmOrder[ARCH_History]&0x08)
			{   *(double *)app=theMaterials[mpm[p]->MatID()]->GetHistory(4,mpm[p]->GetHistoryPtr(0));
				app+=sizeof(double);
			}
		}
		
		// ------- concentration and gradients convert to wt fraction units using csat for this material
		// for pore pressure it is poroelasticity
        if(mpmOrder[ARCH_Concentration]=='Y')
		{	double csat=mpm[p]->GetConcSaturation();
			
			*(double *)app=mpm[p]->pConcentration*csat;
            app+=sizeof(double);
			
			if(mpm[p]->pDiffusion!=NULL)
			{	*(double *)app=mpm[p]->pDiffusion[gGRADx]*csat;
				app+=sizeof(double);
				
				*(double *)app=mpm[p]->pDiffusion[gGRADy]*csat;
				app+=sizeof(double);
				
				if(threeD)
				{	*(double *)app=mpm[p]->pDiffusion[gGRADz]*csat;
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
		
        // ------- total heat energy (Legacy units J)
        // energies in material point based on energy per unit mass
        if(mpmOrder[ARCH_HeatEnergy]=='Y')
        {   *(double *)app = UnitsController::Scaling(1.e-9)*mpm[p]->mp*mpm[p]->GetHeatEnergy();
            app+=sizeof(double);
        }
		
		// ------- element crossings since last archive - now cumulative
        if(mpmOrder[ARCH_ElementCrossings]=='Y')
        {	*(int *)app=mpm[p]->GetElementCrossings();
			app+=sizeof(int);
		}
		
		// ------- initial rotation angle
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
		
		// Normal vector for damaged, softening materials
        if(mpmOrder[ARCH_DamageNormal]=='Y')
		{	Vector dnorm = theMaterials[mpm[p]->MatID()]->GetDamageNormal(mpm[p],threeD);
			if(threeD)
			{	*(double *)app=dnorm.x;
				app+=sizeof(double);
				*(double *)app=dnorm.y;
				app+=sizeof(double);
				*(double *)app=dnorm.z;
				app+=sizeof(double);
			}
			else
			{	*(double *)app=dnorm.x;
				app+=sizeof(double);
				*(double *)app=dnorm.y;
				app+=sizeof(double);
			}
		}

		// Particle spin momentum (Legacy Units J-sec)
		// zero unless tracking particle spin
		if(mpmOrder[ARCH_SpinMomentum]=='Y')
		{
			Vector Lp = MakeVector(0.,0.,0.);
			double Lscale = 1.;
			if(threeD)
			{	*(double *)app=Lscale*Lp.x;
				app+=sizeof(double);
				*(double *)app=Lscale*Lp.y;
				app+=sizeof(double);
				*(double *)app=Lscale*Lp.z;
				app+=sizeof(double);
			}
			else
			{	*(double *)app=Lscale*Lp.z;
				app+=sizeof(double);
			}
		}
		
		// Particle spin velocity (Legacy units 1/sec)
		// zero unless using affine MPM methods
		if(mpmOrder[ARCH_SpinVelocity]=='Y')
		{
			// angular velocity
			Vector wp = MakeVector(0.,0.,0.);
			
			// add to archive
			if(threeD)
			{	*(double *)app=wp.x;
				app+=sizeof(double);
				*(double *)app=wp.y;
				app+=sizeof(double);
				*(double *)app=wp.z;
				app+=sizeof(double);
			}
			else
			{	*(double *)app=wp.z;
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
			
			// softening material damage normal
			if(mpmOrder[ARCH_DamageNormal]=='Y')
			{	app+=Reverse(app,sizeof(double));
				app+=Reverse(app,sizeof(double));
 				if(threeD)
					app+=Reverse(app,sizeof(double));
            }

			// particle angular momentum
			if(mpmOrder[ARCH_SpinMomentum]=='Y')
			{	app+=Reverse(app,sizeof(double));
 				if(threeD)
				{	app+=Reverse(app,sizeof(double));
					app+=Reverse(app,sizeof(double));
				}
            }
			
			// particle angular velocity
			if(mpmOrder[ARCH_SpinMomentum]=='Y')
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
		catch(CommonException& err)
		{   // give up on hopefully temporary file problem
			cout << "# File error - check disk for amount of free space" << endl;
			cout << "# " << err.Message() << endl;
			cout << "# Will skip this file and try to continue" << endl;
			afile.close();
			delete [] aptr;
			return;
		}
    }
    
	// clear material point record buffer
	delete [] aptr;
    
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
	catch(CommonException& err)
	{   // give up on hopefully temporary file problem
		cout << "# " << err.Message() << endl;
		cout << "# Will leave file open and try to continue" << endl;
	}
}

// Archive global results if it is time
void ArchiveData::GlobalArchive(double atime)
{
	// possible if archiving decohesion, but not anything else
	if(globalFile==NULL) return;
	
	// clear previous ones
	lastArchived.clear();
	
	// each global quantity
	GlobalQuantity *nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
	    nextGlobal=nextGlobal->AppendQuantity(lastArchived);
    
	// time (Legacy units ms)
	char fline[1000],numStr[100];
	sprintf(fline,"%g",UnitsController::Scaling(1000.)*atime);
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
	
    catch(CommonException& err)
	{   // divert to standard output and try to continue
		cout << "# File error - check disk for amount of free space" << endl;
		cout << "# " << err.Message() << endl;
		cout << "# Data: " << fline << endl;
		if(global.is_open()) global.close();
	}
}

// Archive global results if it is time
void ArchiveData::Decohesion(double atime,MPMBase *mptr,double alpha,double beta,double gamma,
							 		double GI,double GII1,double GII2,double decohesionCode)
{
	// write in the main results file instead
	if(decohesionFile==NULL)
	{
#pragma omp critical (output)
		{	cout << "# Decohesion: t=" << atime*UnitsController::Scaling(1000.);
			cout << " x=(" << mptr->pos.x << "," << mptr->pos.y << "," << mptr->pos.z << ")";
			cout << " angles=(" << alpha << "," << beta << "," << gamma << ")";
			cout << "  (GI,GII1,GII2,Gtot)=(" << GI << "," << GII1 << "," << GII2 << "," << GI+GII1+GII2 << ")" << endl;
		}
		return;
	}
	
#pragma omp critical (output)
	{	// check for first in this mode
		int mode=0;
		int code=(int)(100*decohesionCode+0.5);
		while(decohesionModes[mode]>0 && code!=decohesionModes[mode]) mode++;
	
		// if not found, then new mode (better be 10 or fewer options in any softening material)
		if(decohesionModes[mode]==0)
		{	decohesionModes[mode] = code;
			decohesionModes[mode+1] = 0;
			cout << "# Decohesion initiation t=" << atime*UnitsController::Scaling(1000.);
			cout << " mode=" << decohesionCode << endl;
		}
	}
	
	// time (Legacy units ms)
	char fline[1000],numStr[200];
	sprintf(fline,"%g",UnitsController::Scaling(1000.)*atime);
	
	// material ID
	sprintf(numStr,"\t%d\t%.2lf",mptr->MatID()+1,decohesionCode);
	strcat(fline,numStr);
	
	// position
	sprintf(numStr,"\t%g\t%g\t%g",mptr->pos.x,mptr->pos.y,mptr->pos.z);
	strcat(fline,numStr);

	// angles
	sprintf(numStr,"\t%g\t%g\t%g",alpha,beta,gamma);
	strcat(fline,numStr);

	// energies
	sprintf(numStr,"\t%g\t%g\t%g\t%g",GI,GII1,GII2,GI+GII1+GII2);
	strcat(fline,numStr);

#pragma omp critical (decohesion)
	{	// append to global results file
		ofstream global;
		try
		{	global.open(decohesionFile,ios::out | ios::app);
			if(!global.is_open())
				FileError("File error opening decohesion results",decohesionFile,"ArchiveData::Decohesion");
			global << fline << endl;
			if(global.bad())
				FileError("File error writing decohesion results",decohesionFile,"ArchiveData::Decohesion");
			global.close();
			if(global.bad())
				FileError("File error closing decohesion results",decohesionFile,"ArchiveData::Decohesion");
		}
	
		catch(CommonException& err)
		{   // divert to standard output and try to continue
			cout << "# File error - check disk for amount of free space" << endl;
			cout << "# " << err.Message() << endl;
			cout << "# Data: " << fline << endl;
			if(global.is_open()) global.close();
		}
	}
}

// Archive the results if it is time
void ArchiveData::ArchiveVTKFile(double atime,vector< int > quantity,vector< int > quantitySize,
											vector< char * > quantityName,vector< int > qparam,double **vtk,int onemat)
{
    char fname[300],fline[300];
	
    // get relative path name to the file
	if(onemat<0)
		GetFilePathNum(fname,"%s%s_%d.vtk",fmobj->mstep);
	else
	{	char matname[20];
		sprintf(matname,"%d",onemat);
		strcpy(fline,"%s%s_mat_");
		strcat(fline,matname);
		strcat(fline,"_%d.vtk");
		GetFilePathNum(fname,fline,fmobj->mstep);
	}
    
    // open the file
	ofstream afile;
	afile.open(fname, ios::out);
	if(!afile.is_open())
        FileError("Cannot open a vtk archive file",fname,"ArchiveData::ArchiveVTKFile");
	
    // required header line
	afile << "# vtk DataFile Version 4.0" << endl;
	
	// title (Legacy time units ms)
	sprintf(fline,"step:%d time:%15.7e %s",fmobj->mstep,atime*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
	afile << fline << endl;
	
	// header
	afile << "ASCII" << endl;
	
	// The points
	int ptx,pty,ptz;
	mpmgrid.GetGridPoints(&ptx,&pty,&ptz);
	if(mpmgrid.IsStructuredEqualElementsGrid())
	{	// regular grid ewith equal element sizes
		afile << "DATASET STRUCTURED_POINTS" << endl;
		
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
		
		Vector csz = mpmgrid.GetCellSize();
		afile << "SPACING "  << csz.x << " " << csz.y;
		if(fmobj->IsThreeD())
			afile << " " << csz.z << endl;
		else
			afile << " " << csz.x << endl;
	}
	
	// title (Legacy time units ms)
	afile << "FIELD FieldData 2" << endl;
	afile << "TIME 1 1 double" << endl;
	afile << atime*UnitsController::Scaling(1.e3) << endl;
	afile << "STEP 1 1 int" << endl;
	afile << fmobj->mstep << endl;
	
	// the data
	afile << "POINT_DATA " << nnodes << endl;
	
	// export selected data
	int i,offset=0;
	unsigned int q;
	double *vtkquant=NULL;
	
    // contact force special case
    int archiveStepInterval=1;
    if(GetDoingArchiveContact())
    {   archiveStepInterval=fmobj->mstep-lastArchiveContactStep;
        lastArchiveContactStep=fmobj->mstep;
		for(int im=0;im<maxMaterialFields;im++) ZeroVector(&contactForce[im]);
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
			{	// Non buffered quantitities (don't respect material selection)
				case VTK_NUMBERPOINTS:
                    // number of points (including rigid contact and mirrored fields)
                    afile << nd[i]->NumberParticles() << endl;
                    break;
				
				case VTK_TEMPERATURE:
					afile << nd[i]->gCond.gTValue << endl;
					break;
				
				case VTK_RIGIDCONTACTFORCES:
				{	// contact force (Legacy units N)
					// average over steps since last archive
					double scale = -1./(double)archiveStepInterval;
					Vector fcontact;
					nd[i]->AddGetContactForce(true,contactForce,scale,&fcontact);
					afile << fcontact.x << " " << fcontact.y << " " << fcontact.z << endl;
					break;
				}
                
                case VTK_VOLUMEGRADIENT:
                {   Vector grad;
                    nd[i]->GetMatVolumeGradient(qparam[q],&grad);			//  qparam[q] is material number
					afile << grad.x << " " << grad.y << " " << grad.z << endl;
                    break;
                }
				
				case VTK_BCFORCES:
					// currently not implement (or documented)
					if(nd[i]->fixedDirection&XYZ_SKEWED_DIRECTION)
					{	//Vector fbc = nd[i]->GetCMatFtot();
						//afile << fbc.x << " " << fbc.y << " " << fbc.z << endl;
						afile << "0. 0. 0." << endl;
					}
					else
						afile << "0. 0. 0." << endl;
					break;
				
				// scalars
				case VTK_MASS:
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
				
				// vectors
				case VTK_VELOCITY:
				case VTK_DISPLACEMENT:
					// Displacement (Legacy units mm)
					if(vtk==NULL) break;
					afile << vtkquant[offset] << " " << vtkquant[offset+1] << " " << vtkquant[offset+2] << endl;
					break;
				
				// tensors
				case VTK_PLASTICSTRAIN:
				case VTK_STRESS:
				case VTK_STRAIN:
				case VTK_TOTALSTRAIN:
					// stress Legacy units MPa, strains are absolute
					if(vtk==NULL) break;
					afile << vtkquant[offset] << " " << vtkquant[offset+3] << " " << vtkquant[offset+4] << endl;
					afile << vtkquant[offset+3] << " " << vtkquant[offset+1] << " " << vtkquant[offset+5] << endl;
					afile << vtkquant[offset+4] << " " << vtkquant[offset+5] << " " << vtkquant[offset+2] << endl;
					break;
				
				// tensor in different order
				case VTK_DEFGRAD:
					if(vtk==NULL) break;
					afile << vtkquant[offset] << " " << vtkquant[offset+1] << " " << vtkquant[offset+2] << endl;
					afile << vtkquant[offset+3] << " " << vtkquant[offset+4] << " " << vtkquant[offset+5] << endl;
					afile << vtkquant[offset+6] << " " << vtkquant[offset+7] << " " << vtkquant[offset+8] << endl;
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
    GetFilePathNum(fname,"%s%s_History_%d.txt",fmobj->mstep);
    
    // open the file
	ofstream afile;
	afile.open(fname, ios::out);
	if(!afile.is_open())
        FileError("Cannot open a particle history archive file",fname,"ArchiveData::ArchiveHistoryFile");
	
    // header line
	afile << "Particle History Data File" << endl;
	
	// title
	sprintf(fline,"step:%d time:%15.7e %s",fmobj->mstep,atime*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
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
        char *hptr = mpm[p]->GetHistoryPtr(0);
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
{	nextArchTime-=archTimes[archBlock];
	if(firstGlobal!=NULL && globalTime>=0.) nextGlobalTime-=globalTime;
}

// report a file error to some file
// throws CommonException()
void ArchiveData::FileError(const char *msg,const char *filename,const char *method)
{	char errNo[50];
	sprintf(errNo,"%d",errno);
	char *errMsg=new (std::nothrow) char[strlen(msg)+strlen(filename)+strlen(errNo)+15];
	if(errMsg !=NULL)
	{	strcpy(errMsg,msg);
		strcat(errMsg," (file: ");
		strcat(errMsg,filename);
		strcat(errMsg,", #");
		strcat(errMsg,errNo);
		strcat(errMsg,").");
		throw CommonException(errMsg,method);
	}
	else
		throw CommonException("File error",method);
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
			FileError("File error opening log file",logFile,"ArchiveData::ClearLogFile");
		logstream.close();
		if(logstream.bad())
			FileError("File error closing log file",logFile,"ArchiveData::ClearLogFile");
		logStartTime=fmobj->CPUTime();
	}
	
    catch(CommonException& err)
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
	
    catch(CommonException& err)
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
bool ArchiveData::WillArchive(void)
{ 	if(mtime+timestep<nextArchTime) return false;
	return true;
}

// Record size after it is calculated
int ArchiveData::GetRecordSize(void) { return recSize; }
 
// set archive orders
void ArchiveData::SetMPMOrder(const char *xData) { strcpy(mpmOrder,xData); }
void ArchiveData::SetMPMOrderByte(int byteNum,char setting) { mpmOrder[byteNum] = setting; }
char ArchiveData::GetMPMOrderByte(int byteNum) { return mpmOrder[byteNum]; }
void ArchiveData::SetCrackOrder(const char *xData) { strcpy(crackOrder,xData); }
void ArchiveData::SetCrackOrderByte(int byteNum,char setting) { crackOrder[byteNum] = setting; }
char ArchiveData::GetCrackOrderByte(int byteNum) { return crackOrder[byteNum]; }
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

// get last archived value
double ArchiveData::GetLastArchived(int qIndex)
{	// if no archived yet say false or invalid
	if(qIndex<0 || qIndex>=lastArchived.size()) return 0.;
	return lastArchived[qIndex];
}

// Propgation Counter
void ArchiveData::IncrementPropagationCounter(void) { propgationCounter++; }
void ArchiveData::SetMaxiumPropagations(int maxp)
{	int blocks = (int)maxProps.size();
	if(blocks>0) maxProps[blocks-1] = maxp;
}

// archive times
// allows settings, so increase size by one to save time, start time and max props
double *ArchiveData::GetArchTimePtr(void)
{	archTimes.push_back(0.);
	firstArchTimes.push_back(0.);
	maxProps.push_back(0);
	int blocks = (int)archTimes.size();
	return &archTimes[blocks-1];
}

// first archive times
// error if none ready
// if first block, create new first block
double *ArchiveData::GetFirstArchTimePtr(void)
{	int blocks = (int)archTimes.size();
	if(blocks==0)
		return NULL;
	else if(blocks==1)
	{	// create initial block with time larger then first archive time
		// but don't know that time yet - it is set at the end as flagged by initial time step < 0
		archTimes.push_back(archTimes[0]);
		archTimes[0] = -1.;
		firstArchTimes.push_back(0.);
		maxProps.push_back(maxProps[0]);
		maxProps[0] = 0.;
		blocks++;
	}
	return &firstArchTimes[blocks-1];
}

// global time pointer
double *ArchiveData::GetGlobalTimePtr(void) { return &globalTime; }

// for contact forces
Vector *ArchiveData::GetLastContactForcePtr(void) { return contactForce; }


