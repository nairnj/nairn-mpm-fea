/********************************************************************************
    CommonArchiveData.hpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on Jan 13, 2006.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include <fstream>

#include "System/CommonArchiveData.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Elements/ElementBase.hpp"

/********************************************************************************
	CommonArchiveData: Constructors and Destructor
********************************************************************************/

CommonArchiveData::CommonArchiveData()
{
	inputDir=NULL;			// input directory
	archiveRoot=NULL;		// root file name (no spaces and does not end in '.')
	archiveParent=NULL;		// archive parent folder (not ending in '/') or empty for same folder
	archiveMesh=FALSE;		// list mesh in output file (FALSE) or in files (TRUE)
	forceUnique=FALSE;		// when TRUE forces unique folder to be created for archiving
}

CommonArchiveData::~CommonArchiveData()
{
}

/********************************************************************************
	CommonArchiveData: methods
********************************************************************************/

/*******************************************************************
	Output file mesh - list or send to a file
*******************************************************************/

// Nodal Points
void CommonArchiveData::ArchiveNodalPoints(int np)
{
	char fname[256];
	int i;
	ofstream outfile;
	
	if(archiveMesh)
	{	sprintf(fname,"%s%s_Nodes.txt",inputDir,archiveRoot);
		outfile.open(fname);
		if(outfile)
    	{	sprintf(fname,"%s_Nodes.txt",archiveRoot);
			cout << "File: " << fname << endl << endl;
			for(i=1;i<=nnodes;i++)
				nd[i]->PrintNodalPoint(outfile);
			outfile.close();
			return;
		}
	}
	
	// list in output file (by request or if file error)
    if(np==AXI_SYM)
        cout << " Node         r               z               t" << endl;
    else
        cout << " Node         x               y               z" << endl;
	cout << "--------------------------------------------------------" << endl;
	for(i=1;i<=nnodes;i++)
		nd[i]->PrintNodalPoint(cout);
	cout << endl;
}

// Elements
void CommonArchiveData::ArchiveElements(int np)
{
	char fname[256];
	int i;
	ofstream outfile;
	
	if(archiveMesh)
	{	sprintf(fname,"%s%s_Elems.txt",inputDir,archiveRoot);
		outfile.open(fname);
		if(outfile)
    	{	sprintf(fname,"%s_Elems.txt",archiveRoot);
			cout << "File: " << fname << endl << endl;
			for(i=0;i<nelems;i++)
				theElements[i]->PrintElement(outfile,np);
			outfile.close();
			return;
		}
	}
	
	// list in output file (by request or if file error)
#ifdef MPM_CODE
    cout << "  No. ID    Nd1   Nd2   Nd3   Nd4   Nd5   Nd6   Nd7   Nd8\n"
   	<< "-----------------------------------------------------------\n";
#else
    cout << "  No. ID Mat  Ang(deg) Thick(mm)  Nd1   Nd2   Nd3   Nd4   Nd5   Nd6   Nd7   Nd8\n"
   	<< "---------------------------------------------------------------------------------\n";
#endif
	for(i=0;i<nelems;i++)
        theElements[i]->PrintElement(cout,np);
    cout << endl;
}

/*******************************************************************
	CommonArchiveData: Accessors
*******************************************************************/

// Set archive mesh option
void CommonArchiveData::SetArchiveMesh(bool newSetting) { archiveMesh=newSetting; }

// Set path to directory of the input file
void CommonArchiveData::SetInputDirPath(const char *cmdFile)
{
	// delete old
	if(inputDir!=NULL) delete [] inputDir;
	
	// create new string and remove file name (but keep last '/')
	//  (note: may end up empty if cmdFile is relative path to file in working directory)
	// second half of string will have file name of the input file only
	inputDir=new char[strlen(cmdFile)+2];
	strcpy(inputDir,cmdFile);
    int i=strlen(inputDir)-1;
	inputDir[i+2]=0;			// new trailere on second half of the tring
    while(i>=0)
	{	inputDir[i+1]=inputDir[i];		// shift ith character of file name to the right
    	if(inputDir[i]=='/')
		{	inputDir[i+1]=0;
			break;
		}
        i--;
    }
	if(i<0) inputDir[0]=0;
	
	// string from inputDir is directory path for input file
	// string from &inputDir[strlen(inputDir)+1] is input file name only
}

// attach inputDir to possible relative path name (caller must delete)
char *CommonArchiveData::ExpandInputPath(const char *partialName)
{
	char *path;
	
	// just the name is returned if full path name already
	if(inputDir==NULL || partialName[0]=='/')
	{	path=new char[strlen(partialName)+1];
		strcpy(path,partialName);
	}
	
	// start with inputDir
	else
	{	path=new char[strlen(inputDir)+strlen(partialName)+1];
		strcpy(path,inputDir);
		strcat(path,partialName);
	}
	
	return path;
}

// Set archiveRoot and make sure does not end in period and remove spaces
// Find parent folder (without the terminal /) or empty string if no parent
void CommonArchiveData::SetArchiveRoot(char *newRoot,bool makeUnique)
{
	if(archiveRoot!=NULL) delete [] archiveRoot;
	archiveRoot=new char[strlen(newRoot)+5];		// 5 saves room for unique folder 1/-999/
	strcpy(archiveRoot,newRoot);
	unsigned i=strlen(archiveRoot);
	if(archiveRoot[i-1]=='.') archiveRoot[i-1]=0;
	
	// spaces can cause problems
	for(i=0;i<strlen(archiveRoot);i++)
	{	if(archiveRoot[i]==' ') archiveRoot[i]='_';
	}
	
	// back up to find parent folder in this relative path
	if(archiveParent!=NULL) delete [] archiveParent;
	archiveParent=new char[strlen(archiveRoot)+5];	// 5 saves room for unique folder /1-/999
	strcpy(archiveParent,archiveRoot);
    for(i=strlen(archiveParent);i>0;i--)
    {   if(archiveParent[i]=='/')
		{   archiveParent[i]=0;
			break;
		}
    }
	if(i==0) archiveParent[0]=0;		// empty parent
	forceUnique=makeUnique;	
}
char *CommonArchiveData::GetArchiveRoot(void) { return archiveRoot; }



