/********************************************************************************
    CommonArchiveData.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 13, 2006.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"

#include <fstream>

#include "System/CommonArchiveData.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Elements/ElementBase.hpp"

#ifdef WINDOWS_EXE
char CommonArchiveData::folderDelim = '\\';
#else
char CommonArchiveData::folderDelim = '/';
#endif

#pragma mark CommonArchiveData: Constructors and Destructor

// Construtor
CommonArchiveData::CommonArchiveData()
{
	// path to folder containing the input command file
	// For windows dos command, user should input dos style with backslashes
	inputDir=NULL;

	// path to folder that is starting point for creating archive files
	// For windows dos command will use backslashes
	outputDir=NULL;

	// root file name only (no spaces and does not end in '.')
	// archiveRoot uses unix slashes and archiveDosRoot uses backslashes
	archiveRoot=NULL;
	archiveDosRoot = NULL;

	// archive parent folder (not ending in '/' or backslash) or empty for same folder
	// For windows dos command, uses backslashes
	archiveParent=NULL;

	archiveMesh=false;		// list mesh in output file (FALSE) or in files (TRUE)
	forceUnique=false;		// when TRUE forces unique folder to be created for archiving
}

// destructor
CommonArchiveData::~CommonArchiveData()
{
}

#pragma mark CommonArchiveData: methods

// Nodal Points - list or send to a file
void CommonArchiveData::ArchiveNodalPoints(int np)
{
	char fname[256];
	int i;
	ofstream outfile;
	
	if(archiveMesh)
	{	GetFilePath(fname,"%s%s_Nodes.txt");
		outfile.open(fname);
		if(outfile)
    	{	sprintf(fname,"%s_Nodes.txt",archiveRoot);		// name for output mpm file
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

// Elements - list or send to a file
void CommonArchiveData::ArchiveElements(int np)
{
	char fname[256];
	int i;
	ofstream outfile;
	
	if(archiveMesh)
	{	GetFilePath(fname,"%s%s_Elems.txt");
		outfile.open(fname);
		if(outfile)
    	{	sprintf(fname,"%s_Elems.txt",archiveRoot);			// name for output mpm file
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
    cout << "  No. ID Mat  Ang(deg) Thick(mm)  Nd1   Nd2   Nd3   Nd4   Nd5   Nd6   Nd7   Nd8   Nd9\n"
   	<< "---------------------------------------------------------------------------------------\n";
#endif
	for(i=0;i<nelems;i++)
        theElements[i]->PrintElement(cout,np);
    cout << endl;
}

#pragma mark CommonArchiveData: Accessors

// Set archive mesh option
void CommonArchiveData::SetArchiveMesh(bool newSetting) { archiveMesh=newSetting; }

// Set path to directory of the input file inputDir, but it is only used to copy
//     It is only used when copying inpput commands file to the archive folder
// cmdFile is from command arg and will be relative or full path (dos path for dos command)
// Also set outputDir - it is inputDir or empty if useWorkingDir is true
//     It is used as starting point for archiving any files
// throws std::bad_alloc
void CommonArchiveData::SetInputDirPath(const char *cmdFile,bool useWorkingDir)
{
	// delete old
	if(inputDir!=NULL) delete [] inputDir;
	
	// create new string and remove file name (but keep last '/' or backslash for Windows)
	//  (note: may end up empty if cmdFile is relative path to file in working directory)
	// second half of string will have file name of the input file only
	inputDir=new char[strlen(cmdFile)+2];
	strcpy(inputDir,cmdFile);
    int i=(int)strlen(inputDir)-1;
	inputDir[i+2]=0;			// new trailer on second half of the string
    while(i>=0)
	{	inputDir[i+1]=inputDir[i];		// shift ith character of file name to the right
		if(inputDir[i]==folderDelim)
		{	inputDir[i+1]=0;
			break;
		}
        i--;
    }
	if(i<0) inputDir[0]=0;
	
	// string from inputDir is directory path for input file (or empty if none)
	// string from &inputDir[strlen(inputDir)+1] is input file name only
	// In Windows exe, expect to have back slashes if has a director path
	
	// use outputDir for saving files; it will either be same as inputDir or empty
	//	for current working directory
	if(useWorkingDir)
	{	// empty to save relative to working directory
		outputDir = new char[1];
		outputDir[0] = 0;
	}
	else
	{	// write to same folder as the input file
		outputDir = inputDir;
	}
}

// Get path for file using outputDir. This will be relative to the input file
//  unless keyed to workingDirectory, and then will be relative to that directory instead
// Used when taking input from a file such as BMP files
// throws std::bad_alloc
char *CommonArchiveData::ExpandOutputPath(const char *partialName)
{
	char *path;
	
	// if full path, return it
#ifdef WINDOWS_EXE
	if(partialName[0] == '/')
	{	// convert to dos full path assume c disk
		path=new char[strlen(partialName)+3];
		strcpy(path, "C:");
		strcat(path, partialName);
		MakeDOSPath(path);
		return path;
	}
	else if(partialName[1] == ':' && partialName[2] == '\\')
	{	// assume dos full path
		path=new char[strlen(partialName)+1];
		strcpy(path,partialName);
		return path;
	}
#endif

	// if no outputDir or unix full path, return the name
	if(outputDir==NULL || partialName[0]=='/')
	{	path=new char[strlen(partialName)+1];
		strcpy(path,partialName);
		return path;
	}
	
	// start with outputDir
	path=new char[strlen(outputDir)+strlen(partialName)+1];
	strcpy(path,outputDir);
	strcat(path,partialName);
#ifdef WINDOWS_EXE
	MakeDOSPath(path);
#endif
	return path;
}

// Set archiveRoot and make sure does not end in period and remove spaces
// Find parent folder (without the terminal / or backslash) or empty string if no parent
// User can enter using slashes or back slashes
// throws std::bad_alloc
bool CommonArchiveData::SetArchiveRoot(char *newRoot,bool makeUnique)
{
	// second set not allowed
	if(archiveRoot!=NULL) return false;
	
	archiveRoot=new char[strlen(newRoot)+5];		// 5 saves room for unique folder /1 - /999
	strcpy(archiveRoot,newRoot);
	unsigned i=(unsigned)strlen(archiveRoot);
	if(archiveRoot[i-1]=='.') archiveRoot[i-1]=0;
	
	// spaces can cause problems and default to slashes
	for(i=0;i<strlen(archiveRoot);i++)
	{	if(archiveRoot[i]==' ') archiveRoot[i]='_';
		if(archiveRoot[i]=='\\') archiveRoot[i]='/';
	}
	
	// back up to find parent folder in this relative path
	if(archiveParent!=NULL) delete [] archiveParent;
	archiveParent=new char[strlen(archiveRoot)+5];	// 5 saves room for unique folder /1 - /999
	strcpy(archiveParent,archiveRoot);
    for(i=(int)strlen(archiveParent);i>0;i--)
    {   if(archiveParent[i]=='/')
		{   archiveParent[i]=0;
			break;
		}
    }
	if(i==0) archiveParent[0]=0;		// empty parent
	forceUnique=makeUnique;

#ifdef WINDOWS_EXE
	MakeDOSPath(archiveParent);
	if (archiveDosRoot != NULL) delete[] archiveRoot;
	archiveDosRoot=new char[strlen(archiveRoot)+1];
	strcpy(archiveDosRoot, archiveRoot);
	MakeDOSPath(archiveDosRoot);
#endif

	return true;
}

char *CommonArchiveData::GetArchiveRoot(void) { return archiveRoot; }

// build archiving file path (fname better be long enough)
void CommonArchiveData::GetFilePath(char *fname, const char *prefix)
{	// build name for unix
#ifdef WINDOWS_EXE
	sprintf(fname, prefix, outputDir, archiveDosRoot);
#else
	sprintf(fname, prefix, outputDir, archiveRoot);
#endif
}

// build archiving file path with a number (fname better be long enough)
void CommonArchiveData::GetFilePathNum(char *fname, const char *prefix,int pnum)
{	// build name for unix
#ifdef WINDOWS_EXE
	sprintf(fname, prefix, outputDir, archiveDosRoot,pnum);
#else
	sprintf(fname, prefix, outputDir, archiveRoot,pnum);
#endif
}

// build arelative rchiving file path with a number (fname better be long enough)
void CommonArchiveData::GetRelativeFilePathNum(char *fname, const char *prefix,int pnum)
{	// build name for unix
#ifdef WINDOWS_EXE
	sprintf(fname, prefix, archiveDosRoot,pnum);
#else
	sprintf(fname, prefix, archiveRoot,pnum);
#endif
}

