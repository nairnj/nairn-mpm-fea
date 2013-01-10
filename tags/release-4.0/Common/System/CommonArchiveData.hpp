/********************************************************************************
    CommonArchiveData.hpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on Jan 15, 2006
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _COMMONARCHIVEDATA_

#define _COMMONARCHIVEDATA_

class CommonArchiveData
{
    public:
	
        //  Constructors and Destructor
        CommonArchiveData();
		virtual ~CommonArchiveData();
		
		// abstract methods
	
		// methods
		void ArchiveNodalPoints(int);
		void ArchiveElements(int);

		// accessors
		bool SetArchiveRoot(char *,bool);
		char *GetArchiveRoot(void);
		void SetArchiveMesh(bool);
		void SetInputDirPath(const char *,bool);
		char *ExpandOutputPath(const char *);
		
	protected:
		char *inputDir,*archiveRoot,*globalFile,*archiveParent;
		char *outputDir;
		bool archiveMesh,forceUnique;
};

#endif
