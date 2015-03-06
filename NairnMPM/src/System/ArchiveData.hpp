/********************************************************************************
    ArchiveData.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    

	Dependencies
		CommonArchiveData.hpp
********************************************************************************/

#ifndef _ARCHIVEDATA_

#define _ARCHIVEDATA_

#include "System/CommonArchiveData.hpp"

class BoundaryCondition;

// Archiving both points and cracks
#define ARCH_ByteOrder 0
#define ARCH_Defaults 1

// Archiving options for material points
enum { ARCH_Velocity=2,ARCH_Stress,ARCH_Strain,ARCH_PlasticStrain,
            ARCH_OldOrigPosition,ARCH_WorkEnergy,ARCH_DeltaTemp,ARCH_PlasticEnergy,
            ARCH_ver2Empty, ARCH_ShearComponents, ARCH_StrainEnergy,
            ARCH_History, ARCH_Concentration,ARCH_HeatEnergy,ARCH_ElementCrossings,
            ARCH_RotStrain, ARCH_MAXMPMITEMS };

// Archiving options for crack segments
enum { ARCH_JIntegral=2,ARCH_StressIntensity,ARCH_BalanceResults,ARCH_MAXCRACKITEMS };

// VTK archiving options
enum { VTK_MASS=0, VTK_VELOCITY, VTK_STRESS, VTK_STRAIN, VTK_DISPLACEMENT, VTK_PLASTICSTRAIN,
		VTK_TEMPERATURE, VTK_CONCENTRATION, VTK_WORKENERGY, VTK_PLASTICENERGY,VTK_MATERIAL,
		VTK_RIGIDCONTACTFORCES, VTK_TOTALSTRAIN, VTK_PRESSURE, VTK_EQUIVSTRESS, VTK_RELDELTAV,
        VTK_EQUIVSTRAIN, VTK_HEATENERGY, VTK_BCFORCES, VTK_VOLUMEGRADIENT, VTK_NUMBERPOINTS };

#define HEADER_LENGTH 64

class ArchiveData : public CommonArchiveData
{
    public:
		bool threeD;
		
        //  Constructors and Destructor
		ArchiveData();
		
		// methods
		void ArchiveVelocityBCs(BoundaryCondition *);
		bool MakeArchiveFolder(void);
		bool BeginArchives(bool);
		void ArchiveResults(double);
		void ArchiveVTKFile(double,vector< int >,vector< int >,vector< char * >,vector< int >,double **);
		void ArchiveHistoryFile(double,vector< int >);
		void FileError(const char *,const char *,const char *);
		
		// log file methods
#ifdef LOG_PROGRESS
		void ClearLogFile(void);
		void WriteLogFile(const char *,double *);
#endif
		
		// accessors
		int WillArchiveJK(bool);
		int WillArchive(void);
		void ForceArchiving(void);
		int GetRecordSize(void);
		void SetMPMOrder(const char *);
		void SetCrackOrder(const char *);
		bool PointArchive(int);
		bool CrackArchive(int);
		void SetDoingArchiveContact(bool);
		bool GetDoingArchiveContact(void);
		int GetArchiveContactStepInterval(void);
		Vector GetLastContactForce(void);
		void SetLastContactForce(Vector);
		bool PassedLastArchived(int,double);
		void IncrementPropagationCounter(void);
		void SetMaxiumPropagations(int);
		double *GetArchTimePtr(void);
		double *GetFirstArchTimePtr(void);
		double *GetGlobalTimePtr(void);
	
	private:
		int archBlock;
		vector<double> archTimes,firstArchTimes,maxProps;
		int propgationCounter;
		double nextArchTime,nextGlobalTime,globalTime;
	
		char *globalFile;
		int recSize;							// archive record size
		int mpmRecSize;							// particle record size
		int crackRecSize;						// crack particle record size
		char mpmOrder[50];						// flags for archiving of particles
		char crackOrder[20];					// flags for archiving of crack particles
		char archHeader[HEADER_LENGTH+1];		// compiler header information once
		int lastArchiveContactStep;					// last time VTK archive was written
		bool doingArchiveContact;					// true if VTK archiving is active
		float *timeStamp;						// pointer to header location for time
		Vector lastContactForce;				// last contact force archived
		vector<double> lastArchived;			// store last global quantitites save to global archive
#ifdef LOG_PROGRESS
		char *logFile;							// file for tracking progress
		double logStartTime;
#endif
	
		// methods
		void CalcArchiveSize(void);
		void SetArchiveHeader(void);
		void GlobalArchive(double);
		void CreateGlobalFile(void);
};

extern ArchiveData *archiver;

#endif
