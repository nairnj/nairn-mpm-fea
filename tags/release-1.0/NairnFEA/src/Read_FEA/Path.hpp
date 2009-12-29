/********************************************************************************
    PathsController.hpp
    NairnFEA
    
    Created by John Nairn on 6/23/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _PATH_

#define _PATH_

class Keypoint;

enum { FIX_PATH_NODES=0,ROTATE_PATH_NODES,LOAD_PATH_NODES,LOAD_PATH_EDGES };
#define INTERNAL_PATH -1000

class Path : public LinkedObject
{
    public:
		// face=0 (path unmeshed), element face # (meshed once), INTERNAL_PATH (meshed twice)
		// ...  -element face # (meshed area and interface - another area not allowed)
		int intervals,face;
		double ratio;
		char pathID[33];
		int firstMainNode;
		int firstMidsideNode;
		int nodeIncrement;
		int firstElem;
		int elemIncrement;
		bool bcFlags[3];
		Path *subPath[2];
	
        //  Constructors and Destructor
		Path(const char *,int,double);
		
		// methods
		int AddKeypoint(char *);
		void ReorientPath(void);
		void AddEdgeBCsToPath(int,int,double *);
		void FirstKeypointToNode(void);
		void ConvertToRatio(void);
		double BinarySolve(double,double,int,short);
		double Distance(double);
		double NormalSum(double);
		int nodeAtIndex(int);
		const char *VerifyBCOption(int,int);
		
		// Accessors
		int ValidPath(void);
		Keypoint *FirstKeypoint(void);
		Keypoint *LastKeypoint(void);
		void GetFirstAndMidXY(Vector *,Vector *);
		bool IsMeshed(void);
		int SetKeys(Keypoint *,Keypoint *,Keypoint *);
	
	private:
		int numKeypts;
		Keypoint *keys[3];
};

#endif
