/********************************************************************************
	MaterialContactNode.hpp
	nairn-mpm-fea

	Created by John Nairn on 3/15/18.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	Dependencies
		none
********************************************************************************/

#ifndef _MATERIALCONTACTNODE_

#define _MATERIALCONTACTNODE_

class NodalPoint;

class MaterialContactNode
{
	public:
	
		// variables (changed in MPM time step)
		static vector< MaterialContactNode * > materialContactNodes;
	
		// constructors and destructors
		MaterialContactNode(NodalPoint *,MaterialContactNode *);
		virtual ~MaterialContactNode();
		void LinkToNode(void);
	
		// common methods
		void NodalMaterialContact(double,int);
		void AddMaterialContactPoint(int,short);
		void PrepareForLists(void);
		NodalPoint *GetTheNode(void);
	
		// accessors
		void SetPrevNode(MaterialContactNode *);
		MaterialContactNode *GetPrevNode(void);
		vector<int> ParticleLists(int);
	
		// class methods
		static bool ContactOnKnownNodes(double,int);
		static void ReleaseContactNodes(void);
	
	protected:
		NodalPoint *theNode;
		MaterialContactNode *prevNode;
		vector< int > *lists;
};

#endif
