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
class ContactLaw;

class MaterialContactNode
{
	public:
		Vector contactNorm;		// save for VTK archiving
		// list of material contact nodes this time step
		static vector< MaterialContactNode * > materialContactNodes;
	
		// constructors and destructors
		MaterialContactNode(NodalPoint *,MaterialContactNode *);
		virtual ~MaterialContactNode();
		void LinkToNode(void);
	
		// common methods
		void NodalMaterialContact(double,int);
		void AddMaterialContactPoint(int,short);
		void PrepareForLists(void);
		void NodalXPICIncrement(double,int);
	
		// accessors
		NodalPoint *GetTheNode(void);
		void SetPrevNode(MaterialContactNode *);
		MaterialContactNode *GetPrevNode(void);
		vector<int> ParticleLists(int);
		FMPMContact *GetContactInfo(int vfld,int matnum);
	
		// class methods
		static bool ContactOnKnownNodes(double,int);
		static void ReleaseContactNodes(void);
		static void GetMaterialContactData(void);
		static void SetContactInfo(FMPMContact *,double,int,Vector *,Vector *,ContactLaw *,double);
		static void ContactInterfaceInfo(FMPMContact *,double,Vector *);
		static void ContactSetLowMass(FMPMContact *,bool,int);
		static void ContactSetDelPi(FMPMContact *,Vector *,double);

	protected:
		NodalPoint *theNode;
		MaterialContactNode *prevNode;
		vector< int > *lists;
		FMPMContact *mmContact[4];			// better is MAX_FIELD_FOR_CRACKS
};

#endif
