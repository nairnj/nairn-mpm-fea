/********************************************************************************
	CrossedCrack.hpp
	nairn-mpm-fea

	Created by John Nairn on 10/4/2014.
	Copyright (c) 2014 John A. Nairn, All rights reserved.

	Dependencies
		LinkedObject.hpp
********************************************************************************/

#ifndef _CROSSEDCRACK_

#define _CROSSEDCRACK_

class CrackHeader;

class CrossedCrack : public LinkedObject
{
    public:
		CrackHeader *crack;
	
		CrossedCrack(CrackHeader *aCrack);
};

#endif
