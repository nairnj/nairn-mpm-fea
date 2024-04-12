/*********************************************************************
	Reservoir.hpp
	Nairn Research Group MPM Code
	
	Created by John Nairn on 6/14/2021.
	Copyright (c) 2021, All rights reserved.
	
	Dependencies
		none
*********************************************************************/

#ifndef _JAN_RESERVOIR_

#define _JAN_RESERVOIR_

class MPMBase;

class Reservoir
{
	public:
	
		// constructors
		Reservoir(Vector *,Vector *);
		void output(void);

		// methods
		void AddParticle(MPMBase *);
		void DeleteParticle(MPMBase *);
		MPMBase *InjectParticle(Vector *,Vector *,int);
		void ChangeParticleSize(MPMBase *,Vector *);
		long currentCount(void) const;
		vector <double> GetResizings(void);
		void ClearResizings(void);

	protected:
		// accessors
		int IndexForClass(Vector *,int) const;
	
	protected:
		Vector store;
		Vector storeSize;
		vector <Vector> resSizes;
		vector <int> matIDs;
		vector <MPMBase *> firstPt;
		vector <MPMBase *> lastPt;
		vector <double> resizings;
		long count;
};

extern Reservoir *mpmReservoir;

#endif
