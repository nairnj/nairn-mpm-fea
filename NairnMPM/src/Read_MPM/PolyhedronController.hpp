/********************************************************************************
	PolyhedronController.hpp
	nairn-mpm-fea

	Created by John Nairn on 1/7/11.
	Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _POLYHEDRONCONTROLLER_

#define _POLYHEDRONCONTROLLER_

#include "Read_XML/ShapeController.hpp"
#include "System/MPMPrefix.hpp"


class PolyTriangle;

enum { NO_FACES=0,TRICLINIC_POINTS,TRICLINIC_VECTORS,BOX_CORNERS,PYRAMID } ;

class PolyhedronController : public ShapeController
{
	public:
    
		// constructors
        PolyhedronController(int);
		virtual ~PolyhedronController();
	
        // Initialize
        virtual void SetParameter(const char *,const char *);
        virtual void SetProperty(char *,CommonReadHandler *);
        virtual bool FinishParameter(void);
        virtual bool FinishSetup(void);
        virtual bool HasAllParameters(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
    
        // accessors
		virtual bool Is2DShape(void);
		virtual const char *GetShapeName(void);
	
	private:
		int style;
		vector< PolyTriangle * > faces;
		Vector pmin,pmax;
		char order[9];
	
};

#endif

