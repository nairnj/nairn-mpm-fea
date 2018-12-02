/********************************************************************************
	ShellController.hpp
	NairnMPM

	Created by John Nairn on 1/22/2014.
	Copyright (c) 2014 John A. Nairn, All rights reserved.

	Dependencies
		BoxController.hpp, ShapeController.hpp
********************************************************************************/

#ifndef _SHELLCONTROLLER_

#define _SHELLCONTROLLER_

#include "Read_XML/BoxController.hpp"

class Expression;

class ShellController : public BoxController
{
	public:
	
		// initialize
		ShellController(int);
		virtual void SetProperty(const char *,char *,CommonReadHandler *);
		virtual bool FinishSetup(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
    
		// accessors
		virtual bool Is2DShape(void);
		virtual const char *GetShapeName(void);
	
	protected:
		Expression *radiusFunction;
		Expression *thicknessFunction;
		static double varHeight;
};

#endif
