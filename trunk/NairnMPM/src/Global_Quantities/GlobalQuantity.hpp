/********************************************************************************
    GlobalQuantity.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Jan 12 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Dependencies
		none
********************************************************************************/

#ifndef _GLOBALQUANTITY_

#define _GLOBALQUANTITY_

// possible global averages
enum { UNKNOWN_QUANTITY,AVG_SXX,AVG_SYY,AVG_SXY,AVG_SZZ,AVG_SXZ,AVG_SYZ,
			AVG_EXXE,AVG_EYYE,AVG_EXYE,AVG_EZZE,AVG_EXZE,AVG_EYZE,
			AVG_EXXP,AVG_EYYP,AVG_EXYP,AVG_EZZP,AVG_EXZP,AVG_EYZP,
			AVG_EXX,AVG_EYY,AVG_EXY,AVG_EZZ,AVG_EXZ,AVG_EYZ,
			KINE_ENERGY,WORK_ENERGY,STRAIN_ENERGY,PLAS_ENERGY,
			AVG_VELX,AVG_VELY,AVG_VELZ,AVG_TEMP,
			AVG_DISPX,AVG_DISPY,AVG_DISPZ,WTFRACT_CONC,STEP_NUMBER,CPU_TIME,ELAPSED_TIME,
			GRID_ALPHA,INTERFACE_ENERGY,HISTORY_VARIABLE,TOT_FCONX,TOT_FCONY,TOT_FCONZ,
			HEAT_ENERGY,GRID_KINE_ENERGY,TOT_REACTX,TOT_REACTY,TOT_REACTZ,
            ENTROPY_ENERGY,INTERNAL_ENERGY,HELMHOLZ_ENERGY,PARTICLE_ALPHA,
			AVG_FXX,AVG_FXY,AVG_FXZ,AVG_FYX,AVG_FYY,AVG_FYZ,AVG_FZX,AVG_FZY,AVG_FZZ };

class GlobalQuantity
{
    public:
        // constructors and destructors
        GlobalQuantity();
		GlobalQuantity(char *,int);
    
        // methods
		GlobalQuantity *AppendName(char *);
		GlobalQuantity *AppendColor(char *);
		GlobalQuantity *AppendQuantity(vector<double> &);
	
		// accessors
		bool IncludeThisMaterial(int);
		GlobalQuantity *GetNextGlobal(void);
		void SetNextGlobal(GlobalQuantity *);
		bool IsSameQuantity(int,int,int);
		
		// class methods
		static int DecodeGlobalQuantity(char *,int *);
	
	private:
		int whichMat;
		int quantity;
		int subcode;
		int colorID;
		GlobalQuantity *nextGlobal;
		char *name;
};

extern GlobalQuantity *firstGlobal;

#endif
