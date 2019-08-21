/********************************************************************************
	Atomic.hpp
	NairnMPMFEA
 
 	Created by John Nairn on 8/7/18.
 	Copyright (c) 2018 John A. Nairn, All rights reserved.
	
	Dependencies
 		none
********************************************************************************/

#ifndef _ATOMIC_

#define _ATOMIC_

// Anything >= OP_PLUS cannot end an expression
enum { ATOM_NUMBER,ATOM_VARIABLE,FUNCTION_NAME,
		OP_CLOSE_GROUP,OP_PLUS,OP_MINUS,OP_TIMES,OP_DIV,OP_POW,OP_OPEN_GROUP,OP_DIVIDE_GROUP };

// supported functions (set string to for debugging)
enum { SIN_FXN=0,COS_FXN,TAN_FXN,ASIN_FXN,ACOS_FXN,ATAN_FXN,SINH_FXN,COSH_FXN,TANH_FXN,
	LOG_FXN,LOG10_FXN,ABS_FXN,INT_FXN,SQRT_FXN,SIGN_FXN,EXP_FXN,RAND_FXN,ERF_FXN,ERFC_FXN,
	COSRAMP_FXN,RAMP_FXN,BOX_FXN,SINBOX_FXN,SGN_FXN,TRI_FXN};

class Atomic
{
	public:
	
		// contructors
		Atomic(char);
		Atomic(double);
		Atomic(char *,int);
		Atomic(int,int);
		virtual ~Atomic();
		void TransferAtom(Atomic *);
		Atomic *GetCopy(bool) const;
	
		// methods
	
		// accessors
		Atomic *GetNextAtom(void) const;
		void SetNextAtom(Atomic *);
		int GetCode(void) const;
		void SetCode(int);
		double GetValue(void) const;
		void SetValue(double);
		char *GetVarName(void) const;
		int GetVarID(void) const;
		void SetReadOnlyVarName(bool);
		int GetFunctionCode(void) const;
		const char *GetFunctionName(void) const;
		int FindFunctionCode(char *) const;
		void Describe(void) const;
	
	private:
		int code;
		int fxnCode;
		double value;
		char *varname;
		Atomic *nextAtom;
		bool readOnlyVarName;
	
};

#endif
