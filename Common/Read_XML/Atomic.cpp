/********************************************************************************
	Atomic.hpp
	NairnMPMFEA
 
 	Created by John Nairn on 8/7/18.
 	Copyright (c) 2018 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/Atomic.hpp"
#include "Read_XML/Expression.hpp"

static char opChar[] = {' ',' ',' ',')','+','-','*','/','^','(',','};
static const char *fxnName[] = {"sin","cos","tan","asin","acos","atan","sinh","cosh",
						"tanh","log","log10","abs","int","sqrt","sign","exp","rand",
						"erf","erfc","cosramp","ramp","box","sinbox","sgn","tri","cdfinv"};

#pragma mark CONTRUCTORS AND DESTRUCTORS

// create with an operator
Atomic::Atomic(char opChar)
{
	switch(opChar)
	{	case '+':
			code = OP_PLUS;
			break;
		case '-':
			code = OP_MINUS;
			break;
		case '*':
			code = OP_TIMES;
			break;
		case '/':
			code = OP_DIV;
			break;
		case '^':
			code = OP_POW;
			break;
		case '(':
			code = OP_OPEN_GROUP;
			break;
		case ')':
			code = OP_CLOSE_GROUP;
			break;
		case ',':
			code = OP_DIVIDE_GROUP;
			break;
		default:
			break;
	}
	varname = NULL;
	nextAtom = NULL;
}

// create with a number
Atomic::Atomic(double nval)
{
	code = ATOM_NUMBER;
	value = nval;
	varname = NULL;
	nextAtom = NULL;
	readOnlyVarName = false;
}

// create with a valid variable name or function name
// Invalid function will set fxnCode to -1
Atomic::Atomic(char *subString,int sCode)
{
	code = sCode;
	readOnlyVarName = false;

	// variable is string, function reduced to code
	if(sCode==ATOM_VARIABLE)
	{	varname = subString;
	}
	else
	{	fxnCode = FindFunctionCode(subString);
		varname = NULL;
	}
	nextAtom = NULL;
}

// create with a valid functionc
Atomic::Atomic(int sCode,int validCode)
{
	if(sCode==FUNCTION_NAME)
	{	code = sCode;
		fxnCode = validCode;
	}
	else
		code = validCode;
	varname = NULL;
	nextAtom = NULL;
	readOnlyVarName = false;
}

// Destructor
Atomic::~Atomic()
{
	if(varname!=NULL && !readOnlyVarName)
		delete [] varname;
}

// Copy source atom to this atom, does not link them up
void Atomic::TransferAtom(Atomic *source)
{
	// use the same code
	code = source->GetCode();
	
	switch(code)
	{	case ATOM_NUMBER:
			value = source->GetValue();
			break;
		case ATOM_VARIABLE:
		{	char *sourceVarname = source->GetVarName();
			
			// if same then done
			if(varname!=NULL)
			{	if(strcmp(varname,sourceVarname)==0) break;
				delete [] varname;
			}
			
			// make copy, should never be needed during evaluations
			varname = new char[strlen(sourceVarname)+1];
			strcpy(varname,sourceVarname);
			break;
		}
		case FUNCTION_NAME:
			fxnCode = source->GetFunctionCode();
			break;
		default:
			// code is enough - an operature
			break;
	}
}

// Copy this atomis and return pointer
Atomic *Atomic::GetCopy(bool roVarName) const
{
	Atomic *copy;
	
	switch(code)
	{	case ATOM_NUMBER:
			copy = new Atomic(value);
			break;
		case ATOM_VARIABLE:
		{	if(roVarName)
			{	copy = new Atomic(varname,ATOM_VARIABLE);
				copy->SetReadOnlyVarName(true);
			}
			else
			{	// need to copy varname because contructor doesn't
				char *copyVarName = new char[strlen(varname)+1];
				strcpy(copyVarName,varname);
				copy = new Atomic(copyVarName,ATOM_VARIABLE);
			}
			break;
		}
		case FUNCTION_NAME:
			copy = new Atomic(FUNCTION_NAME,fxnCode);
			break;
		default:
			copy = new Atomic(OP_PLUS,code);
			break;
	}
	
	return copy;
}


#pragma mark ACCESSORS

// set and get next atomic
Atomic *Atomic::GetNextAtom(void) const { return nextAtom; }
void Atomic::SetNextAtom(Atomic *na) { nextAtom = na; }

// set and get code
void Atomic::SetCode(int ac) { code = ac; }
int Atomic::GetCode(void) const { return code; }

// set and get code
void Atomic::SetValue(double vv) { value = vv; }
double Atomic::GetValue(void) const { return value; }

// get function code
int Atomic::GetFunctionCode(void) const { return fxnCode; }
const char *Atomic::GetFunctionName(void) const { return fxnName[fxnCode]; }

// get variable name
char *Atomic::GetVarName(void) const { return varname; }

// first letter as number based on ASCII value from 'A' up
int Atomic::GetVarID(void) const { return (int)varname[0]-65; }

// make read only (to no delete when atom deleted)
void Atomic::SetReadOnlyVarName(bool setting) { readOnlyVarName = setting; }

// get code for function (or -1 if not value)
int Atomic::FindFunctionCode(char *fxn) const
{
	if(strcmp(fxn,"sin")==0) return SIN_FXN;
	if(strcmp(fxn,"cos")==0) return COS_FXN;
	if(strcmp(fxn,"tan")==0) return TAN_FXN;
	if(strcmp(fxn,"asin")==0) return ASIN_FXN;
	if(strcmp(fxn,"acos")==0) return ACOS_FXN;
	if(strcmp(fxn,"atan")==0) return ATAN_FXN;
	if(strcmp(fxn,"sinh")==0) return SINH_FXN;
	if(strcmp(fxn,"cosh")==0) return COSH_FXN;
	if(strcmp(fxn,"tanh")==0) return TANH_FXN;
	if(strcmp(fxn,"log")==0) return LOG_FXN;
	if(strcmp(fxn,"ln")==0) return LOG_FXN;
	if(strcmp(fxn,"log10")==0) return LOG10_FXN;
	if(strcmp(fxn,"abs")==0) return ABS_FXN;
	if(strcmp(fxn,"int")==0) return INT_FXN;
	if(strcmp(fxn,"sqrt")==0) return SQRT_FXN;
	if(strcmp(fxn,"sign")==0) return SIGN_FXN;
	if(strcmp(fxn,"exp")==0) return EXP_FXN;
	if(strcmp(fxn,"rand")==0) return RAND_FXN;
	if(strcmp(fxn,"erf")==0) return ERF_FXN;
	if(strcmp(fxn,"erfc")==0) return ERFC_FXN;
	
	if(strcmp(fxn,"cosramp")==0) return COSRAMP_FXN;
	if(strcmp(fxn,"ramp")==0) return RAMP_FXN;
	if(strcmp(fxn,"box")==0) return BOX_FXN;
	if(strcmp(fxn,"sinbox")==0) return SINBOX_FXN;
	if(strcmp(fxn,"sgn")==0) return SGN_FXN;
	if(strcmp(fxn,"tri") == 0) return TRI_FXN;
    if(strcmp(fxn,"cdfinv") == 0) return NORMALCDFINVERSE_FXN;
    if(strcmp(fxn,"mod")==0) return MOD_FXN;

	return -1;
}

// For dubbugging
void Atomic::Describe(void) const
{
	switch(code)
	{	case ATOM_NUMBER:
			cout << "... Number: " << value << endl;
			break;
		case ATOM_VARIABLE:
			cout << "... Variable: " << varname << endl;
			break;
		case FUNCTION_NAME:
			cout << "... Function " << fxnName[fxnCode] << "(Code: " << fxnCode << ")" << endl;
			break;
		default:
			cout << "... Operator " << opChar[code] << " (Code: " << code << ")" << endl;
			break;
	}
}
