/********************************************************************************
	Expression.hpp
	NairnMPMFEA
 
 	Created by John Nairn on 8/7/18.
 	Copyright (c) 2018 John A. Nairn, All rights reserved.
	
	Dependencies
 		none
********************************************************************************/

#ifndef _JAN_EXPRESSION_

#define _JAN_EXPRESSION_

#define MAX_EXPRESSIONS 12

class Atomic;

typedef struct {
	unsigned int location;
	unsigned int length;
} ExprRange;

class Expression
{
	public:
	
		// contructors
		Expression(void);
		Expression(const char *);
		virtual ~Expression();
	
		// const methods
		double EvaluateFunction(unordered_map<string, double>) const;
		double EvaluateFunction(double *) const;
		double EvaluateTokens(Atomic *,bool) const;
		void DoFunction(Atomic *) const;
		void OperateTerms(Atomic *firstOpAtom,int,int) const;
		Atomic *GetFunctionArg(Atomic *,Atomic *,double &) const;
		double XYZTValue(Vector *,double) const;
		double TValue(double) const;
	
		// non-const methods
		void TokenizeExpr(void);
		void AddToken(Atomic *);
		void AddSubstringToken(ExprRange,bool,bool);
		void ValidateCodeOrder(int,int) const;
	
		// accessors
		void SetString(const char *);
		const char *GetString(void) const;
		void Describe(void) const;
		void Describe(Atomic *) const;
	
		// class methods
		static Expression *CreateExpression(char *,const char *);
		static bool CreateFunction(char *&,int);
		static void DeleteFunction(int);
		static double FunctionValue(int,double,double,double,double,double,double);
	
		static Expression fxn[MAX_EXPRESSIONS];

	private:
		char *exprStr;
		Atomic *firstAtom;
		Atomic *currentAtom;
		bool exprHasGroups;
		int numAtoms;
};

#endif
