/********************************************************************************

mathexpr.cpp version 2.0

Copyright Yann OLLIVIER 1997-2000

This software may be freely distributed as is, including this whole notice.

It may be modified, but any modification should include this whole 
notice as is and a description of the changes made.

This software may not be sold.

This software or any modified version of it may befreely used in free 
programs. The program should include a copy of this whole notice.
If you want to use it in a program you sell, contact me

This software comes with absolutely no warranty.

Modifications and notes added by John Nairn, April 2007

Downloaded from http://www.yann-ollivier.org/mathlib/mathexpr.html but
	there is little documentation there. Here are some notes:

1. Define varibles as RVar objects such as

	RVar xvar("x",&x);
	RVar yvar("y",&y);

	Can have multiletter names and no checking on valid characters in names

2. Define function in string (can have any extra spaces)

	ROperation op(eqn,2,vararray2);
	
	eqn in char * to equation, number is number of possible variables
	in array of RVar objects at the end. The expression can have following functions
	
	sin(), cos(), tan(), exp(), log(), abs(), sqrt(), asin(), atan()
		acos(), ln()
	
	allows tg() for tan(), arcsin() for asin(), arccos() for acos(),
		arctg() for atan
	
	I added int(), sign(), Sinh(), Cosh(), and Tanh()
	I changed log() to be base 10

	angles are in radians
	
	exponentials in expression must use capital E as in 1.2E3
	I added l.c. e as well
	
	operators - + * / ^ #
	
	# is nth root such as 2#x for square root of x

3. Get value by loading variables and calling op.Val()

4. Other methods

	op.Expr(char varPrefix) - print the function
	op.HasError() - 1 or 0 if error or not
	op.ContainVar() - 1 or 0 if has variable or not
	op.NMembers() - number of members
	op.Substitute() - replace variable with expression
	
Changes by John Nairn, April 2007

1. Exponential notation can use lowercase e as in 1.0e3
2. Added int() which is actually the floor() function so int(-2.4) is -3
3. Added sign() which is 1 for positive and 0 for 0 or negative
4. Removed Diff() because will not use and new functions not implemented
5. Removed NthMember() (only needed by Diff())
6. Made log(x) base 10 while ln(x) is natural log
7. Added Sinh(), Cosh(), and Tanh() (upper case needed to avoid conflicts with sin, etc.)
8. Revised Expr() to remove embedded spaces and option to convert to
	format used in my parser (e.g. '#' before variables and log and log10 for log functions)

********************************************************************************/

#include "mathexpr.hpp"

#pragma mark USAGE UTILITILES

// if change this number, change number of NULLs in open[] to 1 less
#define MAX_FUNCTIONS 9

double xvalue,yvalue,zvalue,dvalue,thetavalue;
int numVars = 8;
PRVar vararray[8];
PROperation op=NULL;
PROperation opex[8]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

// create function of x,y,z (or R=x and Z=y in axisymmetric) or D and T for polar coordinates from x and y
// and dt (used in rigid velocities)
// deletes input string and replaced with formatted expression
// called should DeleteFunction() when done using it
bool CreateFunction(char *&eqn)
{
	// create variables
	vararray[0]=new RVar("x",&xvalue);
	vararray[1]=new RVar("y",&yvalue);
	vararray[3]=new RVar("Z",&yvalue);
	vararray[2]=new RVar("R",&xvalue);
	vararray[4]=new RVar("D",&dvalue);
	vararray[5]=new RVar("T",&thetavalue);
	vararray[6]=new RVar("z",&zvalue);
	vararray[7]=new RVar("dt",&zvalue);
	
	// create the function
	op=new ROperation(eqn,numVars,vararray);
	
	if(op->HasError())
	{	DeleteFunction();
		return false;
	}
	
	// return interpreted expression
	delete [] eqn;
	eqn=op->Expr('#');
	return true;
}

// create extra function of x,y,Z (or r=x and z=y) or D and T for polar coordinates from x and y
// only allowed if CreateFunction(expr) was already called and not deleted
// caller should DeleteFunction(i) when done using it
// can be 1, 2, or 3
bool CreateFunction(char *&eqn,int i)
{
	// redirect if first one
	if(op==NULL && i==1) return CreateFunction(eqn);
	
	// fails if base not there or index bad
	if(op==NULL || i>MAX_FUNCTIONS || i<2) return false;
	
	// create the function
	opex[i-2]=new ROperation(eqn,numVars,vararray);
	
	if(opex[i-2]->HasError())
	{	DeleteFunction(i);
		return false;
	}
	
	// return interpreted expression
	delete [] eqn;
	eqn=opex[i-2]->Expr('#');
	return true;
}
	
// delete objects allocated for the function
void DeleteFunction(void)
{
	if(op==NULL) return;
	int i;
	for(i=0;i<numVars;i++)
		delete vararray[i];
	delete op;
	op=NULL;
}

// delete extra function, but if i<0 delete all extra functions
void DeleteFunction(int i)
{	// if < 0 delete all
	if(i<0)
	{	for(int j=1;j<=MAX_FUNCTIONS;j++)
			DeleteFunction(j);
		return;
	}
	
	// if 1 delete root one with variables
	if(i==1)
		DeleteFunction();
	
	// delete extra ones
	else if(i>=2 && i<=MAX_FUNCTIONS && opex[i-2]!=NULL)
	{	delete opex[i-2];
		opex[i-2]=NULL;
	}
}

// get function value for given x and y and origin
double FunctionValue(int i,double x,double y,double z,double xorig,double yorig,double zorig)
{
	xvalue=x;
	yvalue=y;
	zvalue=z;
	
	// polar coordinates
	x-=xorig;
	y-=yorig;
	dvalue=sqrt(x*x+y*y);
	if(dvalue==0)
		thetavalue=0.;
	else if(x>=0.)
		thetavalue=asin(y/dvalue);
	else
		thetavalue=PI_CONSTANT-asin(y/dvalue);
	
	if(i==1)
		return op->Val();
	else
		return opex[i-2]->Val();
}

#pragma mark STRING UTILITIES

// return substring from charcter i1 to i2 (inclusive) of s in new string
char* MidStr(const char *s,int i1,int i2)
{
	if(i1<0 || i2>=(int)strlen(s) || i1>i2)
	{	char* cp = new char[1];
		cp[0] = '\0';
		return cp;
	}
	char *s1=new char[i2-i1+2];
	int i;
	for(i=i1;i<=i2;i++) s1[i-i1]=s[i];
	s1[i2-i1+1]=0;
	return s1;
}

// copy string and return pointer to new one
char* CopyStr(const char *s)
{	char *s1=new char[strlen(s)+1];
	char *s12=s1;
	const char *s2=s;
	while((*s12++=*s2++));
	return s1;
}

// Insert character at position n in string s
void InsStr(char*&s,int n,char c)
{	if(n<0||n>(int)strlen(s)) return;
	char *s1=new char[strlen(s)+2];
	int i;
	for(i=0;i<n;i++) s1[i]=s[i];
	s1[n]=c;
	for(i=n+1;s[i-1];i++) s1[i]=s[i-1];
	s1[i]=0;
	delete [] s;
	s=s1;
}

// see if two strings are the same
signed char EqStr(const char *s,const char *s2)
{	if(strlen(s)!=strlen(s2)) return 0;
	int i;
	for(i=0;s[i];i++)
	{	if(s[i]!=s2[i]) return 0;
	}
	return 1;
}

signed char CompStr(const char *s,int n,const char *s2)
{	if(n<0 || n>=(int)strlen(s) || n+(int)strlen(s2)>(int)strlen(s)) return 0;
	int i;
	for(i=0;s2[i];i++)
	{	if(s[i+n]!=s2[i])
			return 0;
	}
	return 1;
}

// Delete character n in string s and replace it with new string
void DelStr(char*&s,int n)
{	char *s1=new char[strlen(s)];
	int i;
	for(i=0;i<n;i++) s1[i]=s[i];
	for(i=n;s[i+1];i++) s1[i]=s[i+1];
	s1[i]=0;
	delete [] s;
	s=s1;
}

#pragma mark RVar Class

RVar::RVar(const RVar & rvarp)
{	if(this==&rvarp) return;
	pval=rvarp.pval;
	name=CopyStr(rvarp.name);
}

RVar::RVar(const char *namep,double *pvalp)
{	pval=pvalp;
	name=CopyStr(namep);
}

RVar::~RVar()
{	if(name!=NULL) delete[] name;
}

#pragma mark Function Class

RFunction::RFunction()
{
	type=-1;
	name=new char[1];
	name[0]=0;
	nvars=0;
	ppvar=NULL;
	pfuncval=NULL;
	op=ErrVal;
	buf=NULL;
}

RFunction::RFunction(double ((*pfuncvalp)(double)))
{
	type=0;
	pfuncval=pfuncvalp;
	name=new char[1];
	name[0]=0;
	nvars=1;
	ppvar=NULL;
	op=ErrVal;
	buf=NULL;
}

RFunction::RFunction(const RFunction& rfunc)
{
	if(this==&rfunc) return;
	type=rfunc.type;
	op=rfunc.op;
	pfuncval=rfunc.pfuncval;
	name=CopyStr(rfunc.name);
	nvars=rfunc.nvars;
	if(rfunc.ppvar!=NULL&&nvars)
    {	ppvar=new PRVar[nvars];
		int i;
		for(i=0;i<nvars;i++) ppvar[i]=rfunc.ppvar[i];
		buf=new double[nvars];
	}
	else
	{	ppvar=NULL;
		buf=NULL;
	}
}

RFunction::RFunction(const ROperation& opp,RVar* pvarp):op(opp)
{
	type=1;
	name=new char[1];
	name[0]=0;
	nvars=1;
	ppvar=new PRVar[1];
	ppvar[0]=pvarp;
	buf=new double[1];
}

RFunction::RFunction(const ROperation& opp, int nvarsp,RVar**ppvarp):op(opp)
{
	type=1;
	name=new char[1];
	name[0]=0;
	nvars=nvarsp;
	if(nvars)
    {	ppvar=new PRVar[nvars];
		int i;
		for(i=0;i<nvars;i++) ppvar[i]=ppvarp[i];
		buf=new double[nvars];
	}
	else
	{	ppvar=NULL;
		buf=NULL;
	}
}

RFunction::~RFunction()
{
	if(name!=NULL) delete [] name;
	if(ppvar!=NULL) delete [] ppvar;
	if(buf!=NULL) delete [] buf;
}

RFunction& RFunction::operator=(const RFunction& rfunc)
{
	if(this==&rfunc) return *this;
	type=rfunc.type;
	op=rfunc.op;
	pfuncval=rfunc.pfuncval;
	delete [] name;
	name=CopyStr(rfunc.name);
	if(ppvar!=NULL) delete [] ppvar;
	ppvar=NULL;
	if(buf!=NULL) delete [] buf;
	buf=NULL;
	nvars=rfunc.nvars;
	if(type==1&&nvars)
    {	ppvar=new PRVar[nvars];
		buf=new double[nvars];
		int i;
		for(i=0;i<nvars;i++) ppvar[i]=rfunc.ppvar[i];
	}
	return *this;
}

void RFunction::SetName(const char*s)
{	if(name!=NULL) delete [] name;
	name=CopyStr(s);
}

double RFunction::Val(double x) const
{
	if(type==-1||nvars>=2) return ErrVal;
	if(type==0) return (*pfuncval)(x);
	double xb=*(*ppvar)->pval,y;
	*(*ppvar)->pval=x;  // Warning : could cause trouble if this value is used in a parallel process
	y=op.Val();
	*(*ppvar)->pval=xb;
	return y;
}

double RFunction::Val(double*pv) const
{
	if(type==-1) return ErrVal;
	if(type==0) return (*pfuncval)(*pv);
	double y;
	int i;
	for(i=0;i<nvars;i++)
    {	buf[i]=*ppvar[i]->pval;
		// Warning : could cause trouble if this value is used in a parallel process
		*ppvar[i]->pval=pv[i];
	}
	y=op.Val();
	for(i=0;i<nvars;i++) *ppvar[i]->pval=buf[i];
	return y;
}

#pragma mark ROperation Class

ROperation::ROperation()
{	op=ErrOp;
	mmb1=NULL;
	mmb2=NULL;
	ValC=ErrVal;
	pvar=NULL;
	pvarval=NULL;
	pfunc=NULL;
	containfuncflag=0;
	pinstr=NULL;
	pvals=NULL;
	ppile=NULL;
	pfuncpile=NULL;
	BuildCode();
}

ROperation::~ROperation()
{	Destroy();
}

ROperation::ROperation(const ROperation&ROp)
{	op=ROp.op;
	pvar=ROp.pvar;
	pvarval=ROp.pvarval;
	ValC=ROp.ValC;
	pfunc=ROp.pfunc;
	containfuncflag=0;
	pinstr=NULL;
	pvals=NULL;
	ppile=NULL;
	pfuncpile=NULL;
	if(ROp.mmb1!=NULL) 
		mmb1=new ROperation(*(ROp.mmb1));
	else
		mmb1=NULL;
	if(ROp.mmb2!=NULL)
		mmb2=new ROperation(*(ROp.mmb2));
	else
		mmb2=NULL;
	BuildCode();
}

ROperation::ROperation(double x)
{	if(x==ErrVal)
	{	op=ErrOp;
		mmb1=NULL;
		mmb2=NULL;
		ValC=ErrVal;
	}
	else if(x>=0)
	{	op=Num;
		mmb1=NULL;
		mmb2=NULL;
		ValC=x;
	}
	else
	{	op=Opp;
		mmb1=NULL;
		mmb2=new ROperation(-x);
		ValC=ErrVal;
	}
	pvar=NULL;
	pvarval=NULL;
	pfunc=NULL;
	containfuncflag=0;
	pinstr=NULL;
	pvals=NULL;
	ppile=NULL;
	pfuncpile=NULL;
	BuildCode();
}

ROperation::ROperation(const RVar&varp)
{	op=Var;
	mmb1=NULL;
	mmb2=NULL;
	ValC=ErrVal;
	pvar=&varp;
	pvarval=varp.pval;
	containfuncflag=0;
	pfunc=NULL;
	pinstr=NULL;
	pvals=NULL;
	ppile=NULL;
	pfuncpile=NULL;
	BuildCode();
}

ROperation& ROperation::operator=(const ROperation& ROp)
{	if(this==&ROp) return *this;
	Destroy();
	op=ROp.op;
	pvar=ROp.pvar;
	pvarval=ROp.pvarval;
	ValC=ROp.ValC;
	pfunc=ROp.pfunc;
	containfuncflag=0;
	pinstr=NULL;
	pvals=NULL;
	ppile=NULL;
	pfuncpile=NULL;
	if(ROp.mmb1!=NULL)
		mmb1=new ROperation(*(ROp.mmb1));
	else
		mmb1=NULL;
	if(ROp.mmb2!=NULL)
		mmb2=new ROperation(*(ROp.mmb2));
	else
		mmb2=NULL;
	BuildCode();
	return *this;
}

int operator==(const ROperation& op,const double v)
{	return(op.op==Num&&op.ValC==v);
}

int operator==(const ROperation& op1,const ROperation& op2)
{
	if(op1.op!=op2.op) return 0;
	if(op1.op==Var) return(*(op1.pvar)==*(op2.pvar));
	if(op1.op==Fun) return(op1.pfunc==op2.pfunc); // *op1.pfunc==*op2.pfunc could imply infinite loops in cases of self-dependence
	if(op1.op==Num) return(op1.ValC==op2.ValC);
	if(op1.mmb1==NULL&&op2.mmb1!=NULL) return 0;
	if(op1.mmb2==NULL&&op2.mmb2!=NULL) return 0;
	if(op2.mmb1==NULL&&op1.mmb1!=NULL) return 0;
	if(op2.mmb2==NULL&&op1.mmb2!=NULL) return 0;
	return(((op1.mmb1==NULL&&op2.mmb1==NULL)||(*(op1.mmb1)==*(op2.mmb1)))&&
		((op1.mmb2==NULL&&op2.mmb2==NULL)||(*(op1.mmb2)==*(op2.mmb2))));
}

int operator!=(const ROperation& op1,const ROperation& op2)
{
	if(op1.op!=op2.op)return 1;
	if(op1.op==Var)return(op1.pvar!=op2.pvar);
	if(op1.op==Fun)return(!(op1.pfunc==op2.pfunc)); // *op1.pfunc==*op2.pfunc could imply infinite loops in cases of self-dependence
	if(op1.op==Num)return(op1.ValC!=op2.ValC);
	if(op1.mmb1==NULL&&op2.mmb1!=NULL)return 1;
	if(op1.mmb2==NULL&&op2.mmb2!=NULL)return 1;
	if(op2.mmb1==NULL&&op1.mmb1!=NULL)return 1;
	if(op2.mmb2==NULL&&op1.mmb2!=NULL)return 1;
	return(((op1.mmb1!=NULL||op2.mmb1!=NULL)&&(*(op1.mmb1)!=*(op2.mmb1)))||
		((op1.mmb2!=NULL||op2.mmb2!=NULL)&&(*(op1.mmb2)!=*(op2.mmb2))));
}

ROperation ROperation::operator+() const
{	return *this;
}

ROperation ROperation::operator-() const
{	if(op==Num)return -ValC;
	ROperation resultat;
	if(op==Opp)
		resultat=*mmb2;
	else
	{	resultat.op=Opp;
		resultat.mmb2=new ROperation(*this);
	}
	return resultat;
}

ROperation operator,(const ROperation& op1,const ROperation& op2)
{	ROperation resultat;
	resultat.op=Juxt;
	resultat.mmb1=new ROperation(op1);
	resultat.mmb2=new ROperation(op2);
	return resultat;
}

ROperation operator+(const ROperation& op1,const ROperation& op2)
{
	if(op1.op==Num&&op2.op==Num) return op1.ValC+op2.ValC;
	if(op1==0.) return op2;
	if(op2==0.) return op1;
	if(op1.op==Opp) return op2-*(op1.mmb2);
	if(op2.op==Opp) return op1-*(op2.mmb2);
	ROperation resultat;
	resultat.op=Add;
	resultat.mmb1=new ROperation(op1);
	resultat.mmb2=new ROperation(op2);
	return resultat;
}

ROperation operator-(const ROperation& op1,const ROperation& op2)
{
	if(op1.op==Num&&op2.op==Num) return op1.ValC-op2.ValC;
	if(op1==0.) return -op2;
	if(op2==0.) return op1;
	if(op1.op==Opp) return -(op2+*(op1.mmb2));
	if(op2.op==Opp) return op1+*(op2.mmb2);
	ROperation resultat;
	resultat.op=Sub;
	resultat.mmb1=new ROperation(op1);
	resultat.mmb2=new ROperation(op2);
	return resultat;
}

ROperation operator*(const ROperation& op1,const ROperation& op2)
{
	if(op1.op==Num&&op2.op==Num) return op1.ValC*op2.ValC;
	if(op1==0.||op2==0.) return 0.;
	if(op1==1.) return op2;
	if(op2==1.) return op1;
	if(op1.op==Opp) return -(*(op1.mmb2)*op2);
	if(op2.op==Opp) return -(op1**(op2.mmb2));
	ROperation resultat;
	resultat.op=Mult;
	resultat.mmb1=new ROperation(op1);
	resultat.mmb2=new ROperation(op2);
	return resultat;
}

ROperation operator/(const ROperation& op1,const ROperation& op2)
{
	if(op1.op==Num&&op2.op==Num) return (op2.ValC?op1.ValC/op2.ValC:ErrVal);
	if(op1==0.0) return 0.;
	if(op2==1.) return op1;
	if(op2==0.) return ErrVal;
	if(op1.op==Opp) return -(*(op1.mmb2)/op2);
	if(op2.op==Opp) return -(op1/(*(op2.mmb2)));
	ROperation resultat;
	resultat.op=Div;
	resultat.mmb1=new ROperation(op1);
	resultat.mmb2=new ROperation(op2);
	return resultat;
}

ROperation operator^(const ROperation& op1,const ROperation& op2)
{
	if(op1==0.) return 0.;
	if(op2==0.) return 1.;
	if(op2==1.) return op1;
	ROperation resultat;
	resultat.op=Pow;
	resultat.mmb1=new ROperation(op1);
	resultat.mmb2=new ROperation(op2);
	return resultat;
}               
                
ROperation sqrt(const ROperation& op)
{	ROperation rop;
	rop.op=Sqrt;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation abs(const ROperation& op)
{	ROperation rop;
	rop.op=Abs;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation floor(const ROperation& op)
{	ROperation rop;
	rop.op=IntFun;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation sin(const ROperation& op)
{	ROperation rop;
	rop.op=Sin;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation cos(const ROperation& op)
{	ROperation rop;
	rop.op=Cos;
	rop.mmb2=new ROperation(op);
	return rop;

}
ROperation tan(const ROperation& op)
{	ROperation rop;
	rop.op=Tg;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation log(const ROperation& op)
{	ROperation rop;
	rop.op=Ln;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation log10(const ROperation& op)
{	ROperation rop;
	rop.op=Log;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation exp(const ROperation& op)
{	ROperation rop;
	rop.op=Exp;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation acos(const ROperation& op)
{	ROperation rop;
	rop.op=Acos;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation asin(const ROperation& op)
{	ROperation rop;
	rop.op=Asin;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation atan(const ROperation& op)
{	ROperation rop;
	rop.op=Atan;
	rop.mmb2=new ROperation(op);
	return rop;
}

ROperation ApplyOperator(int n,ROperation**pops,ROperation (*func)(const ROperation&,const ROperation&))
{
	if(n<=0) return ErrVal;
	if(n==1) return *pops[0];
	if(n==2) return (*func)(*pops[0],*pops[1]);
	return (*func)(*pops[0],ApplyOperator(n-1,pops+1,func));
}

ROperation RFunction::operator()(const ROperation& op)
{
	ROperation op2;
	op2.op=Fun;
	op2.pfunc=this;
	op2.mmb2=new ROperation(op);
	return op2;
}

//Auxiliary string functions

// suppress spaces
void SupprSpaces(char*&s)//Deletes the old string
{	int i;
	for(i=0;s[i];i++)
	{	if(s[i]==' '||s[i]=='\t'||s[i]=='\n')
			DelStr(s,i--);
	}
}

// numeric or decimal point
signed char IsNumeric(char c)
{	if(c!='0'&&c!='1'&&c!='2'&&c!='3'&&c!='4'
		&&c!='5'&&c!='6'&&c!='7'&&c!='8'&&c!='9'&&c!='.') return 0;
	return 1;
}

// is strain all numeric
signed char IsTNumeric(char *s)
{	int i;
	for(i=0;i<(int)strlen(s);i++)
	{	if(!IsNumeric(s[i]))
			return 0;
	}
	return 1;
}

// Searches from n+1 for ')' matching an open bracket
int SearchCorOpenbracket(char*s,int n)
{
	if(n>=(int)strlen(s)-1) return -1;
	int i,c=1;
	for(i=n+1;s[i];i++)
	{	if(s[i]=='(')
			c++;
		else if(s[i]==')')
			c--;
		if(!c) return i;
	}
	return -1;
}

int SearchCorClosebracket(char*s,int n)  //Searchs the corresponding bracket of a closing bracket
{	if(n<1) return -1;
	int i,c=1;
	for(i=n-1;i>=0;i--)
	{	if(s[i]==')')
			c++;
		else if(s[i]=='(')
			c--;
		if(!c) return i;
	};
	return -1;
}

int SearchOperator(char *s,ROperator op)
{
	char opc;
	switch(op)
	{	case ErrOp:
		case Num:
		case Var:
			return -1;
		case Juxt:
			opc=',';
			break;
		case Add:
			opc='+';
			break;
		case Sub:
			opc='-';
			break;
		case Mult:
			opc='*';
			break;
		case Div:
			opc='/';
			break;
		case Pow:
			opc='^';
			break;
		case NthRoot:
			opc='#';
			break;
		case E10:
			opc='E';
			break;        
		case e10:
			opc='e';
			break;        
		default:
			return -1;
	}
	int i;
	for(i=(int)strlen(s)-1;i>=0;i--)
	{	if(s[i]==opc && (op!=Sub || (i && s[i-1]==')'))) return i;
		if(s[i]==')')
		{	i=SearchCorClosebracket(s,i);
			if(i==-1) return -1;
		}
	}
	return -1;
}

// Remove leading and trailing blanks and/or leading and trailing brackets
void SimplifyStr(char*&s)
{
	if(!strlen(s)) return;
	
	char *s1=s,*s2=s+strlen(s);
	signed char ind=0;
	
	// check for "(text)"
	if(s1[0]=='('&&SearchCorOpenbracket(s1,0)==s2-s1-1)
	{	s1++;
		s2--;
		ind=1;
	}
	
	// check for empty string
	if(s1==s2)
	{	delete [] s;
		s=new char[1]; // ISO C++ forbids initialization in array new
		s[0]=0;
		return;
	}
	
	// skip leading blanks
	if(s1[0]==' ')
	{	ind=1;
		while(s1[0]==' '&&s1<s2) s1++;
	}
	
	// check for empty string
	if(s1==s2)
	{	delete [] s;
		s=new char[1]; // ISO C++ forbids initialization in array new
		s[0]=0;
		return;
	}
	
	// skip trailing balanks
	if(*(s2-1)==' ')
	{	ind=1;
		while(s2>s1&&*(s2-1)==' ') s2--;
	}
	*s2=0;
	
	// replace with remaining string
	s1=CopyStr(s1);
	delete [] s;
	s=s1;
	
	// call again might be needed
	if(ind) SimplifyStr(s);
}

int max(int a, int b) {return (a>b?a:b); }

// see if string is a varaible (which may be mutlicharacters)
// finds the longest one if start matches another (e.g. x and xloc finds xloc)
// returns length of the matched variable
int IsVar(const char *s,int n,int nvar,PRVar *ppvar)
{
	if(n<0||n>(int)strlen(s)) return 0;
	int i;
	int l=0;
	for(i=0;i<nvar;i++)
	{	if(CompStr(s,n,(*(ppvar+i))->name))
			l=max(l,strlen((*(ppvar+i))->name));
	}
	return l;
}

int IsFunction(const char*s,int n)
{
	if(CompStr(s,n,"sin") || CompStr(s,n,"cos") || CompStr(s,n,"exp")
		|| CompStr(s,n,"tan") || CompStr(s,n,"log") || CompStr(s,n,"atg")
		|| CompStr(s,n,"abs") || CompStr(s,n,"int")) return 3;
	if(CompStr(s,n,"tg") || CompStr(s,n,"ln")) return 2;
	if(CompStr(s,n,"sqrt") || CompStr(s,n,"asin") || CompStr(s,n,"atan") ||
		CompStr(s,n,"acos") || CompStr(s,n,"sign") || CompStr(s,n,"Sinh") ||
		CompStr(s,n,"Cosh") || CompStr(s,n,"Tanh")) return 4;
	if(CompStr(s,n,"arcsin") || CompStr(s,n,"arccos") || CompStr(s,n,"arctan")) return 6;
	if(CompStr(s,n,"arctg")) return 5;
	return 0;
}

int IsFunction(const char*s,int n,int nfunc,PRFunction*ppfunc)
	//Not recognized if a user-defined function is eg "sine" ie begins like
	//a standard function
	//IF PATCHED TO DO OTHERWISE, SHOULD BE PATCHED TOGETHER WITH THE
	//PARSER BELOW which treats standard functions before user-defined ones
{
	int l=IsFunction(s,n);
	if(l) return l;
	int i;
	l=0;
	for(i=0;i<nfunc;i++)
	{	if(CompStr(s,n,ppfunc[i]->name))
			l=max(l,strlen(ppfunc[i]->name));
	}
	return l;
}

signed char IsFunction(ROperator op)
{	return (op==Exp || op==Abs || op==Sin || op==Cos || op==Tg || op==Ln ||
		op==Atan || op==Asin || op==Acos || op==Atan || op==Sqrt || op==Opp
		|| op==IntFun || op==Sign || op==Log || op==Sinh || op==Cosh || op==Tanh);
}

void IsolateVars(char*&s,int nvar,PRVar*ppvar,int nfunc,PRFunction*ppfunc)//Deletes the old string
{
	int i,j;
	i=0;
	for(i=0;s[i];i++)
    {	if(s[i]=='(')
		{	i=SearchCorOpenbracket(s,i);
			if(i==-1)return;
			continue;
		}
		if(((j=IsVar(s,i,nvar,ppvar))>IsFunction(s,i,nfunc,ppfunc)) ||
				((CompStr(s,i,"pi") || CompStr(s,i,"PI") || CompStr(s,i,"Pi")) && (j=2)))
		{	InsStr(s,i,'(');
			InsStr(s,i+j+1,')');
			i+=j+1;
			continue;
		}
		if(IsFunction(s,i,nfunc,ppfunc))
		{	i+=IsFunction(s,i,nfunc,ppfunc)-1;
			if(!s[i]) return;
			continue;
		}
	}
}

void IsolateNumbers(char*&s,int nvar,RVar**ppvar,int nfunc,RFunction**ppfunc)//Deletes the old string
{
	int i,i2=0,ind=0,t1,t2;
	for(i=0;s[i];i++)
    {	if(ind && !IsNumeric(s[i]))
		{	ind=0;
			InsStr(s,i2,'(');
			i++;
			InsStr(s,i,')');
			continue;
		}
		t1=IsVar(s,i,nvar,ppvar);
		t2=IsFunction(s,i,nfunc,ppfunc);
		if(t1 || t2)
		{	i+=max(t1,t2)-1;
			continue;
		}
		if(s[i]=='(')
		{	i=SearchCorOpenbracket(s,i);
			if(i==-1) return;
			continue;
		}
		if(!ind && IsNumeric(s[i]))
		{	i2=i;
			ind=1;
		}
	}
	if(ind) InsStr(s,i2,'(');
	i++;
	InsStr(s,i,')');
}

ROperation::ROperation(char*sp,int nvar,PRVar*ppvarp,int nfuncp,PRFunction*ppfuncp)
{
	ValC=ErrVal;
	mmb1=NULL;
	mmb2=NULL;
	pvar=NULL;
	op=ErrOp;
	pvarval=NULL;
	containfuncflag=0;
	pfunc=NULL;
	pinstr=NULL;
	pvals=NULL;
	ppile=NULL;
	pfuncpile=NULL;
	int i,j,k,l;
	signed char flag=1;
	char *s=CopyStr(sp),*s1=NULL,*s2=NULL;
	
	SimplifyStr(s);
	if(!s[0] || !strcmp(s,"Error")) goto fin;
	while(s[0]==':' || s[0]==';')
	{	s1=CopyStr(s+1);
		delete [] s;
		s=s1;
		s1=NULL;
		SimplifyStr(s);
		if(!s[0] || !strcmp(s,"Error")) goto fin;
	}
	
	// if s is decimenal number, store it and go to fin
	if(IsTNumeric(s))
	{	op=Num;
		ValC=atof(s);
		mmb1=NULL;
		mmb2=NULL;
		goto fin;
	}
	
	// if matches the word pi, log it and go to fin
	if(EqStr(s,"pi") || EqStr(s,"PI") || EqStr(s,"Pi"))
    {	op=Num;
		ValC=3.141592653589793238462643383279L;
		mmb1=NULL;
		mmb2=NULL;
		goto fin;
	}
	
	// if start of s matches a variable, find it and go to fin
	if(IsFunction(s,0,nfuncp,ppfuncp)<IsVar(s,0,nvar,ppvarp))
	{	for(i=0;i<nvar;i++)
		{	if(EqStr(s,(*(ppvarp+i))->name))
			{	pvar=ppvarp[i];
				pvarval=pvar->pval;
				op=Var;
				mmb1=NULL;
				mmb2=NULL;
				goto fin;
			}
		}
	}
	
	// mark functions
	for(k=0;s[k];k++)
	{	// skip bracketed section
    	if(s[k]=='(')
		{	k=SearchCorOpenbracket(s,k);
			if(k==-1) break;
			continue;
		}
		
		// look for function name longer than a variable name
		// insert ';' after the function name if brackets sin;(x)
		// but insert ":" is not brackets such as sin:x
		if((l=IsFunction(s,k,nfuncp,ppfuncp)) && l>=IsVar(s,k,nvar,ppvarp))
		{	i=k+l;
			while(s[i]==' ') i++;
			if(s[i]=='(')
			{	j=SearchCorOpenbracket(s,i);
				if(j!=-1)
				{	InsStr(s,i,';');
					k=j+1;
				}
				else
					break;
			}
			else if(s[i]!=':' && s[i]!=';')
			{	InsStr(s,i,':');
				k=i;
			}
		}
	}

	// locate and bracket numbers and variable if needed
	IsolateNumbers(s,nvar,ppvarp,nfuncp,ppfuncp);
	if(nvar) IsolateVars(s,nvar,ppvarp,nfuncp,ppfuncp);
	
	// remove empty space
	SupprSpaces(s);
	
	// search possible operators
	i=SearchOperator(s,Juxt);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=Juxt;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	i=SearchOperator(s,Add);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=Add;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	i=SearchOperator(s,Sub);
	if(i!=-1)
    {	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=Sub;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	if(s[0]=='-')
	{	s2=MidStr(s,1,strlen(s)-1);
		op=Opp;
		mmb1=NULL;
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	for(i=0;s[i];i++)
	{	if(s[i]=='(')
		{	i=SearchCorOpenbracket(s,i);
			if(i==-1) break;
			continue;
		}
		if(IsFunction(s,i,nfuncp,ppfuncp))
		{	k=i+IsFunction(s,i,nfuncp,ppfuncp);
			while(s[k]==' ') k++;
			if(s[k]==';')
			{	//	s=DelStr(s,k);
				j=k;
				while(s[j]!='(') j++;
				j=SearchCorOpenbracket(s,j);
				if(j!=-1)
				{	InsStr(s,j,')');
					InsStr(s,i,'(');
					i=j+2;
				}
			}
			else if(s[k]==':')
			{	//	s=DelStr(s,k);
				for(j=k;s[j];j++)
				{	if(s[j]=='(')
					{	j=SearchCorOpenbracket(s,j);
						break;
					}
				}
				if(j==-1) break;
				for(j++;s[j];j++)
				{	if(s[j]=='(')
					{	j=SearchCorOpenbracket(s,j);
						if(j==-1)
						{	flag=0;
							break;
						}
						continue;
					}
					if(IsFunction(s,j,nfuncp,ppfuncp)) break;
				}
				if(flag==0)
				{	flag=1;
					break;
				}
				while(j>i&&s[j-1]!=')') j--;
				if(j<=i+1) break;
				InsStr(s,i,'(');
				InsStr(s,j+1,')');
				i=j+1;
			}
		}
	}
	for(i=0;s[i]&&s[i+1];i++)
	{	if(s[i]==')'&&s[i+1]=='(')
			InsStr(s,++i,'*');
	}
	
	// locate functions
	if(s[0]=='(' && SearchCorOpenbracket(s,0)==(int)strlen(s)-1)
	{	if(CompStr(s,1,"exp"))
		{	op=Exp;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"abs"))
		{	op=Abs;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"int"))
		{	op=IntFun;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"sign"))
		{	op=Sign;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"sin"))
		{	op=Sin;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"cos"))
		{	op=Cos;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"tan"))
		{	op=Tg;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"log"))
		{	op=Log;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"atg"))
		{	op=Atan;
			s2=MidStr(s,4,strlen(s)-2);
		}
		else if(CompStr(s,1,"tg"))
		{	op=Tg;
			s2=MidStr(s,3,strlen(s)-2);
		}
		else if(CompStr(s,1,"ln"))
		{	op=Ln;
			s2=MidStr(s,3,strlen(s)-2);
		}
		else if(CompStr(s,1,"asin"))
		{	op=Asin;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"acos"))
		{	op=Acos;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"atan"))
		{	op=Atan;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"sqrt"))
		{	op=Sqrt;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"Sinh"))
		{	op=Sinh;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"Cosh"))
		{	op=Cosh;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"Tanh"))
		{	op=Tanh;
			s2=MidStr(s,5,strlen(s)-2);
		}
		else if(CompStr(s,1,"arcsin"))
		{	op=Asin;
			s2=MidStr(s,7,strlen(s)-2);
		}
		else if(CompStr(s,1,"arccos"))
		{	op=Acos;
			s2=MidStr(s,7,strlen(s)-2);
		}
		else if(CompStr(s,1,"arctan"))
		{	op=Atan;
			s2=MidStr(s,7,strlen(s)-2);
		}
		else if(CompStr(s,1,"arctg"))
		{	op=Atan;
			s2=MidStr(s,6,strlen(s)-2);
		}
		else
		{	// look for user function
			for(i=-1,k=0,j=0;j<nfuncp;j++)
			{	if(CompStr(s,1,ppfuncp[j]->name)&&k<(int)strlen(ppfuncp[j]->name))
				{	k=strlen(ppfuncp[j]->name);
					i=j;
				}
			}
			if(i>-1)
			{	op=Fun;
				s2=MidStr(s,strlen(ppfuncp[i]->name)+1,strlen(s)-2);
				pfunc=ppfuncp[i];
			}
		}
		mmb1=NULL;
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		if(op==Fun)
		{	if(mmb2->NMembers()!=pfunc->nvars)
			{	op=ErrOp;
				mmb1=NULL;
				mmb2=NULL;
				goto fin;
			}
		}
		goto fin;
	}
	
	i=SearchOperator(s,Mult);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=Mult;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	i=SearchOperator(s,Div);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=Div;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	i=SearchOperator(s,Pow);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=Pow;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	i=SearchOperator(s,NthRoot);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		if(i==0||s[i-1]!=')')
		{	op=Sqrt;
			mmb1=NULL;
		}
		else
		{	op=NthRoot;
			mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		}
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	i=SearchOperator(s,E10);
	if(i==-1) i=SearchOperator(s,e10);
	if(i!=-1)
	{	s1=MidStr(s,0,i-1);
		s2=MidStr(s,i+1,strlen(s)-1);
		op=E10;
		mmb1=new ROperation(s1,nvar,ppvarp,nfuncp,ppfuncp);
		mmb2=new ROperation(s2,nvar,ppvarp,nfuncp,ppfuncp);
		goto fin;
	}
	op=ErrOp;
	mmb1=NULL;
	mmb2=NULL;
	
fin:
	BuildCode();
	delete [] s;
	if(s1!=NULL) delete [] s1;
	if(s2!=NULL) delete []s2;
}

void ROperation::Destroy()
{
	if(mmb1!=NULL && mmb2!=NULL && mmb1!=mmb2)
	{	delete mmb1;
		delete mmb2;
		mmb1=NULL;
		mmb2=NULL;
	}
	else if(mmb1!=NULL)
	{	delete mmb1;
		mmb1=NULL;
	}
	else if(mmb2!=NULL)
	{	delete mmb2;
		mmb2=NULL;
	}
	if(pinstr!=NULL)
	{	delete [] pinstr;
		pinstr=NULL;
	}
	if(pvals!=NULL)
	{	if(op==ErrOp || op==Num) delete pvals[0];
		delete [] pvals;
		pvals=NULL;
	}
	if(ppile!=NULL)
	{	delete [] ppile;
		ppile=NULL;
	}
	if(pfuncpile!=NULL)
	{	delete [] pfuncpile;
		pfuncpile=NULL;
	}
}

int operator==(const RVar& var1,const RVar& var2)
{	return (var1.pval==var2.pval&&EqStr(var1.name,var2.name));
}

int operator==(const RFunction& f1,const RFunction& f2)
{	if(f1.type!=f2.type) return 0;
	if(f1.type==-1) return 1; // Nonfunction==nonfunction
	if(f1.type==0) return (f1.pfuncval==f2.pfuncval&&EqStr(f1.name,f2.name));
	if(f1.op!=f2.op) return 0;
	if(!EqStr(f1.name,f2.name)) return 0;
	if(f1.nvars!=f2.nvars) return 0;
	int i;
	for(i=0;i<f1.nvars;i++)
	{	if(!(*f1.ppvar[i]==*f2.ppvar[i])) return 0;
	}
	return 1;
}

signed char ROperation::ContainVar(const RVar& varp) const
{	if(op==Var)
	{	if(EqStr(pvar->name,varp.name)&&pvar->pval==varp.pval)
			return 1;
		else
			return 0;
	}
	if(mmb1!=NULL&&mmb1->ContainVar(varp)) return 1;
	if(mmb2!=NULL&&mmb2->ContainVar(varp)) return 1;
	return 0;
}

signed char ROperation::ContainFuncNoRec(const RFunction& func) const // No recursive test on subfunctions
{	if(op==Fun)
	{	if(*pfunc==func)
			return 1;
		else
			return 0;
	}
	if(mmb1!=NULL&&mmb1->ContainFuncNoRec(func)) return 1;
	if(mmb2!=NULL&&mmb2->ContainFuncNoRec(func)) return 1;
	return 0;
}

signed char ROperation::ContainFunc(const RFunction& func) const // Recursive test on subfunctions
{
	if(containfuncflag) return 0;
	if(op==Fun&&*pfunc==func) return 1;
	containfuncflag=1;
	if(op==Fun)
	{	if(pfunc->op.ContainFunc(func))
		{	containfuncflag=0;
			return 1;
		}
	}
	if(mmb1!=NULL&&mmb1->ContainFunc(func))
	{	containfuncflag=0;
		return 1;
	}
	if(mmb2!=NULL&&mmb2->ContainFunc(func))
	{	containfuncflag=0;
		return 1;
	}
	containfuncflag=0;
	return 0;
}

signed char ROperation::HasError(const ROperation *pop) const
{
	if(op==ErrOp) return 1;
	if(op==Fun && pfunc->type==1 && pfunc->op==*(pop==NULL?this:pop)) return 1;
	if(op==Fun && pfunc->type==1 && pfunc->op.HasError((pop==NULL?this:pop))) return 1;
	if(mmb1!=NULL && mmb1->HasError((pop==NULL?this:pop))) return 1;
	if(mmb2!=NULL && mmb2->HasError((pop==NULL?this:pop))) return 1;
	if(op==Fun&&pfunc->type==-1) return 1;
	return 0;
}

int ROperation::NMembers() const //Number of members for an operation like a,b,c...
{
	if(op==Fun) return(pfunc->type==1?pfunc->op.NMembers():pfunc->type==0?1:0);
	if(op!=Juxt)
		return 1;
	else if(mmb2==NULL)
		return 0;
	else
		return 1+mmb2->NMembers();
}

ROperation ROperation::Substitute(const RVar& var,const ROperation& rop) const // Replaces variable var with expression rop
{
	if(!ContainVar(var)) return *this;
	if(op==Var) return rop;
	ROperation r;
	r.op=op;
	r.pvar=pvar;
	r.pvarval=pvarval;
	r.ValC=ValC;
	r.pfunc=pfunc;
	if(mmb1!=NULL)
		r.mmb1=new ROperation(mmb1->Substitute(var,rop));
	else
		r.mmb1=NULL;
	if(mmb2!=NULL)
		r.mmb2=new ROperation(mmb2->Substitute(var,rop));
	else
		r.mmb2=NULL;
	return r;
}

char* ValToStr(double x)
{
	char *s=new char[30];
	if(x==(double)3.141592653589793238462643383279L)
		sprintf(s,"pi");
	else
		sprintf(s,"%.16G",x);
	return s;
}

// Format expression without any embedded spaces
// varPrefix=='#' to format into my string library expression format
//		otherwise format as mathexpr compatible format
char* ROperation::Expr(char varPrefix) const
{
	char *s=NULL,*s1=NULL,*s2=NULL;
	int n=10;
	signed char f=0,g=0;
	if(op==Fun)
	{	if(strlen(pfunc->name)>4)
			n+=strlen(pfunc->name)-4;
	}
	if(mmb1!=NULL)
	{	s1=mmb1->Expr(varPrefix);
		n+=strlen(s1);
		f=IsFunction(mmb1->op);
	}
	if(mmb2!=NULL)
	{	s2=mmb2->Expr(varPrefix);
		n+=strlen(s2);
		g=IsFunction(mmb2->op);
	}
	s=new char[n];
	
	switch(op)
	{	case Num:
			return ValToStr(ValC);
		case Var:
			// precede variable name which character, unless it is blank character
			if(varPrefix=='#')
				sprintf(s,"#%s",pvar->name);
			else
				return CopyStr(pvar->name);
			break;
		case Juxt:
			sprintf(s,"%s,%s",s1,s2);
			break;
		case Add:
			f=f||(mmb1->op==Juxt);
			g=g||(mmb2->op==Juxt);
			if(f&&g)
				sprintf(s,"(%s)+(%s)",s1,s2);
			else if(f)
				sprintf(s,"(%s)+%s",s1,s2);
			else if(g)
				sprintf(s,"%s+(%s)",s1,s2);
			else
				sprintf(s,"%s+%s",s1,s2);
			break;
		case Sub:
			f=f||(mmb1->op==Juxt);
			g=g||(mmb2->op==Juxt||mmb2->op==Add||mmb2->op==Sub);
			if(f&&g)
				sprintf(s,"(%s)-(%s)",s1,s2);
			else if(f)
				sprintf(s,"(%s)-%s",s1,s2);
			else if(g)
				sprintf(s,"%s-(%s)",s1,s2);
			else
				sprintf(s,"%s-%s",s1,s2);
			break;
		case Opp:
			if(mmb2->op==Add||mmb2->op==Sub||mmb2->op==Juxt)
				sprintf(s,"-(%s)",s2);
			else
				sprintf(s,"-%s",s2);
			break;
		case Mult:
			f=f||(mmb1->op==Juxt||mmb1->op==Add||mmb1->op==Sub||mmb1->op==Opp||mmb1->op==Div);
			g=g||(mmb2->op==Juxt||mmb2->op==Add||mmb2->op==Sub||mmb2->op==Opp);
			if(f&&g)
				sprintf(s,"(%s)*(%s)",s1,s2);
			else if(f)
				sprintf(s,"(%s)*%s",s1,s2);
			else if(g)
				sprintf(s,"%s*(%s)",s1,s2);
			else
				sprintf(s,"%s*%s",s1,s2);
			break;
		case Div:
			f=f||(mmb1->op==Juxt||mmb1->op==Add||mmb1->op==Sub||mmb1->op==Opp||mmb1->op==Div);
			g=g||(mmb2->op==Juxt||mmb2->op==Add||mmb2->op==Sub||mmb2->op==Opp||mmb2->op==Mult||mmb2->op==Div);
			if(f&&g)sprintf(s,"(%s)/(%s)",s1,s2);else
			if(f)sprintf(s,"(%s)/%s",s1,s2);else
			if(g)sprintf(s,"%s/(%s)",s1,s2);else
			sprintf(s,"%s/%s",s1,s2);
			break;
		case Pow:
			f=(mmb1->op!=Num&&mmb1->op!=Var);
			g=(mmb2->op!=Num&&mmb2->op!=Var);
			if(f&&g)sprintf(s,"(%s)^(%s)",s1,s2);else
			if(f)sprintf(s,"(%s)^%s",s1,s2);else
			if(g)sprintf(s,"%s^(%s)",s1,s2);else
			sprintf(s,"%s^%s",s1,s2);
			break;
		case Sqrt:
			sprintf(s,"sqrt(%s)",s2);
			break;
		case NthRoot:
			f=(mmb1->op!=Num&&mmb1->op!=Var);
			g=(mmb2->op!=Num&&mmb2->op!=Var);
			if(f&&g)
				sprintf(s,"(%s)#(%s)",s1,s2);
			else if(f)
				sprintf(s,"(%s)#%s",s1,s2);
			else if(g)
				sprintf(s,"%s#(%s)",s1,s2);
			else
				sprintf(s,"%s#%s",s1,s2);
			break;
		case E10:
			f=(mmb1->op!=Num&&mmb1->op!=Var);
			g=(mmb2->op!=Num&&mmb2->op!=Var);
			if(f&&g)sprintf(s,"(%s)E(%s)",s1,s2);else
			if(f)sprintf(s,"(%s)E%s",s1,s2);else
			if(g)sprintf(s,"%sE(%s)",s1,s2);else
			sprintf(s,"%sE%s",s1,s2);
			break;
		case Ln:
			if(varPrefix=='#')
				sprintf(s,"log(%s)",s2);
			else
				sprintf(s,"ln(%s)",s2);
			break;
		case Log:
			if(varPrefix=='#')
				printf(s,"log10(%s)",s2);
			else
				sprintf(s,"log(%s)",s2);
			break;
		case Exp:
			sprintf(s,"exp(%s)",s2);
			break;
		case Sin:
			sprintf(s,"sin(%s)",s2);
			break;
		case Cos:
			sprintf(s,"cos(%s)",s2);
			break;
		case Tg:
			sprintf(s,"tan(%s)",s2);
			break;
		case Atan:
			sprintf(s,"atan(%s)",s2);
			break;
		case Asin:
			sprintf(s,"asin(%s)",s2);
			break;
		case Acos:
			sprintf(s,"acos(%s)",s2);
			break;
		case Abs:
			sprintf(s,"abs(%s)",s2);
			break;
		case IntFun:
			sprintf(s,"int(%s)",s2);
			break;
		case Sign:
			sprintf(s,"sign(%s)",s2);
			break;
		case Sinh:
			sprintf(s,"sinh(%s)",s2);
			break;
		case Cosh:
			sprintf(s,"cosh(%s)",s2);
			break;
		case Tanh:
			sprintf(s,"tanh(%s)",s2);
			break;
		case Fun:
			sprintf(s,"%s(%s)",pfunc->name,s2);
			break;
		default:
			return CopyStr("Error");
	}
  
	if(s1!=NULL) delete[] s1;
	if(s2!=NULL) delete[] s2;
	return s;
}

const double sqrtmaxfloat=sqrt(DBL_MAX);
const double sqrtminfloat=sqrt(DBL_MIN);
const double inveps=.1/DBL_EPSILON;

void  Addition(double*&p)
{	if(*p==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*(--p)=ErrVal;
		return;
	}
	if(*(--p)==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*p=ErrVal;
		return;
	}
	*p+=(*(p+1));
}

void  Soustraction(double*&p)
{	if(*p==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*(--p)=ErrVal;
		return;
	}
	if(*(--p)==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*p=ErrVal;
		return;
	}
	*p-=(*(p+1));
}

void  Multiplication(double*&p)
{	if(fabsl(*p)<sqrtminfloat)
	{	*--p=0;
		return;
	}
	if(*p==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*(--p)=ErrVal;
		return;
	}
	if(fabsl(*(--p))<sqrtminfloat)
	{	*p=0;
		return;
	}
	if(*p==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*p=ErrVal;
		return;
	}
	*p*=(*(p+1));
}

void  Division(double*&p)
{	if(fabsl(*p)<sqrtminfloat || *p==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*(--p)=ErrVal;
		return;
	}
	if(fabsl(*(--p))<sqrtminfloat)
		*p=0;
	else if(*p==ErrVal || fabsl(*p)>sqrtmaxfloat)
	{	*p=ErrVal;
		return;
	}
	*p/=(*(p+1));
}

void  Puissance(double*&p)
{	double v2=*p--,v1=*p;
	if(!v1)
	{	*p=0;
		return;
	}
	if(v2==ErrVal || v1==ErrVal || fabsl(v2*logl(fabsl(v1)))>DBL_MAX_EXP)
	{	*p=ErrVal;
		return;
	}
	*p=((v1>0 || !fmodl(v2,1)) ? powl(v1,v2) : ErrVal);
}

void  RacineN(double*&p)
{	double v2=*p--,v1=*p;
	if(v1==ErrVal || v2==ErrVal || !v1 || v2*logl(fabsl(v1))<DBL_MIN_EXP)
	{	*p=ErrVal;
		return;
	}
	if(v2>=0)
	{	*p=powl(v2,1/v1);
		return;
	}
	*p=((fabsl(fmodl(v1,2))==1) ? -powl(-v2,1/v1) : ErrVal);
}

void  Puiss10(double*&p)
{	if(fabsl(*p)<sqrtminfloat)
	{	*(--p)=0;
		return;
	}
	if(*p==ErrVal || fabsl(*p)>DBL_MAX_10_EXP)
	{	*(--p)=ErrVal;
		return;
	}
	if(fabsl(*(--p))<sqrtminfloat)
		*p=0;
	else if(*p==ErrVal||fabsl(*p)>sqrtmaxfloat)
	{	*p=ErrVal;
		return;
	}
	*p*=pow10l(*(p+1));
}

void  ArcTangente2(double*&p)
{	if(*p==ErrVal || fabsl(*p)>inveps)
	{	*(--p)=ErrVal;
		return;
	}
	if(*(--p)==ErrVal || fabsl(*p)>inveps)
	{	*p=ErrVal;
		return;
	}
	*p=(*p||*(p+1) ? atan2(*p,*(p+1)) : ErrVal);
}

void  NextVal(double*&) {}

void  RFunc(double*&) {}

void  JuxtF(double*&) {}

void  Absolu(double*&p) { *p=((*p==ErrVal) ? ErrVal : fabsl(*p)); }

void  Intpart(double*&p) { *p=((*p==ErrVal) ? ErrVal : floor(*p)); }

void  SignCalc(double*&p)
{	if(*p==ErrVal)
		*p=ErrVal;
	else
		*p= (*p>0.) ? 1. : 0.;
}

void  SinhCalc(double*&p) { *p=((*p==ErrVal) ? ErrVal : (expl(*p)-expl(-*p))/2.); }

void  CoshCalc(double*&p) { *p=((*p==ErrVal) ? ErrVal : (expl(*p)+expl(-*p))/2.); }

void  TanhCalc(double*&p) { *p=((*p==ErrVal) ? ErrVal : (1.-expl(-2.*(*p)))/(1.+expl(-2.*(*p)))); }

void  Oppose(double*&p) { *p=((*p==ErrVal) ? ErrVal : -*p); }

void  ArcSinus(double*&p)
{	*p=((*p==ErrVal || fabsl(*p)>1) ? ErrVal : asinl(*p));
}

void  ArcCosinus(double*&p)
{	*p=((*p==ErrVal||fabsl(*p)>1) ? ErrVal : acosl(*p));
}

void  ArcTangente(double*&p)
{	*p=((*p==ErrVal) ? ErrVal : atanl(*p));
}

void  Logarithme(double*&p)
{	*p=((*p==ErrVal||*p<=0) ? ErrVal : logl(*p));
}

void  LogarithmeTen(double*&p)
{	*p=((*p==ErrVal||*p<=0) ? ErrVal : log10(*p));
}

void  Exponentielle(double*&p)
{	*p=((*p==ErrVal||*p>DBL_MAX_EXP) ? ErrVal : expl(*p));
}

void  Sinus(double*&p)
{	*p=((*p==ErrVal||fabsl(*p)>inveps) ? ErrVal : sinl(*p));
}

void  Tangente(double*&p)
{	*p=((*p==ErrVal||fabsl(*p)>inveps) ? ErrVal : tanl(*p));
}

void  Cosinus(double*&p)
{	*p=((*p==ErrVal||fabsl(*p)>inveps) ? ErrVal : cosl(*p));
}

void  Racine(double*&p)
{	*p=((*p==ErrVal||*p>sqrtmaxfloat||*p<0) ? ErrVal : sqrtl(*p));
}

void FonctionError(double*&p) { *p=ErrVal; }

inline void ApplyRFunc(PRFunction rf,double*&p)
{	p-=rf->nvars-1;
	*p=rf->Val(p);
}

double ROperation::Val() const
{
	pfoncld *p1=pinstr;
	double **p2=pvals,*p3=ppile-1;
	PRFunction *p4=pfuncpile;
	for(;*p1!=NULL;p1++)
    {	if(*p1==&NextVal)
			*(++p3)=**(p2++);
		else if(*p1==&RFunc)
			ApplyRFunc(*(p4++),p3);
		else
			(**p1)(p3);
	}
	return *p3;
}

void BCDouble(pfoncld*&pf,pfoncld*pf1,pfoncld*pf2,
			double**&pv,double**pv1,double**pv2,
			double*&pp,double*pp1,double*pp2,
			RFunction**&prf,RFunction**prf1,RFunction**prf2,
			pfoncld f)
{
	pfoncld*pf3,*pf4=pf1;long n1,n2;
	for(n1=0;*pf4!=NULL;pf4++,n1++);for(n2=0,pf4=pf2;*pf4!=NULL;pf4++,n2++);
	pf=new pfoncld[n1+n2+2];
	for(pf3=pf,pf4=pf1;*pf4!=NULL;pf3++,pf4++)*pf3=*pf4;
	for(pf4=pf2;*pf4!=NULL;pf3++,pf4++)*pf3=*pf4;
	*pf3++=f;*pf3=NULL;//delete[]pf1,pf2;
	double**pv3,**pv4=pv1;
	for(n1=0;*pv4!=NULL;pv4++,n1++);for(n2=0,pv4=pv2;*pv4!=NULL;pv4++,n2++);
	pv=new double*[n1+n2+1];
	for(pv3=pv,pv4=pv1;*pv4!=NULL;pv3++,pv4++)*pv3=*pv4;
	for(pv4=pv2;*pv4!=NULL;pv3++,pv4++)*pv3=*pv4;
	*pv3=NULL;//delete[]pv1,pv2;
	double*pp3,*pp4=pp1;
	for(n1=0;*pp4!=ErrVal;pp4++,n1++);for(n2=0,pp4=pp2;*pp4!=ErrVal;pp4++,n2++);
	pp=new double[n1+n2+1];  // Really need to add and not to take max(n1,n2) in case of Juxt operator
	for(pp3=pp,pp4=pp1;*pp4!=ErrVal;pp3++,pp4++)*pp3=0;
	for(pp4=pp2;*pp4!=ErrVal;pp3++,pp4++)*pp3=0;
	*pp3=ErrVal;//delete[]pp1,pp2;
	PRFunction*prf3,*prf4=prf1;
	for(n1=0;*prf4!=NULL;prf4++,n1++);for(n2=0,prf4=prf2;*prf4!=NULL;prf4++,n2++);
	prf=new PRFunction[n1+n2+1];
	for(prf3=prf,prf4=prf1;*prf4!=NULL;prf3++,prf4++)*prf3=*prf4;
	for(prf4=prf2;*prf4!=NULL;prf3++,prf4++)*prf3=*prf4;
	*prf3=NULL;//delete[]prf1,prf2;
}

void BCSimple(pfoncld*&pf,pfoncld*pf1,double**&pv,double**pv1,
				double*&pp,double*pp1,RFunction**&prf,RFunction**prf1,pfoncld f)
{
	pfoncld*pf3,*pf4=pf1;
	long n;
	for(n=0;*pf4!=NULL;pf4++,n++);
	pf=new pfoncld[n+2];
	for(pf4=pf1,pf3=pf;*pf4!=NULL;pf3++,pf4++) *pf3=*pf4;
	*pf3++=f;
	*pf3=NULL;		//delete[]pf1;
	double**pv3,**pv4=pv1;
	for(n=0;*pv4!=NULL;pv4++,n++);
	pv=new double*[n+1];
	for(pv3=pv,pv4=pv1;*pv4!=NULL;pv3++,pv4++) *pv3=*pv4;
	*pv3=NULL;		//delete[]pv1;
	double*pp3,*pp4=pp1;
	for(n=0;*pp4!=ErrVal;pp4++,n++);
	pp=new double[n+1];
	for(pp3=pp,pp4=pp1;*pp4!=ErrVal;pp3++,pp4++) *pp3=0;
	*pp3=ErrVal;	//delete[]pp1;
	RFunction**prf3,**prf4=prf1;
	for(n=0;*prf4!=NULL;prf4++,n++);
	prf=new RFunction*[n+1];
	for(prf3=prf,prf4=prf1;*prf4!=NULL;prf3++,prf4++) *prf3=*prf4;
	*prf3=NULL;	//delete[]prf1;
}

void BCFun(pfoncld*&pf,pfoncld*pf1,double**&pv,double**pv1,
			double*&pp,double*pp1,RFunction**&prf,RFunction**prf1,PRFunction rf)
{
	pfoncld*pf3,*pf4=pf1;long n;
	for(n=0;*pf4!=NULL;pf4++,n++);
	pf=new pfoncld[n+2];
	for(pf4=pf1,pf3=pf;*pf4!=NULL;pf3++,pf4++)*pf3=*pf4;
	*pf3++=&RFunc;*pf3=NULL;//delete[]pf1;
	double**pv3,**pv4=pv1;
	for(n=0;*pv4!=NULL;pv4++,n++);
	pv=new double*[n+1];
	for(pv3=pv,pv4=pv1;*pv4!=NULL;pv3++,pv4++)*pv3=*pv4;
	*pv3=NULL;//delete[]pv1;
	double*pp3,*pp4=pp1;
	for(n=0;*pp4!=ErrVal;pp4++,n++);
	pp=new double[n+1];
	for(pp3=pp,pp4=pp1;*pp4!=ErrVal;pp3++,pp4++)*pp3=0;
	*pp3=ErrVal;//delete[]pp1;
	PRFunction*prf3,*prf4=prf1;
	for(n=0;*prf4!=NULL;prf4++,n++);
	prf=new PRFunction[n+2];
	for(prf4=prf1,prf3=prf;*prf4!=NULL;prf3++,prf4++)*prf3=*prf4;
	*prf3++=rf;*prf3=NULL;//delete[]pf1;
}

void ROperation::BuildCode()
{
	//  if(mmb1!=NULL)mmb1->BuildCode();if(mmb2!=NULL)mmb2->BuildCode();
	if(pinstr!=NULL)
	{	delete [] pinstr;
		pinstr=NULL;
	}
	if(pvals!=NULL)
	{	delete [] pvals;	//does not delete pvals[0] in case it was to be deleted... (no way to know)
		pvals=NULL;
	}
	if(ppile!=NULL)
	{	delete [] ppile;
		ppile=NULL;
	}
	if(pfuncpile!=NULL)
	{	delete [] pfuncpile;
		pfuncpile=NULL;
	}
	
	switch(op)
	{	case ErrOp:
			pinstr=new pfoncld[2];
			pinstr[0]=&NextVal;
			pinstr[1]=NULL;
			pvals=new double*[2];
			pvals[0]=new double(ErrVal);
			pvals[1]=NULL;
			ppile=new double[2];
			ppile[0]=0;ppile[1]=ErrVal;
			pfuncpile=new RFunction*[1];
			pfuncpile[0]=NULL;
			break;
		case Num:
			pinstr=new pfoncld[2];pinstr[0]=&NextVal;pinstr[1]=NULL;
			pvals=new double*[2];pvals[0]=new double(ValC);pvals[1]=NULL;
			ppile=new double[2];ppile[0]=0;ppile[1]=ErrVal;
			pfuncpile=new RFunction*[1];pfuncpile[0]=NULL;
			break;
		case Var:
			pinstr=new pfoncld[2];pinstr[0]=&NextVal;pinstr[1]=NULL;
			pvals=new double*[2];pvals[0]=pvarval;pvals[1]=NULL;
			ppile=new double[2];ppile[0]=0;ppile[1]=ErrVal;
			pfuncpile=new RFunction*[1];pfuncpile[0]=NULL;
			break;
		case Juxt:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,
					pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&JuxtF);
			break;
		case Add:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&Addition);
			break;
		case Sub:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&Soustraction);
			break;
		case Mult:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&Multiplication);
			break;
		case Div:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&Division);
			break;
		case Pow:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&Puissance);
			break;
		case NthRoot:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&RacineN);
			break;
		case E10:
			BCDouble(pinstr,mmb1->pinstr,mmb2->pinstr,
				pvals,mmb1->pvals,mmb2->pvals,ppile,mmb1->ppile,mmb2->ppile,pfuncpile,mmb1->pfuncpile,mmb2->pfuncpile,&Puiss10);
			break;
		case Opp:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Oppose);
			break;
		case Sin:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Sinus);
			break;
		case Sqrt:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Racine);
			break;
		case Ln:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Logarithme);
			break;
		case Log:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&LogarithmeTen);
			break;
		case Exp:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Exponentielle);
			break;
		case Cos:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Cosinus);
			break;
		case Tg:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Tangente);
			break;
		case Atan:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,(mmb2->NMembers()>1?&ArcTangente2:&ArcTangente));
			break;
		case Asin:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&ArcSinus);
			break;
		case Acos:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&ArcCosinus);
			break;
		case Abs:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Absolu);
			break;
		case IntFun:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&Intpart);
			break;
		case Sign:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&SignCalc);
			break;
		case Sinh:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&SinhCalc);
			break;
		case Cosh:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&CoshCalc);
			break;
		case Tanh:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&TanhCalc);
			break;
		case Fun:
			BCFun(pinstr,mmb2->pinstr,pvals,mmb2->pvals,ppile,
				mmb2->ppile,pfuncpile,mmb2->pfuncpile,pfunc);
			break;
		default:
			BCSimple(pinstr,mmb2->pinstr,pvals,mmb2->pvals,
				ppile,mmb2->ppile,pfuncpile,mmb2->pfuncpile,&FonctionError);
			break;
	}
}
