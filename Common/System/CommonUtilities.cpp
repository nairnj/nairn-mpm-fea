/********************************************************************************
    CommonUtilities.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 1/12/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include <stdarg.h>

#define T_SQRT3    1.73205080756887729352744634151   // sqrt(3)

// local globals
static int section=1;

// local globals for random numbers
static unsigned long randomMax;
static double randomScale;
static double halfRandomBoxSize;

#pragma mark Miscellaneous Functions

/*********************************************************************
    Start next section in results window
*********************************************************************/

void PrintSection(const char *text)
{
    char sline[200];
	size_t ssize=200;
	snprintf(sline,ssize,"*****%3d. %s\n\n",section,text);
    cout << sline;
	section++;
}

/************************************************************
	See if two double precision numbers are sufficiently
		equal (currently 100 ppm)
	If both are less than 1e-12, they are assumed equal
************************************************************/

#define MAX_DIFF 1.0e-16
#define MAX_REL_DIFF 1.0e-7

bool DbleEqual(double A, double B)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	double diff = fabs(A - B);
	if (diff <= MAX_DIFF)
		return true;

	// Check relative difference
	A = fabs(A);
	B = fabs(B);
	double largest = (B > A) ? B : A;
	if (diff <= largest * MAX_REL_DIFF)
		return true;

	// assuming unequal
	return false;
}

/************************************************************
    Reverse block of bytes to convert between Intel
	processor and Mac processor
************************************************************/

int Reverse(char *nptr,int nlen)
{
    int i,imax;
    char *rptr,temp;
    
    rptr=nptr+nlen-1;
    imax=nlen>>1;
    for(i=0;i<imax;i++)
    {	temp=*nptr;
        *nptr++=*rptr;
        *rptr--=temp;
    }
    return nlen;
}

/************************************************************
	Factor integer num into prime factors (not including 1)
************************************************************/

void PrimeFactors(int num,vector<int> &factors)
{	int maxFactor = (int)(sqrt((double)num)+.001);
	for(int i=2;i<=maxFactor;i++)
	{	// does i go into num evenly?
		if((num % i) == 0)
		{	factors.push_back(i);
			PrimeFactors(num/i, factors);
			return;
		}
	}
	// if did not find any factors, it is a prime
	factors.push_back(num);
}

#pragma mark Vector Functions

// return vector from components
Vector MakeVector(double x,double y,double z)
{	Vector v;
	v.x=x;
	v.y=y;
	v.z=z;
	return v;
}

// return vector from components all times scale factor
Vector MakeScaledVector(double x,double y,double z,double scale)
{	Vector v;
	v.x=x*scale;
	v.y=y*scale;
	v.z=z*scale;
	return v;
}

// return vector from difference of two vectors
Vector SetDiffVectors(const Vector *a,const Vector *b)
{	Vector v;
	v.x=a->x-b->x;
	v.y=a->y-b->y;
	v.z=a->z-b->z;
	return v;
}

// return vector from sum of two vectors
Vector SetSumVectors(const Vector *a,const Vector *b)
{	Vector v;
	v.x=a->x+b->x;
	v.y=a->y+b->y;
	v.z=a->z+b->z;
	return v;
}

// return vector from components all times scale factor
Vector SetScaledVector(const Vector *a,const double scale)
{	Vector v;
	v.x=a->x*scale;
	v.y=a->y*scale;
	v.z=a->z*scale;
	return v;
}

// zero a vector and return pointer to zeroed v
Vector *ZeroVector(Vector *v)
{	v->x=v->y=v->z=0.;
	return v;
}

// copy vector and return pointer to new one
Vector *CopyVector(Vector *v,const Vector *oldv)
{	v->x=oldv->x;
	v->y=oldv->y;
	v->z=oldv->z;
	return v;
}

// scale a vector and return pointer to changed vector
Vector *ScaleVector(Vector *v,const double scale)
{	v->x*=scale;
	v->y*=scale;
	v->z*=scale;
	return v;
}

// scale a vector and copy to new one, then return pointer to new v
Vector *CopyScaleVector(Vector *v,const Vector *oldv,const double scale)
{	v->x=scale*oldv->x;
	v->y=scale*oldv->y;
	v->z=scale*oldv->z;
	return v;
}

// add vector toadd to vector v and return pointer to changed v
Vector *AddVector(Vector *v,const Vector *toadd)
{	v->x+=toadd->x;
	v->y+=toadd->y;
	v->z+=toadd->z;
	return v;
}

// sub vector toadd from vector v and return pointer to changed v
Vector *SubVector(Vector *v,const Vector *toadd)
{	v->x-=toadd->x;
	v->y-=toadd->y;
	v->z-=toadd->z;
	return v;
}

// add vector toadd scaled by scale to vector v and return pointer to changed v
Vector *AddScaledVector(Vector *v,const Vector *toadd,const double scale)
{	v->x+=toadd->x*scale;
	v->y+=toadd->y*scale;
	v->z+=toadd->z*scale;
	return v;
}

// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotVectors(const Vector *v1,const Vector *v2)
{	return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z ;
}

// Dot prodiuct of two 2D vectors
double DotVectors2D(const Vector *v1,const Vector *v2)
{	return v1->x*v2->x + v1->y*v2->y ;
}

// Vector cross product of two vectors in a third vector
Vector *CrossProduct(Vector *cp,const Vector *a,const Vector *b)
{	cp->x = a->y*b->z - b->y*a->z ;
	cp->y = b->x*a->z - a->x*b->z ;
	cp->z = a->x*b->y - b->x*a->y ;
	return cp;
}

// Z component of vector cross product of two vectors in x-y plane
double CrossProduct2D(const Vector *a,const Vector *b)
{	return a->x*b->y - b->x*a->y ;
}

// Print vector to cout when debugging
void PrintVector(const char *label,const Vector *v)
{	if(strlen(label)>0) cout << label;
	cout << "(" << v->x << ", " << v->y << ", " << v->z << ") ";
}

// Given piecewise pairs in xpts and ypts, find y(x)
// Assumes xpts and ypts are the same length and xpts is monotonically increasing
// Returns 1 if arrays are empty or single value if arrays have only 1 value
// Points outside range use linear extrapolation
double PiecewiseInterpolate(double x,vector<double> xpts,vector<double> ypts)
{
    // empty returns 1
    if(xpts.size()==0) return 1;
    
    // single value is constant
    if(xpts.size()==1) return ypts[0];
    
    // get interval for x value
    int intvl = (int)xpts.size()-1;
    for(int i=1;i<xpts.size();i++)
    {   // found the interval
        if(x<=xpts[i])
        {   intvl = i;
            break;
        }
    }
    
    // Here x<=xpts[1], then intvl = 1
    //   x>xpts[size-1], then intvl = size-1
    //   otherwise xpts[intvl-1] < x <= xpts[intvl]
    // Note: assums xpts[intvl]>xpts[intvl-1]
    double fract = (xpts[intvl]-x)/(xpts[intvl]-xpts[intvl-1]);
    return fract*ypts[intvl-1] + (1.-fract)*ypts[intvl];
}

#pragma mark Tensor Functions

// return tensor from components
Tensor MakeTensor(double xx,double yy,double zz,double yz,double xz,double xy)
{	Tensor t;
	t.xx = xx;
	t.yy = yy;
	t.zz = zz;
#ifdef MPM_CODE
	t.yz = yz;
	t.xz = xz;
#endif
	t.xy = xy;
	return t;
}

// return tensor from components
Tensor MakeTensor2D(double xx,double yy,double zz,double xy)
{	Tensor t;
	t.xx = xx;
	t.yy = yy;
	t.zz = zz;
#ifdef MPM_CODE
	t.yz = 0.;
	t.xz = 0.;
#endif
	t.xy = xy;
	return t;
}

// zero a vector and return pointer to zeroed v
Tensor *ZeroTensor(Tensor *t)
{	t->xx=0.;
	t->yy=0.;
	t->zz=0.;
#ifdef MPM_CODE
	t->yz=0.;
	t->xz=0.;
#endif
	t->xy=0.;
	return t;
}

// add tensor toadd to tensor t and return pointer to changed t
Tensor *AddTensor(Tensor *t,Tensor *toadd)
{	t->xx+=toadd->xx;
	t->yy+=toadd->yy;
	t->zz+=toadd->zz;
#ifdef MPM_CODE
	t->yz+=toadd->yz;
	t->xz+=toadd->xz;
#endif
	t->xy+=toadd->xy;
	return t;
}

// subtract tensor tosub from tensor t and return pointer to changed t
Tensor *SubTensor(Tensor *t,Tensor *tosub)
{	t->xx-=tosub->xx;
	t->yy-=tosub->yy;
	t->zz-=tosub->zz;
#ifdef MPM_CODE
	t->yz-=tosub->yz;
	t->xz-=tosub->xz;
#endif
	t->xy-=tosub->xy;
	return t;
}

// scale tensor t and return pointer to changed t
Tensor *ScaleTensor(Tensor *t,double scale)
{	t->xx*=scale;
	t->yy*=scale;
	t->zz*=scale;
#ifdef MPM_CODE
	t->yz*=scale;
	t->xz*=scale;
#endif
	t->xy*=scale;
	return t;
}

// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotTensors2D(const Tensor *t1,const Tensor *t2)
{	return t1->xx*t2->xx + t1->yy*t2->yy + t1->zz*t2->zz + t1->xy*t2->xy;
}

#ifdef MPM_CODE
// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotTensors(const Tensor *t1,const Tensor *t2)
{	return t1->xx*t2->xx + t1->yy*t2->yy + t1->zz*t2->zz
	+ t1->yz*t2->yz + t1->xz*t2->xz + t1->xy*t2->xy;
}

// Find Eigenvalues of symmetric tensor in Voight form
// stress is true for stress vector or false for strain vector
// is2D is true for 2D tensor or false for 3D
Vector TensorEigenvalues(Tensor *t,bool stress,bool is2D)
{
	Vector lam;
	
	if(is2D)
	{   // solving x^2 + bx + c = 0 (see Numerical Recipes in C, page 156)
		double b = -(t->xx+t->yy);
		double c = stress ? t->xx*t->yy - t->xy*t->xy : t->xx*t->yy - 0.25*t->xy*t->xy ;
		double arg = b*b-4.*c;
		if(arg<0.)
		{	// assuming here all matrices are positive definite, which means
			// a negative value should be interpreted as zero
			lam.x = -0.5*b;
		}
		else
		{	arg = sqrt(arg);
			lam.x = b>0 ? -0.5*(b+arg) : -0.5*(b-arg) ;
		}
		lam.y = lam.x==0. ? 0. : c/lam.x;
		lam.z = t->zz;
	}
	else
	{	double mm, c1, c0;
		
		// Determine coefficients of characteristic poynomial. We write
		//       | a   d   f  |
		//  m =  | d   b   e  |
		//       | f   e   c  |
		double fde,dd,ee,ff;		// products of elements
		if(stress)
		{	fde = t->xy*t->yz*t->xz;
			dd = t->xy*t->xy;
			ee = t->yz*t->yz;
			ff = t->xz*t->xz;
		}
		else
		{	fde = 0.125*t->xy*t->yz*t->xz;
			dd = 0.25*t->xy*t->xy;
			ee = 0.25*t->yz*t->yz;
			ff = 0.25*t->xz*t->xz;
		}
		mm  = t->xx + t->yy + t->zz;
		c1 = (t->xx*t->yy + t->xx*t->zz + t->yy*t->zz)
					- (dd + ee + ff);				// a*b + a*c + b*c - d^2 - e^2 - f^2
		c0 = t->zz*dd + t->xx*ee + t->yy*ff - t->xx*t->yy*t->zz
					- 2.0*fde;						// c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
		
		double p, sqrt_p, q, c, s, phi;
		p = mm*mm - 3.0*c1;
		q = mm*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
		sqrt_p = sqrt(fabs(p));
		
		phi = 27.0 * ( 0.25*c1*c1*(p - c1) + c0*(q + 27.0/4.0*c0));
		phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
		
		c = sqrt_p*cos(phi);
		s = (1.0/T_SQRT3)*sqrt_p*sin(phi);
		
		lam.y  = (1.0/3.0)*(mm - c);
		lam.z  = lam.y + s;
		lam.x  = lam.y + c;
		lam.y -= s;
	}
	
	return lam;
}

// convert tensor to matrix (stress is true for straight convertion of false to assume strain tensor)
Matrix3 TensorToMatrix(Tensor *t,bool stress)
{	return stress ?
			Matrix3(t->xx,t->xy,t->xz,t->xy,t->yy,t->yz,t->xz,t->yz,t->zz) :
			Matrix3(t->xx,0.5*t->xy,0.5*t->xz,0.5*t->xy,t->yy,0.5*t->yz,0.5*t->xz,0.5*t->yz,t->zz) ;
}

// convert tensor to matrix (stress is true for straight convertion of false to assume strain tensor)
Matrix3 TensorToMatrix2D(Tensor *t,bool stress)
{	return stress ?
			Matrix3(t->xx,t->xy,t->xy,t->yy,t->zz) :
			Matrix3(t->xx,0.5*t->xy,0.5*t->xy,t->yy,t->zz) ;
}

#endif

// get element of contracted tensor by ID
double Tensor_i(Tensor *t,int contractedID)
{
	switch(contractedID)
	{	case XX:
			return t->xx;
		case YY:
			return t->yy;
		case ZZ:
			return t->zz;
#ifdef MPM_CODE
		case YZ:
			return t->yz;
		case XZ:
			return t->xz;
#endif
		case XY:
			return t->xy;
		default:
			return 0.;
	}
}

// get element of tensor by row and colume
// For FEA can only get xx or 1,1, xy or 1,2 or 2,1, yy or 2,2, anbd zz or 3,3
double Tensor_ij(Tensor *t,int row,int col)
{
	switch(row)
	{	case 1:
			switch(col)
			{	case 1:
					return t->xx;
				case 2:
					return t->xy;
#ifdef MPM_CODE
				case 3:
					return t->xz;
#endif
				default:
					return 0.;
			}
		case 2:
			switch(col)
			{	case 1:
					return t->xy;
				case 2:
					return t->yy;
#ifdef MPM_CODE
				case 3:
					return t->yz;
#endif
				default:
					return 0.;
			}
		case 3:
			switch(col)
			{
#ifdef MPM_CODE
				case 1:
					return t->xz;
				case 2:
					return t->yz;
#endif
				case 3:
					return t->zz;
				default:
					return 0.;
			}
		default:
			return 0.;
	}
}

// Print tensor to cout when debuggingin form of a matrix
void PrintTensor(const char *label,Tensor *t)
{
	int i,lead=(int)strlen(label);
	for(i=0;i<lead;i++) cout << " ";
#ifdef MPM_CODE
	cout << "(" << t->xx << ", " << t->xy << "," << t->xz << ")" << endl;
	cout << label << "(" << t->xy << ", " << t->yy << "," << t->yz << ")" << endl;
#else
	cout << "(" << t->xx << ", " << t->xy << ", 0.0)" << endl;
	cout << label << "(" << t->xy << ", " << t->yy << ", 0.0)" << endl;
#endif
	for(i=0;i<lead;i++) cout << " ";
#ifdef MPM_CODE
	cout << "(" << t->xz << ", " << t->yz << "," << t->zz << ")" << endl;
#else
	cout << "( 0.0, 0.0," << t->zz << ")" << endl;
#endif
}

#ifdef MPM_CODE

TensorAntisym *ZeroTensorAntisym(TensorAntisym *a)
{	a->xy=a->xz=a->yz=0.;
	return a;
}

#endif

#pragma mark Rect Functions

// true if (x,y) is in the Rect
bool XYInRect(double x,double y,Rect *rect)
{   if(x<rect->xmin) return false;
    if(x>rect->xmax) return false;
    if(y<rect->ymin) return false;
    if(y>rect->ymax) return false;
    return true;
}

// true if (x,y) of vectir is in the Rect
bool PtInRect(Vector *pt,Rect *rect)
{   if(pt->x<rect->xmin) return false;
    if(pt->x>rect->xmax) return false;
    if(pt->y<rect->ymin) return false;
    if(pt->y>rect->ymax) return false;
    return true;
}

// area of the rect
double RectArea(Rect *rect)
{   return (rect->xmax-rect->xmin)*(rect->ymax-rect->ymin);
}

#pragma mark String Functions

// convert character to lower case
unsigned char ConvertToLower(unsigned char ch)
{	if (ch >= 'A' && ch <= 'Z')
		ch = 'a' + (ch - 'A');
	return ch;
}

// convert path with / to DOS path with back slashes
// return pointer to input string
char *MakeDOSPath(char *unixPath)
{	for (unsigned int i = 0; i < strlen(unixPath); i++)
	{	if (unixPath[i] == '/')
			unixPath[i] = '\\';
	}
	return unixPath;
}

// case insensitive compare
// return 0 if match, return <0 or >0 if first case insentive mismatch character
//      has lower value in s1 or in s2, respectively
int CIstrcmp(const char *s1, const char *s2)
{
	unsigned char *us1 = (unsigned char *)s1;
	unsigned char *us2 = (unsigned char *)s2;
	
	while(ConvertToLower(*us1) == ConvertToLower(*us2++))
	{	// they are equal, so still matching
		
		// if matching at the end, then found a match
		if(*us1 == '\0') return (0);
		
		// otherwise on to next character
		us1++;
	}
	
	// found mismatch
	us2--;
	return (int)(ConvertToLower(*us1) - ConvertToLower(*us2));
}

// file extension of last component in the path (or empty string if none)
void GetFileExtension(const char *fileName,char *ext,int maxLength)
{
	int endLoc = (int)strlen(fileName)-1;
	int dotLoc = endLoc;
	while(dotLoc>=0 && fileName[dotLoc]!='.' && fileName[dotLoc]!='/') dotLoc--;
	
	// found an extension
	ext[0] = 0;
	if(dotLoc>=0)
	{	// empty if hit folder divider
		if(fileName[dotLoc]=='/') return;
		
		// first character after the period
		dotLoc++;
		int ei = 0;
		
		// last one gets the zero terminator
		while(ei<maxLength && dotLoc<=endLoc)
			ext[ei++]=fileName[dotLoc++];
		ext[ei] = 0;
	}
}

#pragma mark Debugging Functions

// output doubles preceded by labels
// fmt contains N labels delimited by spaces or commas (multiple delim = one)
// variable number arguments must be exactly N but label starting in % is
// output text alone without referencing a variable
// Example: dout("a b c",a,b,c) ==> a=4.5 b=3 c=145 (i.e., equal signs are provided)
void dout(const char* fmt...)
{
	va_list args;
	va_start(args, fmt);
	char label[100];
 
	cout << "# ";
	int length=0;
	bool hadLabel=false;
	while(true)
	{	if(*fmt==' ' || *fmt==',' || *fmt == '\0')
		{	// requires a label
			if(length>0)
			{	if(label[0]!='%')
				{	double d = va_arg(args, double);
					if(hadLabel) cout << ",";
					label[length]=0;
					cout << label << "=" << d;
				}
				else
				{	if(hadLabel) cout << " ";
					label[length]=0;
					cout << (label+1);
				}
				length = 0;
				hadLabel = true;
			}
			if(*fmt == '\0') break;
		}
		else
			label[length++] = *fmt;
		++fmt;
	}
	cout << endl;
 
	va_end(args);
}

// Same as dout, but in a critical (output) block
void doutCritical(const char* fmt...)
{
	va_list args;
	va_start(args, fmt);
	char label[100];

#pragma omp critical (output)
	{	cout << "# ";
		int length=0;
		bool hadLabel=false;
		while(true)
		{	if(*fmt==' ' || *fmt==',' || *fmt == '\0')
			{	// requires a label
				if(length>0)
				{	if(label[0]!='%')
					{	double d = va_arg(args, double);
						if(hadLabel) cout << ",";
						label[length]=0;
						cout << label << "=" << d;
					}
					else
					{	if(hadLabel) cout << " ";
						label[length]=0;
						cout << (label+1);
					}
					length = 0;
					hadLabel = true;
				}
				if(*fmt == '\0') break;
			}
			else
				label[length++] = *fmt;
			++fmt;
		}
		cout << endl;
	}
 
	va_end(args);
}

// Same as dout(), but formatted for Mathmatica assignments
// Cannot use %label option
void doutMM(const char* fmt...)
{
	va_list args;
	va_start(args, fmt);
	char label[100];
 
	cout << "# {";
	int length=0;
	bool hadLabel=false;
	while(true)
	{	if(*fmt==' ' || *fmt==',' || *fmt == '\0')
		{	// requires a label
			if(length>0)
			{	double d = va_arg(args, double);
				if(hadLabel) cout << ",";
				label[length]=0;
				cout << label << "->" << d;
				length = 0;
				hadLabel = true;
			}
			if(*fmt == '\0') break;
		}
		else
			label[length++] = *fmt;
		++fmt;
	}
	cout << "}" << endl;
 
	va_end(args);
}

// Same as doutMM(), but in critical block
// Cannot use %label option
void doutMMCritical(const char* fmt...)
{
	va_list args;
	va_start(args, fmt);
	char label[100];

#pragma omp critical (output)
	{	cout << "# {";
		int length=0;
		bool hadLabel=false;
		while(true)
		{	if(*fmt==' ' || *fmt==',' || *fmt == '\0')
			{	// requires a label
				if(length>0)
				{	double d = va_arg(args, double);
					if(hadLabel) cout << ",";
					label[length]=0;
					cout << label << "->" << d;
					length = 0;
					hadLabel = true;
				}
				if(*fmt == '\0') break;
			}
			else
				label[length++] = *fmt;
			++fmt;
		}
		cout << "}" << endl;
	}
 
	va_end(args);
}


#pragma mark Other Functions

// This method is called at start up to initialize parameters for generation
// of random numbers
// Mac rand() o return [1,RAND_MAX-1] with RAND_MAX=2,147,483,647
// Windows rand() returns [0,RAND_MAX] with RAND)MAX=32,767
// Linux(carbon) rand() returns [0,RAND_MAX] with RAND)MAX=2,147,483,647
void InitRandom(unsigned int seed)
{
	// seed random number generator
	if(seed>0)
		srand((unsigned int)seed);
	else
		srand((unsigned int)time(NULL));
	
	if(RAND_MAX<64000)
	{	// This catches Windows currently with RAND_MAX=32767
		// code here extends to larger numbers
		randomMax = 1073741825;
	}
	else
	{	// Other OS's seem to have RAND_MAX=2,147,483,647
		randomMax = RAND_MAX;
	}
	
	// OS-independnent method will get random number from
	// 1 to randomMax-1
	// return value will scale by 1/(randomMax-1) and then subtract 0.5/(randomMax-1)
	randomScale = 1./(double)(randomMax-1);
	halfRandomBoxSize = 0.5*randomScale;

//#define TEST_RANDOM_RANGE
#ifdef TEST_RANDOM
	// divided into nbox
	long nbox=20;
	long nb[20];
	long numr=1000000000;
	for(int i=0;i<nbox;i++) nb[i]=0;
	
	cout << endl;
	cout << "Testing system that has randomMax=" << randomMax << endl;
	cout << "Generating " << numr << " random numbers using Random()" << endl;
	
	for(long nr=0;nr<numr;nr++)
	{	int boxnum = (int)((double)nbox*Random());
		if(boxnum<0 || boxnum>=nbox)
		{	cout << "Random box out of range" << endl;
			break;
		}
		
		// increment count
		nb[boxnum]++;
	}
	
	// Summary
	double bsize = 1./(double)nbox;
	cout << "Expected rate = " << bsize << endl;
	double bmin = 0.;
	for(int i=0;i<nbox;i++)
	{	cout <<"Rate between " << bmin << " and " << bmin+bsize
					<< " = " << (double)nb[i]/(double)numr << endl;
		bmin += bsize;
	}
	cout << endl;

#endif
#ifdef TEST_RANDOM_LONG
	// divided into nbox
	long nb[100];
	long numr=1000000000;
	long rmin=7,rmax=28,nbox=rmax-rmin+1;
	for(int i=0;i<nbox;i++) nb[i]=0;
	
	cout << endl;
	cout << "Testing system that has randomMax=" << randomMax << endl;
	cout << "Generating " << numr << " random numbers using RandomLong()" << endl;
	
	for(long nr=0;nr<numr;nr++)
	{	long boxnum = RandomLong(rmin,rmax)-rmin;
		if(boxnum<0 || boxnum>=nbox)
		{	cout << "Random box out of range" << endl;
			break;
		}
		
		// increment count
		nb[boxnum]++;
	}
	
	// Summary
	double bsize = 1./(double)nbox;
	cout << "Expected rate = " << bsize << endl;
	long bmin = rmin;
	for(int i=0;i<nbox;i++)
	{	cout <<"Rate for " << bmin << " = " << (double)nb[i]/(double)numr << endl;
		bmin++;
	}
	cout << endl;
#endif
#ifdef TEST_RANDOM_RANGE
	// divided into nbox
	long nbox=20;
	long nb[100];
	long numr=1000000000;
	double rmin=-.3;
	double rmax=.7;
	double bsize=(rmax-rmin)/(double)nbox;
	for(int i=0;i<nbox+1;i++) nb[i]=0;
	
	cout << endl;
	cout << "Testing system that has randomMax=" << randomMax << endl;
	cout << "Generating " << numr << " random numbers using RandomRange()" << endl;
	
	for(long nr=0;nr<numr;nr++)
	{	double r = RandomRange(rmin,rmax);
		int boxnum = (int)((r-rmin)/bsize);
		if(boxnum<0 || boxnum>nbox)
		{	cout << "Random box " << boxnum << " for r = " << r << " is out of range " << endl;
			break;
		}
		
		// increment count
		nb[boxnum]++;
	}
	
	// Summary
	cout << "Expected rate = " << bsize << endl;
	double bmin = rmin;
	for(int i=0;i<nbox;i++)
	{	cout <<"Rate between " << bmin << " and " << bmin+bsize
					<< " = " << (double)nb[i]/(double)numr << endl;
		bmin += bsize;
	}
	cout <<"Rate equal to " << rmax << " = " << (double)nb[nbox]/(double)numr << endl;
	cout << endl;
#endif
}

// This method return random long on the interval [rmin,rmax] with
// endpoint included (rmax must be greater than rmin)
long RandomLong(long rmin,long rmax)
{	// get long from 0 to span-1
	long span = rmax-rmin+1;
	long zerom1 = (long)((double)span*Random());
	
	// adjust to rmin
	return rmin+zerom1;
}

// This method return random double on the interval [rmin,rmax] with
// endpoints included (rmax must be greater than rmin)
double RandomRange(double rmin,double rmax)
{	double range = rmax-rmin;
	double returnRange = 1.-2.*halfRandomBoxSize;
	return rmin + Random()*range/returnRange;
}

// This method return random double 0 < r < 1 independent of the OS
// and it is probablity that randum number is r +/- halfRandomBoxSize
// The minimum value is halfRandomBoxSize and maximum value is 1.-halfRandomBoxSize
double Random(void)
{
	// For OS with low RAND_MAX, combine two random numbers to get more precision
	if(RAND_MAX<64000)
	{	// Two 15 bit numbers 0 to 32767, second one shift 15 bits left
		unsigned long r1 = (unsigned long)rand()%(unsigned long)32768;
		unsigned long r2 = ((unsigned long)rand()%(unsigned long)32768) << 15;
		// long number from 1 to randomMax-1
		unsigned long rn = r1+r2+1;
		return randomScale*rn - halfRandomBoxSize;
	}

	// For other OS's, use RAND_MAX (it is usually 2147483647)
	// get long from 1 to RAND_MAX-1 by ignoring 0 and RAND_MAX if they occur
	unsigned long rn = 0;
	while(rn==0 || rn==randomMax) rn = (unsigned long)rand();
	return randomScale*rn - halfRandomBoxSize;
}

// The cummulative dist function is p(x,mean,s) = 0.5(1+erf((x-mean)/(s*sqrt(2))))
// Given p(d,0,1), return d and set (x-mean)/s = d or x = mean + s*d
// Or relative to mean: value/mean = 1 + s*d/mean = 1 + d*(CV)
//      where (CV) is coefficient of variation.
// More info found at http://www.johndcook.com/blog/normal_cdf_inverse/
double NormalCDFInverse(double p)
{	if(p <= 0.0 || p >= 1.0)
	{   // invalid numbers return zero or get mean value
		return 0.;
	}
	else if(p < 0.5)
	{	// F^-1(p) = - G^-1(p)
		return -RationalApproximation( sqrt(-2.0*log(p)) );
	}
	else if (p < 1.0)
	{	// F^-1(p) = G^-1(1-p)
		return RationalApproximation( sqrt(-2.0*log(1-p)) );
	}
	else
		return 1.;
}

// Abramowitz and Stegun formula 26.2.23.
// The absolute value of the error should be less than 4.5 e-4.
// More info found at http://www.johndcook.com/blog/normal_cdf_inverse/
double RationalApproximation(double t)
{	double c[] = {2.515517, 0.802853, 0.010328};
	double d[] = {1.432788, 0.189269, 0.001308};
	return t - ((c[2]*t + c[1])*t + c[0]) /
	(((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

// Return smooth step function that goes from 0 (when x=xmin) to 1 (when x=xmax)
//   as well is 0 for x<xmin and 1 for x>xmax
// order can be 5 or 3 and anything else is 1
// deriv is 1 for derivative, anything else for function
double SmoothStepRange(double x,double xmin,double xmax,int order,int deriv)
{
	double range = xmax-xmin;
	double xi = (x-xmin)/range;
	double stepVal = SmoothStep(xi,order,deriv);
	if(deriv==1) stepVal /= range;
	return stepVal;
}

// Return smooth step function that goes from 0 (when x=0) to 1 (when x=1).
//   as well is 0 for x<0 and 1 for x>1
// order can be 5 or 3 and anything else is order 1 (i.e. h(x)=x)
// deriv is 1 for derivative, anything else for function
double SmoothStep(double x,int order,int deriv)
{
	// check for out of range
	if(x<0.)
		return 0.;
	else if(x>1.)
		return deriv==1 ? 0. : 1. ;
	
	// order 5, 3, anything else is linear
	if(order==5)
	{
		return deriv==1 ? 30.*x*x*(1.-x)*(1.-x) :
			x*x*x*(x*(6.*x-15.)+10.) ;
	}
	else if(order==3)
	{
		return deriv==1 ? 6*x*(1-x) :
				x*x*(3.-2.*x);
	}
	
	// default to first order
	return deriv==1 ? 1 : x ;
}

// stable calculation of real roots to quadratic equation a*x^2 + b*x + c = 0
// Roots returned in r1 and r2
// a cannot be zero, but it is not checked for here
// return true on real roots or false if no real roots
bool RealQuadraticRoots(double a,double b,double c,double &r1,double &r2)
{
	// find terms that must be positive for real roots
	double arg = b*b - 4.*a*c;
	if(arg<0.)
	{	if(fabs(arg)<1.e-6*b*b)
		{	// if b^2 is close for 4.*a*c (difference very small relative to b^2)
			// accept as equal to zero and get two identical real roots
			r1 = -0.5*b/a;
			r2 = r1;
			return true;
		}
		else
		{	// no real roots
			return false;
		}
	}
	
	// find q = -0.5*(b + sgn(b) sqrt(b*b-4*a*c)
	double q = b>0. ? -0.5*(b+sqrt(arg)) : -0.5*(b-sqrt(arg));
	
	// roots are q/a and c/q
	r1 = q/a;
	r2 = c/q;
	
	// found real roots
	return true;
}


#ifdef MPM_CODE
/****************************************************************
 *  Functions for an intersect of plane and  unit cube
 ****************************************************************/

static double ci[5]={-1.,1.,1.,-1.,-1.};

// Find intersected edges
// Prechecked to have intersection and intersection is not a corner
// On input v = 0.5*(dx,dy,dz)
Vector IntersectEdges(int i,int j,Vector *v,Vector *n)
{
	// corner vectors (0,1,2,3 on +dx face, 4,5,6,7 on -dx face)
	Vector Vertex_i = (i<4) ? MakeVector(v->x,ci[i]*v->y,-ci[i+1]*v->z) :
						MakeVector(-v->x,ci[i-4]*v->y,-ci[i-3]*v->z);
	Vector Edge_ij = (j<4) ? MakeVector(v->x,ci[j]*v->y,-ci[j+1]*v->z) :
						MakeVector(-v->x,ci[j-4]*v->y,-ci[j-3]*v->z);
	
	// Fined edge vector difference of coordinates
	Edge_ij.x -=Vertex_i.x;
	Edge_ij.y -=Vertex_i.y;
	Edge_ij.z -=Vertex_i.z;
	
	// Start calculating lambda (assume has an intersection
	double lambda = -DotVectors(&Vertex_i,n)/DotVectors(&Edge_ij,n);
	
	// if not in [0,1] then it doesn't intersect
	if(lambda>1.0 || lambda<0.0)
	{	cout << "# failed to find intersection " << i << "," << j ;
		PrintVector(", c=",&Vertex_i);
		PrintVector(", e=",&Edge_ij);
		cout << ", lambda=" << lambda << endl;
	}
	
	// Find point and return it
	AddScaledVector(&Vertex_i,&Edge_ij,lambda);
	return Vertex_i;
}

// Find are of parallelogram by vectors (v-c) and (-v-c) in plane with normal n
double PGramArea(Vector *c,Vector *v,Vector *n)
{
	// get sides of parallelogram (v-c and -v-c)
	Vector v1 = *v;
	SubVector(&v1,c);
	Vector v2 = SetScaledVector(v,-1.);
	SubVector(&v2,c);
	
	// get cross product
	Vector cp;
	CrossProduct(&cp,&v1,&v2);
	
	// return area of parallogram enclosed by v1 and v2
	// assume in plane with normal n
	return fabs(DotVectors(&cp,n));
}

// get the area of interection between plane through center of box (dx,dy,dz) and the box
double AreaOverVolume3D(Vector *norm,double dx,double dy,double dz)
{
	// For notes see JANOSU-14-44+
	
	// examine corners on face with n=(1,0,0)
	double cx = -dx*norm->x;
	int sgn[5],firstZero=-1;;
	
	// look at corners in order (dx,-dy,-dz),(dx,dy,-dz),(dx,dy,dz),(dx,-dy,dz)
	for(int i=0;i<4;i++)
	{	double cyz = ci[i]*dy*norm->y - ci[i+1]*dz*norm->z;
		
		// find sign (or zero) of n.xc = cyz - cx
		if(DbleEqual(cx,cyz))
		{	// n.xc is zero
			if(firstZero>=0)
			{	// Case 1: found two intersections. Found  diagonal plane, but which diagonal?
				if(i==1)
				{	// edge 0-1
					return sqrt(dx*dx+dz*dz)/(dx*dz);
				}
				else if(i==2)
				{	// diagonal 0-2 or edge 1-2, respectively
					double area = firstZero==0 ? sqrt(dy*dy+dz*dz)/(dy*dz) : sqrt(dx*dx+dy*dy)/(dx*dy);
					return area;
				}
				// corner 3 with 0, 1, or 2
				if(firstZero==0)
				{	// edge 0-3
					return sqrt(dx*dx+dy*dy)/(dx*dy);
				}
				else if(firstZero==1)
				{	// diagonal 1-3
					return sqrt(dy*dy+dz*dz)/(dy*dz);
				}
				// edge 2-3
				return sqrt(dx*dx+dz*dz)/(dx*dz);
			}
			
			// first one
			firstZero = i;
			sgn[i]=0;
		}
		else if(cyz>cx)
			sgn[i]=1;
		else
			sgn[i]=-1;
	}
	
	// store the size
	Vector v1,v2,sz = MakeVector(0.5*dx,0.5*dy,0.5*dz);
	sgn[4]=sgn[0];			// to help some algorithms
	
	// count postive values
	int numPositive=0;
	for(int i=0;i<4;i++)
	{	if(sgn[i]>0) numPositive++;
	}
	
	// Was there one intersection?
	if(firstZero>=0)
	{	// The one corner
		Vector c = MakeVector(0.5*dx,0.5*ci[firstZero]*dy,-0.5*ci[firstZero+1]*dz);
		
		if(numPositive==0 || numPositive==3)
		{	// Case 2: one corner intersection, zero edge intersections
			//cout << "# Case 2: one corner intersection, zero edge intersections" << endl;
			// find any edge except the one ending in the 0 corner
			int c1 = firstZero==0 ? 3 : firstZero-1;
			v1 = IntersectEdges(c1,c1+4,&sz,norm);
		}
		
		else
		{	// Case 3: one corner intersection, one edge intersections
			// a sign change - find it
			int i;
			for(i=0;i<4;i++)
			{	if(sgn[i]*sgn[i+1]<0)
				{	int c1 = i==3 ? 0 : i+1;
					v1 = IntersectEdges(i,c1,&sz,norm);
					break;
				}
			}
			if(i==4) cout << "# failed to find intersected edge" << endl;
		}
				
		// area of parallelogram
		return PGramArea(&c,&v1,norm)/(dx*dy*dz);
	}
	
	// No corner intersections
	
	if(numPositive==4 || numPositive==0)
	{	// Case 4: no intersections or all nodes on +dx face same side of plane
		// It is parallelgram intersecting edges opn the for dx edges; pick any two
		Vector c = IntersectEdges(0,4,&sz,norm);
		v1 = IntersectEdges(1,5,&sz,norm);
		return PGramArea(&c,&v1,norm)/(dx*dy*dz);
	}
	
	else if(numPositive==1 || numPositive==3)
	{	// Case 5: two intersections on adjacent edges will intersect with a hexagon
		// find the corner with differing sign
		int c1 = -1;
		for(int i=0;i<4;i++)
		{	if((numPositive==1 && sgn[i]==1) || (numPositive==3 && sgn[i]==-1))
			{	c1 = i;
				break;
			}
		}
		if(c1<0) cout << "# cound not find the one corner" << endl;
		
		// find two instersection on +dx plane
		int c2 = c1==3 ? 0 : c1+1;				// corner after c1
		v1 = IntersectEdges(c1,c2,&sz,norm);
		c2 = c1==0 ? 3 : c1-1;					// corner before c1
		v2 = IntersectEdges(c1,c2,&sz,norm);
		Vector v3 = IntersectEdges(c2,c2+4,&sz,norm);	// one edge to -dx face
		
		// first parallolgram
		double area = PGramArea(&v2,&v1,norm);
		
		// second parallogram v2-v3 and -v1-v3
		SubVector(&v2,&v3);
		ScaleVector(&v1,-1.);
		SubVector(&v1,&v3);
		Vector cp;
		CrossProduct(&cp,&v1,&v2);
		area += fabs(DotVectors(&cp,norm));
		return area/(dx*dy*dz);
	}
	
	// Case 6: intersects 2 opposite edges: parallelogram  with corners on along dy or dz edges
	if(sgn[0]==sgn[1])
	{	v1 = IntersectEdges(1,2,&sz,norm);
		v2 = IntersectEdges(3,0,&sz,norm);
	}
	else
	{	v1 = IntersectEdges(0,1,&sz,norm);
		v2 = IntersectEdges(2,3,&sz,norm);
	}
	return PGramArea(&v1,&v2,norm)/(dx*dy*dz);
}
#endif

#ifdef MPM_CODE
/****************************************************************
 *  LambertFunctions or ProductLog branches -1 and 0
 ****************************************************************/

/* specfunc/lambert.c
 *
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman
 * Author:  G. Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#define MAX_ITERS 10
#define DBL_EPSILON 2.2204460492503131e-16

static const double c[12] =
{   -1.0,
     2.331643981597124203363536062168,
    -1.812187885639363490240191647568,
     1.936631114492359755363277457668,
    -2.353551201881614516821543561516,
     3.066858901050631912893148922704,
    -4.175335600258177138854984177460,
     5.858023729874774148815053846119,
    -8.401032217523977370984161688514,
     12.250753501314460424,
    -18.100697012472442755,
     27.029044799010561650
};

//  Halley iteration https://en.wikipedia.org/wiki/Halley%27s_method
static double halley_iteration(double x,double w_initial,unsigned int max_iters)
{
    double w = w_initial;
    unsigned int i;

    for(i=0; i<max_iters; i++)
    {   double tol;
        const double e = exp(w);
        const double p = w + 1.0;
        double t = w*e - x;

        if (w > 0)
            t = (t/p)/e;                    // Newton iteration
        else
            t /= e*p - 0.5*(p + 1.0)*t/p;   // Halley iteration

        w -= t;

        tol = 10 * DBL_EPSILON * fmax(fabs(w), 1.0/(fabs(p)*e));

        if(fabs(t) < tol) return w;
    }

    // if too many iterations
    return 0;
}

// series which appears for q near zero;
//  only the argument is different for the different branches
static double series_eval(double r)
{
    const double t_8 = c[8] + r*(c[9] + r*(c[10] + r*c[11]));
    const double t_5 = c[5] + r*(c[6] + r*(c[7]  + r*t_8));
    const double t_1 = c[1] + r*(c[2] + r*(c[3]  + r*(c[4] + r*t_5)));
    return c[0] + r*t_1;
}

// Branch zero of Lambert function for real x (ProductLog(0,x) in Mathematica)
double gsl_sf_lambert_W0(double x)
{
    const double q = x + exp(-1);

    if(x == 0.0)
    {   return 0.0;
    }
    else if(q <= 0.0)
    {   // Strictly speaking, <0 this is an error. Allowed here in case roundoff
        return -1.0;
    }
    else if(q < 1.0e-03)
    {   // series near -1/e in sqrt(q)
        const double r = sqrt(q);
        return series_eval(r);
    }

    // all other values
    double w;

    if (x < 1.0)
    {   // obtain initial approximation from series near x=0;
        const double p = sqrt(2.0 * exp(1) * q);
        w = -1.0 + p*(1.0 + p*(-1.0/3.0 + p*11.0/72.0));
    }
    else
    {   // obtain initial approximation from rough asymptotic
        w = log(x);
        if(x > 3.0) w -= log(w);
    }
    
    return halley_iteration(x, w, MAX_ITERS);
}

// Branh -1 of Lambert function for real x (ProductLog(-1,x) in Mathematica)
// only defined for -1/e < x < 0
double gsl_sf_lambert_Wm1(double x)
{
    if(x > 0.0)
        return gsl_sf_lambert_W0(x);
     else if(x == 0.0)
        return 0.0;
    
    const double q = x + exp(-1);
    double w;

    if (q < 0.0)
    {   // As in the W0 branch above, return some reasonable answer anyway.
        return -1.0;
    }

    if(x < -1.0e-6)
    {   // Obtain initial approximation from series about q = 0,
        // as long as we're not very close to x = 0.
        const double r = -sqrt(q);
        w = series_eval(r);
        if(q < 3.0e-3)
        {   // this approximation is good enough
            return w;
        }
    }
    else
    {   // Obtain initial approximation from asymptotic near zero.
        const double L_1 = log(-x);
        const double L_2 = log(-L_1);
        w = L_1 - L_2 + L_2/L_1;
    }

    return halley_iteration(x, w, MAX_ITERS);
}
#endif

/*
	Lanczos approximation to Gamma[x] (for x real here)
	https://en.wikipedia.org/wiki/Lanczos_approximation
 
	The code uses reflection (thus the if-else structure) is necessary, even though
	it may look strange, as it allows to extend the approximation to values of z where
	z < 0.5, where the Lanczos method is not valid.
*/

static const double pval[8] =
{	676.5203681218851
	,-1259.1392167224028
	,771.32342877765313
	,-176.61502916214059
	,12.507343278686905
	,-0.13857109526572012
	,9.9843695780195716e-6
	,1.5056327351493116e-7
};

double gamma_fxn(double z)
{	if(z<0.5)
	{	// # Reflection formula
		return PI_CONSTANT/(sin(PI_CONSTANT*z)*gamma_fxn(1.-z));
	}
	z -= 1.;
	double x = 0.99999999999980993;
	for(int i=0;i<8;i++)
		x += pval[i]/(z+i+1.);
	double t = z+7.5;		// z+len(pval)-0.5
	double y = sqrt(2.*PI_CONSTANT)*pow(t,z+0.5)*exp(-t)*x;
	return y;
}
