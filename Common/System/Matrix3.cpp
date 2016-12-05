/*************************************************************************************
	Matrix3.cpp
	nairn-mpm-fea

	Created by John Nairn on 12/8/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
 
 A class for holding 3X3 matrices, with some special handling for 2D versions 
    of these matrices (i.e., with xz, yz, zx, and zy terms zero)
 
 Creation
    Matrix() - zero
    Matrix(xx,xy,xz,yx,yy,yz,zx,zy,zz) - full matrix
    Matrix(xx,xy,yz,yy,zz) - 2D matrix

 Methods
    Zero() - sets all elements to zero
    Transpose() - returns new matrix that is transpose of M
    Exponential(kmax) - find exp(M) using kmax terms in the Taylor series expansion
	Scale(double) - multiply all elements by scaling factor
	Inverse() - returns new matrix with inverse of M
 
 Operators
    m1=m2 - matrix copy (built in method)
    m1*m2 or m1*=m2 - matrix multiplication
    m1+m2 or m1+=m2 - matrix addition
    m1-m2 or m1-=m2 - matrix subtraction
    cout << m1 - print matrix elements in one-line string
    m(i,j) or m(i,j)=x - get or set element i,j (0 based or 0 to 2 for each)
        (warning: for efficiency, i and j and not checked to be in bounds)
        (warning: setting 20, 21, 02, or 12 of 2D matrix will not mark it
                    as now a 3D matrix)
 
 Accessors
    set(double) - set all elements to number (and set not a 2D matrix)
		(warning: always marks the matrix as no longer 2D, unles set to 0.)
    set(double a[][3]) - set all elements to those in the array a
		(warning: always marks the matrix as no longer 2D)
    get(double a[][3]) - copy all elements to the array a
    bool getIs2D() - is this a 2D matrix?
    setIs2D(bool) - set 2D marker for this matrix (and zero out elements
                    20, 21, 02, or 12 of the matrix)
    double trace() - trace of the matrix
    double determinant() - determinant of the matrix
    double second_invariant() - second invarient of the matrix
    characteristics(double &c0,double &c1,double &c2) - get the coefficients
        for the characteristic polynomial defined by
        c(z) = det(zI-M) = z^3 - c0 - c1 z - c2 z^2

*************************************************************************************/

#include "stdafx.h"
#include "Matrix3.hpp"
#include <float.h>

// eigenanalysis
inline void dsyev2(double A, double B, double C, double *rt1, double *rt2,
                   double *cs, double *sn);
int dsyevq3(double A[3][3], double Q[3][3], double *w);
void dsytrd3(double A[3][3], double Q[3][3], double *d, double *e);

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define M_SQRT3    1.73205080756887729352744634151   // sqrt(3)

#pragma mark Matrix3:Constructors and Destructor

// empty matrix
Matrix3::Matrix3()
{	
	Zero();
}

// define matrix in all elements specified
Matrix3::Matrix3(double xx,double xy,double xz,
				 double yx,double yy,double yz,
				 double zx,double zy,double zz)
{
	m[0][0] = xx;
	m[0][1] = xy;
	m[0][2] = xz;
	m[1][0] = yx;
	m[1][1] = yy;
	m[1][2] = yz;
	m[2][0] = zx;
	m[2][1] = zy;
	m[2][2] = zz;
	is2D = FALSE;
}

// define matrix in 2D mechanics with xz, yz, zx, and zy equal to zero
Matrix3::Matrix3(double xx,double xy,double yx,double yy,double zz)
{
	m[0][0] = xx;
	m[0][1] = xy;
	m[0][2] = 0.;
	m[1][0] = yx;
	m[1][1] = yy;
	m[1][2] = 0.;
	m[2][0] = 0.;
	m[2][1] = 0.;
	m[2][2] = zz;
	is2D = TRUE;
}

// define matrix as outer product of two vectors (x1,y1,z1) and (x2,y2,z2)
Matrix3::Matrix3(double x1,double y1,double z1,double x2,double y2,double z2)
{
	m[0][0] = x1*x2;
	m[0][1] = x1*y2;
	m[0][2] = x1*z2;
	m[1][0] = y1*x2;
	m[1][1] = y1*y2;
	m[1][2] = y1*z2;
	m[2][0] = z1*x2;
	m[2][1] = z1*y2;
	m[2][2] = z1*z2;
	is2D = FALSE;
}

#pragma mark Matrix3:methods

// zero the matrix
void Matrix3::Zero(void)
{	set(0.);
	is2D = TRUE;
}

// get the transpose
Matrix3 Matrix3::Transpose(void) const
{   Matrix3 mT(m[0][0],m[1][0],m[2][0],
               m[0][1],m[1][1],m[2][1],
               m[0][2],m[1][2],m[2][2]);
    mT.setIs2D(is2D);
    return mT;
}

// Form the triple product R.M.R^T in one step (R need not be a rotation matrix)
Matrix3 Matrix3::RMRT(Matrix3 &R) const
{	if(is2D && R.getIs2D())
	{	Matrix3 mT(R(0,0)*R(0,0)*m[0][0] + R(0,1)*R(0,0)*m[1][0] + R(0,0)*R(0,1)*m[0][1] + R(0,1)*R(0,1)*m[1][1],
				   R(0,0)*R(1,0)*m[0][0] + R(0,1)*R(1,0)*m[1][0] + R(0,0)*R(1,1)*m[0][1] + R(0,1)*R(1,1)*m[1][1],
				   R(1,0)*R(0,0)*m[0][0] + R(1,1)*R(0,0)*m[1][0] + R(1,0)*R(0,1)*m[0][1] + R(1,1)*R(0,1)*m[1][1],
				   R(1,0)*R(1,0)*m[0][0] + R(1,1)*R(1,0)*m[1][0] + R(1,0)*R(1,1)*m[0][1] + R(1,1)*R(1,1)*m[1][1],
				   R(2,2)*R(2,2)*m[2][2]);
		return mT;
	}
	else
	{	Matrix3 mT;
		int i,j;
		for(i=0;i<3;i++)
		{	for(j=0;j<3;j++)
			{	double mij = m[i][j];
				mT(0,0) += R(0,i)*R(0,j)*mij;
				mT(0,1) += R(0,i)*R(1,j)*mij;
				mT(0,2) += R(0,i)*R(2,j)*mij;
				mT(1,0) += R(1,i)*R(0,j)*mij;
				mT(1,1) += R(1,i)*R(1,j)*mij;
				mT(1,2) += R(1,i)*R(2,j)*mij;
				mT(2,0) += R(2,i)*R(0,j)*mij;
				mT(2,1) += R(2,i)*R(1,j)*mij;
				mT(2,2) += R(2,i)*R(2,j)*mij;
			}
		}
		mT.setIs2D(false);
		return mT;
	}
}

// Form the triple product R.M.R^T in one step (R need not be a rotation matrix)
Matrix3 Matrix3::RTMR(Matrix3 &R) const
{	if(is2D && R.getIs2D())
	{	Matrix3 mT(R(0,0)*R(0,0)*m[0][0] + R(1,0)*R(0,0)*m[1][0] + R(0,0)*R(1,0)*m[0][1] + R(1,0)*R(1,0)*m[1][1],
				   R(0,0)*R(0,1)*m[0][0] + R(1,0)*R(0,1)*m[1][0] + R(0,0)*R(1,1)*m[0][1] + R(1,0)*R(1,1)*m[1][1],
				   R(0,1)*R(0,0)*m[0][0] + R(1,1)*R(0,0)*m[1][0] + R(0,1)*R(1,0)*m[0][1] + R(1,1)*R(1,0)*m[1][1],
				   R(0,1)*R(0,1)*m[0][0] + R(1,1)*R(0,1)*m[1][0] + R(0,1)*R(1,1)*m[0][1] + R(1,1)*R(1,1)*m[1][1],
				   R(2,2)*R(2,2)*m[2][2]);
		return mT;
	}
	else
	{	Matrix3 mT;
		int i,j;
		for(i=0;i<3;i++)
		{	for(j=0;j<3;j++)
			{	double mij = m[i][j];
				mT(0,0) += R(i,0)*R(j,0)*mij;
				mT(0,1) += R(i,0)*R(j,1)*mij;
				mT(0,2) += R(i,0)*R(j,2)*mij;
				mT(1,0) += R(i,1)*R(j,0)*mij;
				mT(1,1) += R(i,1)*R(j,1)*mij;
				mT(1,2) += R(i,1)*R(j,2)*mij;
				mT(2,0) += R(i,2)*R(j,0)*mij;
				mT(2,1) += R(i,2)*R(j,1)*mij;
				mT(2,2) += R(i,2)*R(j,2)*mij;
			}
		}
		mT.setIs2D(false);
		return mT;
	}
}

// exponential of matrix to kmax terms
Matrix3 Matrix3::Exponential(int kmax) const
{	
	if(is2D)
	{	// first term
		if(kmax==1)
		{	return Matrix3(1. + m[0][0], m[0][1],
						   m[1][0], 1. + m[1][1], 1. + m[2][2]);
		}
		
		// kmax is 2 or higher
        int k;
		double c0 = m[0][1]*m[1][0] - m[0][0]*m[1][1];		// -det(A)
		double c1 = m[0][0] + m[1][1];						// Tr(A)
        double beta0 = 0., beta1 = 1.;
		double alpha0 = 1., alpha1 = 1.;
		double betaz = m[2][2], ezz = 1. + betaz;
        double factor, temp;
		
        for(k=2;k<=kmax;k++)
        {	// update beta using beta(k,0) = (1/k) c0 beta(k-1,1)
            //               and beta(k,1) = (1/k) (c1 beta(k-1,1) + beta(k-1,0))
            factor = 1/(double)k;
            temp = beta1;
            beta1 = factor*(c1*temp + beta0);
            beta0 = factor*c0*temp;
            betaz *= factor*m[2][2];
            alpha0 += beta0;
            alpha1 += beta1;
            ezz += betaz;
        }
			
		// return alpha0*I + alpha1*m
		return Matrix3(alpha0 + alpha1*m[0][0], alpha1*m[0][1],
						   alpha1*m[1][0], alpha0 + alpha1*m[1][1], ezz);
	}
	
	// done if only 1 term
	if(kmax==1)
	{	return Matrix3(1. + m[0][0], m[0][1], m[0][2],
						m[1][0], 1. + m[1][1], m[1][2],
						m[2][0], m[2][1], 1. + m[2][2]);
	}
	
	// get square of this matrix
	Matrix3 m2 = *this;
	m2 *= m2;
	
	// just two terms
	if(kmax==2)
	{	return Matrix3(1. + m[0][0] + 0.5*m2(0,0), m[0][1] + 0.5*m2(0,1), m[0][2] + 0.5*m2(0,2),
					   m[1][0] + 0.5*m2(1,0), 1. + m[1][1] + 0.5*m2(1,1), m[1][2] + 0.5*m2(1,2),
					   m[2][0] + 0.5*m2(2,0), m[2][1] + 0.5*m2(2,1), 1. + m[2][2] + 0.5*m2(2,2));
	}
	
	// kmax is 3 or more
    int k;
	double c0,c1,c2,factor,temp;
	characteristics(c0,c1,c2);
    double beta0 = 0.,beta1 = 0.,beta2 = 0.5;
	double alpha0 = 1.,alpha1 = 1.,alpha2 = 0.5;
    
    for(k=3;k<=kmax;k++)
    {	// update beta using beta(k,0) = (1/k) c0 beta(k-1,2)
        //			     and beta(k,1) = (1/k) (c1 beta(k-1,2) + beta(k-1,0))
        //			     and beta(k,2) = (1/k) (c2 beta(k-1,2) + beta(k-1,1))
        factor = 1./(double)k;
        temp = beta2;
        beta2 = factor*(c2*temp + beta1);
        beta1 = factor*(c1*temp + beta0);
        beta0 = factor*c0*temp;
        alpha0 += beta0;
        alpha1 += beta1;
        alpha2 += beta2;
    }

	// return alpha0*I + alpha1*m + alpha2*m2
	return Matrix3(alpha0 + alpha1*m[0][0] + alpha2*m2(0,0), alpha1*m[0][1] + alpha2*m2(0,1), alpha1*m[0][2] + alpha2*m2(0,2),
				   alpha1*m[1][0] + alpha2*m2(1,0), alpha0 + alpha1*m[1][1] + alpha2*m2(1,1), alpha1*m[1][2] + alpha2*m2(1,2),
				   alpha1*m[2][0] + alpha2*m2(2,0), alpha1*m[2][1] + alpha2*m2(2,1), alpha0 + alpha1*m[2][2] + alpha2*m2(2,2));
}

// matrix times a vector
Vector Matrix3::Times(Vector *v) const
{   Vector mv;
	mv.x = m[0][0]*v->x + m[1][0]*v->y + m[2][0]*v->z;
    mv.y = m[0][1]*v->x + m[1][1]*v->y + m[2][1]*v->z;
	mv.z = m[0][2]*v->x + m[1][2]*v->y + m[2][2]*v->z;
    return mv;
}

// scale all elements by factor
void Matrix3::Scale(double factor)
{
	m[0][0] *= factor;
	m[0][1] *= factor;
	m[1][0] *= factor;
	m[1][1] *= factor;
	m[2][2] *= factor;
	if(!is2D)
	{	m[0][2] *= factor;
		m[1][2] *= factor;
		m[2][0] *= factor;
		m[2][1] *= factor;
	}
}

// scale x-y elements, but not z direction
void Matrix3::Scale2D(double factor)
{
	m[0][0] *= factor;
	m[0][1] *= factor;
	m[1][0] *= factor;
	m[1][1] *= factor;
}

// inner product term by term of this matrix
double Matrix3::DotProduct(void)
{	double dot = m[0][0]*m[0][0] + m[1][1]*m[1][1] + m[2][2]*m[2][2]
					+ m[0][1]*m[0][1] + m[1][0]*m[1][0];
	if(!is2D)
	{	dot += m[0][2]*m[0][2] + m[2][0]*m[2][0]
					+ m[1][2]*m[1][2] + m[2][1]*m[2][1];
	}
	return dot;
}

// Return inverse in a new matrix3
Matrix3 Matrix3::Inverse(void) const
{	// special case for 2D
	if(is2D)
	{	double subdet = 1./(m[0][0]*m[1][1] - m[1][0]*m[0][1]);
		Matrix3 inv(m[1][1]*subdet,-m[0][1]*subdet,-m[1][0]*subdet,m[0][0]*subdet,1./m[2][2]);
		return inv;
	}
	
	// general 3D
	Matrix3 inv3(m[1][1]*m[2][2]-m[2][1]*m[1][2], -(m[0][1]*m[2][2]-m[2][1]*m[0][2]), m[1][0]*m[1][2]-m[1][1]*m[0][2],
				 -(m[1][0]*m[2][2]-m[2][0]*m[1][2]), m[0][0]*m[2][2]-m[2][0]*m[0][2], -(m[0][0]*m[1][2]-m[1][0]*m[0][2]),
				   m[1][0]*m[2][1]-m[2][0]*m[1][1], -(m[0][0]*m[2][1]-m[2][0]*m[0][1]), m[0][0]*m[1][1]-m[1][0]*m[0][1]);
	inv3.Scale(1./determinant());
	return inv3;
}

// Find Eigenvalues of positive definite matrix return in a vector
// Assumes all are real and positive (any other will fail)
Vector Matrix3::Eigenvalues(void) const
{
    Vector lam;
    
    if(is2D)
    {   // solving x^2 + bx + c = 0 (see Numerical Recipes in C, page 156)
        double b = -(m[0][0]+m[1][1]);
        double c = m[0][0]*m[1][1] - m[1][0]*m[0][1];
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
        lam.z = m[2][2];
    }
    else
    {	double mm, c1, c0;
		
		// Determine coefficients of characteristic poynomial. We write
		//       | a   d   f  |
		//  m =  | d*  b   e  |
		//       | f*  e*  c  |
		double de = m[0][1] * m[1][2];                                  // d * e
		double dd = SQR(m[0][1]);                                       // d^2
		double ee = SQR(m[1][2]);                                       // e^2
		double ff = SQR(m[0][2]);                                       // f^2
		mm  = m[0][0] + m[1][1] + m[2][2];
		c1 = (m[0][0]*m[1][1] + m[0][0]*m[2][2] + m[1][1]*m[2][2]) 
				- (dd + ee + ff);										// a*b + a*c + b*c - d^2 - e^2 - f^2
		c0 = m[2][2]*dd + m[0][0]*ee + m[1][1]*ff - m[0][0]*m[1][1]*m[2][2]
				- 2.0 * m[0][2]*de;										// c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
		
		double p, sqrt_p, q, c, s, phi;
		p = SQR(mm) - 3.0*c1;
		q = mm*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
		sqrt_p = sqrt(fabs(p));
		
		phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
		phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
		
		c = sqrt_p*cos(phi);
		s = (1.0/M_SQRT3)*sqrt_p*sin(phi);
		
		lam.y  = (1.0/3.0)*(mm - c);
		lam.z  = lam.y + s;
		lam.x  = lam.y + c;
		lam.y -= s;
    }
    
    return lam;
}

// Warning: This method only works for symmetric, Hermitian matrices
// Find Eigenvectors and return three orthogonal vectors in columns of the matrix
// Input is Vector with eigenvalues of this matrix (not needed for 3D matrices)
// No error checking (that is symmetric and no error in algorithm)
Matrix3 Matrix3::Eigenvectors(Vector &Eigenvals) const
{
    Matrix3 Eigenvecs;
	
    if(is2D)
	{	double cs,sn;
		dsyev2(m[0][0],m[0][1],m[1][1],&(Eigenvals.x),&(Eigenvals.y),&cs,&sn);
		Eigenvecs.set(cs,-sn,sn,cs,1.);
    }
    
    else
	{	// set to 3D matrix
		double w[3];
		double Q[3][3];
		double A[3][3];
		A[0][0] = m[0][0];
		A[0][1] = m[0][1];
		A[0][2] = m[0][2];
		A[1][1] = m[1][1];
		A[1][2] = m[1][2];
		A[2][2] = m[2][2];
		dsyevq3(A, Q, w);
		Eigenvecs.set(Q);
		Eigenvals.x = w[0];
		Eigenvals.y = w[1];
		Eigenvals.z = w[2];
	}

    return Eigenvecs;
}

// Polar decomponsition of F through right stretch matrix
//      F = RU = (RQ) Lambda QT
// The target matrix is assumed to be F
// Function returns U and optionally R = F U^-1 (if pointer is not NULL)
//		and optionally (lam1,lam2,lam3) in stretches (if not NULL)
// It does not get Q, but if needed, they are eigenvectors of the
//		returned U matrix
Matrix3 Matrix3::RightDecompose(Matrix3 *R,Vector *stretches) const
{
    if(is2D)
    {   // 2D has simple formulae for R = ((Fsum,Fdif),(-Fdif,Fsum))
        double Fsum = m[0][0]+m[1][1];
        double Fdif = m[0][1]-m[1][0];
        double denom = sqrt(Fsum*Fsum+Fdif*Fdif);
        Fsum /= denom;
        Fdif /= denom;
        
        // U is R^T * F
        Matrix3 U(Fsum*m[0][0]-Fdif*m[1][0],Fsum*m[0][1]-Fdif*m[1][1],
                  Fdif*m[0][0]+Fsum*m[1][0],Fdif*m[0][1]+Fsum*m[1][1],m[2][2]);
        
        // if R pointer not NULL, return it too
        if(R!=NULL)
        {   R->set(Fsum,Fdif,-Fdif,Fsum,1.);
        }
        
		// Return Eigenvalues of U if asked
		if(stretches!=NULL)
		{	*stretches = U.Eigenvalues();
		}
		
        return U;
    }
    
    // rest is for 3D matrix
    
    // Get C=F^TF and C^2
    Matrix3 C = Transpose()*(*this);
    Matrix3 C2 = C*C;
    
    // Eigenvalues of C are lamda^2
    Vector Eigenvals = C.Eigenvalues();
    double lam1 = sqrt(Eigenvals.x);
    double lam2 = sqrt(Eigenvals.y);
    double lam3 = sqrt(Eigenvals.z);
    
    // invariants of U
    double i1 = lam1+lam2+lam3;
    double i2 = lam1*lam2+lam1*lam3+lam2*lam3;
    double i3 = lam1*lam2*lam3;
    
    // set coefficients
    double d1 = 1./(i1*i2-i3);
    double c2 = -d1;                    // coefficient of C2
    double c1 = (i1*i1-i2)*d1;          // coefficient of C
    double cI = i1*i3*d1;               // coefficient of I
    
    // Get U = d1*(-C^2 + (i1*i1-i2)*C + i1*i3*I) (symmetric)
	double U01 = c2*C2(0,1)+c1*C(0,1);
	double U02 = c2*C2(0,2)+c1*C(0,2);
	double U12 = c2*C2(1,2)+c1*C(1,2);
    Matrix3 U(c2*C2(0,0)+c1*C(0,0)+cI, U01,                     U02,
			  U01,                     c2*C2(1,1)+c1*C(1,1)+cI, U12,
			  U02,                     U12,                     c2*C2(2,2)+c1*C(2,2)+cI);
    
    // if R pointer not NULL, find R too
    if(R!=NULL)
    {   c1 = 1/i3;                      // coefficient of C
        double cU = -i1*c1;             // coefficient of U
        cI = i2*c1;                     // coefficient of I
        // Get Uinv = (1/i3)*(C - i1*U + i2*I) (symmetric)
		U01 = c1*C(0,1)+cU*U(0,1);
		U02 = c1*C(0,2)+cU*U(0,2);
		U12 = c1*C(1,2)+cU*U(1,2);
        Matrix3 Uinv(c1*C(0,0)+cU*U(0,0)+cI, U01,                    U02,
                     U01,                    c1*C(1,1)+cU*U(1,1)+cI, U12,
                     U02,                    U12,                    c1*C(2,2)+cU*U(2,2)+cI);
        
        // R = FU^-1
        *R = (*this)*Uinv;
    }
	
	// Return Eigenvalues of U if asked
	if(stretches!=NULL)
	{	stretches->x = lam1;
		stretches->y = lam2;
		stretches->z = lam3;
	}
    
    return U;
}

// Polar decomponsition of F through left stretch matrix
//      F = VR = Q Lambda (QTR)
// The target matrix is assumed to be F
// Function returns V and optionally R = V^-1 F (if pointer is not NULL)
//		and optionally (lam1,lam2,lam3) in stretches (if not NULL)
// It does not get Q, but if needed, they are eigenvectors of the
//		returned V matrix
Matrix3 Matrix3::LeftDecompose(Matrix3 *R,Vector *stretches) const
{
    if(is2D)
    {   // 2D has simple formulae for R = ((Fsum,Fdif),(-Fdif,Fsum))
        double Fsum = m[0][0]+m[1][1];
        double Fdif = m[0][1]-m[1][0];
        double denom = sqrt(Fsum*Fsum+Fdif*Fdif);
        Fsum /= denom;
        Fdif /= denom;
        
        // V is F* R*T
        Matrix3 V(m[0][0]*Fsum+m[0][1]*Fdif,-m[0][0]*Fdif+m[0][1]*Fsum,
				  m[1][0]*Fsum+m[1][1]*Fdif,-m[1][0]*Fdif+m[1][1]*Fsum,m[2][2]);
        
        // if R pointer not NULL, return it too
        if(R!=NULL)
        {   R->set(Fsum,Fdif,-Fdif,Fsum,1.);
        }
        
		// Return Eigenvalues of V if asked
		if(stretches!=NULL)
		{	*stretches = V.Eigenvalues();
		}
		
        return V;
    }
    
    // rest is for 3D matrix
    
    // Get B=FF^T and B^2
    Matrix3 B = (*this)*Transpose();
    Matrix3 B2 = B*B;
    
    // Eigenvalues of B are lamda^2
    Vector Eigenvals = B.Eigenvalues();
    double lam1 = sqrt(Eigenvals.x);
    double lam2 = sqrt(Eigenvals.y);
    double lam3 = sqrt(Eigenvals.z);
    
    // invariants of V
    double i1 = lam1+lam2+lam3;
    double i2 = lam1*lam2+lam1*lam3+lam2*lam3;
    double i3 = lam1*lam2*lam3;
    
    // set coefficients
    double d1 = 1./(i1*i2-i3);
    double c2 = -d1;                    // coefficient of B2
    double c1 = (i1*i1-i2)*d1;          // coefficient of B
    double cI = i1*i3*d1;               // coefficient of I
    
    // Get V = (1/d1)*(-B^2 + (i1*i1-i2)*B + i1*i3*I)
    Matrix3 V(c2*B2(0,0)+c1*B(0,0)+cI, c2*B2(0,1)+c1*B(0,1),    c2*B2(0,2)+c1*B(0,2),
			  c2*B2(1,0)+c1*B(1,0),    c2*B2(1,1)+c1*B(1,1)+cI, c2*B2(1,2)+c1*B(1,2),
			  c2*B2(2,0)+c1*B(2,0),    c2*B2(2,1)+c1*B(2,1),    c2*B2(2,2)+c1*B(2,2)+cI);
    
    // if R pointer not NULL, find R too
    if(R!=NULL)
    {   c1 = 1/i3;                      // coefficient of B
        double cV = -i1*c1;             // coefficient of V
        cI = i2*c1;                     // coefficient of I
        // Get Vinv = (1/i3)*(B - i1*V + i2*I)
        Matrix3 Vinv(c1*B(0,0)+cV*V(0,0)+cI, c1*B(0,1)+cV*V(0,1),     c1*B(0,2)+cV*V(0,2),
                     c1*B(1,0)+cV*V(1,0),    c1*B(1,1)+cV*V(1,1)+cI,  c1*B(1,2)+cV*V(1,2),
                     c1*B(2,0)+cV*V(2,0),    c1*B(2,1)+cV*V(2,1),     c1*B(2,2)+cV*V(2,2)+cI);
        
        // R = V^-1 F
        *R = Vinv*(*this);
    }
	
	// Return Eigenvalues of U if asked
	if(stretches!=NULL)
	{	stretches->x = lam1;
		stretches->y = lam2;
		stretches->z = lam3;
	}
    
    return V;
}

/*
// Get Rotation Matrix using target matrix F
// Eigenvals are for matrix X and matrix U calculated by RightStretch()
// Not currently used - verify before use
Matrix3 Matrix3::Rotation(const Vector &Eigenvals,const Matrix3 &U) const
{
    // Matrix m = F here;
    
    double U1[3][3];
    
    double lam1 = sqrt(Eigenvals.x);
    double lam2 = sqrt(Eigenvals.y);
    double lam3 = sqrt(Eigenvals.z);
    double i1=lam1+lam2+lam3;
    double i2=lam1*lam2+lam1*lam3+lam2*lam3;
    double i3=lam1*lam2*lam3;
    
    //   Calculation of inv U
    //   U1=1/i3*(C-i1*U+i2*eye(3));
    U1[0][0] = 1/i3*(m[0][0]*m[0][0]+m[1][0]*m[1][0]+m[2][0]*m[2][0]-i1*U(0,0)+i2);
    U1[0][1] = 1/i3*(m[0][0]*m[0][1]+m[1][0]*m[1][1]+m[2][0]*m[2][1]-i1*U(0,1));
    U1[0][2] = 1/i3*(m[0][0]*m[0][2]+m[1][0]*m[1][2]+m[2][0]*m[2][2]-i1*U(0,2));
    
    U1[1][0] = 1/i3*(m[0][1]*m[0][0]+m[1][1]*m[1][0]+m[2][1]*m[2][0]-i1*U(1,0));
    U1[1][1] = 1/i3*(m[0][1]*m[0][1]+m[1][1]*m[1][1]+m[2][1]*m[2][1]-i1*U(1,1)+i2);
    U1[1][2] = 1/i3*(m[0][1]*m[0][2]+m[1][1]*m[1][2]+m[2][1]*m[2][2]-i1*U(1,2));
    
    U1[2][0] = 1/i3*(m[0][2]*m[0][0]+m[1][2]*m[1][0]+m[2][2]*m[2][0]-i1*U(2,0));
    U1[2][1] = 1/i3*(m[0][2]*m[0][1]+m[1][2]*m[1][1]+m[2][2]*m[2][1]-i1*U(2,1));
    U1[2][2] = 1/i3*(m[0][2]*m[0][2]+m[1][2]*m[1][2]+m[2][2]*m[2][2]-i1*U(2,2)+i2);
    
    Matrix3 Ui;
    Ui.set(U1);
    cout << Ui*U << endl;
    
    // R = F Uinv
    Matrix3 rotate(m[0][0]*U1[0][0]+m[0][1]*U1[1][0]+m[0][2]*U1[2][0],
                     m[0][0]*U1[0][1]+m[0][1]*U1[1][1]+m[0][2]*U1[2][1],
                     m[0][0]*U1[0][2]+m[0][1]*U1[1][2]+m[0][2]*U1[2][2],
                     m[1][0]*U1[0][0]+m[1][1]*U1[1][0]+m[1][2]*U1[2][0],
                     m[1][0]*U1[0][1]+m[1][1]*U1[1][1]+m[1][2]*U1[2][1],
                     m[1][0]*U1[0][2]+m[1][1]*U1[1][2]+m[1][2]*U1[2][2],
                     m[2][0]*U1[0][0]+m[2][1]*U1[1][0]+m[2][2]*U1[2][0],
                     m[2][0]*U1[0][1]+m[2][1]*U1[1][1]+m[2][2]*U1[2][1],
                     m[2][0]*U1[0][2]+m[2][1]*U1[1][2]+m[2][2]*U1[2][2]);
    
    return rotate;
}
*/

// Get 6X6 rotation matrix for rotation of stress tensor (called R_{\sigma} in my notes)
// To get inverse, call this method for transpose (and not just transpose the 6X6)
void Matrix3::GetRStress(double r[][6]) const
{
	r[0][0] = m[0][0]*m[0][0];
	r[0][1] = m[0][1]*m[0][1];
	r[0][2] = m[0][2]*m[0][2];
	r[0][3] = 2.*m[0][1]*m[0][2];
	r[0][4] = 2.*m[0][0]*m[0][2];
	r[0][5] = 2.*m[0][0]*m[0][1];
	
	r[1][0] = m[1][0]*m[1][0];
	r[1][1] = m[1][1]*m[1][1];
	r[1][2] = m[1][2]*m[1][2];
	r[1][3] = 2.*m[1][1]*m[1][2];
	r[1][4] = 2.*m[1][0]*m[1][2];
	r[1][5] = 2.*m[1][1]*m[1][0];

	r[2][0] = m[2][0]*m[2][0];
	r[2][1] = m[2][1]*m[2][1];
	r[2][2] = m[2][2]*m[2][2];
	r[2][3] = 2.*m[2][2]*m[2][1];
	r[2][4] = 2.*m[2][2]*m[2][0];
	r[2][5] = 2.*m[2][0]*m[2][1];
	
	r[3][0] = m[1][0]*m[2][0];
	r[3][1] = m[1][1]*m[2][1];
	r[3][2] = m[2][2]*m[1][2];
	r[3][3] = m[1][2]*m[2][1] + m[1][1]*m[2][2];
	r[3][4] = m[1][2]*m[2][0] + m[1][0]*m[2][2];
	r[3][5] = m[1][1]*m[2][0] + m[1][0]*m[2][1];

	r[4][0] = m[0][0]*m[2][0];
	r[4][1] = m[0][1]*m[2][1];
	r[4][2] = m[2][2]*m[0][2];
	r[4][3] = m[0][2]*m[2][1] + m[0][1]*m[2][2];
	r[4][4] = m[0][2]*m[2][0] + m[0][0]*m[2][2];
	r[4][5] = m[0][1]*m[2][0] + m[0][0]*m[2][1];
	
	r[5][0] = m[0][0]*m[1][0];
	r[5][1] = m[1][1]*m[0][1];
	r[5][2] = m[0][2]*m[1][2];
	r[5][3] = m[0][2]*m[1][1] + m[0][1]*m[1][2];
	r[5][4] = m[0][2]*m[1][0] + m[0][0]*m[1][2];
	r[5][5] = m[0][1]*m[1][0] + m[0][0]*m[1][1];
}

// Get 6X6 rotation matrix for rotation of strain tensor (called R_{\varepsilon} in my notes)
// To get inverse, call this method for transpose (and not just transpose the 6X6)
// Note since Rstrain = (Rstress^-1)^T, take RStress terms, interchange order of all m terms
//    to get Rstress^-1 then interchange r terms to get the transpose
void Matrix3::GetRStrain(double r[][6]) const
{
	r[0][0] = m[0][0]*m[0][0];
	r[1][0] = m[1][0]*m[1][0];
	r[2][0] = m[2][0]*m[2][0];
	r[3][0] = 2.*m[1][0]*m[2][0];
	r[4][0] = 2.*m[0][0]*m[2][0];
	r[5][0] = 2.*m[0][0]*m[1][0];
	
	r[0][1] = m[0][1]*m[0][1];
	r[1][1] = m[1][1]*m[1][1];
	r[2][1] = m[2][1]*m[2][1];
	r[3][1] = 2.*m[1][1]*m[2][1];
	r[4][1] = 2.*m[0][1]*m[2][1];
	r[5][1] = 2.*m[1][1]*m[0][1];
	
	r[0][2] = m[0][2]*m[0][2];
	r[1][2] = m[1][2]*m[1][2];
	r[2][2] = m[2][2]*m[2][2];
	r[3][2] = 2.*m[2][2]*m[1][2];
	r[4][2] = 2.*m[2][2]*m[0][2];
	r[5][2] = 2.*m[0][2]*m[1][2];
	
	r[0][3] = m[0][1]*m[0][2];
	r[1][3] = m[1][1]*m[1][2];
	r[2][3] = m[2][2]*m[2][1];
	r[3][3] = m[1][2]*m[2][1] + m[1][1]*m[2][2];
	r[4][3] = m[2][1]*m[0][2] + m[0][1]*m[2][2];
	r[5][3] = m[1][1]*m[0][2] + m[0][1]*m[1][2];
	
	r[0][4] = m[0][0]*m[0][2];
	r[1][4] = m[1][0]*m[1][2];
	r[2][4] = m[2][2]*m[2][0];
	r[3][4] = m[2][0]*m[1][2] + m[1][0]*m[2][2];
	r[4][4] = m[0][2]*m[2][0] + m[0][0]*m[2][2];
	r[5][4] = m[1][0]*m[0][2] + m[0][0]*m[1][2];
	
	r[0][5] = m[0][0]*m[0][1];
	r[1][5] = m[1][1]*m[1][0];
	r[2][5] = m[2][0]*m[2][1];
	r[3][5] = m[2][0]*m[1][1] + m[1][0]*m[2][1];
	r[4][5] = m[2][0]*m[0][1] + m[0][0]*m[2][1];
	r[5][5] = m[0][1]*m[1][0] + m[0][0]*m[1][1];
}

#pragma mark Matrix3:operators

// +=
Matrix3 &Matrix3::operator+=(const Matrix3 &rhs)
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			m[i][j] += rhs(i,j);
	}
	if(!rhs.getIs2D()) is2D = FALSE;
	return *this;
}

// +
const Matrix3 Matrix3::operator+(const Matrix3 &rhs) const
{	Matrix3 result = *this;
	result += rhs;
	return result;
}

// -=
Matrix3 &Matrix3::operator-=(const Matrix3 &rhs)
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			m[i][j] -= rhs(i,j);
	}
	if(!rhs.getIs2D()) is2D = FALSE;
	return *this;
}

// -
const Matrix3 Matrix3::operator-(const Matrix3 &rhs) const
{	Matrix3 result = *this;
	result -= rhs;
	return result;
}

// *=
Matrix3 &Matrix3::operator*=(const Matrix3 &rhs)
{	double cm[3][3],rm[3][3];
	get(cm);
	rhs.get(rm);
	if(is2D && rhs.getIs2D())
	{	set(0.);
		m[0][0] = cm[0][0]*rm[0][0] + cm[0][1]*rm[1][0];
		m[0][1] = cm[0][0]*rm[0][1] + cm[0][1]*rm[1][1];
		m[1][0] = cm[1][0]*rm[0][0] + cm[1][1]*rm[1][0];
		m[1][1] = cm[1][0]*rm[0][1] + cm[1][1]*rm[1][1];
		m[2][2] = cm[2][2]*rm[2][2];
	}
	else
	{	int i,j,k;
		for(i=0;i<3;i++)
		{	for(j=0;j<3;j++)
			{	m[i][j] = 0.;
				for(k=0;k<3;k++)
					m[i][j] += cm[i][k]*rm[k][j];
			}
		}
		is2D = FALSE;
	}
	return *this;
}

// *
const Matrix3 Matrix3::operator*(const Matrix3 &rhs) const
{	Matrix3 result = *this;
	result *= rhs;
	return result;
}

// printing
ostream &operator<<(ostream &os, const Matrix3 &mat)
{	
	os << "[(" << mat(0,0) << "," << mat(0,1) << "," << mat(0,2) << "),("
			<< mat(1,0) << "," << mat(1,1) << "," << mat(1,2) << "),("
			<< mat(2,0) << "," << mat(2,1) << "," << mat(2,2) << ")]" ;
	return os;
}

// elements (0 based) (no checking of limits)
double &Matrix3::operator()(unsigned row,unsigned col)
{	return m[row][col];
}

// elements (0 based) of const object (no checking of limits)
const double &Matrix3::operator()(unsigned row,unsigned col) const
{	return m[row][col];
}


#pragma mark Matrix3:accessors

// all element to a constant (and set to not 2D)
void Matrix3::set(double constant)
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			m[i][j] = constant;
	}
	is2D = constant==0. ? TRUE : FALSE ;
}

// set one element (0 based indices)
void Matrix3::set(int row,int col,double value) { m[row][col] = value; }

// copy this object to an array
void Matrix3::get(double c[][3]) const
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			c[i][j] = m[i][j];
	}
}

// get one element (0 based indices)
void Matrix3::get(int row,int col,double *value) const { *value = m[row][col]; }

// get one element using contraction notation of 0 to 5 for xx,yy,zz,yz,xz,xy, 6 to 8 for zy,zx,yx
// applying scalling if an off-axis value (e.g., 2 to get engineering strain)
// uses top half of the matrix
double Matrix3::get(int tid,double scale) const
{
	switch(tid)
	{	case 0:
			return m[0][0];
		case 1:
			return m[1][1];
		case 2:
			return m[2][2];
		case 3:
			return scale*m[1][2];
		case 4:
			return scale*m[0][2];
		case 5:
			return scale*m[0][1];
		case 6:
			return scale*m[2][1];
		case 7:
			return scale*m[2][0];
		case 8:
			return scale*m[1][0];
		default:
			break;
	}
	return 0.;
}

// set this object to contents of an array and makes it 3D
void Matrix3::set(double c[][3])
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			m[i][j] = c[i][j];
	}
	is2D = FALSE;
}

// set this object to five elements and make it 2D
void Matrix3::set(double c0,double c1,double c2,double c3,double c4)
{	m[0][0] = c0;
    m[0][1] = c1;
    m[1][0] = c2;
    m[1][1] = c3;
    m[2][2] = c4;
    m[0][2] = 0.;
    m[1][2] = 0.;
    m[2][0] = 0.;
    m[2][1] = 0.;
    is2D = TRUE;
}

// set this object to nine elements and make it 3D
void Matrix3::set(double c00,double c01,double c02,double c10,double c11,double c12,double c20,double c21,double c22)
{	m[0][0] = c00;
    m[0][1] = c01;
    m[0][2] = c02;
    m[1][0] = c10;
    m[1][1] = c11;
    m[1][2] = c12;
    m[2][0] = c20;
    m[2][1] = c21;
    m[2][2] = c22;
    is2D = FALSE;
}

// swap two columns
void Matrix3::SwapColumns(int i,int j)
{
	double temp = m[0][i];
	m[0][i] = m[0][j];
	m[0][j] = temp;
	temp = m[1][i];
	m[1][i] = m[1][j];
	m[1][j] = temp;
	temp = m[2][i];
	m[2][i] = m[2][j];
	m[2][j] = temp;
	if(i==2 || j==2) is2D = false;
}

// is this a 2D matrix?
bool Matrix3::getIs2D(void) const { return is2D; }
void Matrix3::setIs2D(bool new2D)
{   is2D = new2D;
    if(is2D)
    {   m[0][2] = 0.;
        m[1][2] = 0.;
        m[2][0] = 0.;
        m[2][1] = 0.;
    }
}

// set to other matrix
void Matrix3::set(Matrix3 mat)
{
	m[0][0] = mat(0,0);
	m[0][1] = mat(0,1);
	m[1][0] = mat(1,0);
	m[1][1] = mat(1,1);
	m[2][2] = mat(2,2);
	if(mat.getIs2D())
	{	if(!is2D)
		{	is2D = true;
			m[0][2] = 0.;
			m[1][2] = 0.;
			m[2][0] = 0.;
			m[2][1] = 0.;
		}
	}
	else
	{	is2D = false;
		m[0][2] = mat(0,2);
		m[1][2] = mat(1,2);
		m[2][0] = mat(2,0);
		m[2][1] = mat(2,1);
	}
}

// trace
double Matrix3::trace(void) const { return m[0][0] + m[1][1] + m[2][2]; }

// determinant
double Matrix3::determinant(void) const
{	if(is2D)
	{	return m[2][2]*(m[0][0]*m[1][1]-m[1][0]*m[0][1]);
	}
	return m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
			-m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2])
			+m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);
}

// second invariant
double Matrix3::second_invariant(void) const
{	double ii = m[0][0]*(m[1][1] + m[2][2]) + m[1][1]*m[2][2] - m[0][1]*m[1][0];
	if(!is2D) ii -= (m[0][2]*m[2][0] + m[1][2]*m[2][1]);
	return ii;
}

// Characteristic equation c(z) = det(zI-M) = z^3 - c0 - c1 z - c2 z^2
void Matrix3::characteristics(double &c0,double &c1,double &c2) const
{
	c2 = trace();
	c1 = -second_invariant();
	c0 = determinant();
}

#pragma mark Matrix3::Class Methods

// return an identity matrix
Matrix3 Matrix3::Identity(void)
{	Matrix3 im(1.,0.,0.,1.,1.);
    return im;
}

// return an rotation matrix for ccw rotation about Z axis for angle (in radians)
Matrix3 Matrix3::CCWZRotation(double ang)
{	double cosw = cos(ang);
	double sinw = sin(ang);
	Matrix3 rm(cosw,-sinw,sinw,cosw,1.);
    return rm;
}

// Return a rotation matrix from Euler angles (in radians or in degrees)
// Follows method in Arfken, pg 179, which is counter clockwise Z, Y', Z''
Matrix3 Matrix3::Rotation3D(double alpha,double beta,double gamma,bool indegrees)
{	if(indegrees)
	{	double convert = PI_CONSTANT/180.;
		alpha *= convert;
		beta *= convert;
		gamma *= convert;
	}
	double ca= cos(alpha),sa=sin(alpha);
	double cb= cos(beta),sb=sin(beta);
	double cg= cos(gamma),sg=sin(gamma);
	Matrix3 rm(  cg*cb*ca-sg*sa,   cg*cb*sa+sg*ca,   -cg*sb,
			    -sg*cb*ca-cg*sa,  -sg*cb*sa+cg*ca,    sg*sb,
			          sb*ca,           sb*sa,           cb   );
    return rm;
}


#pragma mark Eigenanalysis

// ----------------------------------------------------------------------------
// Numerical diagonalization of 2x2 and 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Calculates the eigensystem of a real symmetric 2x2 matrix
//    [ A  B ]
//    [ B  C ]
// in the form
//    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
//    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
// On call, rt1 and rt2 are the eigen values. They are reordered if needed
// such that rt1 >= rt2. The orthongal eigenvectors are
// e1 = (cs, sn), e2 = (-sn, cs)
// ----------------------------------------------------------------------------
inline void dsyev2(double A, double B, double C, double *rt1, double *rt2,
                   double *cs, double *sn)
{
	double df = A - C;
	double rt = sqrt(df*df + 4.0*B*B);
	double t;
	
	// order the eigenvalues
	if(*rt2 > *rt1)
	{	double temp = *rt1;
		*rt1 = *rt2;
		*rt2 = temp;
	}
	
	// Calculate eigenvectors
	if (df > 0.0)
		*cs = df + rt;
	else
		*cs = df - rt;
	
	if (fabs(*cs) > 2.0*fabs(B))
	{	t   = -2.0 * B / *cs;
		*sn = 1.0 / sqrt(1.0 + t*t);
		*cs = t * (*sn);
	}
	else if (fabs(B) == 0.0)
	{	*cs = 1.0;
		*sn = 0.0;
	}
	else
	{	t   = -0.5 * (*cs) / B;
		*cs = 1.0 / sqrt(1.0 + t*t);
		*sn = t * (*cs);
	}
	
	if (df > 0.0)
	{	t   = *cs;
		*cs = -(*sn);
		*sn = t;
	}
}


int dsyevq3(double A[3][3], double Q[3][3], double *w)
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ----------------------------------------------------------------------------
{
	const int n = 3;
	double e[3];                   // The third element is used only as temporary workspace
	double g, r, p, f, b, s, c, t; // Intermediate storage
	int nIter;
	int m;
	
	// Transform A to real tridiagonal form by the Householder method
	dsytrd3(A, Q, w, e);
	
	// Calculate eigensystem of the remaining real symmetric tridiagonal matrix
	// with the QL method
	//
	// Loop over all off-diagonal elements
	for (int l=0; l < n-1; l++)
	{
		nIter = 0;
		while (1)
		{
			// Check for convergence and exit iteration loop if off-diagonal
			// element e(l) is zero
			for (m=l; m <= n-2; m++)
			{
				g = fabs(w[m])+fabs(w[m+1]);
				if (fabs(e[m]) + g == g)
					break;
			}
			if (m == l)
				break;
			
			if (nIter++ >= 30)
				return -1;
			
			// Calculate g = d_m - k
			g = (w[l+1] - w[l]) / (e[l] + e[l]);
			r = sqrt(SQR(g) + 1.0);
			if (g > 0)
				g = w[m] - w[l] + e[l]/(g + r);
			else
				g = w[m] - w[l] + e[l]/(g - r);
			
			s = c = 1.0;
			p = 0.0;
			for (int i=m-1; i >= l; i--)
			{
				f = s * e[i];
				b = c * e[i];
				if (fabs(f) > fabs(g))
				{
					c      = g / f;
					r      = sqrt(SQR(c) + 1.0);
					e[i+1] = f * r;
					c     *= (s = 1.0/r);
				}
				else
				{
					s      = f / g;
					r      = sqrt(SQR(s) + 1.0);
					e[i+1] = g * r;
					s     *= (c = 1.0/r);
				}
				
				g = w[i+1] - p;
				r = (w[i] - g)*s + 2.0*c*b;
				p = s * r;
				w[i+1] = g + p;
				g = c*r - b;
				
				// Form eigenvectors
				for (int k=0; k < n; k++)
				{
					t = Q[k][i+1];
					Q[k][i+1] = s*Q[k][i] + c*t;
					Q[k][i]   = c*Q[k][i] - s*t;
				}
			}
			w[l] -= p;
			e[l]  = g;
			e[m]  = 0.0;
		}
	}
	
	return 0;
}

void dsytrd3(double A[3][3], double Q[3][3], double *d, double *e)
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
	const int n = 3;
	double u[n], q[n];
	double omega, f;
	double K, h, g;
	
	// Initialize Q to the identitity matrix
	for (int i=0; i < n; i++)
	{
		Q[i][i] = 1.0;
		for (int j=0; j < i; j++)
			Q[i][j] = Q[j][i] = 0.0;
	}
	
	// Bring first row and column to the desired form
	h = SQR(A[0][1]) + SQR(A[0][2]);
	if (A[0][1] > 0)
		g = -sqrt(h);
	else
		g = sqrt(h);
	e[0] = g;
	f    = g * A[0][1];
	u[1] = A[0][1] - g;
	u[2] = A[0][2];
	
	omega = h - f;
	if (omega > 0.0)
	{
		omega = 1.0 / omega;
		K     = 0.0;
		for (int i=1; i < n; i++)
		{
			f    = A[1][i] * u[1] + A[i][2] * u[2];
			q[i] = omega * f;                  // p
			K   += u[i] * f;                   // u* A u
		}
		K *= 0.5 * SQR(omega);
		
		for (int i=1; i < n; i++)
			q[i] = q[i] - K * u[i];
		
		d[0] = A[0][0];
		d[1] = A[1][1] - 2.0*q[1]*u[1];
		d[2] = A[2][2] - 2.0*q[2]*u[2];
		
		// Store inverse Householder transformation in Q
		for (int j=1; j < n; j++)
		{
			f = omega * u[j];
			for (int i=1; i < n; i++)
				Q[i][j] = Q[i][j] - f*u[i];
		}
		
		// Calculate updated A[1][2] and store it in e[1]
		e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
	}
	else
	{
		for (int i=0; i < n; i++)
			d[i] = A[i][i];
		e[1] = A[1][2];
	}
}




