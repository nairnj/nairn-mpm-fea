/*************************************************************************************
	Matrix3.cpp
	NairnFEA and NairnMPM

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

#include "Matrix3.hpp"

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

// Find Eigenvalues and return real and imaginary parts in separate arrays
// return TRUE is positive definite or FALSE if not
// Input double pointers must have room for three values
bool Matrix3::Eigenvalues(double *lam) const
{
    if(is2D)
    {   double b = -(m[0][0]+m[1][1]);
        double c = m[0][0]*m[1][1] - m[1][0]*m[0][1];
        double arg = b*b-4.*c;
        if(arg<0)
        {   if(fabs(arg/b)>1.e-12) return FALSE;
            arg = 0.;
        }
        arg = sqrt(arg);
        lam[0] = 0.5*(-b+arg);
        lam[1] = 0.5*(-b-arg);
        lam[2] = m[2][2];
        return TRUE;
    }
    
    
    // coefficients for c(z) = det(zI-M) = z^3 - c0 - c1 z - c2 z^2
    double c0,c1,c2;
    characteristics(c0,c1,c2);
    
    return FALSE;
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

// copy this object to an array
void Matrix3::get(double c[][3]) const
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			c[i][j] = m[i][j];
	}
}

// set this object to contents of an array
void Matrix3::set(double c[][3])
{	int i,j;
	for(i=0;i<3;i++)
	{	for(j=0;j<3;j++)
			m[i][j] = c[i][j];
	}
	is2D = FALSE;
}

// is this a 2D matrix?
bool Matrix3::getIs2D(void) const { return is2D; }
void Matrix3::setIs2D(bool new2D)
{   is2D = new2D;
    m[0][2] = 0.;
    m[1][2] = 0.;
    m[2][0] = 0.;
    m[2][1] = 0.;
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

// coefficients for c(z) = det(zI-M) = z^3 - c0 - c1 z - c2 z^2
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



