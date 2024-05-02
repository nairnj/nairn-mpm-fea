/*************************************************************************************
	Matrix4.cpp
	nairn-mpm-fea

	Created by John Nairn on 3/27/2018.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
*************************************************************************************/
 
#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Matrix4.hpp"
#include "Matrix3.hpp"

#pragma mark Matrix4:Constructors and Destructor

// empty matrix
Matrix4::Matrix4()
{
	set(0.);
}

// define matrix in all elements specified
Matrix4::Matrix4(double xx,double xy,double xz,double xw,
				 double yx,double yy,double yz,double yw,
				 double zx,double zy,double zz,double zw,
				 double wx,double wy,double wz,double ww)
{
	m[0][0] = xx;
	m[0][1] = xy;
	m[0][2] = xz;
	m[0][3] = xw;
	m[1][0] = yx;
	m[1][1] = yy;
	m[1][2] = yz;
	m[1][3] = yw;
	m[2][0] = zx;
	m[2][1] = zy;
	m[2][2] = zz;
	m[2][3] = zw;
	m[3][0] = wx;
	m[3][1] = wy;
	m[3][2] = wz;
	m[3][3] = ww;
}

// define symmetric matrix
Matrix4::Matrix4(double xx,double xy,double xz,double xw,
				 		   double yy,double yz,double yw,
				 					 double zz,double zw,
				 							   double ww)
{
	m[0][0] = xx;
	m[0][1] = xy;
	m[0][2] = xz;
	m[0][3] = xw;
	m[1][0] = xy;
	m[1][1] = yy;
	m[1][2] = yz;
	m[1][3] = yw;
	m[2][0] = xz;
	m[2][1] = yz;
	m[2][2] = zz;
	m[2][3] = zw;
	m[3][0] = xw;
	m[3][1] = yw;
	m[3][2] = zw;
	m[3][3] = ww;
}

// define diagonal matrix
Matrix4::Matrix4(double xx,double yy,double zz,double ww)
{
	set(0.);
	m[0][0] = xx;
	m[1][1] = yy;
	m[2][2] = zz;
	m[3][3] = ww;
}

#pragma mark Matrix4:methods

// zero the matrix
void Matrix4::Zero(void)
{	set(0.);
}

// get the transpose
Matrix4 Matrix4::Transpose(void) const
{   Matrix4 mT(m[0][0],m[1][0],m[2][0],m[3][0],
			   m[0][1],m[1][1],m[2][1],m[3][1],
			   m[0][2],m[1][2],m[2][2],m[3][2],
			   m[0][3],m[1][3],m[2][3],m[3][3]);
	return mT;
}

// Return inverse in a new matrix4
Matrix4 Matrix4::Inverse(void) const
{	Matrix4 inv4(cofactor(0,0), -cofactor(1,0),  cofactor(2,0), -cofactor(3,0),
				-cofactor(0,1),  cofactor(1,1), -cofactor(2,1),  cofactor(3,1),
				 cofactor(0,2), -cofactor(1,2),  cofactor(2,2), -cofactor(3,2),
				-cofactor(0,3),  cofactor(1,3), -cofactor(2,3),  cofactor(3,3));
	inv4.Scale(1./determinant());
	return inv4;
}

// scale all elements by factor
void Matrix4::Scale(double factor)
{
	m[0][0] *= factor;
	m[0][1] *= factor;
	m[0][2] *= factor;
	m[0][3] *= factor;
	m[1][0] *= factor;
	m[1][1] *= factor;
	m[1][2] *= factor;
	m[1][3] *= factor;
	m[2][0] *= factor;
	m[2][1] *= factor;
	m[2][2] *= factor;
	m[2][3] *= factor;
	m[3][0] *= factor;
	m[3][1] *= factor;
	m[3][2] *= factor;
	m[3][3] *= factor;
}

// matrix times vector v, result return in p
// both must be dimension at least 4
void Matrix4::Times(double *v,double *p) const
{
	p[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2] + m[0][3]*v[3];
	p[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2] + m[1][3]*v[3];
	p[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2] + m[2][3]*v[3];
	p[3] = m[3][0]*v[0] + m[3][1]*v[1] + m[3][2]*v[2] + m[3][3]*v[3];
}

#pragma mark Matrix4:operators

// *=
Matrix4 &Matrix4::operator*=(const Matrix4 &rhs)
{
	double cm[4][4],rm[4][4];
	get(cm);
	rhs.get(rm);
	int i,j,k;
	for(i=0;i<4;i++)
	{	for(j=0;j<4;j++)
		{	m[i][j] = 0.;
			for(k=0;k<4;k++)
				m[i][j] += cm[i][k]*rm[k][j];
		}
	}
	return *this;
}

// *
const Matrix4 Matrix4::operator*(const Matrix4 &rhs) const
{	Matrix4 result = *this;
	result *= rhs;
	return result;
}

// printing
ostream &operator<<(ostream &os, const Matrix4 &mat)
{
	os << "[(" << mat(0,0) << "," << mat(0,1) << "," << mat(0,2) << "," << mat(0,3)<< "),("
	<< mat(1,0) << "," << mat(1,1) << "," << mat(1,2) << "," << mat(1,3) << "),("
	<< mat(2,0) << "," << mat(2,1) << "," << mat(2,2) << "," << mat(2,3) << "),("
	<< mat(3,0) << "," << mat(3,1) << "," << mat(3,2) << "," << mat(3,3) << ")]" ;
	return os;
}

// elements (0 based) (no checking of limits)
double &Matrix4::operator()(unsigned row,unsigned col)
{	return m[row][col];
}

// elements (0 based) of const object (no checking of limits)
const double &Matrix4::operator()(unsigned row,unsigned col) const
{	return m[row][col];
}

#pragma mark Matrix4:accessors

// all elements to a constant
void Matrix4::set(double constant)
{	int i,j;
	for(i=0;i<4;i++)
	{	for(j=0;j<4;j++)
			m[i][j] = constant;
	}
}

// copy this object to an array
void Matrix4::get(double c[][4]) const
{	int i,j;
	for(i=0;i<4;i++)
	{	for(j=0;j<4;j++)
			c[i][j] = m[i][j];
	}
}

// set one element (0 based indices)
void Matrix4::set(int row,int col,double value) { m[row][col] = value; }

// cofactor determinant (0 based indices)
double Matrix4::cofactor(int row,int col) const
{
	double c[3][3];
	int crow=0,ccol;
	
	for(int i=0;i<4;i++)
	{	if(i==row) continue;
		ccol = 0;
		for(int j=0;j<4;j++)
		{	if(j==col) continue;
			c[crow][ccol] = m[i][j];
			ccol++;
		}
		crow++;
	}

	Matrix3 cof;
	cof.set(c);
	return cof.determinant();
}

// determinant
double Matrix4::determinant(void) const
{	return m[0][0]*cofactor(0,0) - m[1][0]*cofactor(1,0)
				+ m[2][0]*cofactor(2,0) - m[3][0]*cofactor(3,0);
}

#pragma mark Matrix4::Class Methods

// return an identity matrix
Matrix4 Matrix4::Identity(void)
{	Matrix4 im(1.,1.,1.,1.);
	return im;
}



