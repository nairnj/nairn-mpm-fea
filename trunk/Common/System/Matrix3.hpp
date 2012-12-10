/********************************************************************************
	Matrix3.hpp
	NairnFEA and NairnMPM
 
	Created by John Nairn on 12/8/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#ifndef _MATRIX3_
#define _MATRIX3_

class Matrix3
{
	public:
	
		//  Constructors and Destructor
		Matrix3();
		Matrix3(double,double,double,double,double,double,
						double,double,double);
		Matrix3(double,double,double,double,double);
	
		// printing
		friend ostream &operator<<(ostream &os, const Matrix3 &);
	
		// methods
		void Zero();
		Matrix3 Exponential(int);
	
		// operators
		Matrix3 &operator=(const Matrix3 &);
		Matrix3 &operator+=(const Matrix3 &);
		const Matrix3 operator+(const Matrix3 &) const;
		Matrix3 &operator-=(const Matrix3 &);
		const Matrix3 operator-(const Matrix3 &) const;
		Matrix3 &operator*=(const Matrix3 &);
		const Matrix3 operator*(const Matrix3 &) const;
		double &operator()(unsigned,unsigned);
		const double &operator()(unsigned row,unsigned col) const;
	
		// accessors
		void get(double c[3][3]) const;
		void set(double);
		void set(int,int,double);
		void set(double c[3][3]);
		bool getIs2D(void) const;
		double determinant(void) const;
		double second_invariant(void) const;
		double trace(void) const;
		void characteristics(double &,double &,double &);
	
	private:
		bool is2D;
		double m[3][3];
};

#endif
