/********************************************************************************
	Matrix4.hpp
	nairn-mpm-fea

	Created by John Nairn on 3/28/2018.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
********************************************************************************/

#ifndef _MATRIX4_
#define _MATRIX4_

class Matrix4
{
	public:
	
		//  Constructors and Destructor
		Matrix4();
		Matrix4(double,double,double,double,
				double,double,double,double,
				double,double,double,double,
				double,double,double,double);
		Matrix4(double,double,double,double,
					   double,double,double,
							  double,double,
									 double);
		Matrix4(double,double,double,double);
	
		// operators
		Matrix4 &operator*=(const Matrix4 &);
		const Matrix4 operator*(const Matrix4 &) const;
		friend ostream &operator<<(ostream &os, const Matrix4 &);
		double &operator()(unsigned,unsigned);
		const double &operator()(unsigned row,unsigned col) const;
	
		// methods
		void Zero();
		Matrix4 Transpose(void) const;
		Matrix4 Inverse(void) const;
		void Scale(double);
		void Times(double *v,double *p) const;

		// accessors
		void set(double);
		void set(int,int,double);
		void get(double c[][4]) const;
		double cofactor(int,int) const;
		double determinant(void) const;
	
		// class methods
		static Matrix4 Identity(void);
	
	private:
		double m[4][4];
};

#endif
