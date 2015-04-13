/********************************************************************************
	Matrix3.hpp
	nairn-mpm-fea
 
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
        Matrix3(double,double,double,double,double,double);
	
		// printing
		friend ostream &operator<<(ostream &os, const Matrix3 &);
	
		// methods
		void Zero();
    
        // read only methods
        Matrix3 Transpose(void) const;
		Matrix3 RMRT(Matrix3 &) const;
		Matrix3 RTMR(Matrix3 &) const;
		Matrix3 Exponential(int) const;
		Vector Times(Vector *) const;
		void Scale(double);
        void Scale2D(double);
		Matrix3 Inverse(void) const;
        Vector Eigenvalues(void) const;
        Matrix3 RightDecompose(Matrix3 *,Vector *) const;
		Matrix3 LeftDecompose(Matrix3 *,Vector *) const;
		Matrix3 Eigenvectors(Vector &) const;
	
		// not currently used
        //Matrix3 RightStretch(const Vector &) const;
        //Matrix3 Rotation(const Vector &,const Matrix3 &) const;
	
		// operators
		Matrix3 &operator+=(const Matrix3 &);
		const Matrix3 operator+(const Matrix3 &) const;
		Matrix3 &operator-=(const Matrix3 &);
		const Matrix3 operator-(const Matrix3 &) const;
		Matrix3 &operator*=(const Matrix3 &);
		const Matrix3 operator*(const Matrix3 &) const;
		double &operator()(unsigned,unsigned);
		const double &operator()(unsigned row,unsigned col) const;
	
		// accessors
		void set(double);
        void set(int,int,double);
		void set(double c[][3]);
        void set(double,double,double,double,double);
        void set(double,double,double,double,double,double,double,double,double);
        void setIs2D(bool);
    
        // read only accessors
        void get(double c[][3]) const;
        void get(int,int,double *) const;
		double get(int,double) const;
        bool getIs2D(void) const;
        double determinant(void) const;
        double second_invariant(void) const;
        double trace(void) const;
        void characteristics(double &,double &,double &) const;
    
        // class methods
        static Matrix3 Identity(void);
        static Matrix3 CCWZRotation(double);
		static Matrix3 Rotation3D(double,double,double,bool);
	
	private:
		bool is2D;
		double m[3][3];
};

#endif
