/********************************************************************************
    Spline2D.hpp
    nairn-mpm-fea
 
    Created by John Nairn on 9/7/23.
    Copyright (c) 2023 John A. Nairn, All rights reserved.
********************************************************************************/

#ifndef _JAN_SPLINE2D_
#define _JAN_SPLINE2D_

enum { SPLINE_NOERR=0,SPLINE_XERR,SPLINE_YERR,
    SPLINE_TOOFEWKNOTS,SPLINE_INVALID_DERIVATIVE,
    SPLINE_NOTPARAMETRIC,SPLINE_NOFIT,SPLINE_OUTOFRANGE
};

class Spline2D
{
    public:
    
        //  Constructors and Destructor
        Spline2D();
        virtual ~Spline2D();
        void ClearSplineFit(void);
     
        // methods
        int XYHermiteSpline(vector<double>,vector<double>,double,bool,double,bool);
        bool HermiteInterpolant(vector<double>,vector<double>,double,double);
        int XYFreeSpline(vector<double>,vector<double>);
        bool FreeInterpolant(vector<double> &,vector<double> &);
        int XYBSplineCurve(vector<double> &,vector<double> &,int,int);
        int BSplineCurve(vector<double> &,vector<double> &,int,int);
        int XYBSplinePoint(double,int,vector<double> &);
        double BSplinePoint(double,int,int *);

    protected:
        vector<double> bx;
        vector<double> by;
        vector<double> knots;
    
        // baase methods
        double Barjp(int,int,int,vector<double>);
        bool HermiteInterpolant(vector<double>,vector<double>,vector<double> &,double,double);
        bool FreeInterpolant(vector<double> &,vector<double> &,vector<double> &);
        double BSplineValue(int,int,int,double);
        double BSplineDeriv(int,int,int,double);
        bool LinSolve(vector<vector<double>> &,vector<double> &);
        bool Ludcmp(vector<vector<double>> &,vector<int> &,double *);
        bool Lubksb(vector<vector<double>> &,vector<int> &,vector<double> &);

};

// Create and return Spline2D object
Spline2D *ParametricHermiteSpline(vector<double>,vector<double>,double,bool,double,bool);
Spline2D *HermiteSpline(vector<double>,vector<double>,double,double);
Spline2D *ParametricFreeSpline(vector<double>,vector<double>);
Spline2D *FreeSpline(vector<double>,vector<double>);

// for testing
bool SplineTesting(void);


#endif /* Spline2D_hpp */
