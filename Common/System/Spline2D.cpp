/********************************************************************************
    Spline2D.cpp
    nairn-mpm-fea
 
    Documentation
 
    1. Interpolation with strictly increasing x data for y(x)
            HermiteSpline()
            FreeSpline()

        After getting spline, plot with
            BSplineCurve,BSplinePoint

    2. Interpolation of arbitrary x and y data by parametric x(t) and y(t)
            ParametricHermiteSpline()
            ParametricFreeSpline()
 
        After getting spline, plot with
            XYBSplineCurve,XYBSplinePoint
 
    Future
        Least squares interpolant
        Natural interpolant
        Make linear solver static function
 
    Created by John Nairn on 9/7/23.
    Copyright (c) 2023 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Spline2D.hpp"

#pragma mark Spline2D:Constructors and Destructor

// empty matrix
Spline2D::Spline2D()
{
}

// clear spline fit
void Spline2D::ClearSplineFit(void)
{
    knots.clear();
    bx.clear();
    by.clear();
}

// destructor (may not be needed)
Spline2D::~Spline2D()
{   ClearSplineFit();
}

#pragma mark Spline2D:Calculate Splines

/**********************************************************************************
    Create a spline interpolate through x and y data. They need not be monotonic
    but should be the same length and have at least 4 (because doing cubic splines)
    The length and number are not checked, so must check before calling
 
    Gets parameter fits x[t] and y[t] where t is picked to be strictly
    increaseing

    dydx1 and dydxn are initial and final slopes for the interpolant
    if x[1]==x[0], provide dxdy1 and set flip1 to true
    if x[n-1]==x[n-2], provide dxdyn and set flipn to true
 
    return 0 if acceptable or error code if not
**********************************************************************************/
int Spline2D::XYHermiteSpline(vector<double> sx,vector<double> sy,double dydx1,bool flip1,double dydxn,bool flipn)
{
    // length based on x data
    unsigned long n = sx.size();
    
    // cardinal knots of length n
    vector<double> cx;
    for(int i=0;i<n;i++) cx.push_back((double)i);
    
    // clear knots (bx and by cleared when each axis is done)
    knots.clear();
    
    // not that a matrix in two interpolations are the same. Could
    // speed up slightly be doing just once
    
    // do x axis using cx and sx
    double dx1 = flip1 ? dydx1*(sy[1]-sy[0]) : sx[1]-sx[0] ;
    double dxn = flipn ? dydxn*(sy[n-1]-sy[n-2]) : sx[n-1]-sx[n-2] ;
    if(!HermiteInterpolant(cx,sx,bx,dx1,dxn))
        return SPLINE_XERR;
    
    // do y axis
    double dy1 = flip1 ? sy[1]-sy[0] : dydx1*(sx[1]-sx[0]) ;
    double dyn = flipn ? sy[n-1]-sy[n-2] : dydxn*(sx[n-1]-sx[n-2]) ;
    if(!HermiteInterpolant(cx,sy,by,dy1,dyn))
        return SPLINE_YERR;
    
    return SPLINE_NOERR;
}

/********************************************************************************
    Create a spline interpolate through x and y data. They need not be monotonic
    but should be the same length and have at least 4 (because doing cubic splines)
    The length and number are not checked, so must check before calling
 
    Gets parameter fits x[t] and y[t] where t is picked to be monotonically
    increaseing
 
    A "Free" interpolant where g''' continuous at x2 and xmi1

    return 0 if acceptable or error code if not
**********************************************************************************/
int Spline2D::XYFreeSpline(vector<double> sx,vector<double> sy)
{
    // length based on x data
    unsigned long n = sx.size();
    
    // cardinal knots of length n
    vector<double> cx;
    for(int i=0;i<n;i++) cx.push_back((double)i);
    
    // clear knots (bx and by cleared when each axis is done)
    knots.clear();
    
    // not that a matrix in two interpolations are the same. Could
    // speed up slightly be doing just once

    // fit x axis
    if(!FreeInterpolant(cx,sx,bx))
        return SPLINE_XERR;
    
    // fit y axis
    if(!FreeInterpolant(cx,sy,by))
        return SPLINE_YERR;
    
    return SPLINE_NOERR;
}

/****************************************************************
    Fit interpolating cubic spline function through y[1] to
        y[n] using x[1] to x[n] to define the knots. The x
        values must be strictly increasing.
 
    Vector x has n values and assumed > 4 (for cubic splines)
    Vector y assumed at least n values
    Class vector knots assumed prefilled for Hermite interpolant
 
    Spline fit output in b
    return true (or false on LinSolve error)
****************************************************************/
bool Spline2D::HermiteInterpolant(vector<double> x,vector<double> y,double f1p,double fnp)
{   knots.clear();          // assume single fit here and not parametric fit
    by.clear();             // bx cleared in call
    return HermiteInterpolant(x,y,bx,f1p,fnp);
}
bool Spline2D::HermiteInterpolant(vector<double> x,vector<double> y,vector<double> &c,double f1p,double fnp)
{
    // get size and clear vector
    int n = (int)(x.size());
    c.clear();
    
    // knots for the spline fit: x0,x0,x0,x0,x1,x2,...xn-2,xn-1,xn-1,xn-1,xn-1
    // length n+6, if not empty assume set already
    if(knots.size()==0)
    {   int i;
        for(i=0;i<3;i++) knots.push_back(x[0]);
        for(i=0;i<n;i++) knots.push_back(x[i]);
        for(i=0;i<3;i++) knots.push_back(x[n-1]);
    }
        
    // tridiagonal system (size n+2)
    vector<vector<double>> a;
    for(int i=0;i<n+2;i++)
    {   vector<double> row;
        row.assign(n+2,0.);
        a.push_back(row);
    }

    // first row (1)
    a[0][0] = BSplineDeriv(1,3,n+2,x[0]);
    a[0][1] = BSplineDeriv(2,3,n+2,x[0]);
    c.push_back(f1p);
    
    // internal rows (2 to n+1)
    for(int i=2;i<=n+1;i++)
    {   a[i-1][i-2] = BSplineValue(i-1,3,n+2,x[i-2]);
        a[i-1][i-1] = BSplineValue(i,3,n+2,x[i-2]);
        a[i-1][i] = BSplineValue(i+1,3,n+2,x[i-2]);
        c.push_back(y[i-2]);            // y is zero based, pushes 0 to n-1
    }
    
    // last row (n+2)
    a[n+1][n] = BSplineDeriv(n+1,3,n+2,x[n-1]);
    a[n+1][n+1] = BSplineDeriv(n+2,3,n+2,x[n-1]);
    c.push_back(fnp);
    
    // solve and then done
    if(!LinSolve(a,c)) return false;
    
    return true;
}

/****************************************************************
    Fit interpolating cubic spline function through y[1] to
        y[n] using x[1] to x[n] to define the knots. The x
        values must be strictly increasing.

    Vector x has n values and assumed > 4 (for cubic splines)
    Vector y assumed at least n values
    Class vector knots assumed prefilled for free interpolant

    Spline fit output in c
****************************************************************/
bool Spline2D::FreeInterpolant(vector<double> &x,vector<double> &y)
{   knots.clear();          // assume single fit here and not parametric fit
    by.clear();             // bx cleared in call
    return FreeInterpolant(x,y,bx);
}
bool Spline2D::FreeInterpolant(vector<double> &x,vector<double> &y,vector<double> &c)
{
    // get size and cloear vector
    int n = (int)(x.size());
    c.clear();
     
    // knots x0,x0,x0,x0,x2,...xn-3,xn-1,xn-1,xn-1,xn-1
    // length n+4, if not empty assume set already
    if(knots.size()==0)
    {   int i;
        for(i=0;i<4;i++) knots.push_back(x[0]);
        // free ends skips x1 and xn-2 (zero based)
        for(i=2;i<n-2;i++) knots.push_back(x[i]);
        for(i=0;i<4;i++) knots.push_back(x[n-1]);
    }

    // matrix (size n)
    vector<vector<double>> ac;
    for(int i=0;i<n;i++)
    {   vector<double> row;
        row.assign(n,0.);
        ac.push_back(row);
    }

    // decompose the controlling matrix
    int i,j,r=4;
    for(i=0;i<n;i++)
    {   // x[i] is at knots[r] thus
        //      f(t) = Sum(j=r-3,r) c[j] B(j,3)(t)
        for(j=r-4;j<r;j++)
            ac[i][j] = BSplineValue(j+1,3,n,x[i]);
        
        // right hand side
        c.push_back(y[i]);
        
        // advance interval
        if(i!=0 && i<n-3) r++;
    }
    
    // solve and then done
    if(!LinSolve(ac,c)) return false;

    return true;
}

#pragma mark Spline2D:Generate spline curves

/*************************************************************************
    Given parametric fits to x(t) and y(t) in bx, by, and knots, get
    (x,y) curve using npts points
 
    derive can only be 0, 1, or 2 for parameteric plots of cubic splines
    if error, return interger code and x and y will be empty
 
    Assumed parametric fits found first by
        XYHermiteSpline() or XYFreeSpline()
**************************************************************************/
int Spline2D::XYBSplineCurve(vector<double> &x,vector<double> &y,int deriv,int npts)
{
    // clear if needed
    x.clear();
    y.clear();
    
    // need bx and by
    if(bx.size()==0 || by.size()==0)
    {   // incomplete parametric fit
        return SPLINE_NOTPARAMETRIC;
    }

    // e.g., interpolant of N points with cubic splines will
    //  have n = N+2, knots length N+6, then d = 3
    int n=(int)(bx.size());
    int d=(int)(knots.size())-n-1;
    if(d<1)
    {   // The vector of knots is too small
        return SPLINE_TOOFEWKNOTS;
    }
    if(deriv>=d)
    {   // The derivative is too high for input order spline fit
        return SPLINE_INVALID_DERIVATIVE;
    }
    if(deriv<0 || deriv>2)
    {   // only 0, 1, 2 supports
        return SPLINE_INVALID_DERIVATIVE;
    }

    // get spacing
    double xval,yval,t,tmin=knots[d],tmax=knots[n];
    double spacing=(tmax-tmin)/(double)(npts-1);
    
    // parametrically take t from tmin to tmax
    int i,ii,j=d+1;
    t = tmin;
    for(i=0;i<npts;i++)
    {   if(t>knots[n]) t = knots[n];
        while(t>knots[j]) j++;
        
        // get x(t) always
        xval = 0.;
        for(ii=j-d;ii<=j;ii++)
            xval += Barjp(ii,1,d,bx)*BSplineValue(ii,d,n,t);
        
        if(deriv==0)
        {   // get y(t)
            yval = 0.;
            for(ii=j-d;ii<=j;ii++)
                yval += Barjp(ii,1,d,by)*BSplineValue(ii,d,n,t);
        }
        else
        {   // d^(deriv)y/dt^(deriv) (t)
            double dyddt = 0.;
            for(ii=j-d+deriv;ii<=j;ii++)
                dyddt += Barjp(ii,deriv+1,d,by)*BSplineValue(ii,d-deriv,n,t);
            
            // d^(deriv)x/dt^(deriv) (t)
            double dxddt = 0.;
            for(ii=j-d+deriv;ii<=j;ii++)
                dxddt += Barjp(ii,deriv+1,d,bx)*BSplineValue(ii,d-deriv,n,t);
            
            if(deriv==1)
            {   yval = dyddt/dxddt;
            }
            else
            {   // dy/dt (t)
                double dydt = 0.;
                for(ii=j-d+1;ii<=j;ii++)
                    dydt += Barjp(ii,2,d,by)*BSplineValue(ii,d-1,n,t);
                
                // dx/dt (t)
                double dxdt = 0.;
                for(ii=j-d+1;ii<=j;ii++)
                    dxdt += Barjp(ii,2,d,bx)*BSplineValue(ii,d-1,n,t);
                
                // get second derivative
                yval = dyddt/dxdt - dxddt*dydt/(dxdt*dxdt);
            }
        }

        // push to vectors
        x.push_back(xval);
        y.push_back(yval);
        
        // next t
        t+=spacing;
    }
    
    return SPLINE_NOERR;
}

/*************************************************************************
    Given interpolant to y[1] to y[n] in bx and knots, get
    (x,y) curve using npts points or a derivative
 
    deriv only up to spline order-1 (e.g., 2 for cubic)
        I think up to order works, by just step functions
    on error, return integer code and x and y are empty
 
    Assumed interpolant fits found first by
        HermiteSpline() or FreeSpline()
**************************************************************************/
int Spline2D::BSplineCurve(vector<double> &x,vector<double> &y,int deriv,int npts)
{
    // clear if needed
    x.clear();
    y.clear();
    
    // need bx and by
    if(bx.size()==0)
    {   // incomplete or missing fit
        return SPLINE_NOFIT;
    }

    // e.g., interpolant of N points with cubic splines will
    //  have n = N+2, knots length N+6, then d = 3
    int n=(int)(bx.size());
    int d=(int)(knots.size())-n-1;
    if(d<1)
    {   // The vector of knots is too small
        return SPLINE_TOOFEWKNOTS;
    }
    if(deriv>=d)
    {   // The derivative is too high for input order spline fit
        return SPLINE_INVALID_DERIVATIVE;
    }

    // get spacing
    double t,tmin=knots[d],tmax=knots[n];
    double spacing=(tmax-tmin)/(double)(npts-1);
    
    // parametrically take t from tmin to tmax
    int i,ii,j=d+1;
    t = tmin;
    for(i=0;i<npts;i++)
    {   if(t>knots[n]) t = knots[n];
        while(t>knots[j]) j++;
        
        // d^(deriv)y/dt^(deriv) (t)
        double yval = 0.;
        for(ii=j-d+deriv;ii<=j;ii++)
            yval += Barjp(ii,deriv+1,d,bx)*BSplineValue(ii,d-deriv,n,t);

        // push to vectors
        x.push_back(t);
        y.push_back(yval);
        
        // next t
        t+=spacing;
    }

    return SPLINE_NOERR;
}

/*************************************************************************
    Given parametric fits to x(t) and y(t) in bx, by, and knots, get
    (x,y) at point t from 0 to n-1 where n is number of points. The
    returned point will be in the t^(th) interval of the interpolant
    numbered from 0
 
    deriv can only be 0, 1, or 2 for parameteric plots of cubic splines
 
    If works, (x,y) values in x[0] and x[1]
    if error, return interger code and x will be empty
    (note t<0 or t> last point is out of range error
 
    Assumed parametric fits found first by
        XYHermiteSpline() or XYFreeSpline()
**************************************************************************/
int Spline2D::XYBSplinePoint(double t,int deriv,vector<double> &x)
{
    // clear if needed
    x.clear();
    
    // need bx and by
    if(bx.size()==0 || by.size()==0)
    {   // incomplete parametric fit
        return SPLINE_NOTPARAMETRIC;
    }

    // e.g., interpolant of N points with cubic splines will
    //  have n = N+2, knots length N+6, then d = 3
    int n=(int)(bx.size());
    int d=(int)(knots.size())-n-1;
    if(d<1)
    {   // The vector of knots is too small
        return SPLINE_TOOFEWKNOTS;
    }
    if(deriv>=d)
    {   // The derivative is too high for input order spline fit
        return SPLINE_INVALID_DERIVATIVE;
    }
    if(deriv<0 || deriv>2)
    {   // only 0, 1, 2 supports
        return SPLINE_INVALID_DERIVATIVE;
    }

    if(t<knots[d] || t>knots[n]) return SPLINE_OUTOFRANGE;
    int i,j=d+1;
    while(t>knots[j]) j++;
    
    // get x(t) always
    double yval,xval = 0;
    for(i=j-d;i<=j;i++)
        xval += Barjp(i,1,d,bx)*BSplineValue(i,d,n,t);

    if(deriv==0)
    {   // get y(t)
        yval = 0.;
        for(i=j-d;i<=j;i++)
            yval += Barjp(i,1,d,by)*BSplineValue(i,d,n,t);
    }
    else
    {   // d^(deriv)y/dt^(deriv) (t)
        double dyddt = 0.;
        for(i=j-d+deriv;i<=j;i++)
            dyddt += Barjp(i,deriv+1,d,by)*BSplineValue(i,d-deriv,n,t);
        
        // d^(deriv)x/dt^(deriv) (t)
        double dxddt = 0.;
        for(i=j-d+deriv;i<=j;i++)
            dxddt += Barjp(i,deriv+1,d,bx)*BSplineValue(i,d-deriv,n,t);
        
        if(deriv==1)
        {   yval = dyddt/dxddt;
        }
        else
        {   // dy/dt (t)
            double dydt = 0.;
            for(i=j-d+1;i<=j;i++)
                dydt += Barjp(i,2,d,by)*BSplineValue(i,d-1,n,t);
            
            // dx/dt (t)
            double dxdt = 0.;
            for(i=j-d+1;i<=j;i++)
                dxdt += Barjp(i,2,d,bx)*BSplineValue(i,d-1,n,t);
            
            // get second derivative
            yval = dyddt/dxdt - dxddt*dydt/(dxdt*dxdt);
        }
    }
    
    // push to vectors
    x.push_back(xval);
    x.push_back(yval);
    
    return SPLINE_NOERR;
}

/*************************************************************************
    Find deriv(th) derivative of the curve at t of the current spline
    fit in bx and knots
 
    If fit trom strictly increasing x, t is x value and derivative is
    wrt to x
    
    If parametric fit, t is from 0 to n-1 (number of points) and returns
    x value only. Use XYSplinePoint() to get point on a parametric
    plot interpolant
 
    Returns 0 if t outside knot span
 
    if error, error code is set, otherwise set to SPLINE_NOERR
**************************************************************************/
double Spline2D::BSplinePoint(double t,int deriv,int *err)
{
    // clear error
    *err = 0;
    
    // e.g., interpolant of N points with cubic splines will
    //  have n = N+2, knots length N+6, then d = 3
    int n=(int)(bx.size());
    int d=(int)(knots.size())-n-1;
    if(d<1)
    {   // The vector of knots is too small
        *err = SPLINE_TOOFEWKNOTS;
        return 0;
    }
    if(deriv>=d)
    {   // The derivative is too high for input order spline fit
        *err = SPLINE_INVALID_DERIVATIVE;
        return 0;
    }
        
    // find interval
    double f = 0.;
    if(t<knots[d] || t>knots[n]) return f;
    int i,j=d+1;
    while(t>knots[j]) j++;
    for(i=j-d+deriv;i<=j;i++)
        f += Barjp(i,deriv+1,d,bx)*BSplineValue(i,d-deriv,n,t);
    
    return f;
}

#pragma mark Spline2D:Internal Methods

/*************************************************************************
    Find coefficient for expansion of derivatives of BSpline curve
    See deBoor, "Package for calculating with B-Splines," page 445
**************************************************************************/
double Spline2D::Barjp(int r,int jp,int d,vector<double> a)
{
    if(jp==1) return a[r-1];
    int j=jp-1;
    if(knots[r+d-j]>knots[r-1])
        return (double(d+1-j))*(Barjp(r,j,d,a)-Barjp(r-1,j,d,a))/(knots[r+d-j]-knots[r-1]);
    else
        return 0.;
}


/*************************************************************************
    Find BSpline basis function B(i,d)(t) by recursion.
    The knots are in class variables
    
    d is polynomial order (d=1 linear, 2 quadratic, etc.)
    n is needed to deal with B(i,0)(t) in last interval
**************************************************************************/
double Spline2D::BSplineValue(int i,int d,int n,double t)
{   double value=0.;
    if(d==0)
    {   if(t>=knots[i-1] && t<knots[i])
            value=1.;
        else if(i==n && t==knots[i])
            value=1.;
    }
    else
    {   if(knots[i-1+d]>knots[i-1])
            value += (t-knots[i-1])*BSplineValue(i,d-1,n,t)/(knots[i-1+d]-knots[i-1]);
        if(knots[i+d]>knots[i])
            value += (knots[i+d]-t)*BSplineValue(i+1,d-1,n,t)/(knots[i+d]-knots[i]);
    }
    return value;
}

/*************************************************************************
    Find first derivative of B(i,d)(t)
    The knots are in class variables
    Note that i is 1 based and knots are stored 0 based
 
    d is polynomial order (d=1, linear, 2 quadratic, etc.)
    n only needed to pass to BSplineValue (and then only when d=0)
**************************************************************************/
double Spline2D::BSplineDeriv(int i,int d,int n,double t)
{
    if(knots[i-1]==knots[i+d])
        return 0.;
    else if(knots[i-1]<knots[i-1+d] && knots[i]==knots[i+d])
        return BSplineValue(i,d-1,n,t)/(knots[i-1+d]-knots[i-1]);
    else if(knots[i-1]==knots[i-1+d] && knots[i]<knots[i+d])
        return -BSplineValue(i+1,d-1,n,t)/(knots[i+d]-knots[i]);
    
    return BSplineValue(i,d-1,n,t)/(knots[i-1+d]-knots[i-1]) - BSplineValue(i+1,d-1,n,t)/(knots[i+d]-knots[i]);
}

#pragma mark Spline2D:Linear Solver

/****************************************************************
    Solve linear system a X = b
    On exit X will replace b and a will be decomposed
    a must be nXn and b must be n (not checked)

    To solve multiple times better to decompose matrix
        and reuse on multiple substitutiones
****************************************************************/
bool Spline2D::LinSolve(vector<vector<double>> &a,vector<double> &b)
{
    // matrix size
    int n = (int)(a.size());
    
    // index vecctor
    vector<int> indx;
    indx.assign(n,0);
    
    // prints a in Mathematica Form
    /*
    cout << "{";
    for(int i=0;i<n;i++)
    {   cout << "{";
        for(int j=0;j<n;j++)
        {   cout << a[i][j];
            if(j<n-1) cout << ",";
        }
        cout << "}";
        if(i<n-1) cout << ",";
    }
    cout << "}" << endl;
    */

    // solve
    double d;
    if(!Ludcmp(a,indx,&d)) return false;
    if(!Lubksb(a,indx,b)) return false;
    
    return true;
}

/****************************************************************
    LU decomposition of a square matrix. The matrix must be over
    the range [0...n-1][0...n-1] and indx is [0...n-1] (not checked)

    Decomposed matrix returned in a and indx needed for further
        use of decomposted matrix in indx
    d is + or - 1 depending on row permutations
    
    Return false or true if works
    Numerical recipes in C, pg 43
****************************************************************/
bool Spline2D::Ludcmp(vector<vector<double>> &a,vector<int> &indx,double *d)
{
    int n,i,j,k,imax=0;
    double big,temp,sum,dum;
    
    // matrix size
    n = (int)(a.size());
    
    // find largest element, exit if row all zero
    *d=1.;
    vector<double> vv;
    for(i=0;i<n;i++)
    {   big=0.;
        for(j=0;j<n;j++)
        {   temp=a[i][j];
            if(temp<0.) temp=-temp;
            if(temp>big) big=temp;
        }
        if(big==0.) return false;
        vv.push_back(1./big);
    }
    
    // decomposition
    for(j=0;j<n;j++)
    {   for(i=0;i<j;i++)
        {   sum=a[i][j];
            for(k=0;k<i;k++)
                sum-=a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.;
        for(i=j;i<n;i++)
        {   sum=a[i][j];
            for(k=0;k<j;k++)
                sum-=a[i][k]*a[k][j];
            a[i][j]=sum;
            if(sum>=0.)
                dum=vv[i]*sum;
            else
                dum=-vv[i]*sum;
            if(dum>=big)
            {   big=dum;
                imax=i;
            }
        }
        if(j!=imax)
        {   for(k=0;k<n;k++)
            {   dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d=-(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        //if(a[j][j]==0.) a[j][j]=1.0e-20;
        if(a[j][j]==0.) return false;
        if(j!=n-1)
        {   dum=1./(a[j][j]);
            for(i=j+1;i<n;i++)
                a[i][j]*=dum;
        }
    }
    return true;
}

/****************************************************************
    LU back substitution of a decomposed square matrix.
    The matrix must be over the range [1...n][1...n]
    
    indx and b must be at least n long (their size will not
        be changed)
        
    Solves A X = B where A is undecomposed matrix. Solution
        replaces B.
        
    Return FALSE on error and set error message
    Return TRUE if it works
    
    Numerical recipes in C, pg 44
****************************************************************/
bool Spline2D::Lubksb(vector<vector<double>> &a,vector<int> &indx,vector<double> &b)
{
    int i,ii=-1,ip,j,n;
    double sum;
    
    // matrix size
    n = (int)(a.size());
    if(indx.size()<n) return false;
    if(b.size()<n) return false;
    
    for(i=0;i<n;i++)
    {   ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if(ii>=0)
        {   for(j=ii;j<=i-1;j++)
                sum-=a[i][j]*b[j];
        }
        else if(sum)
        {    ii=i;
        }
        b[i]=sum;
    }
    for(i=n-1;i>=0;i--)
    {   sum=b[i];
        for(j=i+1;j<n;j++)
            sum-=a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
    
    return true;
}

#pragma mark One Step Spline2D creators

// create spline curve with set end derivatives through arbirary curve
// x need not be montonic
// reguire size y >= size x and size x >= 4 (checked)
// Result in Spline2D object with bx and by for Bezier control points
//      for knots in class. Each axis done by cardinal knots 0  to n-1
//      where n is number of points
// returns NULL on error
// caller should delete the Spline2D object
Spline2D *ParametricHermiteSpline(vector<double> sx,vector<double> sy,
                              double dydx1,bool flip1,double dydxn,bool flipn)
{
    unsigned long n = sx.size();
    if(sy.size()<n || n<4)
    {   // must have sy vector equal to long to x
        // must have n at least 4 (because doing cubinc splines)
        return NULL;
    }
    
    // create spline curve
    Spline2D *theSpline = new Spline2D();
    if(theSpline->XYHermiteSpline(sx,sy,dydx1,flip1,dydxn,flipn)!=0)
    {   delete theSpline;
        return NULL;
    }
    
    return theSpline;
}

// create spline curve with set end derivatives through data when
// x is strictly increasing (not checked)
// reguire size y >= size x and size >= 4 (checked)
// Result in Spline2D object with bx as Bezier control points
//      for knots in class.
// returns NULL on error
// caller should delete the Spline2D object
Spline2D *HermiteSpline(vector<double> sx,vector<double> sy,double dydx1,double dydxn)
{
    unsigned long n = sx.size();
    if(sy.size()<n || n<4)
    {   // must have sy vector equal to long to x
        // must have n at least 4 (because doing cubinc splines)
        return NULL;
    }
    
    // create spline curve
    Spline2D *theSpline = new Spline2D();
    if(!theSpline->HermiteInterpolant(sx,sy,dydx1,dydxn))
    {   delete theSpline;
        return NULL;
    }
    
    return theSpline;
}

// create spline curve through arbirary curve with free conditions
// x need not be montonic
// Free means g'''(x)=0 and 2nd and penultimate points
// reguire size y >= size x and size x >= 4 (checked)
// Result is Spline2D object with bx and by for Bezier control points
//      for knots in class. Each axis done by cardinal knots  0  to n-1
//      where n is number of points
// returns NULL on error
// caller should delete the Spline2D object
Spline2D *ParametricFreeSpline(vector<double> sx,vector<double> sy)
{
    unsigned long n = sx.size();
    if(sy.size()<n || n<4)
    {   // must have sy vector equal to long to x
        // must have n at least 4 (because doing cubinc splines)
        return NULL;
    }
    
    // create spline curve
    Spline2D *theSpline = new Spline2D();
    if(theSpline->XYFreeSpline(sx,sy)!=0)
    {   delete theSpline;
        return NULL;
    }
    
    return theSpline;
}

// create spline curve through points when
// x is strictly increasing (not checked)
// Free means g'''(x)=0 and 2 and penultimate points
// reguire size y >= size x and size x >= 4 (checked)
// Result is Spline2D object with bx for Bezier control points
//      for knots in class.
// returns NULL on error
// caller should delete the Spline2D object
Spline2D *FreeSpline(vector<double> sx,vector<double> sy)
{
    unsigned long n = sx.size();
    if(sy.size()<n || n<4)
    {   // must have sy vector equal to long to x
        // must have n at least 4 (because doing cubinc splines)
        return NULL;
    }
    
    // create spline curve
    Spline2D *theSpline = new Spline2D();
    if(!theSpline->FreeInterpolant(sx,sy))
    {   delete theSpline;
        return NULL;
    }
    
    return theSpline;
}

// To test splines
// Call this in many and exit if returns true (it always does)
// Will also need to include Spline2D.hpp in main
bool SplineTesting(void)
{
    // pick some points
    vector<double> x;
    vector<double> y;
    x.push_back(0.);
    y.push_back(0.);
    x.push_back(2.);
    y.push_back(0.05);
    x.push_back(3.);
    y.push_back(0.2);
    x.push_back(4.);
    y.push_back(0.55);
    x.push_back(5.);
    y.push_back(0.8);
    x.push_back(6.);
    y.push_back(0.95);
    x.push_back(7.);
    y.push_back(1.);
    x.push_back(8.5);
    y.push_back(0.8);
    x.push_back(10.);
    y.push_back(1.2);
    cout << "#setLineType\tnone" << endl;
    cout << "#setSymbolType\tsquare" << endl;
    cout << "#setName\tinput" << endl;
    for(int i=0;i<x.size();i++)
    {   cout << x[i] << "," << y[i] << endl;
    }
    cout << endl;
    
    char clr[100];
    vector<double> xp;
    vector<double> yp;
    int npts = 100;
    int deriv = 0;
    for(int stype=1;stype<=2;stype++)
    {   Spline2D *theSpline = NULL;
        if(stype==1)
        {   theSpline = ParametricHermiteSpline(x,y,0.,false,.1,false);
            if(theSpline!=NULL)
                theSpline->XYBSplineCurve(xp,yp,deriv,npts);
            strcpy(clr,"black");
        }
        else if(stype==2)
        {   theSpline = ParametricFreeSpline(x,y);
            if(theSpline!=NULL)
                theSpline->XYBSplineCurve(xp,yp,deriv,npts);
            strcpy(clr,"blue");
        }
        else if(stype==3)
        {   theSpline = HermiteSpline(x,y,0.,0.);
            if(theSpline!=NULL)
                theSpline->BSplineCurve(xp,yp,deriv,npts);
            strcpy(clr,"red");
        }
        else if(stype==4)
        {   theSpline = FreeSpline(x,y);
            if(theSpline!=NULL)
                theSpline->BSplineCurve(xp,yp,deriv,npts);
            strcpy(clr,"green");
        }

        if(theSpline==NULL)
        {   cout << "*** Error *** for stype " << stype << endl;
            continue;
        }
        else
        {   //delete theSpline;
            //theSpline = NULL;
        }
        
        if(xp.size()>0)
        {   cout << "#setColor\t" << clr << endl;
            cout << "#setName\toutput" << stype << endl;
            for(int i=0;i<xp.size();i++)
            {   cout << xp[i] << "," << yp[i] << "\n";
            }
            cout << endl;
        }
    }
    
    return true;
}

