/*********************************************************************
    Utilities.cpp
    Nairn Research Group FEA Code
    
    Created by John Nairn on 10/23/05.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
*********************************************************************/

#include "stdafx.h"

// Return TRUE or FALSE according to num even
int IsEven(int num)
{
    int half=num>>1;
    if(2*half==num) return TRUE;
    return FALSE;
}

/* Given a meshed line with n intervals and ratio between first
	and last element size, calculate r and a where a is length
	of the first element and r is the size ratio between elements.
	
	This the element sizes are a, ar, ar^2 ... ar^(n-1)
	
	Solution process: the sum of all element sizes in dimensionless
		units is 2 where ratio = a/ar^(n-1) = 1/r^(n-1)
	
	But need special case for ratio=1 of r=1 and a=2/n
*/
void GetLineParameters(int n,double ratio,double *r,double *a)
{
	// for more than 1 interval
	if(n>1)
	{	// use equations if ratio is not close to 1
		if(fabs(ratio-1.)>.001)
		{	*r=pow(1./ratio,1./((double)(n-1)));
			*a=2.*(1.-*r)/(1.-*r/ratio);
		}
		else
		{	*r=1.;
			*a=2./(double)n;
		}
	}
	
	// trival case for 1 interval (ratio can be supported here)
	else
	{	*r=1.;
		*a=2.;
	}
}

// for line with n intervals and ratio r between adjacent element
//  sizes, return ratio of first element to last element size
double GetRatio(int n,double r)
{
	if(DbleEqual(r,1.)) return (double)1.;
	return (n>1) ? 1./pow(r,(double)(n-1)) : (double)1. ;
}



