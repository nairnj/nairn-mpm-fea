/********************************************************************************
    Gauss.cpp
    NairnFEA
    
    Created by John Nairn on Mar 13 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

// prototype
int gelbnd(double **,int,int,double *,double *,int);

/***********************************************************
    Subroutine to solve a linear system defined by
            ax=r
    where a is band diagonal and symmetric stored in
    band storage format such that aij=a[j-i+1][i] and
    storage includes only elements with
            1 <= j-i+1 <= nband
    where nband is bandwidth of the matrix

    Input parameters
        a: matrix in band symmetric storage
        n: order of linear system
        nband: bandwidth (<=n)
        r: right hand side of linear system
        work: work vector of length >=n
        iflag: 0 to reduce matrix or 1 if already reduced

    Output parameters
        a: forward reduced matrix
        r: solution to linear system
        return value: 0 if no errors
                -1 if reduction is getting too severe
                1 if diagonal or pivot is 0

    Method
        Gaussian elimination with no pivoting which means
    diagonal elements better be OK
***********************************************************/

int gelbnd(double **a,int n,int nband,double *r,double *work,int iflag)
{
    int i,j,j2,k,ind,imin1,iplus1,jend;
    double temp,temp2,factor,scale,sum;
    int ierr=0;
    
    // Check for trival unit band length problem
    if(nband==1)
    {	for(i=1;i<=n;i++)
        {   if(a[1][i]==0.) return 1;
            r[i]/=a[1][i];
        }
        return ierr;
    }

    // If iflag=0 then do matrix reduction
    if(iflag==0)
    {	// Make copy of the diagonal elements
        for(i=1;i<=n;i++) work[i]=a[1][i];

        // Forward reduction of coefficient matrix
        for(i=1;i<=n-1;i++)
        {   imin1=i-1;
            iplus1=i+1;
            jend=fmin(nband,n-imin1);

            // Check for over reduced or bad matrix 
            temp=a[1][i];
            if(temp==0.)
                return 1;
            else
            {	temp2=temp/work[i];
                if(temp2<0.) temp2=-temp2;
                if(temp2<1e-12) ierr=-1;
            }
            
            for(j=2;j<=jend;j++)
            {	factor=a[j][i];
                if(factor!=0.)
                {   scale=factor/temp;
                    k=j+imin1;
                    ind=0;
                    for(j2=j;j2<=jend;j2++)
                    {	ind=ind+1;
                        factor=a[j2][i];
                        if(factor!=0.)
                            a[ind][k]-=scale*factor;
                    }
                    a[j][i]=scale;
                }
            }
        }
    }
	   
    // Forward reduction of the right side of the equation
    for(i=1;i<=n;i++)
    {	imin1=i-1;
        jend=fmin(nband,n-imin1);
        for(j=2;j<=jend;j++)
        {   factor=a[j][i];
            if(factor!=0.)
            {	ind=imin1+j;
                r[ind]-=factor*r[i];
            }
        }
        r[i]/=a[1][i];
    }

    /* Solve for unknowns by back substitution and put the
                    answer back into r */
    for(i=n-1;i>=1;i--)
    {	imin1=i-1;
        jend=fmin(nband,n-imin1);
        sum=0.;
        for(j=2;j<=jend;j++)
        {   factor=a[j][i];
            if(factor!=0.)
                sum+=factor*r[imin1+j];
        }
        r[i]-=sum;
    }
    
    return ierr;
}



