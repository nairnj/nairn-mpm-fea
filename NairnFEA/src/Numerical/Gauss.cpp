/********************************************************************************
    Gauss.cpp
    NairnFEA
    
    Created by John Nairn on Mar 13 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

// prototype
int gelbnd(double **,int,int,double *,double *,int);

/*  Subroutine to solve a linear system defined by
            ax=r
    where a is band diagonal and symmetric stored in
    band storage format such that aij=a[i][j-i+1] and
    storage includes only elements with
            1 <= j-i+1 <= nband
    where nband is bandwidth of the matrix. In other words
            a[i][k] = ai,i+k-1

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
                 1 if diagonal is 0

    Method
        Gaussian elimination with no pivoting which means
        diagonal elements better be OK
 
	Some personal notes found at JANOSU-4-197 and JANOSU-7-1
*/

int gelbnd(double **a,int n,int nband,double *r,double *work,int iflag)
{
    int i,j,j2,k,imin1,jend,ind;
    int ierr=0;
    double Uii,Uim,Uik,UikOverUii,Lki,yi,sum;
    
    // Check for trival unit band length problem
    if(nband==1)
    {	for(i=1;i<=n;i++)
        {   if(a[i][1]==0.) return 1;
            r[i]/=a[i][1];
        }
        return ierr;
    }

    // If iflag=0 then do matrix reduction
	// Reducing matrix a into a = LU where
	//   L is lower diagonal matrix (with ones on the diagonal)
	//   U is upper diagonal matrix
	// Because a is symmetric, only need to find U and it replaces
	//	the a elements as it is calculated
    if(iflag==0)
    {   
    	// Make copy of the diagonal elements
        for(i=1;i<=n;i++) work[i] = a[i][1];
        
        // Forward reduction of coefficient matrix in remaining n-1 rows
        for(i=1; i<=n-1; i++)
		{	// Check for over reduced or bad matrix
            Uii = a[i][1];                    // Now Uii since row i has been reduced
			if(Uii == 0.)
			{	// singular matrix - give up with an error
				return 1;
			}
			else if(work[i]!=0.)
			{	// compare to initial diagonal element if possible
				// if small, set warning flag, but continue
				if(fabs(Uii/work[i])<1.e-12) ierr=-1;
			}
			
			// We are reducing rows k = i+1 to n, but because this matrix is banded
			//     we can stop when aik = 0 or k-i<=nband-1 or k<=i+nband-1
			//     we index j = 2 to nband, then k = j+i-1 to jend+i-1
			//     except j is truncated near end of matrix
			
			imin1 = i-1;                      // = i-1
            jend = fmin(nband,n-imin1);       // last non-zero element in this row
			
            // Loop over rows that need to be reduced
			// Changes only row k  and aik, input is aik, Uii
#pragma omp parallel for private(k,ind,j2,Uim,Uik,UikOverUii)
            for(j=2; j<=jend; j++)
			{	// aik = a[i][k-i+1] = a[i][j]
            	Uik = a[i][j];			// Now Uik since row i has been reduced, but not converted to Lki yet
				
				// if zero, nothing to do, otherwise do a reduction
                if(Uik != 0.)
				{	UikOverUii = Uik/Uii;
                    k = j+imin1;				// reducing row k
					
					// in row k, need to reduce elements m = k to n to
					//     akm = akm - (Uik/Uii)*Uim
					// but because matrix is banded, we can only need m-k<=nband-1 or m<=k+nband-1
					// Let j2 run from j to jend and counter ind from 1 to jend-j+1, then m = k+ind-1 = k+j2-j
					//	   akm = a[k][m-k+1] = a[k][ind]
					//     aim = a[i][m-i+1] = a[i][j2]
                    ind = 0;
					for(j2=j; j2<=jend; j2++)
                    {	ind++;						// m-k+1
                        Uim = a[i][j2];				// Now Ui,j2-i+1 = Uim since row i has been reduced
                        if(Uim != 0.)
                            a[k][ind] -= UikOverUii*Uim;
                    }
                }
            }
			
            // separate loop to avoid parallel conflicts in above loop
#pragma omp parallel for
            for(j=2; j<=jend; j++)
			{	// convert element in previous row to Lki
				a[i][j] /= Uii;		// Now Uik/Uii = Lki
			}
			
			// Elements through a[i+1][] now fully reduced to row i+1 of U (with off-diagonals to Lki = Uik/Uii)
        }
		
		// As reduced, a[1][i] is diagonal element Uii and a[j][i] is Li+j-1,i
    }
	   
    // Forward reduction of the right side of the equation
	// Solve Ly = r
    for(i=1; i<=n; i++)
    {	imin1 = i-1;
        jend = fmin(nband,n-imin1);
        yi = r[i];
		
		// Let k = i+j-1, then k from i+1 to end
#pragma omp parallel for private(Lki)
        for(j=2; j<=jend; j++)
        {   Lki = a[i][j];                  // Lki or loop Li+1,i to end
            if(Lki!=0.)
			{	r[imin1+j] -= Lki*yi;		// changing r[i+1] to r[i-1+jend] in parallel
            }
        }
		
		// r[i+1] has been converted to yi
		
		// No longer need yi in r[i], so scale now by diagonal r[i] = yi/Uii
        // which is used in back substitution to account for Lki being in aik
        r[i] /= a[i][1];
    }

    // Solve for unknowns by back substitution and put the answer back into r
	// Solve Ux = y
	// First xn = yn/Unn = rn from above forward reduction
	// Find rest xi = ri - Sum_{k=i+1,n} Lki xk
    for(i=n-1; i>=1; i--)
    {	imin1 = i-1;
        jend = fmin(nband,n-imin1);
		
		// Let k=i+j-1 and loop k=i+1 to end
        sum = 0.;
#pragma omp parallel for reduction(+:sum) private(Lki)
        for(j=2; j<=jend; j++)
        {   Lki = a[i][j];                  // From decomposition is Lki
            if(Lki!=0.)
                sum += Lki*r[imin1+j];		// summand Lki*xk now in r[i+1] to r[i-1+jend] in parallel
        }
		
		// subtract sum to get xi
        r[i] -= sum;
    }
    
    return ierr;
}



