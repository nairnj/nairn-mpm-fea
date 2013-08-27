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
    where nband is bandwidth of the matrix. In other words
            a[k][i] = ai,i+k-1

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
 
	Some personal notes found at JANOSU-4-197
***********************************************************/

int gelbnd(double **a,int n,int nband,double *r,double *work,int iflag)
{
    int i;
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
	// Reducing matrix a into a = LU where
	//   L is lower diagonal matrix (with ones on the diagonal)
	//   U is upper diagonal matrix
	// Because a is symmetric, only need to find U and it replaces
	//	the a elements as it is calculated
    if(iflag==0)
    {   
    	// Make copy of the diagonal elements
        for(i=1;i<=n;i++) work[i] = a[1][i];
        
        // Forward reduction of coefficient matrix in remaining n-1 rows
        for(i=1; i<=n-1; i++)
		{	// Check for over reduced or bad matrix
            double Uii = a[1][i];                    // Now Uii since row i has been reduced
			if(Uii == 0.)
			{	// singular matrix - give up with an error
				return 1;
			}
			else if(work[i]!=0.)
			{	// compare to initial diagonal element if possible
				// if small, set warning flag, but continue
				if(fabs(Uii/work[i])<1e-12) ierr=-1;
			}
			
			// We are reducing rows k = i+1 to n, but because this matrix is banded
			//     we can stop when aik = 0 or k-i<=nband-1 or k<=i+nband-1
			//     we index j = 2 to nband, then k = j+i-1 to jend+i-1
			//     except j is truncated near end of matrix
			
			int imin1 = i-1;                      // = i-1
            int jend = fmin(nband,n-imin1);       // last non-zero element in this row
			
            // Loop over rows that need to be reduced
			// Changes only row k  and aik, input is aik, Uii
#pragma omp parallel for
            for(int j=2; j<=jend; j++)
			{	// aik = a[k-i+1][i] = a[j][i]
            	double Uim,Uik = a[j][i];			// Now Uik since row i has been reduced, but not converted to Lki yet
				
				// if zero, nothing to do, otherwise do a reduction
                if(Uik != 0.)
				{	double UikOverUii = Uik/Uii;
                    int k = j+imin1;				// reducing row k
					
					// in row k, need to reduce elements m = k+1 to n to
					//     akm = akm - (Uik/Uii)*Uim
					// but because matrix is banded, we can only need m-k<=nband-1 or m<=k+nband-1
					// Let j2 run from j to jend and counter ind from 1 to jend-j+1, then m = k+ind-1 = k+j2-j
					//	   akm = a[m-k+1][k] = a[ind][k]
					//     aim = a[m-i+1][i] = a[j2][i]
                    int ind = 0;
					for(int j2=j; j2<=jend; j2++)
                    {	ind++;						// m-k+1
                        Uim = a[j2][i];				// Now Ui,j2-i+1 = Uim since row i has been reduced
                        if(Uim != 0.)
                            a[ind][k] -= UikOverUii*Uim;
                    }
					
					// convert element in previous row to Lki
                    //a[j][i] = UikOverUii;		// Now Uik/Uii = Lki
                }
            }
			
#pragma omp parallel for
            for(int j=2; j<=jend; j++)
			{	// convert element in previous row to Lki
				a[j][i] /= Uii;		// Now Uik/Uii = Lki
			}
			
			// Elements through a[][i+1] now fully reduced to row i+1 of U (with off-diagonals to Lki = Uik/Uii)
        }
		
		// As reduced, a[1][i] is diagonal element Uii and a[j][i] is Li+j-1,i
    }
	   
    // Forward reduction of the right side of the equation
	// Solve Ly = r
    for(i=1; i<=n; i++)
    {	int imin1 = i-1;
        int jend = fmin(nband,n-imin1);
		
		// Let k = i+j-1, then k from i+1 to end
#pragma omp parallel for
        for(int j=2; j<=jend; j++)
        {   double Lki = a[j][i];			// Lki or loop Li+1,i to end
            if(Lki!=0.)
			{	r[imin1+j] -= Lki*r[i];		// changing r[i+1] to r[i-1+jend] in parallel
            }
        }
		
		// r[i+1] has been converted to yi
		
		// No longer need yi in r[i], so scale now by diagonal r[i] = yi/Uii
        r[i] /= a[1][i];
    }

    // Solve for unknowns by back substitution and put the answer back into r
	// Solve U x = y
	// First xn = yn/Unn = rn from above forward reduction
	// Find rest xi = ri - Sum_{k=i+1,n} Lki xk
    for(i=n-1; i>=1; i--)
    {	int imin1 = i-1;
        int jend = fmin(nband,n-imin1);
		
		// Let k=i+j-1 and loop k=i+1 to end
        double sum = 0.;
#pragma omp parallel for reduction(+:sum)
        for(int j=2; j<=jend; j++)
        {   double Lki = a[j][i];				// From decomposition is Lki
            if(Lki!=0.)
                sum += Lki*r[imin1+j];			// summand Lki*xk now in r[i+1] to r[i-1+jend] in parallel
        }
		
		// subtract sum to get xi
        r[i] -= sum;
    }
    
    return ierr;
}



