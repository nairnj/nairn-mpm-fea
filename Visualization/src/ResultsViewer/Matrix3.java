/*
 * Matrix3
 * NairnFEAMPMViz
 * 
 * Created 5/23/2015, John A. Nairn
 * Copyright 2015, RSAC Software, All Rights Reserved
 */


public class Matrix3
{
	private boolean is2D;
	private double m[][];
	
	private static final double M_SQRT3=1.73205080756887729352744634151;   // sqrt(3)
	
	// empty and all zeroed (markes as 2D)
	Matrix3()
	{	m = new double[3][3];
		set(0.);
		is2D = true;
	}
	
	// set to copy of the M matrix
	Matrix3(Matrix3 M)
	{	m = new double[3][3];
		boolean is2D = M.getIs2D();
		m[0][0] = M.get(0,0);
		m[0][1] = M.get(0,1);
		m[1][0] = M.get(1,0);
		m[1][1] = M.get(1,1);
		m[2][2] = M.get(2,2);
		if(!is2D)
		{	m[0][2] = M.get(0,2);
			m[1][2] = M.get(1,2);
			m[2][0] = M.get(2,0);
			m[2][1] = M.get(2,1);
		}
	}
	
	// set 2D elements, non-2D to zero and mark as 2D
	Matrix3(double m00,double m01,double m10,double m11,double m22)
	{	m = new double[3][3];
		m[0][0] = m00;
		m[0][1] = m01;
		m[0][2] = 0.;
		m[1][0] = m10;
		m[1][1] = m11;
		m[1][2] = 0.;
		m[2][0] = 0.;
		m[2][1] = 0.;
		m[2][2] = m22;
		is2D = true;
	}
	
	// set all elements and mark as 2D (even if 3D terms are equal to zero)
	Matrix3(double m00,double m01,double m02,double m10,double m11,double m12,double m20,double m21,double m22)
	{	m = new double[3][3];
		m[0][0] = m00;
		m[0][1] = m01;
		m[0][2] = m02;
		m[1][0] = m10;
		m[1][1] = m11;
		m[1][2] = m12;
		m[2][0] = m20;
		m[2][1] = m21;
		m[2][2] = m22;
		is2D = false;
	}
	
	// Return transpose in a new matrix
	public Matrix3 Transpose()
	{	Matrix3 M = new Matrix3(this);
		M.set(0,1,m[1][0]);
		M.set(1,0,m[0][1]);
		if(!is2D)
		{	M.set(0,2,m[2][0]);
			M.set(1,2,m[2][1]);
			M.set(2,0,m[0][2]);
			M.set(2,1,m[1][2]);
		}
		return M;
	}
	
	// Replace this matrix (M) with M-R
	public void subtract(Matrix3 R)
	{
		if(is2D && R.getIs2D())
		{	m[0][0] -= R.get(0,0);
			m[0][1] -= R.get(0,1);
			m[1][0] -= R.get(1,0);
			m[1][1] -= R.get(1,1);
			m[2][2] -= R.get(2,2);
		}
		else
		{	int i,j;
			for(i=0;i<3;i++)
			{	for(j=0;j<3;j++)
					m[i][j] -= R.get(i,j);
			}
			is2D = false;
		}
	}
	
	// Replace this matrix (M) with M*R
	public void times(Matrix3 R)
	{
		if(is2D && R.getIs2D())
		{	double m00 = m[0][0]*R.get(0,0) + m[0][1]*R.get(1,0);
			double m01 = m[0][0]*R.get(0,1) + m[0][1]*R.get(1,1);
			double m10 = m[1][0]*R.get(0,0) + m[1][1]*R.get(1,0);
			double m11 = m[1][0]*R.get(0,1) + m[1][1]*R.get(1,1);
			double m22 = m[2][2]*R.get(2,2);
			set(m00,m01,m10,m11,m22);
		}
		else
		{	double p[][] = new double[3][3];
			int i,j;
			for(i=0;i<3;i++)
			{	for(j=0;j<3;j++)
				{	p[i][j] = m[i][0]*R.get(0,j) + m[i][1]*R.get(1,j)
								+ m[i][2]*R.get(2,j);
				}
			}
			set(p);
		}
	}
	// Form the triple product R.M.R^T in one step (R need not be a rotation matrix)
	public Matrix3 RMRT(Matrix3 R)
	{	if(is2D && R.getIs2D())
		{	Matrix3 mT = new Matrix3(R.get(0,0)*R.get(0,0)*m[0][0] + R.get(0,1)*R.get(0,0)*m[1][0] + R.get(0,0)*R.get(0,1)*m[0][1] + R.get(0,1)*R.get(0,1)*m[1][1],
					   R.get(0,0)*R.get(1,0)*m[0][0] + R.get(0,1)*R.get(1,0)*m[1][0] + R.get(0,0)*R.get(1,1)*m[0][1] + R.get(0,1)*R.get(1,1)*m[1][1],
					   R.get(1,0)*R.get(0,0)*m[0][0] + R.get(1,1)*R.get(0,0)*m[1][0] + R.get(1,0)*R.get(0,1)*m[0][1] + R.get(1,1)*R.get(0,1)*m[1][1],
					   R.get(1,0)*R.get(1,0)*m[0][0] + R.get(1,1)*R.get(1,0)*m[1][0] + R.get(1,0)*R.get(1,1)*m[0][1] + R.get(1,1)*R.get(1,1)*m[1][1],
					   R.get(2,2)*R.get(2,2)*m[2][2]);
			return mT;
		}
		else
		{	Matrix3 mT = new Matrix3();
			int i,j;
			for(i=0;i<3;i++)
			{	for(j=0;j<3;j++)
				{	double mij = m[i][j];
					mT.set(0,0,mT.get(0,0) + R.get(0,i)*R.get(0,j)*mij);
					mT.set(0,1,mT.get(0,1) + R.get(0,i)*R.get(1,j)*mij);
					mT.set(0,2,mT.get(0,2) + R.get(0,i)*R.get(2,j)*mij);
					mT.set(1,0,mT.get(1,0) + R.get(1,i)*R.get(0,j)*mij);
					mT.set(1,1,mT.get(1,1) + R.get(1,i)*R.get(1,j)*mij);
					mT.set(1,2,mT.get(1,2) + R.get(1,i)*R.get(2,j)*mij);
					mT.set(2,0,mT.get(2,0) + R.get(2,i)*R.get(0,j)*mij);
					mT.set(2,1,mT.get(2,1) + R.get(2,i)*R.get(1,j)*mij);
					mT.set(2,2,mT.get(2,2) + R.get(2,i)*R.get(2,j)*mij);
				}
			}
			return mT;
		}
	}
	
	// Find Eigenvalues of positive definite matrix return in a vector
	// Assumes all are real and positive (any other will fail)
	public double[] Eigenvalues()
	{
		double[] lam = new double[3];
	    
	    if(is2D)
	    {   // solving x^2 + bx + c = 0 (see Numerical Recipes in C, page 156)
	        double b = -(m[0][0]+m[1][1]);
	        double c = m[0][0]*m[1][1] - m[1][0]*m[0][1];
	        double arg = b*b-4.*c;
			if(arg<0.)
			{	// assuming here all matrices are positive definite, which means
				// a negative value should be interpreted as zero
				lam[0] = -0.5*b;
			}
			else
			{	arg = Math.sqrt(arg);
				lam[0] = b>0 ? -0.5*(b+arg) : -0.5*(b-arg) ;
			}
	        lam[1] = c/lam[0];
	        lam[2] = m[2][2];
	    }
	    else
	    {	double mm, c1, c0;
			
			// Determine coefficients of characteristic poynomial. We write
			//       | a   d   f  |
			//  m =  | d*  b   e  |
			//       | f*  e*  c  |
			double de = m[0][1] * m[1][2];                                  // d * e
			double dd = m[0][1]*m[0][1];                                    // d^2
			double ee = m[1][2]*m[1][2];                                    // e^2
			double ff = m[0][2]*m[0][2];                                      // f^2
			mm  = m[0][0] + m[1][1] + m[2][2];
			c1 = (m[0][0]*m[1][1] + m[0][0]*m[2][2] + m[1][1]*m[2][2])
					- (dd + ee + ff);										// a*b + a*c + b*c - d^2 - e^2 - f^2
			c0 = m[2][2]*dd + m[0][0]*ee + m[1][1]*ff - m[0][0]*m[1][1]*m[2][2]
					- 2.0 * m[0][2]*de;										// c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
			
			double p, sqrt_p, q, c, s, phi;
			p = mm*mm - 3.0*c1;
			q = mm*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
			sqrt_p = Math.sqrt(Math.abs(p));
			
			phi = 27.0 * ( 0.25*c1*c1*(p - c1) + c0*(q + 27.0/4.0*c0));
			phi = (1.0/3.0) * Math.atan2(Math.sqrt(Math.abs(phi)), q);
			
			c = sqrt_p*Math.cos(phi);
			s = (1.0/M_SQRT3)*sqrt_p*Math.sin(phi);
			
			lam[1]  = (1.0/3.0)*(mm - c);
			lam[2]  = lam[1] + s;
			lam[0]  = lam[1] + c;
			lam[1] -= s;
	    }
	    
	    return lam;
	}
	
	// Find Eigenvvectors of positive definite matrix return in the Eigenvecs matrix (as column vectors)
	// Provide eigenvalues when call in values for 2D, but not needed for 3D
	// Assumes all are real and positive (any other will fail)
	// No change to target matrix
	public Matrix3 Eigenvectors(double[] Eigenvals)
	{
		Matrix3 Eigenvecs = new Matrix3();
		
	    if(is2D)
		{	double [] vals = dsyev2(m[0][0],m[0][1],m[1][1],Eigenvals[0],Eigenvals[1]);
			Eigenvals[0] = vals[0];
			Eigenvals[1] = vals[1];
			Eigenvecs.set(vals[2],-vals[3],vals[3],vals[2],1.);
	    }
	    
	    else
		{	// set to 3D matrix
			double Q[][] = new double[3][3];
			double A[][] = new double[3][3];
			A[0][0] = m[0][0];
			A[0][1] = m[0][1];
			A[0][2] = m[0][2];
			A[1][1] = m[1][1];
			A[1][2] = m[1][2];
			A[2][2] = m[2][2];
			dsyevq3(A, Q, Eigenvals);			// QL - if finds eigenvalues
			Eigenvecs.set(Q);
	    }
	    
	    return Eigenvecs;
	}

	// Polar decomponsition of F through left stretch matrix
	//  F = VR = Q Lambda (QTR)
	//The target matrix is assumed to be F
	//Function returns V and optionally R = V^-1 F (if pointer is not NULL)
	//	and optionally (lam1,lam2,lam3) in stretches (if not NULL)
	//It does not get Q, but if needed, they are eigenvectors of the
	//	returned V matrix
	public Matrix3 LeftDecompose(Matrix3 R)
	{
		if(is2D)
		{   // 2D has simple formulae for R = ((Fsum,Fdif),(-Fdif,Fsum))
			double Fsum = m[0][0]+m[1][1];
			double Fdif = m[0][1]-m[1][0];
			double denom = Math.sqrt(Fsum*Fsum+Fdif*Fdif);
			Fsum /= denom;
			Fdif /= denom;
    
			// V is F* R*T
			Matrix3 V = new Matrix3(m[0][0]*Fsum+m[0][1]*Fdif,-m[0][0]*Fdif+m[0][1]*Fsum,
					m[1][0]*Fsum+m[1][1]*Fdif,-m[1][0]*Fdif+m[1][1]*Fsum,m[2][2]);
   
			// if R pointer not null, return it too
			if(R!=null)
			{   R.set(Fsum,Fdif,-Fdif,Fsum,1.);
			}
			
			return V;
		}
	
	    // rest is for 3D matrix
	    
	    // Get B and B^2
		Matrix3 B = new Matrix3(this);
		B.times(Transpose());
		Matrix3 B2 = new Matrix3(B);
	    B2.times(B);
	    
	    // Eigenvalues of B are lamda^2
	    double[] Eigenvals = B.Eigenvalues();
	    double lam1 = Math.sqrt(Eigenvals[0]);
	    double lam2 = Math.sqrt(Eigenvals[1]);
	    double lam3 = Math.sqrt(Eigenvals[2]);
	    
	    // invariants of V
	    double i1 = lam1+lam2+lam3;
	    double i2 = lam1*lam2+lam1*lam3+lam2*lam3;
	    double i3 = lam1*lam2*lam3;
	    
	    // set coefficients
	    double d1 = 1./(i1*i2-i3);
	    double c2 = -d1;                    // coefficient of B2
	    double c1 = (i1*i1-i2)*d1;          // coefficient of B
	    double cI = i1*i3*d1;               // coefficient of I
	    
	    // Get V = (1/d1)*(-B^2 + (i1*i1-i2)*B + i1*i3*I)
	    Matrix3 V = new Matrix3(c2*B2.get(0,0)+c1*B.get(0,0)+cI, c2*B2.get(0,1)+c1*B.get(0,1),    c2*B2.get(0,2)+c1*B.get(0,2),
				  c2*B2.get(1,0)+c1*B.get(1,0),    c2*B2.get(1,1)+c1*B.get(1,1)+cI, c2*B2.get(1,2)+c1*B.get(1,2),
				  c2*B2.get(2,0)+c1*B.get(2,0),    c2*B2.get(2,1)+c1*B.get(2,1),    c2*B2.get(2,2)+c1*B.get(2,2)+cI);
	    
	    // if R pointer not NULL, find R too
	    if(R!=null)
	    {   c1 = 1/i3;                      // coefficient of B
	        double cV = -i1*c1;             // coefficient of V
	        cI = i2*c1;                     // coefficient of I
	        
	        // Get Vinv = (1/i3)*(B - i1*V + i2*I)
	        double p[][] = new double[3][3];
	        p[0][0] = c1*B.get(0,0)+cV*V.get(0,0)+cI;
	        p[0][1] = c1*B.get(0,1)+cV*V.get(0,1);
	        p[0][2] = c1*B.get(0,2)+cV*V.get(0,2);
	        p[1][0] = c1*B.get(1,0)+cV*V.get(1,0);
	        p[1][1] = c1*B.get(1,1)+cV*V.get(1,1)+cI;
	        p[1][2] = c1*B.get(1,2)+cV*V.get(1,2);
	        p[2][0] = c1*B.get(2,0)+cV*V.get(2,0);
	        p[2][1] = c1*B.get(2,1)+cV*V.get(2,1);
	        p[2][2] = c1*B.get(2,2)+cV*V.get(2,2)+cI;
	     
	        // R = V^-1 F
	        R.set(p);
	        R.times(this);
	    }
		
		return V;
	}

	// set all elements to a constant
	public void set(double constant)
	{	int i,j;
		for(i=0;i<3;i++)
		{	for(j=0;j<3;j++)
				m[i][j] = constant;
		}
		is2D = constant==0. ? true : false ;
	}
	
	// set to 2D matrix with elements
	public void set(double m00,double m01,double m10,double m11,double m22)
	{	m[0][0] = m00;
		m[0][1] = m01;
		m[1][0] = m10;
		m[1][1] = m11;
		m[2][2] = m22;
		m[0][2] = m[2][0] = 0.;
		m[1][2] = m[2][1] = 0.;
		is2D = true;
	}
	
	// set all elements another array
	public void set(double[][] p)
	{	int i,j;
		for(i=0;i<3;i++)
		{	for(j=0;j<3;j++)
				m[i][j] = p[i][j];
		}
		is2D = false ;
	}
	// set one element (0 based indices) (does not handle is2D)
	public void set(int row,int col,double value) { m[row][col] = value; }
	
	// get one element
	public double get(int row,int col) { return m[row][col]; }
	
	// 2D settting
	public boolean getIs2D() { return is2D; }
	public void setIs2D(boolean setting) { is2D = setting; }
	
	// determinant
	public double determinant()
	{	if(is2D)
		{	return m[2][2]*(m[0][0]*m[1][1]-m[1][0]*m[0][1]);
		}
		return m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
			-m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2])
			+m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);
	}

	// ----------------------------------------------------------------------------
	// Calculates the eigensystem of a real symmetric 2x2 matrix
	//	    [ A  B ]
	//	    [ B  C ]
	// in the form
	//	    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
	//	    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
	// On call, rt1 and rt2 are the eigen values.
	// return is {rt1,rt2,cs,sn) and rt1 and rt2 might be reordered
	// such that rt1 >= rt2. The orthongal eigenvectors are
	// e1 = (cs, sn), e2 = (-sn, cs)
	// ----------------------------------------------------------------------------
	private double[] dsyev2(double A, double B, double C, double rt1, double rt2)
	{
		double df = A - C;
		double rt = Math.sqrt(df*df + 4.0*B*B);
		double t;
		
		// order the eigenvalues
		if(rt2 > rt1)
		{	double temp = rt1;
			rt1 = rt2;
			rt2 = temp;
		}
		
		// Calculate eigenvectors
		double cs,sn;
		if (df > 0.0)
			cs = df + rt;
		else
			cs = df - rt;
		
		if (Math.abs(cs) > 2.0*Math.abs(B))
		{	t   = -2.0 * B / cs;
			sn = 1.0 / Math.sqrt(1.0 + t*t);
			cs = t * sn;
		}
		else if (Math.abs(B) == 0.0)
		{	cs = 1.0;
			sn = 0.0;
		}
		else
		{	t   = -0.5 * cs / B;
			cs = 1.0 / Math.sqrt(1.0 + t*t);
			sn = t * cs;
		}
		
		if (df > 0.0)
		{	t   = cs;
			cs = -sn;
			sn = t;
		}
		
		double[] response = new double[4];
		response[0] = rt1;
		response[1] = rt2;
		response[2] = cs;
		response[3] = sn;
		return response;
	}


	private int dsyevq3(double A[][], double Q[][], double w[])
	// ----------------------------------------------------------------------------
	// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
	// matrix A using the QL algorithm with implicit shifts, preceded by a
	// Householder reduction to tridiagonal form.
	// The function accesses only the diagonal and upper triangular parts of A.
	// The access is read-only.
	// ----------------------------------------------------------------------------
	// Parameters:
	//   A: The symmetric input matrix
	//   Q: Storage buffer for eigenvectors
	//   w: Storage buffer for eigenvalues
	// ----------------------------------------------------------------------------
	// Return value:
	//   0: Success
	//  -1: Error (no convergence)
	// ----------------------------------------------------------------------------
	// Dependencies:
	//   dsytrd3()
	// ----------------------------------------------------------------------------
	{
		int n = 3;
		double e[] = new double[3];    // The third element is used only as temporary workspace
		double g, r, p, f, b, s, c, t; // Intermediate storage
		int nIter;
		int m;
		
		// Transform A to real tridiagonal form by the Householder method
		dsytrd3(A, Q, w, e);
		
		// Calculate eigensystem of the remaining real symmetric tridiagonal matrix
		// with the QL method
		//
		// Loop over all off-diagonal elements
		for (int l=0; l < n-1; l++)
		{
			nIter = 0;
			while(true)
			{
				// Check for convergence and exit iteration loop if off-diagonal
				// element e(l) is zero
				for (m=l; m <= n-2; m++)
				{
					g = Math.abs(w[m])+Math.abs(w[m+1]);
					if (Math.abs(e[m]) + g == g)
						break;
				}
				if (m == l)
					break;
				
				if (nIter++ >= 30)
					return -1;
				
				// Calculate g = d_m - k
				g = (w[l+1] - w[l]) / (e[l] + e[l]);
				r = Math.sqrt(g*g + 1.0);
				if (g > 0)
					g = w[m] - w[l] + e[l]/(g + r);
				else
					g = w[m] - w[l] + e[l]/(g - r);
				
				s = c = 1.0;
				p = 0.0;
				for (int i=m-1; i >= l; i--)
				{
					f = s * e[i];
					b = c * e[i];
					if (Math.abs(f) > Math.abs(g))
					{
						c      = g / f;
						r      = Math.sqrt(c*c + 1.0);
						e[i+1] = f * r;
						c     *= (s = 1.0/r);
					}
					else
					{
						s      = f / g;
						r      = Math.sqrt(s*s + 1.0);
						e[i+1] = g * r;
						s     *= (c = 1.0/r);
					}
					
					g = w[i+1] - p;
					r = (w[i] - g)*s + 2.0*c*b;
					p = s * r;
					w[i+1] = g + p;
					g = c*r - b;
					
					// Form eigenvectors
					for (int k=0; k < n; k++)
					{
						t = Q[k][i+1];
						Q[k][i+1] = s*Q[k][i] + c*t;
						Q[k][i]   = c*Q[k][i] - s*t;
					}
				}
				w[l] -= p;
				e[l]  = g;
				e[m]  = 0.0;
			}
		}
		
		return 0;
	}

	private void dsytrd3(double A[][], double Q[][], double d[], double e[])
	// ----------------------------------------------------------------------------
	// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
	// (unitary) Householder transformations:
	//	            [ d[0]  e[0]       ]
	//	    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
	//	            [       e[1]  d[2] ]
	// The function accesses only the diagonal and upper triangular parts of
	// A. The access is read-only.
	// ---------------------------------------------------------------------------
	{
		int n = 3;
		double u[] = new double[n];
		double q[] = new double[n];
		double omega, f;
		double K, h, g;
		
		// Initialize Q to the identitity matrix
		for (int i=0; i < n; i++)
		{
			Q[i][i] = 1.0;
			for (int j=0; j < i; j++)
				Q[i][j] = Q[j][i] = 0.0;
		}
		
		// Bring first row and column to the desired form
		h = A[0][1]*A[0][1] + A[0][2]*A[0][2];
		if (A[0][1] > 0)
			g = -Math.sqrt(h);
		else
			g = Math.sqrt(h);
		e[0] = g;
		f    = g * A[0][1];
		u[1] = A[0][1] - g;
		u[2] = A[0][2];
		
		omega = h - f;
		if (omega > 0.0)
		{
			omega = 1.0 / omega;
			K     = 0.0;
			for (int i=1; i < n; i++)
			{
				f    = A[1][i] * u[1] + A[i][2] * u[2];
				q[i] = omega * f;                  // p
				K   += u[i] * f;                   // u* A u
			}
			K *= 0.5 * omega*omega;
			
			for (int i=1; i < n; i++)
				q[i] = q[i] - K * u[i];
			
			d[0] = A[0][0];
			d[1] = A[1][1] - 2.0*q[1]*u[1];
			d[2] = A[2][2] - 2.0*q[2]*u[2];
			
			// Store inverse Householder transformation in Q
			for (int j=1; j < n; j++)
			{
				f = omega * u[j];
				for (int i=1; i < n; i++)
					Q[i][j] = Q[i][j] - f*u[i];
			}
			
			// Calculate updated A[1][2] and store it in e[1]
			e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
		}
		else
		{
			for (int i=0; i < n; i++)
				d[i] = A[i][i];
			e[1] = A[1][2];
		}
	}

}
