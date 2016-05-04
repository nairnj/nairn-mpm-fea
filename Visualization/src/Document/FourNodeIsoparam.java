/*******************************************************************
	FourNodeIsoparam.java
	NairnFEAMPMViz

	Created by John Nairn on 2/27/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class FourNodeIsoparam extends ElementBase
{
	// Local globals
	private static double[] xii = {-1.,1.,1.,-1.};
	private static double[] eti = {-1.,-1.,1.,1.};

	// initialize
	FourNodeIsoparam(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}

	/* Find dimensionless coordinates analytically if possible
		methods - on input, xi and eta should be initial
		guess in case needed for numerical
		Note: analytical possible if parallelogram
	*/
	public void getXiEta(Point2D.Double xiEta,Point2D.Double pt,NodalPoint[] ndpt) throws Exception
	{
		double xdel,ydel,dxx,dxy,dyy,dyx,det;
		
		// are edges parallel wrt x coordinate?
		xdel=Math.abs(ndpt[2].x-ndpt[1].x);
		ydel=Math.abs(ndpt[3].x-ndpt[0].x);
		if(!ElementBase.DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10)
		{   super.getXiEta(xiEta,pt,ndpt);
			return;
		}
				
		// ara edges parallel wrt y coordinate?
		xdel=Math.abs(ndpt[2].y-ndpt[1].y);
		ydel=Math.abs(ndpt[3].y-ndpt[0].y);
		if(!ElementBase.DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10)
		{   super.getXiEta(xiEta,pt,ndpt);
			return;
		}
		
		// analytical solution for parallelograms
		xdel=4.*pt.x-(ndpt[0].x+ndpt[1].x+ndpt[2].x+ndpt[3].x);
		dxx=ndpt[1].x+ndpt[2].x-ndpt[0].x-ndpt[3].x;
		dxy=ndpt[2].x+ndpt[3].x-ndpt[0].x-ndpt[1].x;
		ydel=4.*pt.y-(ndpt[0].y+ndpt[1].y+ndpt[2].y+ndpt[3].y);
		dyx=ndpt[1].y+ndpt[2].y-ndpt[0].y-ndpt[3].y;
		dyy=ndpt[2].y+ndpt[3].y-ndpt[0].y-ndpt[1].y;
		det=dxx*dyy-dxy*dyx;
		
		xiEta.x=(dyy*xdel-dxy*ydel)/det;
		xiEta.y=(dxx*ydel-dyx*xdel)/det;
	}

	// shape functions only
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] ndpt)
	{	double temp1,temp2;
		int i;
		
		// shape functions
		for(i=0;i<4;i++)
		{	temp1=(1.+xii[i]*xiEta.x);
			temp2=(1.+eti[i]*xiEta.y);
			sfxn[i]=0.25*temp1*temp2;
		}
	}
	
	// general shape function BMatrix evaluation
	public void getShapeBMatrix(Point2D.Double xiEta,double[] xiDeriv,double[] etaDeriv,double[] asbe,NodalPoint[] eNodes,boolean derivs)
	{
		// shape function and derivatives (dimensionless)
		double sfxn[]=new double[9];
		double temp1,temp2;
		int i;
		for(i=0;i<4;i++)
		{	temp1=(1.+xii[i]*xiEta.x)/4.;
			temp2=(1.+eti[i]*xiEta.y)/4.;
			sfxn[i]=4.*temp1*temp2;
			xiDeriv[i]=xii[i]*temp2;
			etaDeriv[i]=eti[i]*temp1;
		}
		if(derivs) return;
		
		// Get B matrix elements in xiDeriv[] and etaDeriv[] and Ni/r in asbe[]
		// Find Jacobian Matrix
		double[][] jac={{0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.}};
		for(i=0;i<4;i++)
		{	jac[1][1]+=xiDeriv[i]*eNodes[i].x;
			jac[1][2]+=xiDeriv[i]*eNodes[i].y;
			jac[2][1]+=etaDeriv[i]*eNodes[i].x;
			jac[2][2]+=etaDeriv[i]*eNodes[i].y;
		}
		double detjac=jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1];
		temp1=jac[1][1]/detjac;
		jac[1][1]=jac[2][2]/detjac;
		jac[2][2]=temp1;
		jac[1][2]=-jac[1][2]/detjac;
		jac[2][1]=-jac[2][1]/detjac;
		
		// For radial position for axisymmetric analyses
		double asr=0.;
		for(i=0;i<4;i++)
			asr=asr+sfxn[i]*eNodes[i].x;
			
		/* Load J(-1)(1,1) dNi/dxi + J(-1)(1,2) dNi/deta into xiDeriv[i]
				J(-1)(2,1) dNi/dxi + J(-1)(2,2) dNi/deta into etaDeriv[i]
				Ni/r into asbe[i] */
		for(i=0;i<4;i++)
		{	temp1=xiDeriv[i];
			temp2=etaDeriv[i];
			xiDeriv[i]=jac[1][1]*temp1 + jac[1][2]*temp2;
			etaDeriv[i]=jac[2][1]*temp1 + jac[2][2]*temp2;
			if(asr!=0.)
				asbe[i]=sfxn[i]/asr;
			else
				asbe[i]=0.;
		}
	}
	
	// Get dimensionless coordinate of node (0 based)
	public Point2D.Double getNodeXiEta(int i,NodalPoint[] eNodes)
	{	switch(i)
		{	case 0:
				return new Point2D.Double(-1.,-1.);
			case 1:
				return new Point2D.Double(1.,-1.);
			case 2:
				return new Point2D.Double(1.,1.);
			case 3:
				return new Point2D.Double(-1.,1.);
			default:
				break;
		}
		return new Point2D.Double(0.,0.);
	}
}
