/*******************************************************************
	IsoTriangle.java
	NairnFEAMPMViz

	Created by John Nairn on 9/13/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class IsoTriangle extends ElementTriangle
{
	// initialize
	IsoTriangle(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	// shape function only
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] eNodes)
	{
		double pxi3=1.-xiEta.x-xiEta.y;
		sfxn[0]=xiEta.x*(2.*xiEta.x-1.);
		sfxn[1]=xiEta.y*(2.*xiEta.y-1.);
		sfxn[2]=pxi3*(2.*pxi3-1.);
		sfxn[3]=4.*xiEta.x*xiEta.y;
		sfxn[4]=4.*xiEta.y*pxi3;
		sfxn[5]=4.*xiEta.x*pxi3;
	}
	
	// general shape function BMatrix evaluation
	public void getShapeBMatrix(Point2D.Double xiEta,double[] xiDeriv,double[] etaDeriv,double[] asbe,NodalPoint[] eNodes,boolean derivs)
	{
		double sfxn[]=new double[6];
		double pxi3=1.-xiEta.x-xiEta.y;
		sfxn[0]=xiEta.x*(2.*xiEta.x-1.);
		sfxn[1]=xiEta.y*(2.*xiEta.y-1.);
		sfxn[2]=pxi3*(2.*pxi3-1.);
		sfxn[3]=4.*xiEta.x*xiEta.y;
		sfxn[4]=4.*xiEta.y*pxi3;
		sfxn[5]=4.*xiEta.x*pxi3;
		
		xiDeriv[0]=4*xiEta.x-1.;
		xiDeriv[1]=0.;
		xiDeriv[2]=-4.*pxi3+1.;
		xiDeriv[3]=4.*xiEta.y;
		xiDeriv[4]=-4.*xiEta.y;
		xiDeriv[5]=4.*(pxi3-xiEta.x);
		etaDeriv[0]=0.;
		etaDeriv[1]=4*xiEta.y-1.;
		etaDeriv[2]=-4.*pxi3+1.;
		etaDeriv[3]=4.*xiEta.x;
		etaDeriv[4]=4.*(pxi3-xiEta.y);
		etaDeriv[5]=-4.*xiEta.x;
		
		if(derivs) return;

		// Find Jacobian Matrix
		int i;
		double[][] jac={{0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.}};
		jac[1][1]=0.;
		jac[1][2]=0.;
		jac[2][1]=0.;
		jac[2][2]=0.;
		for(i=0;i<6;i++)
		{	jac[1][1]+=xiDeriv[i]*eNodes[i].x;
			jac[1][2]+=xiDeriv[i]*eNodes[i].y;
			jac[2][1]+=etaDeriv[i]*eNodes[i].x;
			jac[2][2]+=etaDeriv[i]*eNodes[i].y;
		}
		double detjac=jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1];

		// invert Jacobian
		double temp2,temp1=jac[1][1]/detjac;
		jac[1][1]=jac[2][2]/detjac;
		jac[2][2]=temp1;
		jac[1][2]=-jac[1][2]/detjac;
		jac[2][1]=-jac[2][1]/detjac;
			
		// For radial position for axisymmetric analyses
		double asr=0.;
		for(i=0;i<6;i++)
			asr=asr+sfxn[i]*eNodes[i].x;
				
		/* Load J(-1)(1,1) ∂Ni/∂xi + J(-1)(1,2) ∂Ni/∂eta into xiDeriv[i]
				J(-1)(2,1) ∂Ni/∂xi + J(-1)(2,2) ∂Ni/∂eta into etaDeriv[i]
				Ni/r into asbe[i] */
		for(i=0;i<6;i++)
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
	
	//--------------------------------------------------------------
	// Accessors probably needing overrides
	//--------------------------------------------------------------

	// centroid for this element and coordinates of the nodes
	public Point2D.Double getCentroid() { return new Point2D.Double(0.3333,0.3333); }
	
	// nodes, sides, and mid side nodes (override as needed)
	public int getNumberNodes() { return 6; }
	public boolean hasMidSideNodes() { return true; }
	
	// Get dimensionless coordinate of node (0 based)
	public Point2D.Double getNodeXiEta(int i,NodalPoint[] eNodes)
	{	switch(i)
		{	case 0:
				return new Point2D.Double(1.,0.);
			case 1:
				return new Point2D.Double(0.,1.);
			case 2:
				return new Point2D.Double(0.,0.);
			case 3:
				return new Point2D.Double(.5,.5);
			case 4:
				return new Point2D.Double(0.,.5);
			case 5:
				return new Point2D.Double(.5,0.);
			default:
				break;
		}
		return new Point2D.Double(0.,0.);
	}
}

