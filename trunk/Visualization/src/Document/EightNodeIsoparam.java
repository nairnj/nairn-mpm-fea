/*
 * EightNodeIsoparam.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 9/7/07.
 * Copyright 2007 RSAC Software. All rights reserved.
 */

import java.awt.geom.*;

public class EightNodeIsoparam extends ElementBase
{
	// Local globals
	private static double[] xii = {-1.,1.,1.,-1.};
	private static double[] eti = {-1.,-1.,1.,1.};

	// initialize
	EightNodeIsoparam(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	// shape function only
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] ndpt)
	{	double temp1,temp2;
		int i,ind1,ind2;
		
		// shape function and derivatives
		for(i=0;i<4;i++)
		{	temp1=(1.+xii[i]*xiEta.x)/4.;
			temp2=(1.+eti[i]*xiEta.y)/4.;
			sfxn[i]=4.*temp1*temp2;
		}
		sfxn[4]=(1.-xiEta.x*xiEta.x)*(1.-xiEta.y)/2.;
		sfxn[5]=(1.+xiEta.x)*(1.-xiEta.y*xiEta.y)/2.;
		sfxn[6]=(1.-xiEta.x*xiEta.x)*(1.+xiEta.y)/2.;
		sfxn[7]=(1.-xiEta.x)*(1.-xiEta.y*xiEta.y)/2.;
		for(i=0;i<4;i++)
		{	if(i==0)
				ind1=7;
			else
				ind1=i+3;
			ind2=i+4;
			sfxn[i]-=(sfxn[ind1]+sfxn[ind2])/2.;
		}
	}
	
	// general shape function BMatrix evaluation
	public void getShapeBMatrix(Point2D.Double xiEta,double[] xiDeriv,double[] etaDeriv,double[] asbe,NodalPoint[] eNodes,boolean derivs)
	{	double temp1,temp2;
		int i,ind1,ind2;
		double sfxn[]=new double[9];
		
		// shape function and derivatives
		for(i=0;i<4;i++)
		{	temp1=(1.+xii[i]*xiEta.x)/4.;
			temp2=(1.+eti[i]*xiEta.y)/4.;
			sfxn[i]=4.*temp1*temp2;
			xiDeriv[i]=xii[i]*temp2;
			etaDeriv[i]=eti[i]*temp1;
		}
		sfxn[4]=(1.-xiEta.x*xiEta.x)*(1.-xiEta.y)/2.;
		sfxn[5]=(1.+xiEta.x)*(1.-xiEta.y*xiEta.y)/2.;
		sfxn[6]=(1.-xiEta.x*xiEta.x)*(1.+xiEta.y)/2.;
		sfxn[7]=(1.-xiEta.x)*(1.-xiEta.y*xiEta.y)/2.;
		xiDeriv[4]=-xiEta.x*(1.-xiEta.y);
		etaDeriv[4]=-(1.-xiEta.x*xiEta.x)/2.;
		xiDeriv[5]=(1.-xiEta.y*xiEta.y)/2.;
		etaDeriv[5]=-xiEta.y*(1.+xiEta.x);
		xiDeriv[6]=-xiEta.x*(1.+xiEta.y);
		etaDeriv[6]=(1.-xiEta.x*xiEta.x)/2.;
		xiDeriv[7]=-(1.-xiEta.y*xiEta.y)/2.;
		etaDeriv[7]=-xiEta.y*(1.-xiEta.x);
		for(i=0;i<4;i++)
		{	if(i==0)
				ind1=7;
			else
				ind1=i+3;
			ind2=i+4;
			sfxn[i]-=(sfxn[ind1]+sfxn[ind2])/2.;
			xiDeriv[i]-=(xiDeriv[ind1]+xiDeriv[ind2])/2.;
			etaDeriv[i]-=(etaDeriv[ind1]+etaDeriv[ind2])/2.;
		}
		
		if(derivs) return;
		
		// Get B matrix elements in xiDeriv[] and etaDeriv[] and Ni/r in sfxn[]
		// Find Jacobian Matrix
		double[][] jac={{0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.}};
		jac[1][1]=0.;
		jac[1][2]=0.;
		jac[2][1]=0.;
		jac[2][2]=0.;
		for(i=0;i<8;i++)
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
		for(i=0;i<8;i++)
			asr+=sfxn[i]*eNodes[i].x;
			
		/* Load J(-1)(1,1) ‚àÇNi/‚àÇxi + J(-1)(1,2) ‚àÇNi/‚àÇeta into xiDeriv[i]
				J(-1)(2,1) ‚àÇNi/‚àÇxi + J(-1)(2,2) ‚àÇNi/‚àÇeta into etaDeriv[i]
				Ni/r into asbe[i] */
		for(i=0;i<8;i++)
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
			case 4:
				return new Point2D.Double(0.,-1.);
			case 5:
				return new Point2D.Double(1.,0.);
			case 6:
				return new Point2D.Double(0.,1.);
			case 7:
				return new Point2D.Double(-1.,0.);
			default:
				break;
		}
		return new Point2D.Double(0.,0.);
	}

	//--------------------------------------------------------------
	// Accessors probably needing overrides
	//--------------------------------------------------------------

	// nodes, sides, and mid side nodes (override as needed)
	public int getNumberNodes() { return 8; }
	public boolean hasMidSideNodes() { return true; }
	
}

