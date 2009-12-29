/*******************************************************************
	CSTriangle.java
	NairnFEAMPMViz

	Created by John Nairn on 9/10/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class CSTriangle extends ElementTriangle
{
	// initialize
	CSTriangle(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	// shape function only
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] eNodes)
	{	double ai=eNodes[1].x*eNodes[2].y-eNodes[2].x*eNodes[1].y;
		double aj=eNodes[2].x*eNodes[0].y-eNodes[0].x*eNodes[2].y;
		double ak=eNodes[0].x*eNodes[1].y-eNodes[1].x*eNodes[0].y;
		double twicearea=(ai+aj+ak);
		double bi=eNodes[1].y-eNodes[2].y;
		double bj=eNodes[2].y-eNodes[0].y;
		double bk=eNodes[0].y-eNodes[1].y;
		double ci=eNodes[2].x-eNodes[1].x;
		double cj=eNodes[0].x-eNodes[2].x;
		double ck=eNodes[1].x-eNodes[0].x;
		sfxn[0]=(ai+bi*xiEta.x+ci*xiEta.y)/twicearea;
		sfxn[1]=(aj+bj*xiEta.x+cj*xiEta.y)/twicearea;
		sfxn[2]=(ak+bk*xiEta.x+ck*xiEta.y)/twicearea;
	}
	
	// general shape function BMatrix evaluation
	public void getShapeBMatrix(Point2D.Double xiEta,double[] xiDeriv,double[] etaDeriv,double[] asbe,NodalPoint[] eNodes,boolean derivs)
	{	double ai=eNodes[1].x*eNodes[2].y-eNodes[2].x*eNodes[1].y;
		double aj=eNodes[2].x*eNodes[0].y-eNodes[0].x*eNodes[2].y;
		double ak=eNodes[0].x*eNodes[1].y-eNodes[1].x*eNodes[0].y;
		double twicearea=(ai+aj+ak);
		double bi=eNodes[1].y-eNodes[2].y;
		double bj=eNodes[2].y-eNodes[0].y;
		double bk=eNodes[0].y-eNodes[1].y;
		double ci=eNodes[2].x-eNodes[1].x;
		double cj=eNodes[0].x-eNodes[2].x;
		double ck=eNodes[1].x-eNodes[0].x;
		
		double sfxn[]=new double[4];
		sfxn[0]=(ai+bi*xiEta.x+ci*xiEta.y)/twicearea;
		sfxn[1]=(aj+bj*xiEta.x+cj*xiEta.y)/twicearea;
		sfxn[2]=(ak+bk*xiEta.x+ck*xiEta.y)/twicearea;
		
		// elements of B Matrix
		xiDeriv[0]=bi/twicearea;
		xiDeriv[1]=bj/twicearea;
		xiDeriv[2]=bk/twicearea;
		etaDeriv[0]=ci/twicearea;
		etaDeriv[1]=cj/twicearea;
		etaDeriv[2]=ck/twicearea;
		
		if(derivs) return;
		
		// Ni/r if r≠0
		if(xiEta.x!=0.)
		{	asbe[0]=sfxn[0]/xiEta.x;
			asbe[1]=sfxn[1]/xiEta.x;
			asbe[2]=sfxn[2]/xiEta.x;
		}
		else
		{	asbe[0]=0.;
			asbe[1]=0.;
			asbe[2]=0.;
		}
	}
	
	// find dimensionless location of a point
	public void getXiEta(Point2D.Double xiEta,Point2D.Double pt,NodalPoint[] ndpt) throws Exception
	{	xiEta.x=pt.x;
		xiEta.y=pt.y;
	}

	// find coordinates for given dimensionless position
	public Point2D.Double findCoords(Point2D.Double xiEta,NodalPoint[] ndpt,boolean displaced) throws Exception
	{	double zeta=1.-xiEta.x-xiEta.y;
		if(displaced)
		{	return new Point2D.Double((ndpt[0].x+ndpt[0].dispx)*xiEta.x + (ndpt[1].x+ndpt[1].dispx)*xiEta.y
					+ (ndpt[2].x+ndpt[2].dispx)*zeta,(ndpt[0].y+ndpt[0].dispy)*xiEta.x
					+ (ndpt[1].y+ndpt[1].dispy)*xiEta.y + (ndpt[2].y+ndpt[2].dispy)*zeta);
		}
		else
		{	return new Point2D.Double(ndpt[0].x*xiEta.x + ndpt[1].x*xiEta.y + ndpt[2].x*zeta,
							ndpt[0].y*xiEta.x + ndpt[1].y*xiEta.y + ndpt[2].y*zeta);
		}
	}

	// find value for given dimensionless position
	public double findValueAt(Point2D.Double xiEta) throws Exception
	{	double zeta=1.-xiEta.x-xiEta.y;
		return plotValues[0]*xiEta.x + plotValues[1]*xiEta.y + plotValues[2]*zeta;
	}

	//--------------------------------------------------------------
	// Accessors probably needing overrides
	//--------------------------------------------------------------

	// nodes, sides, and mid side nodes (override as needed)
	public int getNumberNodes() { return 3; }
	
}

