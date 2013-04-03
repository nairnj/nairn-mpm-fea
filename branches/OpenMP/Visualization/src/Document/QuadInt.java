/*******************************************************************
	QuadInt.java
	NairnFEAMPMViz

	Created by John Nairn on 10/4/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class QuadInt extends Interface
{
	// initialize
	QuadInt(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	// general shape function evaluation - does not get derivatives
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] ndpt)
	{
		// shape functions
		sfxn[0]=(xiEta.x*xiEta.x-xiEta.x)/2.;
		sfxn[1]=1.-xiEta.x*xiEta.x;
		sfxn[2]=(xiEta.x*xiEta.x+xiEta.x)/2.;
		sfxn[3]=-sfxn[2];
		sfxn[4]=-sfxn[1];
		sfxn[5]=-sfxn[0];
	}

	public int getNumberNodes() { return 6; }
	public boolean hasMidSideNodes() { return true; }
}
