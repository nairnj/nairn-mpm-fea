/*******************************************************************
	LinearInt.java
	NairnFEAMPMViz

	Created by John Nairn on 10/4/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class LinearInt extends Interface
{
	// initialize
	LinearInt(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}

	// general shape function evaluation - does not get derivatives
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] ndpt)
	{
		// shape functions
		sfxn[0]=(1.-xiEta.x)/2.;
		sfxn[1]=(1.+xiEta.x)/2.;
		sfxn[2]=-sfxn[1];
		sfxn[3]=-sfxn[0];
	}
	
}
