/*******************************************************************
	GridDispBC.java
	NairnFEAMPMViz

	Created by John Nairn on 5/8/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class ParticleBC extends BoundaryCondition
{
	// initialize
	ParticleBC(int partNum,int bcDof,int theID,double theVal,double theArg)
	{	super(partNum,bcDof,theID,theVal,theArg,0.0);
	}

	// draw the boundary conditino
	public void stroke(MeshPlotView pv,ResultsDocument doc)
	{	// exit if not on yet
		if(doc.currentTime()<argument)
		{	if(style==CONSTANT_VALUE || style==LINEAR_VALUE)
				return;
		}
		
		// Path with generic horizontal arrow  ---->
		GeneralPath arrow=makeArrow(pv.xyBounds,1.0f);
		MaterialPoint mpt=doc.mpmPoints.get(node);
		
		// rotate as needed
		AffineTransform transform=new AffineTransform();
		switch(dof)
		{   case X_DIRECTION:
				if(value<0.)
					transform.rotate(Math.PI);
				break;
			case Y_DIRECTION:
				if(value<0.)
					transform.rotate(-Math.PI/2.);
				else
					transform.rotate(Math.PI/2.);
				break;
			default:
				break;
		}
		arrow.transform(transform);
		
		transform=new AffineTransform();
		transform.translate((float)mpt.x,(float)mpt.y);
		arrow.transform(transform);
		
		pv.fillShape(arrow);
	}
}
