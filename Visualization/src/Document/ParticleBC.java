/*******************************************************************
	GridDispBC.java
	NairnFEAMPMViz

	Created by John Nairn on 5/8/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class ParticleBC extends BoundaryCondition
{
	private int face;
	
	// initialize
	ParticleBC(int partNum,int bcDof,int bcFace,int theID,double theVal,double theArg)
	{	super(partNum,bcDof,theID,theVal,theArg,0.0);
		face = bcFace;			// 0 for load or 1-6 for traction
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
			case N_DIRECTION:
				if(face==1)
				{	if(value>=0.)
						transform.rotate(-Math.PI/2.);
					else
						transform.rotate(Math.PI/2.);
				}
				else if(face==2)
				{	if(value<0.)
						transform.rotate(Math.PI);
				}
				else if(face==3)
				{	if(value<0.)
						transform.rotate(-Math.PI/2.);
					else
						transform.rotate(Math.PI/2.);
				}
				else if(face==4)
				{	if(value>=0.)
						transform.rotate(Math.PI);
				}
				break;
			case T1_DIRECTION:
				if(face==4)
				{	if(value>=0.)
						transform.rotate(-Math.PI/2.);
					else
						transform.rotate(Math.PI/2.);
				}
				else if(face==1)
				{	if(value<0.)
						transform.rotate(Math.PI);
				}
				else if(face==2)
				{	if(value<0.)
						transform.rotate(-Math.PI/2.);
					else
						transform.rotate(Math.PI/2.);
				}
				else if(face==3)
				{	if(value>=0.)
						transform.rotate(Math.PI);
				}
				break;
			default:
				break;
		}
		
		// for traction, move to face
		double dx=0.;
		double dy=0.;
		if(face>0)
		{	Point2D.Double mpcell = pv.getMpSize(mpt,doc);
			if(face==1)
				dy = -mpcell.getY();
			else if(face==2)
				dx = mpcell.getX();
			else if(face==3)
				dy = mpcell.getY();
			else if(face==4)
				dx = -mpcell.getX();
		}
		
		// apply transformation
		arrow.transform(transform);
		
		transform=new AffineTransform();
		transform.translate((float)(mpt.x+dx),(float)(mpt.y+dy));
		arrow.transform(transform);
		
		pv.fillShape(arrow);
	}
}
