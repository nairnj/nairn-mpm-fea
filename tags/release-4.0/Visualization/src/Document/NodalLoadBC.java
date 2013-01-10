/*******************************************************************
	NodalLoadBC.java
	NairnFEAMPMViz

	Created by John Nairn on 9/6/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class NodalLoadBC extends GridDispBC
{
	// initialize
	NodalLoadBC(int nodeNum,int bcDof,double theVal,double theAngle)
	{	super(nodeNum,bcDof,FEA_LOADBC,theVal,0.0,theAngle);
	}

	// draw the boundary conditino
	public void stroke(MeshPlotView pv,ResultsDocument doc)
	{	
		// Path with generic horizontal arrow  ---->
		GeneralPath arrow=makeArrow(pv.xyBounds,1.0f);
		NodalPoint nd=doc.nodes.get(node);
		Point2D.Double tip=new Point2D.Double(0.,0.);
		double length=arrow.getBounds2D().getWidth();
		
		// rotate as needed
		double rotationAngle=-angle;
		AffineTransform transform=new AffineTransform();
		switch(dof)
		{   case X_DIRECTION:
				if(value<0.) transform.rotate(Math.PI);
				if(Math.abs(rotationAngle)<=45.)
				{	if(value<0.)
					{	if(nd.x>pv.xyBounds.getCenterX()) tip.x=length;
					}
					else if(nd.x<pv.xyBounds.getCenterX())
						tip.x=-length;
				}
				else if(Math.abs(rotationAngle)>=135.)
				{	if(value>0.)
					{	if(nd.x>pv.xyBounds.getCenterX()) tip.x=-length;
					}
					else if(nd.x<pv.xyBounds.getCenterX())
						tip.x=length;
				}
				else if(rotationAngle<0.)
				{	if(value>0.)
					{	if(nd.y>pv.xyBounds.getCenterY()) tip.x=-length;
					}
					else if(nd.y<pv.xyBounds.getCenterY())
						tip.x=length;
				}
				else
				{	if(value<0.)
					{	if(nd.y>pv.xyBounds.getCenterY()) tip.x=length;
					}
					else if(nd.y<pv.xyBounds.getCenterY())
						tip.x=-length;
				}
				break;
			case Y_DIRECTION:
				if(value<0.)
					transform.rotate(-Math.PI/2.);
				else
					transform.rotate(Math.PI/2.);
				if(Math.abs(rotationAngle)<=45.)
				{	if(value<0.)
					{	if(nd.y>pv.xyBounds.getCenterY()) tip.y=length;
					}
					else if(nd.y<pv.xyBounds.getCenterY())
						tip.y=-length;
				}
				else if(Math.abs(rotationAngle)>=135.)
				{	if(value>0.)
					{	if(nd.y>pv.xyBounds.getCenterY()) tip.y=-length;
					}
					else if(nd.y<pv.xyBounds.getCenterY())
						tip.y=length;
				}
				else if(rotationAngle<0.)
				{	if(value>0.)
					{	if(nd.x<pv.xyBounds.getCenterX()) tip.y=-length;
					}
					else if(nd.x>pv.xyBounds.getCenterX())
						tip.y=length;
				}
				else
				{	if(value<0.)
					{	if(nd.x<pv.xyBounds.getCenterX()) tip.y=length;
					}
					else if(nd.x>pv.xyBounds.getCenterX())
						tip.y=-length;
				}
				break;
			default:
				break;
		}
		double radAngle=rotationAngle*Math.PI/180.;
		transform.rotate(radAngle);
		arrow.transform(transform);
		
		// translate to node (with tip consideration)
		double ndx=nd.x+Math.cos(radAngle)*tip.x-Math.sin(radAngle)*tip.y;
		double ndy=nd.y+Math.sin(radAngle)*tip.x+Math.cos(radAngle)*tip.y;
		transform=new AffineTransform();
		transform.translate((float)ndx,(float)ndy);
		arrow.transform(transform);
		
		pv.fillShape(arrow);
	}
}
