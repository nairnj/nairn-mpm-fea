/*******************************************************************
	GridDispBC.java
	NairnFEAMPMViz

	Created by John Nairn on 5/8/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.Color;
import java.awt.geom.*;

public class GridDispBC extends BoundaryCondition
{
	protected double angle2;
	
	public static float lineWidth=(float)2.0;
	
	// initialize
	GridDispBC(int nodeNum,int bcDof,int theID,double theVal,double theArg,double theAngle,double xyzAngle)
	{	super(nodeNum,bcDof,theID,theVal,theArg,theAngle);
		if(dof==SKEWXY_DIRECTION) dof = XY_SKEWED_INPUT;
		angle2 = xyzAngle;
	}

	// draw the boundary conditino
	public void stroke(MeshPlotView pv,ResultsDocument doc)
	{	// exit if not on yet
		if(doc.currentTime()<argument)
		{	if(style==CONSTANT_VALUE || style==LINEAR_VALUE)
				return;
		}
		
		// FEA rotation not drawn
		if(style==FEA_ROTATEBC) return;
		
		// get size
		double height=Math.max(pv.xyBounds.getWidth(),pv.xyBounds.getHeight())*BC_SIZE*6./7.;
		double width=height/4.;
		double diam=height/6;
		NodalPoint nd=doc.nodes.get(node);
		
		GeneralPath theBC=new GeneralPath();
		if(dof==TEMPERATURE_DIR)
		{	// diamond on the node
			height = 0.6*height;
			width = 0.5*height;
			theBC.moveTo(0.f,0.f);
			theBC.lineTo((float)(-width),(float)(-0.5*height));
			theBC.lineTo((float)0.,(float)(-height));
			theBC.lineTo((float)width,(float)(-0.5*height));
			theBC.lineTo(0.f,0.f);
		}
		else if(dof==CONCENTRATION_DIR)
		{	// diamond on the node
			height = 0.5*height;
			width = height;
			theBC.append(new Ellipse2D.Double(-0.5*width,-height,width,height), false);
		}
		else
		{	theBC.moveTo(0.f,0.f);
			theBC.lineTo((float)(-width),(float)(-height));
			theBC.lineTo((float)width,(float)(-height));
			theBC.lineTo(0.f,0.f);
			theBC.append(new Ellipse2D.Double(-width,-height-diam,diam,diam),false);
			theBC.append(new Ellipse2D.Double(width-diam,-height-diam,diam,diam),false);
		}
		
		double rotationAngle;
		if(style==FEA_DISPBC)
			rotationAngle=-angle;
		else
			rotationAngle = (dof==XY_SKEWED_INPUT) ? -angle : 0.;
		AffineTransform transform=new AffineTransform();
		switch(dof)
		{   case X_DIRECTION:
			case XY_SKEWED_INPUT:
				if(Math.abs(rotationAngle)<=45.)
				{	if(nd.x<pv.xyBounds.getCenterX())
						transform.rotate(-Math.PI/2.);
					else
						transform.rotate(Math.PI/2.);
				}
				else if(Math.abs(rotationAngle)>=135.)
				{	if(nd.x<pv.xyBounds.getCenterX())
						transform.rotate(Math.PI/2.);
					else
						transform.rotate(-Math.PI/2.);
				}
				else if(rotationAngle>0.)
				{	if(nd.y>pv.xyBounds.getCenterY())
						transform.rotate(Math.PI/2.);
					else
						transform.rotate(-Math.PI/2.);
				}
				else
				{	if(nd.y>pv.xyBounds.getCenterY())
						transform.rotate(-Math.PI/2.);
					else
						transform.rotate(Math.PI/2.);
				}
				break;
			case Y_DIRECTION:
				if(Math.abs(rotationAngle)<=45.)
				{	if(nd.y>pv.xyBounds.getCenterY())
						transform.rotate(Math.PI);
				}
				else if(Math.abs(rotationAngle)>=135.)
				{	if(nd.y<pv.xyBounds.getCenterY())
						transform.rotate(Math.PI);
				}
				else if(rotationAngle>0.)
				{	if(nd.x<pv.xyBounds.getCenterX())
						transform.rotate(Math.PI);
				}
				else
				{	if(nd.x>pv.xyBounds.getCenterX())
						transform.rotate(Math.PI);
				}
				break;
			case TEMPERATURE_DIR:
			case CONCENTRATION_DIR:
			{	double xb = (nd.x-pv.xyBounds.getMinX())/pv.xyBounds.getWidth(),xt = 1.-xb;
				double yb = (nd.y-pv.xyBounds.getMinY())/pv.xyBounds.getHeight(),yt = 1.-yb;
				// which is smallest and assume tha edge
				if(xb<=xt && xb<=yb && xb<=yt)
				{	// closest to left edge
					transform.rotate(-Math.PI/2.);
				}
				else if(xt<=xb && xt<=yb && xt<=yt)
				{	// closest to right edge
					transform.rotate(Math.PI/2.);
				}
				else if(yt<=yb && yt<=xb && yt<=xt)
				{	// closest to right edge
					transform.rotate(Math.PI);
				}
				break;
			}
			default:
				break;
		}
		transform.rotate(rotationAngle*Math.PI/180.);
		theBC.transform(transform);
		
		transform=new AffineTransform();
		transform.translate(nd.x,nd.y);
		theBC.transform(transform);
		
		Color bcColor;
		if(dof==TEMPERATURE_DIR)
		{	bcColor=NFMVPrefs.getPrefColor(NFMVPrefs.tempBCColorKey,NFMVPrefs.tempBCColorDef);
			pv.setLineWidth(lineWidth);
		}
		else if(dof==CONCENTRATION_DIR)
		{	bcColor=NFMVPrefs.getPrefColor(NFMVPrefs.concBCColorKey,NFMVPrefs.concBCColorDef);
			pv.setLineWidth(lineWidth);
		}
		else
		{	bcColor=NFMVPrefs.getPrefColor(NFMVPrefs.meshLineColorKey,NFMVPrefs.meshLineColorDef);
			pv.setLineWidth(ElementBase.lineWidth);
		}
		pv.drawColor(bcColor);
		pv.strokeShape(theBC);
	}
	
	// check if this object rotates node of another object
	public boolean rotates(int nodeNum)
	{	return style==FEA_ROTATEBC && node==nodeNum-1;
	}
}
