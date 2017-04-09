/*******************************************************************
	GridDispBC.java
	NairnFEAMPMViz

	Created by John Nairn on 5/8/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.Color;
import java.awt.geom.*;

public class ParticleBC extends BoundaryCondition
{
	private int face;
	private int BCType;
	private float lineWidth=1.0f;
	
	// initialize
	ParticleBC(int partNum,int bcDof,int bcFace,int theID,double theVal,double theArg)
	{	super(partNum,bcDof,theID,theVal,theArg,0.0);
		face = bcFace;			// 0 for load or 1-6 for traction
		BCType=0;				// might change to TEMPERATURE_DIR or CONCENTRATION_DIR
	}

	// draw the boundary conditino
	public void stroke(MeshPlotView pv,ResultsDocument doc)
	{	// exit if not on yet
		if(doc.currentTime()<argument)
		{	if(style==CONSTANT_VALUE || style==LINEAR_VALUE)
				return;
		}
		
		GeneralPath arrow;
		float length;
		if(BCType!=0)
		{	// jagged line
			length=(float)(Math.max(pv.xyBounds.getWidth(),pv.xyBounds.getHeight())*BC_SIZE)*0.6f;
			float awidth=length/6.f;
			
			// Path for spring shape
			arrow=new GeneralPath();
			arrow.moveTo((float)0.,(float)0.);
			arrow.lineTo(0.111*length,awidth);
			arrow.lineTo(0.333*length,-awidth);
			arrow.lineTo(0.556*length,awidth);
			arrow.lineTo(0.778*length,-awidth);
			arrow.lineTo(length,awidth);
		}
		else
		{	// Path with generic horizontal arrow  ---->
			arrow=makeArrow(pv.xyBounds,0.6f);
			length=(float)(arrow.getBounds2D().getWidth());
		}
		
		// get material point
		MaterialPoint mpt=doc.mpmPoints.get(node);
		
		// rotate as needed
		AffineTransform transform=new AffineTransform();
		if(BCType!=0)
		{	// for temperature and concentration, point out from the face
			if(face==1)
				transform.rotate(-Math.PI/2.);
			else if(face==3)
				transform.rotate(Math.PI/2.);
			else if(face==4)
				transform.rotate(Math.PI);
		}
		else
		{	switch(dof)
			{   case X_DIRECTION:
					if(value<0.)
					{	transform.rotate(Math.PI);
						if(face>0 && face!=4)
							transform.translate(-length,0.0f);
					}
					else if(face>0 && face==4)
						transform.translate(-length,0.0f);
					break;
				case Y_DIRECTION:
					if(value<0.)
					{	transform.rotate(-Math.PI/2.0f);
						if(face>0 && face==3)
							transform.translate(-length,0.0f);
					}
					else
					{	transform.rotate(Math.PI/2.);
						if(face>0 && face==1)
							transform.translate(-length,0.0f);
					}
					break;
				case N_DIRECTION:
					if(face==1)
						transform.rotate(-Math.PI/2.0f);
					else if(face==3)
						transform.rotate(Math.PI/2.0f);
					else if(face==4)
						transform.rotate(Math.PI);
					if(value<0.)
					{	transform.rotate(Math.PI);
						transform.translate(-length,0.0f);
					}
					break;
				case T1_DIRECTION:
					if(face==2)
						transform.rotate(Math.PI/2.0f);
					else if(face==3)
						transform.rotate(Math.PI);
					else if(face==4)
						transform.rotate(-Math.PI/2.0f);
					if(value<0.)
						transform.rotate(Math.PI);
					break;
				default:
					break;
			}
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
		
		Color bcColor;
		if(BCType==TEMPERATURE_DIR)
		{	bcColor=NFMVPrefs.getPrefColor(NFMVPrefs.tempBCColorKey,NFMVPrefs.tempBCColorDef);
			pv.setLineWidth(lineWidth);
			pv.drawColor(bcColor);
			pv.strokeShape(arrow);
		}
		else if(BCType==CONCENTRATION_DIR)
		{	bcColor=NFMVPrefs.getPrefColor(NFMVPrefs.concBCColorKey,NFMVPrefs.concBCColorDef);
			pv.setLineWidth(lineWidth);
			pv.drawColor(bcColor);
			pv.strokeShape(arrow);
		}
		else
		{	bcColor=NFMVPrefs.getPrefColor(NFMVPrefs.meshLineColorKey,NFMVPrefs.meshLineColorDef);
			pv.drawColor(bcColor);
			pv.fillShape(arrow);
		}
	}
	
	// set the BC typed (used to known about heat and concentration flux BCs)
	public void setBCType(int newType) { BCType = newType; }
}
