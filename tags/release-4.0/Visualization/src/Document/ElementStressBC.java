/*******************************************************************
	ElementStressBC.java
	NairnFEAMPMViz

	Created by John Nairn on 9/14/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class ElementStressBC extends GridDispBC
{
	public static final int NORMAL_DIRECTION=1;
	public static final int SHEAR_DIRECTION=2;
	private int dir;
	
	// initialize
	ElementStressBC(int elemNnum,int face,String orient,double str1,double str2,double str3)
	{	super(elemNnum,face-1,FEA_ELEMSTRESSBC,str1,str2,str3);
	
		// 0-based elemNum stored in node
		// 0-based face number stored in dof
		// three stresses stored in value,argument,angle
		// 1 normal and 2 shear, stored in dof
		dir=SHEAR_DIRECTION;
		if(orient.length()>0)
		{	if(orient.charAt(0)=='N')
				dir=NORMAL_DIRECTION;
		}
	}
	
	// draw the boundary conditino
	public void stroke(MeshPlotView pv,ResultsDocument doc)
	{	
		// stroke each half of the shape
		for(int half=1;half<=2;half++)
		{	// Path with generic horizontal arrow  ---->
			GeneralPath arrow=makeArrow(pv.xyBounds,0.8f);
			double length=arrow.getBounds2D().getWidth();
		
			// transformation
			Point2D.Double origin=new Point2D.Double();
			Point2D.Double slope=new Point2D.Double();
			ElementBase theElem=doc.elements.get(node);
			theElem.getFaceInfo(dof,half,origin,slope,doc.nodes);
			
			double rotAngle=0.;
			if(ElementBase.DbleEqual(slope.x,0.))
				rotAngle = slope.y>0. ? Math.PI/2. : -Math.PI/2 ;
			else
			{	if(slope.x>0.)
					rotAngle = Math.atan(slope.y/slope.x);
				else
					rotAngle = Math.atan(slope.y/slope.x)-Math.PI;
			}
			
			// value
			double stress = half==1 ? (value+argument)/2. : (argument+angle)/2.;
			
			// rotate as needed
			AffineTransform transform=new AffineTransform();
			double move=0.1;
			switch(dir)
			{	case SHEAR_DIRECTION:
					if(stress>0.)
						move=-0.5;
					else
					{	rotAngle+=Math.PI;
						move=0.5;
					}
					break;
				case NORMAL_DIRECTION:
					if(stress>0.)
						rotAngle-=Math.PI/2.;
					else
					{	rotAngle+=Math.PI/2.;
						move=1.1;
					}
					break;
				default:
					break;
			}
			transform.rotate(rotAngle);
			arrow.transform(transform);
			
			// translate to element face
			transform=new AffineTransform();
			transform.translate((float)(origin.x+slope.x/2.),(float)(origin.y+slope.y/2.));
			if(dir==NORMAL_DIRECTION)
			{	double norm=move*length/Math.sqrt(slope.x*slope.x+slope.y*slope.y);
				transform.translate((float)(slope.y*norm),(float)(-slope.x*norm));
			}
			else if(dir==SHEAR_DIRECTION)
			{	double norm=length/Math.sqrt(slope.x*slope.x+slope.y*slope.y);
				transform.translate((float)(move*slope.x*norm),(float)(move*slope.y*norm));
				transform.translate((float)(0.5*slope.y*norm),(float)(-0.5*slope.x*norm));
			}
			arrow.transform(transform);
			
			// draw it
			pv.fillShape(arrow);
		}
	}
}

