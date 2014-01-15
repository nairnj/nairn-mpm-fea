/*******************************************************************
	ElementTriangle.java
	NairnFEAMPMViz

	Created by John Nairn on 9/10/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.geom.*;
import java.util.*;

public class ElementTriangle extends ElementBase
{
	// initialize
	ElementTriangle(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	// This is called when preparing a mesh plot. It should store sub elements
	// if needed for plot.
	public void allocateSubElements(ArrayList<NodalPoint>nodes,int density,boolean displaced) throws Exception
	{
		// exit if already have correct subelements
		if(subValues.size()==density*density && displaced==subPathsDisplaced) return;
		subPathsDisplaced=displaced;

		// get parameters
		NodalPoint[] ndpt=new NodalPoint[9];
		getNodalPoints(ndpt,nodes);
		double delta=1./(double)density;
		
		// clear all sub arrays
		int i,j,size=density*density;
		subPaths.clear();
		subValues.clear();
		subColors.clear();
		for(i=0;i<size;i++)
		{	subValues.add(new Double(0.));
			subColors.add(Color.red);
		}
		
		// loop over sub elements
		Point2D.Double xiEta=new Point2D.Double();
		Point2D.Double[] pts=new Point2D.Double[4];
		for(i=0;i<density;i++)
		{   // top corner
			xiEta.x=0.;
			xiEta.y=i*delta;
			pts[0]=findCoords(xiEta,ndpt,displaced);
			xiEta.y+=delta;
			pts[3]=findCoords(xiEta,ndpt,displaced);
			xiEta.y-=delta;
			for(j=0;j<density-i;j++)
			{	xiEta.x+=delta;
				pts[1]=findCoords(xiEta,ndpt,displaced);
				xiEta.y+=delta;
				pts[2]=findCoords(xiEta,ndpt,displaced);
				xiEta.y-=delta;
				
				// make two triangle paths
				GeneralPath newPath=new GeneralPath();
				newPath.moveTo((float)pts[0].x,(float)pts[0].y);
				newPath.lineTo((float)pts[1].x,(float)pts[1].y);
				newPath.lineTo((float)pts[3].x,(float)pts[3].y);
				newPath.lineTo((float)pts[0].x,(float)pts[0].y);
				subPaths.add(newPath);
				
				// second triangle unless at the end
				if(j<density-i-1)
				{	GeneralPath extraPath=new GeneralPath();
					extraPath.moveTo((float)pts[1].x,(float)pts[1].y);
					extraPath.lineTo((float)pts[2].x,(float)pts[2].y);
					extraPath.lineTo((float)pts[3].x,(float)pts[3].y);
					extraPath.lineTo((float)pts[1].x,(float)pts[1].y);
					subPaths.add(extraPath);
				}
				
				// reuse right side points
				pts[3]=pts[2];
				pts[0]=pts[1];
			}
		}
	}
		
	// create sub elements paths
	public void setPlotValues(int density) throws Exception
	{
		if(subValues.size()==0) return;			// can skip if not plotting

		// get parameters
		double delta=1./(double)density;
		double hdelta=delta/2.;
		double qdelta=hdelta/2.;
		
		// loop over sub elements
		int i,j,vindex=0;
		Point2D.Double xiEta=new Point2D.Double();
		for(i=0;i<density;i++)
		{   // find top
			xiEta.x=qdelta;
			xiEta.y=i*delta + qdelta;
			for(j=0;j<density-i;j++)
			{   // get value at (xi,eta)
				subValues.set(vindex,new Double(findValueAt(xiEta)));
				vindex++;

				// second triangle if needed
				if(j<density-i-1)
				{	xiEta.x+=hdelta;
					xiEta.y+=hdelta;
					subValues.set(vindex,new Double(findValueAt(xiEta)));
					xiEta.x-=hdelta;
					xiEta.y-=hdelta;
					vindex++;
				}
				
				// next point
				xiEta.x+=delta;
			}
		}
	}
	
	//--------------------------------------------------------------
	// Accessors probably needing overrides
	//--------------------------------------------------------------

	// nodes, sides, and mid side nodes (override as needed)
	public int getNumberNodes() { return 3; }
	public int getNumberSides() { return 3; }
	
}
