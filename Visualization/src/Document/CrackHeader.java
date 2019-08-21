/*******************************************************************
	CrackHeader.java
	NairnFEAMPMViz

	Created by John Nairn on 3/11/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.util.*;
import java.awt.*;

public class CrackHeader
{
	public ArrayList<CrackSegment> segments;
	
	//---------------------------------------------------------------------
	// initialize
	//---------------------------------------------------------------------
	
	CrackHeader()
	{
		segments=new ArrayList<CrackSegment>(100);
	}
	
	// add crack segment
	public void add(CrackSegment cs)
	{	segments.add(cs);
	}

	//---------------------------------------------------------------------
	// draw
	//---------------------------------------------------------------------
	
	// draw crack data
	public void stroke(MeshPlotView pv,ResultsDocument doc,boolean showCrackPlanes,boolean showCrackSurfaces)
	{
		if(segments.size()==0) return;
		
		if(showCrackPlanes)
		{	Color cplaneColor = NFMVPrefs.getPrefColor(NFMVPrefs.cplaneColorKey,NFMVPrefs.cplaneColorDef);
			pv.drawColor(cplaneColor);
			drawCrack(pv,CrackSegment.PLANE_POS);
		}
		
		if(showCrackSurfaces)
		{	Color csurfaceColor = NFMVPrefs.getPrefColor(NFMVPrefs.caboveColorKey,NFMVPrefs.caboveColorDef);
			pv.drawColor(csurfaceColor);
			drawCrack(pv,CrackSegment.ABOVE_POS);
			csurfaceColor = NFMVPrefs.getPrefColor(NFMVPrefs.cbelowColorKey,NFMVPrefs.cbelowColorDef);
			pv.drawColor(csurfaceColor);
			drawCrack(pv,CrackSegment.BELOW_POS);
			
			// draw traction lines
			pv.setLineWidth(2.f);
			for(int i=0;i<segments.size();i++)
			{	CrackSegment seg=segments.get(i);
				if(seg.tractionMaterial>0)
				{	pv.moveTo(seg.xpos[CrackSegment.ABOVE_POS],seg.ypos[CrackSegment.ABOVE_POS]);
					pv.lineTo(seg.xpos[CrackSegment.BELOW_POS],seg.ypos[CrackSegment.BELOW_POS]);
				}
			}
			pv.setLineWidth(1.f);
		}
	}
	
	// draw on crack line
	public void drawCrack(MeshPlotView pv,int which)
	{	CrackSegment seg=segments.get(0);
		pv.moveTo(seg.xpos[which],seg.ypos[which]);
		
		int i;
		for(i=1;i<segments.size();i++)
		{	seg=segments.get(i);
			pv.lineTo(seg.xpos[which],seg.ypos[which]);
		}
	}
	
	// get all coordinates to a crack surface or crack plane
	// side == PLANE_POS, ABOVE_POX, or BELOW_POC
	public void getSurface(int side,ArrayList<Double>xpts,ArrayList<Double>ypts)
	{	int i;
		CrackSegment seg;
		for(i=0;i<segments.size();i++)
		{	seg=segments.get(i);
			xpts.add(seg.xpos[side]);
			ypts.add(seg.ypos[side]);
		}
	}
	
}
