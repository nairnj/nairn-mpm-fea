
/*******************************************************************
	MeshPlotScroll.java
	NairnFEAMPMViz

	Created by John Nairn on 7/20/2010.
	Copyright 2010 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.Dimension;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.Point;

import javax.swing.*;

public class MeshPlotScroll extends JScrollPane implements ComponentListener
{
	private static final long serialVersionUID = 1L;
	private double scale=1.0;
	private MeshPlotView plotView;

	// initialize
	public MeshPlotScroll(MeshPlotView thePlotView)
	{	super(thePlotView);
		plotView=thePlotView;
		addComponentListener(this);
	}
	
	// change the scale
	public void setScale(double newScale,Point click)
	{	// find current center or use click point
		int hcenter,vcenter,hvalue,vvalue;
		int sBarWidth=getHorizontalScrollBar().getHeight();
		Dimension scrollSize=getSize();
		if(click==null)
		{	hvalue=getHorizontalScrollBar().getValue();
			vvalue=getVerticalScrollBar().getValue();
			
			// scroll bar value + half the visible width is the center at 100%
			hcenter=(int)((hvalue+(scrollSize.width-sBarWidth)/2.)/scale);
			vcenter=(int)((vvalue+(scrollSize.height-sBarWidth)/2.)/scale);
		}
		else
		{	// convert to 100%
			hcenter = (int)((double)click.x/scale);
			vcenter = (int)((double)click.y/scale);
		}
		
		// change scale
		scale=newScale;
		componentResized(null);
		
		// recenter plot on (hcenter,vcenter)
		getHorizontalScrollBar().setMaximum(plotView.getWidth());
		getVerticalScrollBar().setMaximum(plotView.getHeight());
		hvalue=(int)(hcenter*scale-(scrollSize.width-sBarWidth)/2.);
		if(hvalue>0) getHorizontalScrollBar().setValue(hvalue);
		vvalue=(int)(vcenter*scale-(scrollSize.height-sBarWidth)/2.);
		if(vvalue>0) getVerticalScrollBar().setValue(vvalue);
	}
	
	//----------------------------------------------------------------------------
	// Component events
	//----------------------------------------------------------------------------

	// on resize, resive the MeshPlotView
	public void componentResized(ComponentEvent e)
	{	Dimension fullSize=getSize();
		int sbarWidth=getHorizontalScrollBar().getHeight();
		if(sbarWidth<1) sbarWidth=15;
		Dimension newSize=new Dimension((int)(fullSize.width*scale-sbarWidth),(int)(fullSize.height*scale-sbarWidth));
		plotView.setPreferredSize(newSize);
		plotView.setSize(newSize);
	}
	public void	componentHidden(ComponentEvent e) {}
	public void componentMoved(ComponentEvent e) {}
	public void componentShown(ComponentEvent e) {}
	
}

