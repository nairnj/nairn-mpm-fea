
/*******************************************************************
	MeshPlotScroll.java
	NairnFEAMPMViz

	Created by John Nairn on 7/20/2010.
	Copyright 2010 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.Dimension;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

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
	public void setScale(double newScale)
	{	// find current center
		int hvalue=getHorizontalScrollBar().getValue();
		int sBarWidth=getHorizontalScrollBar().getHeight();
		int vvalue=getVerticalScrollBar().getValue();
		Dimension scrollSize=getSize();
		int hcenter=(int)((hvalue+(scrollSize.width-sBarWidth)/2.)/scale);
		int vcenter=(int)((vvalue+(scrollSize.height-sBarWidth)/2.)/scale);
		
		// change scale
		scale=newScale;
		componentResized(null);
		
		// recenter plot
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

