/*
	TimeSelector.java
	NairnFEAMPMViz
	
	This same class is used for time selector in both the results
	control panel and in the movie window. Its location is determined
	by whether movieCtrl==null (in results control panel) or
	movieCtrl!=null (in movie plot window). docCtrl is always the
	document controller.

	Created by John Nairn on Tue Mar 09 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*/

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;
import geditcom.JNFramework.*;

public class TimeSelector extends PlotControl
{
	static final long serialVersionUID=23L;
	
	// variables
	public JSlider select=new JSlider(0,50,0);
	private JLabel timeSelected=new JLabel();
	private int value=-1;
	private MoviePlotWindow movieCtrl;
	private DocViewer docCtrl;
	private String prefix;
	
	// initialize
	TimeSelector(MoviePlotWindow mpc,DocViewer dc)
	{   super(ControlPanel.WIDTH,64,dc);
		setLayout(new GridLayout(2,1));

		// slider bar
		docCtrl=dc;
		movieCtrl=mpc;
		if(movieCtrl==null)
		{	select.setMajorTickSpacing(10);
			select.setPaintTicks(true);
			select.setToolTipText("Select the initial time results to be plotted");
		}
		else
			select.setToolTipText("Select different time results to be plotted");
		add(select);
		
		// the displayed value in a text field
		timeSelected.setHorizontalAlignment(JLabel.CENTER);
		add(timeSelected);
		
		// initialize
		if(movieCtrl!=null)
		{	setEnabled(LoadArchive.PARTICLE_PLOT,0);
			setBorder(BorderFactory.createEmptyBorder());		// of the PlotControl JPanel
			setBackground(Color.lightGray);
			select.setBackground(Color.lightGray);
			timeSelected.setFont(new Font("sanserif",Font.PLAIN,10));
			timeSelected.setText(dc.resDoc.mpmArchives.get(0).getTime() + " "+dc.resDoc.units.timeUnits());
			prefix="";
		}
		else
		{	setEnabled(LoadArchive.NO_PLOT,0);
			timeSelected.setText("Time: 0.0 "+dc.resDoc.units.timeUnits());
			prefix="Time: ";
		}
		select.setFocusable(false);
			
		// the slider listener
		select.addChangeListener(new ChangeListener()
		{   public void stateChanged(ChangeEvent e)
			{	ResultsDocument resDoc=TimeSelector.this.docCtrl.resDoc;
				String theTime=JNUtilities.formatDouble(resDoc.mpmArchives.get(select.getValue()).getTime());
				timeSelected.setText(prefix + theTime + " "+resDoc.units.timeUnits());
				if(!select.getValueIsAdjusting())
				{	if(value!=select.getValue())
					{	value=select.getValue();
						docCtrl.getMovieFrame().changeArchiveIndex(value);
					}
				}
			}
		});
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// update label on rescaling
	public void updateLabel()
	{	ResultsDocument resDoc=TimeSelector.this.docCtrl.resDoc;
		String theTime=JNUtilities.formatDouble(resDoc.mpmArchives.get(select.getValue()).getTime());
		timeSelected.setText(prefix + theTime + " "+resDoc.units.timeUnits());
	}
	
	// called when new file loaded
	public void setEnabled(int selected,int plotComponent)
	{	// only used for MPM plots, except for time plots
		boolean opened = (selected!=LoadArchive.TIME_PLOT) && 
					(selected!=LoadArchive.NO_PLOT) && (docCtrl.resDoc.isMPMAnalysis());
	    if(opened)
		{	if(plotComponent!=PlotQuantity.IMPORTANDPLOTFILE)
			{	select.setValueIsAdjusting(true);
		    	setVisible(true);
		    	select.setEnabled(true);
				select.setMaximum(docCtrl.resDoc.mpmArchives.size()-1);
				int ticks=(docCtrl.resDoc.mpmArchives.size()-1)/10;
				if(ticks<1) ticks=1;
				select.setMajorTickSpacing(ticks);
				//select.setValue(0);
				timeSelected.setEnabled(true);
			}
			else
				setVisible(false);
		}
		else
		{   setVisible(false);
			//select.setEnabled(false);
			//timeSelected.setEnabled(false);
		}
	}
	
	// return selected archive by index
	public int getArchiveIndex() { return select.getValue(); }
	
	// change value and return true or false if this change caused a reload
	// if it did not reload, will need to reload data if component has changed too
	public boolean setArchiveIndex(int newIndex)
	{	if(newIndex!=select.getValue())
		{	select.setValue(newIndex);
			return true;
		}
		else if(value<0)
		{	// only here on first call, so post notification to get it read if in movie controller
			value=select.getValue();
			docCtrl.getMovieFrame().changeArchiveIndex(value);
			return true;
		}
		return false;
	}
	
	// increment archive index if possible
	public boolean incrementArchiveIndex()
	{	int currentIndex=select.getValue();
		if(currentIndex==select.getMaximum()) return false;
		select.setValue(currentIndex+1);
		return true;
	}
	
	// decrement archive index if possible
	public boolean decrementArchiveIndex()
	{	int currentIndex=select.getValue();
		if(currentIndex==select.getMinimum()) return false;
		select.setValue(currentIndex-1);
		return true;
	}

}
