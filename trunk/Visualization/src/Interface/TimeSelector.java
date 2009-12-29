/*******************************************************************
	TimeSelector.java
	NairnFEAMPMViz

	Created by John Nairn on Tue Mar 09 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;
import java.util.*;
import javax.swing.event.*;

public class TimeSelector extends PlotControl
{
	static final long serialVersionUID=23L;
	
	// variables
	private JSlider select=new JSlider(0,50,0);
	private JLabel timeSelected=new JLabel();
	private ArrayList<Double> archiveTimes;
	private int value=-1;
	private MoviePlotWindow movieCtrl;
	
	// initialize
	TimeSelector(MoviePlotWindow mpc,DocViewer dc)
	{   super(ControlPanel.WIDTH,64,dc);
		setLayout(new GridLayout(2,1));

		// slider bar
		movieCtrl=mpc;
		if(movieCtrl==null)
		{	select.setMajorTickSpacing(10);
			select.setPaintTicks(true);
		}
		add(select);
		
		// the displayed value
		timeSelected.setHorizontalAlignment(JLabel.CENTER);
		add(timeSelected);
		
		// initialize
		if(movieCtrl!=null)
		{	setEnabled(LoadArchive.PARTICLE_PLOT);
			setBorder(BorderFactory.createEmptyBorder());		// of the PlotControl JPanel
			setBackground(Color.lightGray);
			select.setBackground(Color.lightGray);
			timeSelected.setFont(new Font("sanserif",Font.PLAIN,10));
			timeSelected.setText(archiveTimes.get(0) + " ms");
		}
		else
		{	setEnabled(LoadArchive.NO_PLOT);
			timeSelected.setText("Time: 0.0 ms");
		}
		select.setFocusable(false);
			
		// the slider listener
		select.addChangeListener(new ChangeListener()
		{   public void stateChanged(ChangeEvent e)
			{	if(movieCtrl!=null)
				{	timeSelected.setText(archiveTimes.get(select.getValue()) + " ms");
					if(!select.getValueIsAdjusting() && value!=select.getValue())
					{	value=select.getValue();
						movieCtrl.changeArchiveIndex(value);
					}
				}
				else
				{	timeSelected.setText("Time: " + archiveTimes.get(select.getValue()) + " ms");
					if(!select.getValueIsAdjusting()) value=select.getValue();
				}
			}
		});
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// called when new file loaded
	public void setEnabled(int selected)
	{	boolean opened = (selected!=LoadArchive.TIME_PLOT) && (selected!=LoadArchive.NO_PLOT) && (docCtrl.resDoc.isMPMAnalysis());
	    if(opened)
		{   archiveTimes=docCtrl.resDoc.archiveTimes;
			select.setEnabled(true);
			select.setMaximum(archiveTimes.size()-1);
			//select.setValue(0);
			timeSelected.setEnabled(true);
		}
		else
		{   select.setEnabled(false);
			timeSelected.setEnabled(false);
		}
	}
	
	// return selected archive by index
	public int getArchiveIndex() { return select.getValue(); }
	
	// change value and return true or false if this change caused a reload
	// if it did not reload, will need to reload data if component has changed too
	public boolean setArchiveIndex(int newIndex)
	{	if(newIndex!=select.getValue())
		{	select.setValue(newIndex);
			// state change will reload the data
			return true;
		}
		else if(value<0)
		{	// only here on first call
			value=select.getValue();
			if(movieCtrl!=null)
				movieCtrl.changeArchiveIndex(value);
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
}
