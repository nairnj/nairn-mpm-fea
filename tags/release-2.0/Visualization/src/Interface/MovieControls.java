/*******************************************************************
	MovieControls.java
	NairnFEAMPMViz

	Created by John Nairn on 3/1/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;
import java.net.*;

public class MovieControls extends JPanel
{
	static final long serialVersionUID=11L;
	
	//----------------------------------------------------------------------------
	// variables and constants
	//----------------------------------------------------------------------------
	
	public static final int HEIGHT=44;

	private JButton rewindMov=new JButton();
	private JButton playMov=new JButton();
	private TimeSelector selectTime=null;
	
	//----------------------------------------------------------------------------
	// initialize
	//----------------------------------------------------------------------------
	
	MovieControls(int width,ResultsDocument gResDoc,MoviePlotWindow movieCtrl)
	{   super();
		setLayout(null);
		setPreferredSize(new Dimension(width,HEIGHT));
		setBackground(Color.lightGray);
		
		int hpos=8;
		if(gResDoc.isMPMAnalysis())
		{	// rewind
			setButtonIcon("Rewind.png",rewindMov);
			rewindMov.setSize(rewindMov.getPreferredSize());
			int centerLoc=(HEIGHT-rewindMov.getHeight())>>1;
			rewindMov.setLocation(hpos,centerLoc);
			rewindMov.setActionCommand("Rewind");
			rewindMov.addActionListener(movieCtrl);
			rewindMov.setFocusPainted(false);
			add(rewindMov);
			
			// play button
			setButtonIcon("Play.png",playMov);
			playMov.setSize(playMov.getPreferredSize());
			centerLoc=(HEIGHT-playMov.getHeight())>>1;
			hpos+=rewindMov.getWidth()+6;
			playMov.setLocation(hpos,centerLoc);
			playMov.setActionCommand("Play");
			playMov.addActionListener(movieCtrl);
			playMov.setFocusPainted(false);
			add(playMov);
			
			LookAndFeel laf=UIManager.getLookAndFeel();
			if(laf!=null)
			{	if(laf.isNativeLookAndFeel())
				{	playMov.setBackground(Color.lightGray);
					rewindMov.setBackground(Color.lightGray);
				}
			}

			// time selector
			selectTime=new TimeSelector(movieCtrl,gResDoc.docCtrl);
			selectTime.setSize(new Dimension(140,HEIGHT));
			hpos+=playMov.getWidth()+12;
			selectTime.setLocation(hpos,0);
			add(selectTime);
			
			hpos+=selectTime.getWidth()+20;
		}
		
		// data panael
		MeshPlotData meshData=new MeshPlotData(movieCtrl.getPlotView(),gResDoc);
		meshData.setSize(new Dimension(200,HEIGHT));
		meshData.setLocation(hpos,0);
		meshData.setBackground(Color.lightGray);
		add(meshData);
		
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// help to set an image of a button
	private void setButtonIcon(String icon,JButton theBtn)
	{	URL btnImage=Main.class.getResource(icon);
		theBtn.setIcon(new ImageIcon(btnImage));
	}
	
	// reset selected  index in mesh view when about to replot and true or false if changed
	public boolean setArchiveIndex(int newIndex) { return selectTime!=null ? selectTime.setArchiveIndex(newIndex) : false; }
	
	public int getArchiveIndex() { return selectTime!=null ? selectTime.getArchiveIndex() : 0; }
	public boolean incrementArchiveIndex() { return selectTime!=null ? selectTime.incrementArchiveIndex() : false; }
	
	// called when movie started or stopped
	public void setPlaying(boolean playing)
	{	if(playing)
			setButtonIcon("Pause.png",playMov);
		else
			setButtonIcon("Play.png",playMov);
	}
}
