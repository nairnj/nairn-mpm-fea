/*******************************************************************
	LoadArchive.java
	NairnFEAMPMViz

	Created by John Nairn on Tue Mar 09 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.event.*;
import javax.swing.*;
import java.net.*;
import java.awt.*;

public class LoadArchive extends PlotControl implements ActionListener
{
	static final long serialVersionUID=7L;
	
	// variables and constants
	private static final int ICON_WIDTH=32;
	
	private ButtonGroup plotGroup=new ButtonGroup();
	private JRadioButton [] plotType;
	private JButton [] rbIcon;
	public static final int NO_PLOT=-1;
	public static final int PARTICLE_PLOT=0;
	public static final int TIME_PLOT=1;
	public static final int MESH_PLOT=2;
	public static final int MESH2D_PLOT=3;
	private int selected;
	
	
	// initilize
	LoadArchive(DocViewer dc)
	{   super(ControlPanel.WIDTH,20,dc);
		setLayout(null);
	
		// plot type buttons
		plotType=new JRadioButton[4];
		rbIcon=new JButton[4];
		int top=6;
		createButtonIcon(PARTICLE_PLOT,"ParticlePlots",top);
		plotType[PARTICLE_PLOT].setSelected(true);
		createButtonIcon(TIME_PLOT,"TimePlots",top);
		createButtonIcon(MESH_PLOT,"MeshPlots",top);
		createButtonIcon(MESH2D_PLOT,"Mesh2DPlots",top);
		
		setEnabled(false);
		setSize(ControlPanel.WIDTH,plotType[MESH2D_PLOT].getY()+plotType[MESH2D_PLOT].getHeight()+12);
	}

	// called when new file is loaded
	public void setEnabled()
	{	if(docCtrl.resDoc.isMPMAnalysis())
		{	plotType[PARTICLE_PLOT].setSelected(true);
			selected=PARTICLE_PLOT;
			plotType[PARTICLE_PLOT].setEnabled(true);
			rbIcon[PARTICLE_PLOT].setEnabled(true);
			plotType[TIME_PLOT].setEnabled(true);
			rbIcon[TIME_PLOT].setEnabled(true);
		}
		else
		{	plotType[MESH_PLOT].setSelected(true);
			selected=MESH_PLOT;
			plotType[PARTICLE_PLOT].setEnabled(false);
			rbIcon[PARTICLE_PLOT].setEnabled(false);
			plotType[TIME_PLOT].setEnabled(false);
			rbIcon[TIME_PLOT].setEnabled(false);
		}
		plotType[MESH_PLOT].setEnabled(true);
		rbIcon[MESH_PLOT].setEnabled(true);
		plotType[MESH2D_PLOT].setEnabled(true);
		rbIcon[MESH2D_PLOT].setEnabled(true);
	}
	
	// help to set an image of a radio button
	private void createButtonIcon(int rbID,String icon,int yloc)
	{
		plotType[rbID]=new JRadioButton();
		int xloc= rbID>0 ? plotType[rbID-1].getX()+plotType[rbID-1].getWidth()+ICON_WIDTH+5 : 5 ;
		plotType[rbID].addActionListener(this);
		plotType[rbID].setActionCommand(icon);
		plotType[rbID].setSize(plotType[rbID].getPreferredSize());
		plotType[rbID].setLocation(xloc,yloc+6);
		plotGroup.add(plotType[rbID]);
		add(plotType[rbID]);
		
		URL btnImage=Main.class.getResource(icon+".png");
		rbIcon[rbID]=new JButton(new ImageIcon(btnImage));
		rbIcon[rbID].addActionListener(this);
		rbIcon[rbID].setActionCommand(icon);
		rbIcon[rbID].setSize(ICON_WIDTH,ICON_WIDTH);
		rbIcon[rbID].setLocation(plotType[rbID].getX()+plotType[rbID].getWidth(),yloc);
		rbIcon[rbID].setFocusPainted(false);
		rbIcon[rbID].setBorderPainted(false);
		rbIcon[rbID].setContentAreaFilled(false);
		add(rbIcon[rbID]);
	}

	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();
		int oldSelected=selected;
	
		if(theCmd.equals("ParticlePlots"))
		{	selected=PARTICLE_PLOT;
		}
		
		else if(theCmd.equals("TimePlots"))
		{	selected=TIME_PLOT;
		}
		
		else if(theCmd.equals("MeshPlots"))
		{	selected=MESH_PLOT;
		}
		
		else if(theCmd.equals("Mesh2DPlots"))
		{	selected=MESH2D_PLOT;
		}

		if(selected!=oldSelected)
		{	plotType[selected].setSelected(true);
			docCtrl.controls.hiliteControls();
		}
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------

	public int getSelected() { return selected; }

	// for connecting lines
	public Rectangle getControlRect()
	{	Point loc=rbIcon[selected].getLocation();
		Point box=getLocation();
		return new Rectangle(box.x+loc.x,box.y,rbIcon[selected].getWidth(),getHeight());
	}
	
}
