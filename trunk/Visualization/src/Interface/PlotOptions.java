/*******************************************************************
	PlotOptions.java
	NairnFEAMPMViz

	Created by John Nairn on 2/9/06.
	Copyright (c) 2006 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;

public class PlotOptions extends PlotControl
{
	static final long serialVersionUID=17L;
	
	// constants and variables
	static final int SHOW_MESH=0;
	static final int SHOW_MESHBCS=1;
	static final int SHOW_NODENUMS=2;
	static final int SHOW_ELEMNUMS=3;
	static final int SHOW_MATPTS=4;
	static final int SHOW_MATPTNUMS=5;
	static final int SHOW_CRACKPLANES=6;
	static final int SHOW_CRACKSURFACES=7;
	static final int SHOW_SQUAREPTS=8;
	static final int SHOW_NODES=9;
	static final int SHOW_DISPLACEDMESH=10;
	static final int TRANSFORM_PTS=11;
	static final int NUM_OPTIONS=12;
	
	// default on
	JCheckBox showPts=new JCheckBox("Show Material Pts",true);
	JCheckBox showMesh=new JCheckBox("Show Mesh",true);
	JCheckBox showCrackSurf=new JCheckBox("Show Crack Surfaces",true);
	JCheckBox showMeshBCs=new JCheckBox("Show BCs",true);
	JCheckBox showDispMesh=new JCheckBox("Show Displaced Mesh",true);
	JCheckBox transformPts=new JCheckBox("Transform Pts",true);
	
	// default off
	JCheckBox showPtNums=new JCheckBox("Show Mat Pt Numbers",false);
	JCheckBox showSquarePts=new JCheckBox("Square Material Pts",false);
	JCheckBox showCrackPlanes=new JCheckBox("Show Crack Planes",false);
	JCheckBox showNodeNums=new JCheckBox("Show Node Numbers",false);
	JCheckBox showElemNums=new JCheckBox("Show Elem Numbers",false);
	JCheckBox showNodes=new JCheckBox("Show Nodes",false);
	
	// initialize
	PlotOptions(DocViewer dc)
	{   super(ControlPanel.WIDTH,120,dc);
		setLayout(new GridLayout(6,2));

		// add check boxes
		add(showPts);
		add(showSquarePts);
		add(showCrackPlanes);
		add(showCrackSurf);
		add(showMesh);
		add(showMeshBCs);
		add(showNodes);
		add(showNodeNums);
		add(showElemNums);
		add(showPtNums);
		add(showDispMesh);
		add(transformPts);
		
		setEnabled(LoadArchive.NO_PLOT);
	}

	// enable or disable check boxes
	public void setEnabled(int plotType)
	{	if(plotType==LoadArchive.PARTICLE_PLOT)
		{	showPts.setEnabled(true);
			showSquarePts.setEnabled(true);
			showCrackPlanes.setEnabled(true);
			showMesh.setEnabled(true);
			showCrackSurf.setEnabled(true);
			showMeshBCs.setEnabled(true);
			showNodes.setEnabled(true);
			showNodeNums.setEnabled(true);
			showElemNums.setEnabled(true);
			showPtNums.setEnabled(true);
			showDispMesh.setEnabled(false);
			transformPts.setEnabled(true);
		}
		else if(plotType==LoadArchive.MESH_PLOT)
		{	showPts.setEnabled(false);
			showSquarePts.setEnabled(false);
			showCrackPlanes.setEnabled(false);
			showMesh.setEnabled(true);
			showCrackSurf.setEnabled(false);
			showMeshBCs.setEnabled(true);
			showNodes.setEnabled(true);
			showNodeNums.setEnabled(true);
			showElemNums.setEnabled(true);
			showPtNums.setEnabled(false);
			showDispMesh.setEnabled(docCtrl.resDoc.isFEAAnalysis());
			transformPts.setEnabled(false);
		}
		else
		{	showPts.setEnabled(false);
			showSquarePts.setEnabled(false);
			showCrackPlanes.setEnabled(false);
			showMesh.setEnabled(false);
			showCrackSurf.setEnabled(false);
			showMeshBCs.setEnabled(false);
			showNodes.setEnabled(false);
			showNodeNums.setEnabled(false);
			showElemNums.setEnabled(false);
			showPtNums.setEnabled(false);
			showDispMesh.setEnabled(false);
			transformPts.setEnabled(false);
		}
	}

	// get list of options for the plot
	public boolean [] getOptions()
	{
		boolean [] flags=new boolean[NUM_OPTIONS];
		flags[SHOW_MESH]=showMesh.isSelected();
		flags[SHOW_MESHBCS]=showMeshBCs.isSelected();
		flags[SHOW_NODENUMS]=showNodeNums.isSelected();
		flags[SHOW_ELEMNUMS]=showElemNums.isSelected();
		flags[SHOW_MATPTS]=showPts.isSelected();
		flags[SHOW_CRACKPLANES]=showCrackPlanes.isSelected();
		flags[SHOW_CRACKSURFACES]=showCrackSurf.isSelected();
		flags[SHOW_SQUAREPTS]=showSquarePts.isSelected();
		flags[SHOW_MATPTNUMS]=showPtNums.isSelected();
		flags[SHOW_NODES]=showNodes.isSelected();
		flags[SHOW_DISPLACEDMESH]=showDispMesh.isSelected();
		flags[TRANSFORM_PTS]=transformPts.isSelected();
		return flags;
	}
}