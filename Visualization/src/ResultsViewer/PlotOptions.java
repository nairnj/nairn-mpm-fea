/*******************************************************************
	PlotOptions.java
	NairnFEAMPMViz

	Created by John Nairn on 2/9/06.
	Copyright (c) 2006 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

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
	static final int CLIP_TO_PARTICLES=12;
	static final int NUM_OPTIONS=13;
	
	// bits to default flags
	static final int ShowMaterialPts = 0x0001;
	static final int ShowMesh = 0x0002;
	static final int ShowCrackSurfaces = 0x0004;
	static final int ShowBCs = 0x0008;
	static final int ShowDisplacedMesh = 0x0010;
	static final int TransformPts = 0x0020;
	static final int SquareMaterialPts = 0x0040;
	static final int ShowMatPtNumbers = 0x0080;
	static final int ShowCrackPlanes = 0x0100;
	static final int ShowNodeNumbers = 0x0200;
	static final int ShowElemNumbers = 0x0400;
	static final int ShowNodes = 0x0800;
	static final int ClipToParticles = 0x1000;
	
	// default on
	JCheckBox showPts=new JCheckBox("Show Material Pts");
	JCheckBox showMesh=new JCheckBox("Show Mesh");
	JCheckBox showCrackSurf=new JCheckBox("Show Crack Surfaces");
	JCheckBox showMeshBCs=new JCheckBox("Show BCs");
	JCheckBox showDispMesh=new JCheckBox("Show Displaced Mesh");
	JCheckBox transformPts=new JCheckBox("Transform Pts");
	JCheckBox showSquarePts=new JCheckBox("Square Material Pts");
	
	// default off
	JCheckBox showPtNums=new JCheckBox("Show Mat Pt Numbers");
	JCheckBox showCrackPlanes=new JCheckBox("Show Crack Planes");
	JCheckBox showNodeNums=new JCheckBox("Show Node Numbers");
	JCheckBox showElemNums=new JCheckBox("Show Elem Numbers");
	JCheckBox showNodes=new JCheckBox("Show Nodes");
	JCheckBox clipParticles=new JCheckBox("Clip To Particles");
	// particle size
	private JLabel sizeSelected=new JLabel("100",JLabel.LEFT);
	public JSlider mpmParticleSize=new JSlider(JSlider.HORIZONTAL,0,200,5);
	int particleSize=100;
	
	PlotOptions(DocViewer dc)
	{   super(ControlPanel.WIDTH,142,dc);
		setLayout(new GridLayout(7,2));

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
		add(clipParticles);
		
		JPanel sizePanel=new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		sizePanel.setLayout(gridbag);
		c.fill = GridBagConstraints.HORIZONTAL;
		
		sizeSelected.setFont(new Font("sanserif",Font.PLAIN,9));
		c.insets=new Insets(3,10,0,0);
		c.gridwidth = 1;
		c.weightx = 0.0;
		gridbag.setConstraints(sizeSelected, c);
		sizePanel.add(sizeSelected);

		particleSize = NFMVPrefs.prefs.getInt(NFMVPrefs.ParticleSizeKey,NFMVPrefs.ParticleSizeDef);
		mpmParticleSize.setValue(particleSize);
		mpmParticleSize.setToolTipText("Scale particle size as percent of cell size (default is 100%)");
		mpmParticleSize.addChangeListener(new ChangeListener()
		{   public void stateChanged(ChangeEvent e)
			{	particleSize=mpmParticleSize.getValue();
				sizeSelected.setText(""+particleSize);
				NFMVPrefs.prefs.putInt(NFMVPrefs.ParticleSizeKey, particleSize);
			}
		});
		c.insets=new Insets(0,0,0,0);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(mpmParticleSize, c);
		sizePanel.add(mpmParticleSize);
		add(sizePanel);
		
		setEnabled(LoadArchive.NO_PLOT);
		sizeSelected.setMinimumSize(sizeSelected.getPreferredSize());
		sizeSelected.setText(""+particleSize);
		
		int optionsFlags = NFMVPrefs.prefs.getInt(NFMVPrefs.PlotOptionsFlagsKey,
									NFMVPrefs.PlotOptionsFlagsDef);
		// bits to default flags
		showMesh.setSelected((optionsFlags&ShowMesh)!=0);
		showMeshBCs.setSelected((optionsFlags&ShowBCs)!=0);
		showNodeNums.setSelected((optionsFlags&ShowNodeNumbers)!=0);
		showElemNums.setSelected((optionsFlags&ShowElemNumbers)!=0);
		showPts.setSelected((optionsFlags&ShowMaterialPts)!=0);
		showCrackPlanes.setSelected((optionsFlags&ShowCrackPlanes)!=0);
		showCrackSurf.setSelected((optionsFlags&ShowCrackSurfaces)!=0);
		showSquarePts.setSelected((optionsFlags&SquareMaterialPts)!=0);
		showPtNums.setSelected((optionsFlags&ShowMatPtNumbers)!=0);
		showNodes.setSelected((optionsFlags&ShowNodes)!=0);
		showDispMesh.setSelected((optionsFlags&ShowDisplacedMesh)!=0);
		transformPts.setSelected((optionsFlags&TransformPts)!=0);
		clipParticles.setSelected((optionsFlags&ClipToParticles)!=0);
		
		showMesh.setFocusable(false);
		showMeshBCs.setFocusable(false);
		showNodeNums.setFocusable(false);
		showElemNums.setFocusable(false);
		showPts.setFocusable(false);
		showCrackPlanes.setFocusable(false);
		showCrackSurf.setFocusable(false);
		showSquarePts.setFocusable(false);
		showPtNums.setFocusable(false);
		showNodes.setFocusable(false);
		showDispMesh.setFocusable(false);
		transformPts.setFocusable(false);
		clipParticles.setFocusable(false);
	}

	// enable or disable check boxes
	public void setEnabled(int plotType)
	{	if(plotType==LoadArchive.PARTICLE_PLOT)
		{	this.setVisible(true);
			showPts.setEnabled(true);
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
			clipParticles.setEnabled(false);
			mpmParticleSize.setEnabled(true);
		}
		else if(plotType==LoadArchive.MESH_PLOT)
		{	this.setVisible(true);
			showPts.setEnabled(false);
			showCrackPlanes.setEnabled(false);
			showMesh.setEnabled(true);
			showCrackSurf.setEnabled(false);
			showMeshBCs.setEnabled(true);
			showNodes.setEnabled(true);
			showNodeNums.setEnabled(true);
			showElemNums.setEnabled(true);
			showPtNums.setEnabled(false);
			if(docCtrl.resDoc.isFEAAnalysis())
			{	showDispMesh.setEnabled(true);
				showSquarePts.setEnabled(false);
				transformPts.setEnabled(false);
				clipParticles.setEnabled(false);
				mpmParticleSize.setEnabled(false);
			}
			else
			{	showDispMesh.setEnabled(false);
				showSquarePts.setEnabled(true);
				transformPts.setEnabled(true);
				clipParticles.setEnabled(true);
				mpmParticleSize.setEnabled(true);
			}
		}
		else
		{	/*
			showPts.setEnabled(false);
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
			clipParticles.setEnabled(false);
			mpmParticleSize.setEnabled(false);
			*/
			this.setVisible(false);
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
		flags[CLIP_TO_PARTICLES]=clipParticles.isSelected();
		
		// save to preferences
		int optionsFlags = 0;
		if(flags[SHOW_MESH]) optionsFlags += ShowMesh;
		if(flags[SHOW_MESHBCS]) optionsFlags += ShowBCs;
		if(flags[SHOW_NODENUMS]) optionsFlags += ShowNodeNumbers;
		if(flags[SHOW_ELEMNUMS]) optionsFlags += ShowElemNumbers;
		if(flags[SHOW_MATPTS]) optionsFlags += ShowMaterialPts;
		if(flags[SHOW_CRACKPLANES]) optionsFlags += ShowCrackPlanes;
		if(flags[SHOW_CRACKSURFACES]) optionsFlags += ShowCrackSurfaces;
		if(flags[SHOW_SQUAREPTS]) optionsFlags += SquareMaterialPts;
		if(flags[SHOW_MATPTNUMS]) optionsFlags += ShowMatPtNumbers;
		if(flags[SHOW_NODES]) optionsFlags += ShowNodes;
		if(flags[SHOW_DISPLACEDMESH]) optionsFlags += ShowDisplacedMesh;
		if(flags[TRANSFORM_PTS]) optionsFlags += TransformPts;
		if(flags[CLIP_TO_PARTICLES]) optionsFlags += ClipToParticles;
		NFMVPrefs.prefs.putInt(NFMVPrefs.PlotOptionsFlagsKey, optionsFlags);

		return flags;
	}
}