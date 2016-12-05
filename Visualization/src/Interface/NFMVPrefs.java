/*
 * NFMVPrefs.java
 * NairnFEAMPMViz Application
 * 
 * Created by John Nairn on 2/17/08
 * Copyright (c) 2008-2010 RSAC Software. All rights reserved
 */

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import javax.swing.*;

import geditcom.JNFramework.*;
import geditcom.plot2d.*;

public class NFMVPrefs extends JNPreferences implements ActionListener
{	
	private static final long serialVersionUID = 1L;

	// Results window split locations
	public static String ResultsTopSplitKey = "Results Top Split";
	public static double ResultsTopSplitDef = 0.35;
	public static String ResultsSplitKey = "Results Split";
	public static double ResultsSplitDef = 0.50;

	// Commands window split location
	public static String CommandsSplitKey = "Commands Split";
	public static double CommandsSplitDef = 0.60;

	// path to NairnMPM binary
	public static String NairnMPMKey = "NairnMPM Path";
	public static String NairnMPMDef = "$(bundle)";

	// path to NairnMPM dtd
	public static String NairnMPMDTDKey = "NairnMPM DTD File Path";
	public static String NairnMPMDTDDef = "$(bundle)";
	public static String NairnMPMValidateKey = "Validate MPM Calculations";
	public static boolean NairnMPMValidateDef = true;

	// path to NairnFEA binary
	public static String NairnFEAKey = "NairnFEA Path";
	public static String NairnFEADef = "$(bundle)";

	// path to NairnFEA DTD
	public static String NairnFEADTDKey = "NairnFEA DTD File Path";
	public static String NairnFEADTDDef = "$(bundle)";
	public static String NairnFEAValidateKey = "Validate MPM Calculations";
	public static boolean NairnFEAValidateDef = true;

	// path to shell command
	public static String ShellKey = "Shell Path";
	public static String ShellDef = "$(windows)";

	// path to workspace folder
	public static String WorkSpaceKey = "Workspace Folder";
	public static String WorkSpaceDef = "";

	// REMOTE_ACCESS
	// remote server
	public static String RemoteServerKey = "Remote Server";
	public static String RemoteServerDef = "";
	public static String RemoteUserKey = "Remote Username";
	public static String RemoteUserDef = "";
	public static String RemoteUserPassKey = "Remote Password";
	public static String RemoteUserPassDef = "";
	// code exec location (local or remote)
	public static String CodeExecLocationKey = "Code Exec Location";
	public static String CodeExecLocationDef = "local";
	public static String RemoteMPMPathKey = "Remote Path to MPM Code";
	public static String RemoteMPMPathDef = "/usr/local/bin/NairnMPM";
	public static String RemoteMPMDTDKey = "Remote Path to MPM DTD";
	public static String RemoteMPMDTDDef = "/usr/local/bin/NairnMPM.dtd";
	public static String RemoteFEAPathKey = "Remote Path to FEA Code";
	public static String RemoteFEAPathDef = "/usr/local/bin/NairnFEA";
	public static String RemoteFEADTDKey = "Remote Path to FEA DTD";
	public static String RemoteFEADTDDef = "/usr/local/bin/NairnFEA.dtd";

	// plot spectrum
	public static String SpectrumKey = "PlotSpectrum";
	public static int SpectrumDef = 1;
	public static String NumContoursKey = "NumberOfPlotContours";
	public static int NumContoursDef = 0;
	public static String NumSubelementsKey = "NumberOfSubelements";
	public static int NumSubelementsDef = 4;

	// plot colors
	public static String backColorKey = "Background Color";
	public static Color backColorDef = new Color(84, 89, 109); // 0.329f,0.349f,0.427f
	public static String meshLineColorKey = "Mesh Line Color";
	public static Color meshLineColorDef = new Color(191, 191, 191); // 0.75f,0.75f,0.75f
	public static String meshNodesColorKey = "Nodal Point Color";
	public static Color meshNodesColorDef = new Color(212, 212, 212); // 0.83f,0.83f,0.83f
	public static String textColorKey = "Text Label Color";
	public static Color textColorDef = new Color(255, 255, 255); // 1.f,1.f,1.f

	private JTextField mpmCodePath = new JTextField();
	private JTextField mpmDTDPath = new JTextField();
	private JTextField feaCodePath = new JTextField();
	private JTextField feaDTDPath = new JTextField();
	private JTextField shellPath = new JTextField();
	private JTextField workPath = new JTextField();
	private JTextField numContours = new JTextField();
	private JTextField numSubelements = new JTextField();
	private JFileChooser chooser = new JFileChooser();

	// prefs
	private JComboBox<String> spectrumBox;

	private JCheckBox validateMPM = new JCheckBox("Validate",NairnMPMValidateDef);
	private JCheckBox validateFEA = new JCheckBox("Validate",NairnFEAValidateDef);

	// REMOTE_ACCESS
	private JTextField remoteServerAddr = new JTextField();
	private JTextField remoteUsername = new JTextField();
	private JPasswordField remoteUserpass = new JPasswordField();
	private JTextField mpmCodePathFld = new JTextField();
	private JTextField mpmDTDPathFld = new JTextField();
	private JTextField feaCodePathFld = new JTextField();
	private JTextField feaDTDPathFld = new JTextField();
	private JRadioButton rdbtnExecLocal = new JRadioButton("Local");
	private JRadioButton rdbtnExecRemote = new JRadioButton("Remote");
	
	private static boolean currentRemoteMode = false;

	// ----------------------------------------------------------------------------
	// constructor
	// ----------------------------------------------------------------------------

	public NFMVPrefs()
	{	super("Preferences");

		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.setFocusable(false);

		// Remove panel
		tabbedPane.addTab("Code", null, codePanel(),"Settings for running calculations");

		// Remote preferences pane
		tabbedPane.addTab("Remote", null, remotePanel(), "Remote server settings");

		// color preferences pane ----------------------------------------
		tabbedPane.addTab("Colors", null, colorPanel(), "Select colors used in plots");

		// FEA and MPM preferences pane
		// -----------------------------------------
		//JPanel panel4 = new JPanel();
		//tabbedPane.addTab("FEA and MPM", null, panel4,"FEA and MPM plot settings");

		getContentPane().add(tabbedPane);
		setSize(600, 480);
		setLocation(60, 30 + JNApplication.menuBarHeight());
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
	}
	
	// create code prferences panel
	public JPanel codePanel()
	{	JPanel panel1 = new JPanel(new GridLayout(6,1));

		panel1.add(filePathPanel("MPM Code", mpmCodePath, prefs.get(
				NairnMPMKey, NairnMPMDef), false, null));
		panel1.add(filePathPanel("MPM DTD File", mpmDTDPath, prefs.get(
				NairnMPMDTDKey, NairnMPMDTDDef), prefs.getBoolean(
				NairnMPMValidateKey, NairnMPMValidateDef), validateMPM));
		panel1.add(filePathPanel("FEA Code", feaCodePath, prefs.get(
				NairnFEAKey, NairnFEADef), false, null));
		panel1.add(filePathPanel("FEA DTD File", feaDTDPath, prefs.get(
				NairnFEADTDKey, NairnFEADTDDef), prefs.getBoolean(
				NairnFEAValidateKey, NairnFEAValidateDef), validateFEA));
		panel1.add(filePathPanel("Shell Command", shellPath, prefs.get(
				ShellKey, ShellDef), false, null));
		panel1.add(filePathPanel("Work Space Directory", workPath, prefs.get(
				WorkSpaceKey, WorkSpaceDef), false, null));

		return panel1;
	}
	
	// REMOTE_ACCESS
	// build panel for remote connections
	public JPanel remotePanel()
	{	JPanel panel3 = new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		panel3.setLayout(gridbag);
		
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets=new Insets(3,6,3,6);
		c.weighty = 0.;

		// Remove server line -------------------------------------------
		String toolTip = "Enter server address (e.g. 'mpm.fsl.orst.edu')";
		remoteLine(panel3,gridbag,c,"Remote Server:",remoteServerAddr,
				prefs.get(RemoteServerKey, RemoteServerDef),RemoteServerKey,toolTip);
		
		// User name and password line -------------------------
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		JLabel label = new JLabel("Username:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);

		c.gridx=1;
		c.weightx = 1.5;
		c.gridwidth = 1;
		remoteUsername.setText(prefs.get(RemoteUserKey, RemoteUserDef));
		remoteUsername.setActionCommand("Username:");
		remoteUsername.addActionListener(this);
		remoteUsername.addFocusListener(new PrefFocusListener(RemoteUserKey));
		remoteUsername.setToolTipText("Enter user name for server '"+remoteServerAddr.getText()+"'");
		gridbag.setConstraints(remoteUsername,c);
		panel3.add(remoteUsername);
		
		c.gridx=2;
		c.weightx = 0.;
		c.gridwidth = 1;
		label = new JLabel("Password:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);
		
		c.gridx=3;
		c.weightx = 1.5;
		c.gridwidth = 1;
		remoteUserpass.setText(prefs.get(RemoteUserPassKey, RemoteUserPassDef));
		remoteUserpass.setActionCommand("Password:");
		remoteUserpass.addActionListener(this);
		remoteUserpass.addFocusListener(new PrefFocusListener(RemoteUserPassKey));
		remoteUserpass.setToolTipText("Enter password for user '"+remoteUsername.getText()+"' on server '"+remoteServerAddr.getText()+"'");
		gridbag.setConstraints(remoteUserpass,c);
		panel3.add(remoteUserpass);
		
		// MPM Code -------------------------------------------
		toolTip = "Enter full path to MPM executeable (can use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"MPM Code Path:",mpmCodePathFld,
				prefs.get(RemoteMPMPathKey, RemoteMPMPathDef),RemoteMPMPathKey,toolTip);
		toolTip = "Enter full path to MPM dtd file (cannot use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"MPM DTD Path:",mpmDTDPathFld,
				prefs.get(RemoteMPMDTDKey, RemoteMPMDTDDef),RemoteMPMDTDKey,toolTip);
		
		// FEA Code -------------------------------------------
		toolTip = "Enter full path to FEA executeable (can use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"FEA Code Path:",feaCodePathFld,
				prefs.get(RemoteFEAPathKey, RemoteFEAPathDef),RemoteFEAPathKey,toolTip);
		toolTip = "Enter full path to FEA dtd file (cannot use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"FEA DTD Path:",feaDTDPathFld,
				prefs.get(RemoteFEADTDKey, RemoteFEADTDDef),RemoteFEADTDKey,toolTip);
		
		// Remove/local buttons
		ButtonGroup execLocation = new ButtonGroup();
		rdbtnExecRemote.addActionListener(this);
		rdbtnExecRemote.setFocusable(false);
		rdbtnExecLocal.addActionListener(this);
		rdbtnExecLocal.setFocusable(false);
		rdbtnExecLocal.setActionCommand("code exec location");
		execLocation.add(rdbtnExecLocal);
		GridBagConstraints gbc_rdbtnExecLocal = new GridBagConstraints();
		gbc_rdbtnExecLocal.fill = GridBagConstraints.BOTH;
		gbc_rdbtnExecLocal.insets = new Insets(18, 20, 3, 5);
		gbc_rdbtnExecLocal.gridx = 0;
		gbc_rdbtnExecLocal.gridy = 8;
		panel3.add(rdbtnExecLocal, gbc_rdbtnExecLocal);
		rdbtnExecRemote.setActionCommand("code exec location");
		execLocation.add(rdbtnExecRemote);
		GridBagConstraints gbc_rdbtnExecRemote = new GridBagConstraints();
		gbc_rdbtnExecRemote.insets = new Insets(0, 20, 0, 5);
		gbc_rdbtnExecRemote.fill = GridBagConstraints.BOTH;
		gbc_rdbtnExecRemote.gridx = 0;
		gbc_rdbtnExecRemote.gridy = 9;
		panel3.add(rdbtnExecRemote, gbc_rdbtnExecRemote);
		
		if(restoreRemoteMode())
			rdbtnExecRemote.setSelected(true);
		else
			rdbtnExecLocal.setSelected(true);
		
		// empty fill space on the bottom
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 4;
		c.weighty = 10;
		c.fill = GridBagConstraints.VERTICAL;
		label = new JLabel(" ");
		gridbag.setConstraints(label,c);
		panel3.add(label);

		return panel3;
	}
	
	// build panel for remote connections
	public JPanel colorPanel()
	{	JPanel panel2 = new JPanel();

		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints pc = new GridBagConstraints(); // for the pane
		panel2.setLayout(gridbag);

		// spectrum box
		GridBagLayout spectrumBag = new GridBagLayout();
		JNBoxPanel spectrumPan = new JNBoxPanel("Color Plot Spectrum",spectrumBag);
		GridBagConstraints c = new GridBagConstraints(); // for the box panel

		JLabel section = new JLabel("Type: ");
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.WEST;
		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 0.0;
		c.insets = new Insets(6, 6, 1, 0); // tlbr
		spectrumBag.setConstraints(section, c);
		spectrumPan.add(section);

		String[] lines = { "Rainbow 1 (Purple to Red)",
				"Rainbow 2 (Blue to Red)", "Rainbow 3 (Purple to Orange)",
				"Grayscale (Black to White)", "Cool Diverging (Blue to Red)" };
		spectrumBox = new JComboBox<String>(lines);
		spectrumBox.setSelectedIndex(ColorPicker.getSpectrumType());
		spectrumBox.addActionListener(this);
		spectrumBox.setFocusable(false);
		spectrumBox.setToolTipText("Select one of the available spectra for coloring the plots.");
		c.insets = new Insets(6, 3, 1, 6); // tlbr
		c.gridx = 1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(spectrumBox, c);
		spectrumPan.add(spectrumBox);

		JLabel section1 = new JLabel("Contours: ");
		c.insets = new Insets(6, 6, 1, 0); // tlbr
		c.gridx = 0;
		c.weightx = 0.0;
		spectrumBag.setConstraints(section1, c);
		spectrumPan.add(section1);

		numContours.setText("" + ColorPicker.getNumberOfContours());
		numContours.setEditable(true);
		numContours.setColumns(5);
		numContours.setToolTipText("Enter 1 for continuous colors or >1 for number of discrete colors and hit return or enter");
		numContours.addActionListener(this);
		numContours.setActionCommand("change contours");
		numContours.addFocusListener(new PrefFocusListener(NumContoursKey));
		c.insets = new Insets(6, 3, 1, 6); // tlbr
		c.gridx = 1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(numContours, c);
		spectrumPan.add(numContours);

		JLabel section2 = new JLabel("Subelements: ");
		c.insets = new Insets(6, 6, 1, 0); // tlbr
		c.gridx = 0;
		c.weightx = 0.0;
		spectrumBag.setConstraints(section2, c);
		spectrumPan.add(section2);

		numSubelements.setText("" + ElementBase.getSubelementDensity());
		numSubelements.setEditable(true);
		numSubelements.setColumns(5);
		numSubelements.setToolTipText("Enter number of subelements (nXn) to use for mesh plots.");
		numSubelements.addActionListener(this);
		numSubelements.setActionCommand("change subelement density");
		numSubelements.addFocusListener(new PrefFocusListener(NumSubelementsKey));
		c.insets = new Insets(6, 3, 1, 6); // tlbr
		c.gridx = 1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(numSubelements, c);
		spectrumPan.add(numSubelements);

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 1.0;

		// add spectrum box to pan
		pc.fill = GridBagConstraints.HORIZONTAL;
		pc.gridx = 0;
		pc.gridwidth = 1;
		pc.weightx = 1.0;
		pc.weighty = 0.;
		pc.insets = new Insets(1, 6, 1, 6); // tlbr
		pc.anchor = GridBagConstraints.CENTER;
		gridbag.setConstraints(spectrumPan, pc);
		panel2.add(spectrumPan);

		// plot element color wells box
		GridBagLayout wellsBag = new GridBagLayout();
		JNBoxPanel wellsPan = new JNBoxPanel("Plot Elements", wellsBag);
		GridBagConstraints ec = new GridBagConstraints(); // for the plot

		JLabel elab1 = new JLabel("Background: ");
		ec.fill = GridBagConstraints.NONE;
		ec.anchor = GridBagConstraints.EAST;
		ec.insets = new Insets(10, 0, 0, 0); // tlbr
		ec.gridx = 0;
		ec.weightx = 0.5;
		ec.gridwidth = 1;
		wellsBag.setConstraints(elab1, ec);
		wellsPan.add(elab1);

		JNColorWell bgWell = new JNColorWell(getPrefColor(backColorKey,backColorDef));
		bgWell.addActionListener(this, "bgColor", "Choose Background Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(bgWell, ec);
		wellsPan.add(bgWell);

		JLabel elab2 = new JLabel("Mesh Lines: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab2, ec);
		wellsPan.add(elab2);

		JNColorWell mlWell = new JNColorWell(getPrefColor(meshLineColorKey,meshLineColorDef));
		mlWell.addActionListener(this, "mlColor", "Choose Mesh Line Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(mlWell, ec);
		wellsPan.add(mlWell);

		JLabel elab3 = new JLabel("Labels: ");
		ec.insets = new Insets(10, 0, 6, 0); // tlbr
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 0;
		ec.weightx = 0.5;
		ec.gridwidth = 1;
		wellsBag.setConstraints(elab3, ec);
		wellsPan.add(elab3);

		JNColorWell labWell = new JNColorWell(getPrefColor(textColorKey,textColorDef));
		labWell.addActionListener(this, "labColor","Choose Color for Text Labels");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(labWell, ec);
		wellsPan.add(labWell);

		JLabel elab4 = new JLabel("Nodal Points: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab4, ec);
		wellsPan.add(elab4);

		JNColorWell npWell = new JNColorWell(getPrefColor(meshNodesColorKey,meshNodesColorDef));
		npWell.addActionListener(this, "npColor", "Choose Nodal Point Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(npWell, ec);
		wellsPan.add(npWell);

		// add plot element color wells box to pane
		gridbag.setConstraints(wellsPan, pc);
		panel2.add(wellsPan);
		
		// empty fill on the bottom
		pc.fill = GridBagConstraints.VERTICAL;
		pc.weighty = 10.;
		JLabel label = new JLabel(" ");
		gridbag.setConstraints(label,pc);
		panel2.add(label);
		
		return panel2;
	}
	
	// line on remote panel with label and text field
	public void remoteLine(JPanel panel3,GridBagLayout gridbag,GridBagConstraints c,String tlab,
			JTextField theFld,String theValue,String aPrefKey,String toolTip)
	{
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		JLabel label = new JLabel(tlab);
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);
		
		c.gridx=1;
		c.weightx = 3.0;
		c.gridwidth = 3;
		theFld.setText(theValue);
		theFld.setActionCommand(tlab);
		theFld.addActionListener(this);
		theFld.addFocusListener(new PrefFocusListener(aPrefKey));
		theFld.setToolTipText(toolTip);
		gridbag.setConstraints(theFld,c);
		panel3.add(theFld);
	}

	// Create panel for file path entry
	public JPanel filePathPanel(String pathName, JTextField pathField,
			String pathText, boolean validPref, JCheckBox validate)
	{	JPanel filePanel = new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		filePanel.setLayout(gridbag);

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 1.0;
		c.insets = new Insets(0, 6, 0, 0); // tlbr
		JLabel nameLabel = new JLabel(pathName + " Path:");
		gridbag.setConstraints(nameLabel, c);
		filePanel.add(nameLabel);

		int rowWidth = 2;
		if (validate != null)
		{	c.gridx++;
			validate.setSelected(validPref);
			validate.setActionCommand(pathName + " Use");
			validate.addActionListener(this);
			validate.setToolTipText("Check to validate input files using the provided DTD file");
			gridbag.setConstraints(validate, c);
			filePanel.add(validate);
			rowWidth++;
		}

		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.EAST;
		c.gridx++;
		c.insets = new Insets(0, 0, 0, 6); // tlbr
		JButton change = new JButton("Change...");
		change.setActionCommand(pathName);
		change.addActionListener(this);
		change.setToolTipText("Click to change the " + pathName + " path");
		change.setFocusable(false);
		gridbag.setConstraints(change, c);
		filePanel.add(change);

		c.anchor = GridBagConstraints.CENTER;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridwidth = rowWidth;
		c.insets = new Insets(0, 6, 0, 9); // tlbr
		pathField.setText(pathText);
		pathField.setEditable(false);
		pathField.setColumns(45);
		pathField.setToolTipText("The current " + pathName
				+ " path; click 'Change...' to change it.");
		gridbag.setConstraints(pathField, c);
		filePanel.add(pathField);

		return filePanel;
	}

	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();

		if(theCmd.equals("MPM Code"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			mpmCodePath.setText(newPath);
			prefs.put(NairnMPMKey,newPath);
		}
		
		else if(theCmd.equals("MPM DTD File Use"))
		{	prefs.putBoolean(NairnMPMValidateKey,validateMPM.isSelected());
		}
		
		else if(theCmd.equals("MPM DTD File"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			mpmDTDPath.setText(newPath);
			prefs.put(NairnMPMDTDKey,newPath);
		}
		
		else if(theCmd.equals("FEA Code"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			feaCodePath.setText(newPath);
			prefs.put(NairnFEAKey,newPath);
		}
		
		else if(theCmd.equals("FEA DTD File Use"))
		{	prefs.putBoolean(NairnFEAValidateKey,validateFEA.isSelected());
		}
		
		else if(theCmd.equals("bgColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(backColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("labColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(textColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("mlColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(meshLineColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("npColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(meshNodesColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("FEA DTD File"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			feaDTDPath.setText(newPath);
			prefs.put(NairnFEADTDKey,newPath);
		}
		
		else if(theCmd.equals("Shell Command"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(windows)")>=0)
				newPath = "$(windows)";
			shellPath.setText(newPath);
			prefs.put(ShellKey,newPath);
		}
		
		else if(theCmd.equals("Work Space Directory"))
		{	chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int result=chooser.showOpenDialog(this);
			String newPath=chooser.getSelectedFile().getPath();
			chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			if(result==JFileChooser.CANCEL_OPTION) return;
			workPath.setText(newPath);
			prefs.put(WorkSpaceKey,newPath);
			setWorkspace(NairnFEAMPMViz.GetAppChooser());
		}
		
		else if(theCmd.equals("change contours"))
		{	try
			{	int newContours=Integer.valueOf(numContours.getText());
				if(newContours<1) newContours=1;
				int oldNum=ColorPicker.getNumberOfContours();
				if(newContours!=oldNum)
				{	prefs.putInt(NumContoursKey,newContours);
					ColorPicker.setNumberOfContours();
				}
			}
			catch(Exception ie)
			{	Toolkit.getDefaultToolkit().beep();
			}
		}

		else if(theCmd.equals("change subelement density"))
		{	try
			{	int newDensity=Integer.valueOf(numSubelements.getText());
				if(newDensity<1) newDensity=1;
				int oldNum=ElementBase.getSubelementDensity();
				if(newDensity!=oldNum)
				{	prefs.putInt(NumSubelementsKey,newDensity);
					ElementBase.setSubelementDensity();
				}
			}
			catch(Exception ie)
			{	Toolkit.getDefaultToolkit().beep();
			}
		}

	// REMOTE_ACCESS
		else if (theCmd.equalsIgnoreCase("Username:"))
		{	prefs.put(RemoteUserKey, remoteUsername.getText());
		}

		else if (theCmd.equalsIgnoreCase("Remote Server:"))
		{	prefs.put(RemoteServerKey, this.remoteServerAddr.getText());
		}

		else if (theCmd.equalsIgnoreCase("Password:"))
		{	prefs.put(RemoteUserPassKey, new String(remoteUserpass.getPassword()));
		}

		else if (theCmd.equalsIgnoreCase("MPM Code Path:"))
		{	prefs.put(RemoteMPMPathKey, mpmCodePathFld.getText());
		}

		else if (theCmd.equalsIgnoreCase("MPM DTD Path:"))
		{	prefs.put(RemoteMPMDTDKey, mpmDTDPathFld.getText());
		}

		else if (theCmd.equalsIgnoreCase("FEA Code Path:"))
		{	prefs.put(RemoteFEAPathKey, feaCodePathFld.getText());
		}

		else if (theCmd.equalsIgnoreCase("FEA DTD Path:"))
		{	prefs.put(RemoteFEADTDKey, feaDTDPathFld.getText());
		}

		else if (theCmd.equalsIgnoreCase("code exec location"))
		{	if (this.rdbtnExecLocal.isSelected())
			{	prefs.put(CodeExecLocationKey, "local");
				currentRemoteMode = false;
			}
			else if (this.rdbtnExecRemote.isSelected())
			{	prefs.put(CodeExecLocationKey, "remote");
				currentRemoteMode = true;
			}
		}
		
		else
		{	JComboBox<?> cb = (JComboBox<?>) e.getSource();
			if(cb == spectrumBox)
			{	int oldType = ColorPicker.getSpectrumType();
				int newType = cb.getSelectedIndex();
				if(newType != oldType)
				{	prefs.putInt(SpectrumKey, newType);
					ColorPicker.setSpectrumType();
				}
			}
			else
				System.out.println("Unrecognized preferences commane: "+theCmd);
		}
	}

	// ----------------------------------------------------------------------------
	// class methods
	// ----------------------------------------------------------------------------

	// create preferences storage object
	public static void initializePrefs()
	{	createPrefs("com/geditcom/NairnFEAMPMViz");
		restoreRemoteMode();
	
		// default close action
		NFMVPrefs.setDefaultQuitCloseLastWindow(false);

		// restore open docs
		NFMVPrefs.setDefaultRestoreOpenDocs(true);
		
		// default settings
		ColorPicker.setNumberOfContours();
		ElementBase.setSubelementDensity();
	}

	// set a chooser to current work space
	public static void setWorkspace(JFileChooser chooser)
	{	String workPath=prefs.get(WorkSpaceKey,WorkSpaceDef);
		if(workPath.length()>0)
		{	try
			{	File workDirectory=new File(workPath);
				chooser.setCurrentDirectory(workDirectory);
			}
			catch(Exception e) {}
		}
	}
	
	// return true it set to run remotely
	// REMOTE_ACCESS
	public static boolean getRemoteMode() { return currentRemoteMode; }
	public static void setRemoteMode(boolean tempRemote) { currentRemoteMode = tempRemote; }
	public static boolean restoreRemoteMode()
	{	if(prefs!=null)
			currentRemoteMode = prefs.get(CodeExecLocationKey, CodeExecLocationDef).equalsIgnoreCase("remote");
		else
			currentRemoteMode = false;
		return currentRemoteMode;
	}
	
	// listenere for fields tied to string preferences
	public class PrefFocusListener implements FocusListener
	{
		private String prefKey;
		
		public PrefFocusListener(String aKey)
		{	prefKey = aKey;
		}
		
		public void focusLost(FocusEvent e)
		{	String newValue;
			if(prefKey.equals(NumContoursKey))
			{	try
				{	int newContours=Integer.valueOf(numContours.getText());
					if(newContours<1) newContours=1;
					int oldNum=ColorPicker.getNumberOfContours();
					if(newContours!=oldNum)
					{	prefs.putInt(NumContoursKey,newContours);
						ColorPicker.setNumberOfContours();
					}
				}
				catch(Exception ie)
				{	Toolkit.getDefaultToolkit().beep();
				}
			}
			else if(prefKey.equals(NumSubelementsKey))
			{	try
				{	int newDensity=Integer.valueOf(numSubelements.getText());
					if(newDensity<1) newDensity=1;
					int oldNum=ElementBase.getSubelementDensity();
					if(newDensity!=oldNum)
					{	prefs.putInt(NumSubelementsKey,newDensity);
						ElementBase.setSubelementDensity();
					}
				}
				catch(Exception ie)
				{	Toolkit.getDefaultToolkit().beep();
				}
			}
			else
			{	if(prefKey.equals(RemoteUserPassKey))
					newValue = new String(((JPasswordField)e.getComponent()).getPassword());
				else
					newValue = ((JTextField)e.getComponent()).getText();
				prefs.put(prefKey, newValue);
			}
		}
		
		public void focusGained(FocusEvent e)
		{
		}
	}

}
