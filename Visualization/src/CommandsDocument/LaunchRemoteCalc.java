/*
 * CmdViewer
 * NairnFEAMPMViz Application
 * 
 * Created 
 */

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import geditcom.JNFramework.*;

public class LaunchRemoteCalc extends JNDialog
{
	private JPanel values=new JPanel();
	private JTextField remoteOutputName = new JTextField();
	private JTextField localSaveFolder = new JTextField();
	private JRadioButton overWrite = new JRadioButton("Overwrite");
	private JRadioButton createUnique = new JRadioButton("Create Unique");
	private JRadioButton clearParent = new JRadioButton("Clear Parent");
	private JComboBox<String> doDownload;
	
	private String remoteFolder;
	private String remoteFileName;
	private JFileChooser chooser = new JFileChooser();
	
	private static final long serialVersionUID = 2839653163815510656L;
	
	public static final int DOWNLOAD_TO_FOLDER = 0;
	public static final int NO_DOWNLOAD = 1;
	public static final int OPEN_ON_SERVER = 2;

	public LaunchRemoteCalc(CmdViewer cmdCtrl,String initName,boolean initUnique,String initOutput,
								int initDown,boolean initClear)
	{	super(cmdCtrl,"Launch Remote Calculations","Select output file on server and local folder for optional download",
			" OK ","Cancel",null);
	
		GridBagLayout gridbag = new GridBagLayout();
		values.setLayout(gridbag);
	
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets=new Insets(0,0,0,0);

		// Remove save file  -------------------------------------------
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		JLabel label = new JLabel("Remote Output:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		values.add(label);
			
		c.gridx=1;
		c.weightx = 10.;
		c.gridwidth = 2;
		remoteOutputName.setText(initName);
		remoteOutputName.setToolTipText("Enter path to output file (using /'s) and must have at least one parent folder");
		gridbag.setConstraints(remoteOutputName,c);
		values.add(remoteOutputName);
		
		// unique check box line -------------------------------------
		c.insets=new Insets(2,0,6,0);
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		label = new JLabel("");
		gridbag.setConstraints(label,c);
		values.add(label);
		
		JPanel panel1 = new JPanel(new GridLayout(1,3));
		ButtonGroup writeOption = new ButtonGroup();
		
		overWrite.setFocusable(false);
		overWrite.setToolTipText("Choose to overwrite files with same name in the parent folder");
		writeOption.add(overWrite);
		panel1.add(overWrite);
		
		createUnique.setFocusable(false);
		createUnique.setToolTipText("Choose to force unique folder on the server within the parent folder");
		writeOption.add(createUnique);
		panel1.add(createUnique);
		
		clearParent.setFocusable(false);
		clearParent.setToolTipText("Choose to delete prior contents of the parent folder before the new calculations");
		writeOption.add(clearParent);
		panel1.add(clearParent);
		
		if(initUnique)
			createUnique.setSelected(true);
		else if(initClear)
			clearParent.setSelected(true);
		else
			overWrite.setSelected(true);
			
		c.gridx=1;
		c.weightx = 10.;
		c.gridwidth = 2;
		gridbag.setConstraints(panel1,c);
		values.add(panel1);
		
		// Download save destination  -------------------------------------------
		c.insets=new Insets(0,0,0,0);
		c.gridx=0;
		c.weightx=0.;
		c.gridwidth = 1;
		label = new JLabel("Local Folder:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		values.add(label);
		
		c.gridx=1;
		c.weightx = 10.;
		c.gridwidth = 2;
		localSaveFolder.setText(initOutput);
		localSaveFolder.setEditable(false);
		localSaveFolder.setToolTipText("Click 'Change...' to choose folder for downloaded output");
		gridbag.setConstraints(localSaveFolder,c);
		values.add(localSaveFolder);
		
		// change button -------------------------------------
		c.insets=new Insets(4,0,0,0);
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		label = new JLabel("");
		gridbag.setConstraints(label,c);
		values.add(label);
		
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.WEST;
		c.gridx=1;
		c.weightx=5.;
		JButton changeBtn = new JButton("Change...");
		changeBtn.setFocusable(false);
		changeBtn.setActionCommand("local Folder");
		changeBtn.addActionListener(this);
		gridbag.setConstraints(changeBtn,c);
		values.add(changeBtn);
		
		c.gridx=2;
		String [] lines={"Download to this Folder","Do not Download","Home Directory on Server"};
		doDownload=new JComboBox<String>(lines);
		gridbag.setConstraints(doDownload,c);
		doDownload.setSelectedIndex(initDown);
		doDownload.setFocusable(false);
		values.add(doDownload);
		
		// add to dialog
		add(values,BorderLayout.CENTER);
		
		// finish up
		setSize(600,210,"  ","   ");
	}
	
	// handle commands
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();

		if(theCmd.equals("local Folder"))
		{	chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			String newPath = localSaveFolder.getText();
			try
			{	File startDirectory=new File(newPath);
				chooser.setCurrentDirectory(startDirectory);
			}
			catch(Exception ee) {}
			int result=chooser.showSaveDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			newPath=chooser.getSelectedFile().getPath();
			localSaveFolder.setText(newPath);
		}
		
		else
			super.actionPerformed(e);
	}
	
	// verify remote path is valid
	public boolean dialogDone(int result)
	{
		String path = remoteOutputName.getText();
		int offset = path.lastIndexOf("/");
		
		if(offset<1 || offset>=path.length()-1)
		{	JOptionPane.showMessageDialog(this,"The remote folder must have at least a folder and a name (e.g., 'folder/name')");
			return false;
		}
		
		remoteFolder = path.substring(0, offset);
		remoteFileName = path.substring(offset+1);
		
		if(getClearContents() && !getMakeUnique())
		{	String message = "Are you sure you want to clear entire entire contents\n of folder '~/"
					+remoteFolder+"' before the new calculations";
			int decide = JOptionPane.showConfirmDialog(this, message, "Delete Folder?", JOptionPane.OK_CANCEL_OPTION);
			if(decide==JOptionPane.CANCEL_OPTION) return false;
		}
		
		if(getDoDownload()!=NO_DOWNLOAD)
		{	File fldr = new File(localSaveFolder.getText());
			if(!fldr.exists() || !fldr.isDirectory())
			{	JOptionPane.showMessageDialog(this,"The local folder must exist and must be a folder.");
				return false;
			}
		}

		return true;
	}
	
	// accessors
	public String getRemoteFolder() { return remoteFolder; }
	public String getRemoteFileName() { return remoteFileName; }
	public boolean getMakeUnique() { return createUnique.isSelected(); }
	public boolean getClearContents() { return clearParent.isSelected(); }
	public String getLocalFolder() { return localSaveFolder.getText(); }
	public int getDoDownload() { return doDownload.getSelectedIndex(); }
}
