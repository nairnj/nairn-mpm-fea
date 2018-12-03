/*
 * ConsolePane.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 14 Feb 2008.
 * Copyright (c) 2008 RSAC Software. All rights reserved.
 */

import java.awt.*;
import java.io.*;
import javax.swing.*;

import geditcom.JNFramework.*;

public class ConsolePane extends JNConsolePane
{
	static final long serialVersionUID=33L;
	
	//-----------------------------------------------------------------
	// Class variables and constants
	//-----------------------------------------------------------------
	
	// file to save the output results and a chooser
	public File outputFile=null;
	private JFileChooser chooser=new JFileChooser();
	public String remoteFilePath = null;
	public String remoteHomePath = null;
	public boolean uniqueOutput = false;
	public int downloadResults = LaunchRemoteCalc.DOWNLOAD_TO_FOLDER;
	public boolean clearPriorContents = false;
	
	//-----------------------------------------------------------------
	// Initialize - scrolling edit field
	//-----------------------------------------------------------------
	
	public ConsolePane(Font theFont)
	{	super();
	
		if(theFont!=null)
			console.setFont(theFont);
		
		// set file filter
		JNFileFilter filter=new JNFileFilter();
		filter.addExtension("mpm");
		filter.addExtension("fea");
		filter.setDescription("MPM or FEA Output Files");
		chooser.setFileFilter(filter);
		
		NFMVPrefs.setWorkspace(chooser);
	}
	
	//-----------------------------------------------------------------
	// Methods
	//-----------------------------------------------------------------
	
	// get file name prior to running an analysis
	// useOutput is provided when initiated by script, otherwise the
	//    user selects output file to save the results now.
	public boolean setOutputPath(File inputFile,String extension,File useOutput)
	{
		// use or choose
		if(useOutput == null)
		{	// set to directory
			if(outputFile==null)
				chooser.setCurrentDirectory(inputFile);
			else
				chooser.setSelectedFile(outputFile);
			
			// REMOTE_ACCESS - customize to get path on server or a local file
			//System.out.println("...debug get file from chooser "+chooser);
			try
			{	int result = chooser.showSaveDialog(this);
				if(result != JFileChooser.APPROVE_OPTION) return false;
				outputFile=chooser.getSelectedFile();
			}
			catch(HeadlessException e)
			{	String msg = "The Java file chooser failed to launch.\n";
				msg = msg+e.getLocalizedMessage();
				JNApplication.appBeep();
				JOptionPane.showMessageDialog(null,msg);
				return false;
			}
		}
		else
			outputFile = useOutput;
		
		// add extension if needed
		if(!outputFile.getName().endsWith(extension))
			outputFile=new File(outputFile.getParent(),outputFile.getName()+"."+extension);
			
		return true;
	}
	
	// when run remotely, just set output path when done
	public boolean setOutputPath(String pathName)
	{	outputFile=new File(pathName);
		return true;
	}
	
	// get the file
	public File getFile() { return outputFile; }
	
	public String getOutputPath()
	{	if(outputFile==null) return null;
		return outputFile.getPath();
	}
	
	public String getOutputName()
	{	if(outputFile==null) return null;
		return outputFile.getName();
	}
	
	// save commands and return saved file (or null)
	public boolean saveOutput(JFrame fileWindow)
	{
		// save to saveFile
		try
		{	FileWriter theFile=new FileWriter(outputFile);
			theFile.write(console.getText());
			theFile.flush();
			theFile.close();
		}
		catch (Exception fe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(fileWindow,"Error writing analysis results: " + fe);
			return false;
		}

		return true;
	}
	
	// store remote file path
	public String getRemoteFilePath() { return remoteFilePath; }
	public void setRemoteFilePath(String path) { remoteFilePath = new String(path); }
	
	// store remote home path (when ussed
	public String getRemoteHomePath() { return remoteHomePath; }
	public void setRemoteHomePath(String path) { remoteHomePath = new String(path); }
	
	// when run remotely, store path to remote file
	// folder/folder/name
	public void getRemoteFilePath(String pathName)
	{	remoteFilePath=new String(pathName);
	}
	
	// last word of text if was a submitted job
	public String processID()
	{	String idText=console.getText();
		int start=Math.max(0,idText.lastIndexOf(' '));
		String theID=idText.substring(start);
		if(theID.length()==0)
			return "";
		else
			return " with process ID: "+theID;
	}
}