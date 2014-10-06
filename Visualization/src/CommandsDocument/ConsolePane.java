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
	
	//-----------------------------------------------------------------
	// Initialize - scrolling edit field
	//-----------------------------------------------------------------
	
	public ConsolePane()
	{	super();
		
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
			{	chooser.setSelectedFile(outputFile);
			}
		
			int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return false;
		
			outputFile=chooser.getSelectedFile();
		}
		else
			outputFile = useOutput;
		
		// add extension if needed
		if(!outputFile.getName().endsWith(extension))
			outputFile=new File(outputFile.getParent(),outputFile.getName()+"."+extension);
			
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