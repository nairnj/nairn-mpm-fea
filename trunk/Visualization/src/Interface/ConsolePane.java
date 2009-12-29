/*******************************************************************
	ConsolePane.java
	NairnFEAMPMViz

	Created by John Nairn on 14 Feb 2008.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.io.*;
import javax.swing.*;

public class ConsolePane extends JScrollPane
{
	static final long serialVersionUID=33L;
	
	//-----------------------------------------------------------------
	// Class variables and constants
	//-----------------------------------------------------------------
	
	// members
	public JTextArea soutPane = new JTextArea( );
	public File outputFile=null;
  
	// command file chooser
	private JFileChooser chooser=new JFileChooser();
	
	//-----------------------------------------------------------------
	// Initialize - scrolling edit field
	//-----------------------------------------------------------------
	
	public ConsolePane()
	{	super();
		setViewportView(soutPane);
		
		// text editing panel
		soutPane.setEditable(false);
		soutPane.setFont(new Font("Monospaced",Font.PLAIN,12));
		
		// always show scroll bar
		setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		
		// set file filter
		NFMVFilter filter=new NFMVFilter();
		filter.addExtension("mpm");
		filter.addExtension("fea");
		filter.setDescription("MPM for FEA Output Files");
		chooser.setFileFilter(filter);
		
		NFMVPrefs.setWorkspace(chooser);
	}
	
	//-----------------------------------------------------------------
	// Methods
	//-----------------------------------------------------------------
	
	// get file name prior to running an analysis
	public boolean setOutputPath(File inputFile,String extension)
	{
		// set to directory
		if(outputFile==null)
			chooser.setCurrentDirectory(inputFile);
		else
		{	chooser.setSelectedFile(outputFile);
		}
		
		int result = chooser.showSaveDialog(this);
		if(result != JFileChooser.APPROVE_OPTION) return false;
		
		outputFile=chooser.getSelectedFile();
		
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
	
	public void clear() { soutPane.setText(""); }
	
	// append line and scroll to the  end
	public void appendLine(String moreText)
	{	soutPane.append(moreText+"\n");
		JScrollBar vertBar=getVerticalScrollBar();
		vertBar.setValue(vertBar.getMaximum());
	}
	
	// save commands and return saved file (or null)
	public boolean saveOutput(JFrame fileWindow)
	{
		// save to saveFile
		try
		{	FileWriter theFile=new FileWriter(outputFile);
			theFile.write(soutPane.getText());
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
	{	String idText=soutPane.getText();
		int start=Math.max(0,idText.lastIndexOf(' '));
		String theID=idText.substring(start);
		if(theID.length()==0)
			return "";
		else
			return " with process ID: "+theID;
	}
}