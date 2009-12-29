/*******************************************************************
	CommandDocument.java
	NairnFEAMPMViz

	Created by John Nairn on 14 Feb 2008.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.io.*;

public class CommandDocument
{
	static final long serialVersionUID=32L;
	
	//---------------------------------------------------------
	// variables and constants
	//---------------------------------------------------------
	
	public File inputFile=null;
	public CmdViewer docCtrl;
	
	//---------------------------------------------------------
	// Initialize
	//---------------------------------------------------------
	
	public CommandDocument()
	{
	}
	
	//-----------------------------------------------------------------
	// Accessor Methods
	//-----------------------------------------------------------------
	
	// set file of the mpm file
	public void setFile(File xmlFile) { inputFile=xmlFile; }
	public File getFile() { return inputFile; }
	
	// file name or untitled
	public String getName()
	{	if(inputFile==null)
			return "untitled";
		else
			return inputFile.getName();
	}
	public String getFullPath()
	{	if(inputFile==null)
			return "";
		else
			return inputFile.getPath();
	}
	
	// set controller
	public void setDocController(CmdViewer dc) { docCtrl=dc; }
	
}
