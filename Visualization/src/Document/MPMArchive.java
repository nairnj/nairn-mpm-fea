/*******************************************************************
	MPMArchive.java
	NairnFEAMPMViz

	Created by John Nairn on 11/13/2022.
	Copyright 2022 RSAC Software. All rights reserved.
*******************************************************************/

import java.io.File;

import geditcom.JNFramework.JNUtilities;

public class MPMArchive
{
	private int mstep;
	private double mtime;
	private File archive;
	
	//---------------------------------------------------------------------
	// initialize
	//---------------------------------------------------------------------
	
	MPMArchive(int astep,double atime,File afile)
	{
		mstep = astep;
		mtime = atime;
		archive = afile;
	}
	
	// get name of the archive file
	public String getName()
	{	return archive.getName();
	}
	
	// get file object
	public File getFile()
	{	return archive;
	}
	
	// get the time
	public double getTime()
	{	return mtime;
	}

	// get the step
	public double getStep()
	{	return mstep;
	}
	
	// Scripting attributes for internal scripts
	public String gcis_getAttribute(String [] atoms,int i,CmdViewer server)
	{
	    String attr = server.grabAtom(atoms,i);
	    
	    if(attr.equals("step"))
			return Integer.toString(mstep);
		
		else if(attr.equals("filename"))
			return archive.getName();
	    
		else if(attr.equals("time"))
			return JNUtilities.formatDouble(mtime);
	    
		else if(attr.equals("class"))
			return "MPMArchive";

		return null;
	}


}
