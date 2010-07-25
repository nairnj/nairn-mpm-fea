/*******************************************************************
	NFMVFilter.java
	NairnFEAMPMViz

	Created by John Nairn on 2/17/08.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.util.*;
import java.io.File;
import javax.swing.filechooser.*;

public class NFMVFilter extends FileFilter
{
	private ArrayList<String> extensions;
	private String description="NairnFEAMPMViz Files";
	
	// initialize
	public NFMVFilter()
	{
		extensions=new ArrayList<String>(10);
	}
	
	// add an extension to the list
	public void addExtension(String extension)
	{
		extensions.add("."+extension);
	}

	// true if directory or ends in an added extennsion
	public boolean accept(File f)
	{	if(f.isDirectory()) return true;
	
		String fileName=f.getName();
		for(int i=0;i<extensions.size();i++)
		{	if(fileName.endsWith(extensions.get(i)))
				return true;
		}
		return false;
	}
	
    //The description of this filter
    public void setDescription(String describe) { description=describe; }
    public String getDescription() { return description; }
}