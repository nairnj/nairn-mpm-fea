/*
 * ISNoneType.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 9 NOV 2022.
 * Copyright (c) 2022 RSAC Software. All rights reserved.
 */


public class ISNoneType
{

	// ----------------------------------------------------------------------------
	// Initialize
	// ----------------------------------------------------------------------------

	public ISNoneType()
	{ 
	}
	
	// ----------------------------------------------------------------------------
	// Methods
	// ----------------------------------------------------------------------------
	
	// Scripting attributes for internal scripts for GEDCOMObject
	public String gcis_getAttribute(String [] atoms,int i,CmdViewer server)
	{
	    String attr = server.grabAtom(atoms,i);
	    
		// class
		if(attr.equals("class"))
			return "none";

		return null;
	}
	
}
