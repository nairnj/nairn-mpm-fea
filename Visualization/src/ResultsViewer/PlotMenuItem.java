/*******************************************************************
	PlotMenuItem.java
	NairnFEAMPMViz

	Created by John Nairn on 2/8/06.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

public class PlotMenuItem
{
	private String item;
	private int tag;
	
	PlotMenuItem(String newItem,int newTag)
	{	item=new String(newItem);
		tag=newTag;
	}

	PlotMenuItem(String newItem)
	{	item=new String(newItem);
		tag=0;
	}
	
	public String toString( ) { return item; }
	
	public int getTag( ) { return tag; }
	
	public void setString(String newItem) { item = newItem; }
}
