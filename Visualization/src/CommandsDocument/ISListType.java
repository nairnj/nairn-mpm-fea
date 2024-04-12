/*
 * ISListType.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 9 NOV 2022.
 * Copyright (c) 2022 RSAC Software. All rights reserved.
 */

import java.util.ArrayList;

public class ISListType
{
	private ArrayList<Object> list;

	// ----------------------------------------------------------------------------
	// Initialize
	// ----------------------------------------------------------------------------

	public ISListType(ArrayList<Object> seed,boolean copyList)
	{ 
		if(seed==null)
			list = new ArrayList<Object>(200);
		else if(copyList)
			list = new ArrayList<Object>(seed);
		else
			list = seed;
	}
	
	public ISListType(ArrayList<Double> seed)
	{ 
		if(seed==null)
			list = new ArrayList<Object>(200);
		else
			list = new ArrayList<Object>(seed);
	}
	
	// ----------------------------------------------------------------------------
	// Methods
	// ----------------------------------------------------------------------------

	// add commands supported by ISListType as developed
	public static boolean supportsCommand(String testCmd)
	{	String cmd = testCmd.toLowerCase();
		
		// GC List and List commands
		if(cmd.equals("addobject")) return true;
		if(cmd.equals("addstring")) return true;
		if(cmd.equals("insertobject")) return true;
		if(cmd.equals("insertstring")) return true;
		if(cmd.equals("get")) return true;
		if(cmd.equals("sort")) return true;
		if(cmd.equals("join")) return true;
		if(cmd.equals("remove")) return true;
		if(cmd.equals("addlist")) return true;
		if(cmd.equals("sublist")) return true;
		if(cmd.equals("pop")) return true;
		
		return false;
	}

	// add object to the list
	public void gcis_addObject(Object obj)
	{	list.add(obj);
	}

	// sort the list
	/*
	- (void)gcis_sort:(NSString *)style
	{
		[list sortUsingSelector:@selector(caseInsensitiveCompare:)];
	}
	*/

	// join the list components
	/*
	- (id)gcis_join:(NSString *)joinString
	{
		return [list componentsJoinedByString:joinString];
	}
	*/

	// insert object, but return NO if out of bounds
	public boolean gcis_insertObject(Object newObj,int index)
	{
		if(index<0 || index>list.size()) return false;
		list.add(index,newObj);
		return true;
	}

	// insert object, but return NO is out of bounds
	public boolean gcis_removeObjectAtIndex(int index)
	{
		if(index<0 || index>list.size()) return false;
		list.remove(index);
		return true;
	}

	// Scripting attributes for internal scripts
	public String gcis_getAttribute(String [] atoms,int i,CmdViewer server)
	{
	    String attr = server.grabAtom(atoms,i);
	    
	    if(attr.equals("count"))
			return Integer.toString(list.size());
		
		// class
		else if(attr.equals("class"))
			return "list";

		return null;
	}

	// ----------------------------------------------------------------------------
	// Accessors
	// ----------------------------------------------------------------------------

	// length of the list
	public int count() { return list.size(); }

	// fetch and object
	public Object objectAtIndex(int index) { return list.get(index); }

	// remove object withnotindex checking
	public void removeObjectAtIndex(int index) { list.remove(index); }

	// pop last object
	public Object getRemoveLastObject()
	{	int index = list.size()-1;
		if(index<0) return null;
		Object last = list.get(index);
		list.remove(index);
		return last;
	}

	// ----------------------------------------------------------------------------
	// Accessors
	// ----------------------------------------------------------------------------

	public ArrayList<Object>getList() { return list; }
	/*
	- (void)setList:(NSMutableArray *)newList
	{	[newList retain];
		[list release];
		list = newList;
	}
	*/

}
