/*
 * ISDictType.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 9 NOV 2022.
 * Copyright (c) 2022 RSAC Software. All rights reserved.
 */

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import geditcom.JNFramework.JNUtilities;

public class ISDictType
{
	private HashMap<String, Object> dict = null;

	// ----------------------------------------------------------------------------
	// Initialize
	// ----------------------------------------------------------------------------

	public ISDictType()
	{ 
		dict = new HashMap<String, Object>(100);
	}
	
	// ----------------------------------------------------------------------------
	// Methods
	// ----------------------------------------------------------------------------
	
	// Scripting attributes for internal scripts for dictionary
	public Object gcis_getObjectAttribute(String attr,CmdViewer server)
	{
		if(attr.contentEquals("allKeys"))
		{	ISListType keys = new ISListType(null);
			Set<String> allKeys = dict.keySet();
			Iterator<String> keyIterator = allKeys.iterator();
			while(keyIterator.hasNext())
				keys.gcis_addObject(keyIterator.next());
			return keys;
		}
		
		else if(attr.contentEquals("allValues"))
		{	ISListType values = new ISListType(null);
			Collection<Object> allValues = dict.values();
			Iterator<Object> valIterator = allValues.iterator();
			while(valIterator.hasNext())
				values.gcis_addObject(valIterator.next());
			return values;
		}
		
		/*
		else if([attr isEqualToString:@"sortedKeys"])
		{	ISListType *list = [ISListType listWithList:[dict allKeys]];
			[list gcis_sort:nil];
			return list;
		}
		
		else if([attr isEqualToString:@"keysSortedByNumericValues"])
		{	NSArray *sortedKeys = [dict keysSortedByValueUsingComparator:^(id obj1, id obj2)
			{	// sort non numbers to the end
				if(![obj1 isKindOfClass:[NSNumber class]])
				{	// first not a number, if second is a number it comes first
					if([obj2 isKindOfClass:[NSNumber class]])
						return (NSComparisonResult)NSOrderedAscending;
					return (NSComparisonResult)NSOrderedSame;
				}
				else if(![obj2 isKindOfClass:[NSNumber class]])
				{	// first number, so it is higher
					return (NSComparisonResult)NSOrderedDescending;
				}
				// both numbers
				return [obj1 compare:obj2];
			}];
			return [ISListType listWithList:sortedKeys];
		}
		
		else if([attr isEqualToString:@"keysSortedByStringValues"])
		{	NSArray *sortedKeys = [dict keysSortedByValueUsingComparator:^(id obj1, id obj2)
			{	// sort non numbers to the end
				NSString *obj1s,*obj2s;
				if(![obj1 isKindOfClass:[NSString class]])
				{	// try to convert to a string
	NS_DURING
					obj1s = [obj1 stringValue];
	NS_HANDLER
					obj1s = nil;
	NS_ENDHANDLER
				}
				if(![obj2 isKindOfClass:[NSString class]])
				{	// try to convert to a string
	NS_DURING
					obj2s = [obj2 stringValue];
	NS_HANDLER
					obj2s = nil;
	NS_ENDHANDLER
				}
				if(obj1s==nil)
				{	// first not a string, if second is a string it comes first
					if(obj2s!=nil)
						return (NSComparisonResult)NSOrderedAscending;
					return (NSComparisonResult)NSOrderedSame;
				}
				else if(obj2s==nil)
				{	// first string, so it is higher
					return (NSComparisonResult)NSOrderedDescending;
				}
				// both strings
				return [obj1s localizedCompare:obj2s];
			}];
			return [ISListType listWithList:sortedKeys];
		}
		*/
		
		return null;
	}

	// Scripting attributes for internal scripts for GEDCOMObject
	public String gcis_getAttribute(String [] atoms,int i,CmdViewer server)
	{
	    String attr = server.grabAtom(atoms,i);
	    
	    if(attr.equals("count"))
			return Integer.toString(dict.size());
		
		// class
		else if(attr.equals("class"))
			return "dictionary";

		return null;
	}

	/*
	// Getting scripting properties for internal scripts for GEDCOMObject
	- (BOOL)gcis_setAttribute:(NSString *)attr value:(NSString *)value sender:(RunServer *)server
	{
		return NO;
	}
	*/
	
	// set object for a key
	public void gcis_setObjectforKey(Object newObj,String key)
	{
		dict.put(key,newObj);
	}

	// get object
	public Object gcis_objectForKey(String key)
	{
		return dict.get(key);
	}

	// remove object
	public void gcis_removeValueForKey(String key)
	{
		dict.remove(key);
	}

	// get integer for key
	// throws exception in an error
	public int gcis_integerForKey(String key) throws Exception
	{	Object arg = dict.get(key);
		if(arg==null)
			throw new Exception("The key "+key+" is missing when integer was expected");
		if(arg.getClass().equals(Double.class))
			return ((Double)arg).intValue();
		else if(!arg.getClass().equals(String.class))
			throw new Exception("The key "+key+" is not a number or a string");
		return JNUtilities.readInt((String)arg);
	}

	// get double for key
	// throws exception in an error
	public double gcis_doubleForKey(String key) throws Exception
	{	Object arg = dict.get(key);
		if(arg==null)
			throw new Exception("The key "+key+" is missing when double was expected");
		if(arg.getClass().equals(Double.class))
			return ((Double)arg).doubleValue();
		else if(!arg.getClass().equals(String.class))
			throw new Exception("The key "+key+" is not a number or a string");
		return JNUtilities.readDouble((String)arg);
	}
	
	// ----------------------------------------------------------------------------
	// Accessors
	// ----------------------------------------------------------------------------

	public HashMap<String, Object>getDictionary() { return dict; }

}
