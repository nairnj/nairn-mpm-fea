/***************************************************************
 	JNNotificationCenter
 	JNFramework and NairnVEAMPMViz
 
 	1. Add notifications by name and object to send to target
 	2. object is associated with the notification (may be null to
 		get all notifications of that name)
 	3. target in JNNotificationListerner
 	4. Be sure to remove added notification if object goes away
 	5. Post as needed and sender should match the object to be heard
 
 	Created by John Nairn on 9/30/10.
	Copyright (c) 2010 RSAC Software. All rights reserved.
*****************************************************************/

//package geditcom.JNFramework;

import java.util.Hashtable;
import java.util.Vector;
import java.util.Set;

public class JNNotificationCenter extends Object
{
	// thread-safe singleton
	private static JNNotificationCenter instance = new JNNotificationCenter();
	
	// prevent other creation
	private JNNotificationCenter()
	{	
	}
	
	// get the instance
	public static JNNotificationCenter getInstance() { return instance; }
	
	// list of objects to send for each notification type
	Hashtable<String, Vector<JNNotificationObject>> registrations = new Hashtable<String,Vector<JNNotificationObject>>();
	
	// Register to send notification of type notificationName posted by object to target
	// If object is null, then send all notifications of type notificationName to target
	public void addNameAndObjectForTarget(String notificationName, Object object,JNNotificationListener target)
	{	Vector<JNNotificationObject> vector = registrations.get(notificationName);
		JNNotificationObject notify=new JNNotificationObject(notificationName,object,target);
		if (vector==null)
		{	vector = new Vector<JNNotificationObject>();
			vector.add(notify);
			registrations.put(notificationName,vector);
		}
		else if(!notificationContains(vector,notify))
		{	vector.add(notify);
		}
		//System.out.println("Added "+notificationName+" for "+object.getClass().getName()+" target "+target.getClass().getName()+", now have "+vector.size());
	}
	
	// return true is this notification contains the same object
	private boolean notificationContains(Vector<JNNotificationObject> vector,JNNotificationObject notify)
	{	for(JNNotificationObject notifyObject : vector)
		{	if (notifyObject.equals(notify))
				return true;
		}
		return false;
	}
	
	// remove all notifications wth any name for a target
	// Objects that register for notifications must remove themselves before they disappear
	public void removeAllForTarget(JNNotificationListener target)
	{	Set<String> keySet=registrations.keySet();
		Object key[]=keySet.toArray();
		int i;
		for(i=keySet.size()-1;i>=0;i--)
		{	removeAllForTarget((String)key[i],target);
		}
	}
	
	// remove all notifications with given name for a target
	// Objects that register for notifications must remove themselves before they disappear
	public void removeAllForTarget(String notificationName, JNNotificationListener target)
	{	Vector<JNNotificationObject> vector = registrations.get(notificationName);
		if(vector==null) return;
		int index;
		for(index=vector.size()-1;index>=0;index--)
		{	JNNotificationObject notifyObject=vector.get(index);
			if(notifyObject.getTarget()==target)
			{	//System.out.println("Removed "+notificationName+" for "+target.getClass().getName());
				vector.remove(index);
			}
		}
		checkEmpty(notificationName,vector);
	}
	
	// remove a notification name if it has no targets
	private void checkEmpty(String notificationName,Vector<JNNotificationObject> vector)
	{	if(vector.size()==0)
		{	vector=null;
			registrations.remove(notificationName);
		}
	}

	// Post notificaiton by sender to any target that wants to hear about it
	// Those that want to hear have registered for this notification name with wither
	//		object==object or object==null. The object is stored in sender member of 
	//		the JNNotificationObject. It will be same as object unless object is null
	//		info is any object with more details about the notification
	public void postNotification(String notificationName,Object sender,Object info)
	{	Vector<JNNotificationObject> vector = registrations.get(notificationName);
		if(vector==null || vector.size()==0) return;
		for (JNNotificationObject notify : vector)
		{	if(notify.getObject()==null || notify.getObject()==sender)
			{	notify.setUserInfo(info);
				notify.setSender(sender);
				//System.out.println("Send "+notificationName+" to "+notify.getTarget().getClass().getName());
				notify.getTarget().receiveNotification(notify);
			}
		}
	}
	
}
