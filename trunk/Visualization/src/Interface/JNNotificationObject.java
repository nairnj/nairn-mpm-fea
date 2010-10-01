/*******************************************************************
	JNNotificationObject.java
	JNFramework and NairnFEAMPMViz

	Created by John Nairn on 30 SEP 2010.
	Copyright (c) 2010 RSAC Software. All rights reserved.
*******************************************************************/

//package geditcom.JNFramework;

public class JNNotificationObject extends Object
{
	private String name;
	private Object object;
	private JNNotificationListener target;
	private Object userInfo;
	private Object sender;
	
	// initialize
	public JNNotificationObject(String notifyName,Object notifyObject,JNNotificationListener notifyTarget)
	{	name=new String(notifyName);
		object=notifyObject;
		target=notifyTarget;
	}
	
	// Accessors
	
	// name of the notification
	public String getName() { return name; }
	
	// object associated with the notification (typcially the sender, but need not be)
	public Object getObject() { return object; }
	
	// object to send this notification when requested
	public JNNotificationListener getTarget() { return target; }
	
	// sender - same as object unless object is null
	public Object getSender() { return sender; }
	public void setSender(Object theSender) { sender=theSender; }
	
	// user info
	public Object getUserInfo() { return userInfo; }
	public void setUserInfo(Object info) { userInfo=info; }
	
	// check if two objects have the same information
	public boolean equals(JNNotificationObject other)
	{	if(other.getName()!=name) return false;
		if(other.getObject()!=object) return false;
		if(other.getTarget()!=target) return false;
		return true;
	}
}
