/*******************************************************************
	Main.java
	NairnFEAMPMViz
	
	static final long serialVersionUID=33L;
		Used 1-27 and 30-35 and 37

	Created by John Nairn on  3/1/07.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import javax.swing.*;

public class Main
{
	public static final int NATIVE_LNF=0;
	public static final int XPLATFORM_LNF=1;
	public static int useLnf=NATIVE_LNF;
	
	//Schedule a job for the event-dispatching thread:
	//creating and showing this application's GUI.
	public static void main(String[] args)
	{	javax.swing.SwingUtilities.invokeLater(new Runnable()
		{	public void run()
			{	createAndShowGUI();
            }
        });
    }
	
	// Create the GUI and show it.  For thread safety,
	// this method should be invoked from the
	// event-dispatching thread.
    private static void createAndShowGUI()
	{	// undocumented anti-aliasing of text in swing
		String os=System.getProperty("os.name").toLowerCase();
		if(os.indexOf("mac")>=0 || os.indexOf("win")>=0)
		{	String antialising = "swing.aatext";
			if (null == System.getProperty(antialising))
				System.setProperty(antialising, "true");
		}
		
		try
		{	if(useLnf==XPLATFORM_LNF)
				UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
			else
				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}
		catch(Exception e)
		{
		}
			
		//Make sure we have nice window decorations.
        JFrame.setDefaultLookAndFeelDecorated(true);

        // Create the application controller
        new NairnFEAMPMViz();

        //Display the window.
        //frame.setVisible(true);
    }
}
