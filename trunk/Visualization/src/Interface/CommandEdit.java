/*******************************************************************
	CommandEdit.java
	NairnFEAMPMViz

	Created by John Nairn on Sat Mar 06 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;

public class CommandEdit extends JScrollPane
{
	static final long serialVersionUID=34L;
	
	//-----------------------------------------------------------------
	// Class variables and constants
	//-----------------------------------------------------------------
	
	// members
	public JTextArea textPane = new JTextArea( );
	private boolean changed=false;
	
	// command file chooser
	private JFileChooser chooser=new JFileChooser();
  
	//-----------------------------------------------------------------
	// Initialize - scrolling edit field
	//-----------------------------------------------------------------
	
	public CommandEdit()
	{	super();
		setViewportView(textPane);
		
		// text editing panel
		//textPane.setEditable(false);
		textPane.setFont(new Font("Monospaced",Font.PLAIN,12));
		
		textPane.getDocument().addDocumentListener(new DocumentListener()
		{   public void insertUpdate(DocumentEvent e) { changed=true; }
			public void removeUpdate(DocumentEvent e) { changed=true; }
			public void changedUpdate(DocumentEvent e) {}
		});
		
		// set file filter
		NFMVFilter filter=new NFMVFilter();
		filter.addExtension("fmcmd");
		filter.addExtension("fcmd");
		filter.addExtension("mcmd");
		filter.addExtension("cmd");
		filter.setDescription("FEA or MPM Input Files");
		chooser.setFileFilter(filter);
		
		NFMVPrefs.setWorkspace(chooser);
	}
	
	//--------------------------------------------------------------------
	// Opening new file starts here
	//--------------------------------------------------------------------
	
	// load text into the command field
	public void setCommands(String fileText)
	{	textPane.setText(fileText);
		textPane.select(0,0);
		changed=false;
	}

	// save commands and return saved file (or null)
	public File saveCommands(File saveFile,boolean getNameFirst,JFrame fileWindow)
	{
		// exit if not changed (i.e. already saved)
		if(!changed && !getNameFirst && saveFile!=null)
			return saveFile;
		
		// get file if need it
		if(getNameFirst || saveFile==null)
		{	int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return null;
			
			saveFile=NairnFEAMPMViz.CheckFileStatus(chooser.getSelectedFile(),fileWindow,"fmcmd");
			if(saveFile==null) return null;
		}
		
		// save to saveFile
		try
		{	FileWriter theFile=new FileWriter(saveFile);
			theFile.write(textPane.getText());
			theFile.flush();
			theFile.close();
		}
		catch (Exception fe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(fileWindow,"Error writing XML input commands: " + fe);
			return null;
		}
		
		changed=false;
		return saveFile;
	}
	
	// save commands and return saved file (or null)
	public File saveCopyOfCommands(File saveFile,JFrame fileWindow)
	{
		// save to saveFile
		try
		{	FileWriter theFile=new FileWriter(saveFile);
			theFile.write(textPane.getText());
			theFile.flush();
			theFile.close();
		}
		catch (Exception fe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(fileWindow,"Error writing XML input  commands: " + fe);
			return null;
		}
		
		return saveFile;
	}
	
	// insert DTD path
	public boolean insertDTD(String dtdPath)
	{
		String initText=textPane.getText();
		int offset=initText.indexOf("<!DOCTYPE"),endOffset=0;
		int docOffset=offset;
		if(offset>0)
		{	offset=initText.indexOf("JANFEAInput",offset+9);
			if(offset>0)
			{	offset=initText.indexOf("SYSTEM",offset+6);
			
				// if filed <!DOCTYPE JANFEAInput SYSEM located quoted path name
				if(offset>0)
				{	char startChar=' ';
					while(offset<initText.length())
					{	startChar=initText.charAt(offset);
						if(startChar=='\'' || startChar=='\"') break;
						if(startChar=='\n' || startChar=='\r')
						{	offset=0;
							break;
						}
						offset++;
					}
					if(offset>0)
					{	offset++;
						endOffset=offset;
						while(endOffset<initText.length())
						{	char endChar=initText.charAt(endOffset);
							if(endChar==startChar) break;
							if(endChar=='\n' || endChar=='\r')
							{	offset=0;
								break;
							}
							endOffset++;
						}
					}
				}
			}
		}
		if(offset>0)
		{	// insert path if needed
			String oldDTD=initText.substring(offset,endOffset);
			if(!oldDTD.equals(dtdPath))
			{	textPane.setText(initText.substring(0,offset)+dtdPath+
							initText.substring(endOffset,initText.length()));
				changed=true;
			}
		}
		else if(docOffset>0)
		{	docOffset+=9;
			textPane.setText(initText.substring(0,docOffset)+" JANFEAInput SYSTEM '"
				+dtdPath+"'"+initText.substring(docOffset,initText.length()));
			changed=true;
		}
		else
		{	// insert <!DOCTYPE
			docOffset=initText.indexOf("?>");
			if(docOffset<=0) return false;
			docOffset+=2;
			textPane.setText(initText.substring(0,docOffset)
						+"\n<!DOCTYPE JANFEAInput SYSTEM '"+dtdPath+"'>"
						+initText.substring(docOffset,initText.length()));
			changed=true;
		}
		
		return true;
	}

	//--------------------------------------------------------------------
	// Accessors
	//--------------------------------------------------------------------
	
	public boolean getChanged() { return changed; }
	
	public boolean isMPMAnalysis()
	{	return textPane.getText().indexOf("<MPMHeader>")>0;
	}

}