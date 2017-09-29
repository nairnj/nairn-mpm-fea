/*
 * RegionsPiece.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 15 Nov 2016.
 * Copyright (c) 2016 RSAC Software. All rights reserved.
 */

import java.util.ArrayList;

public class RegionPiece
{
	private CmdViewer doc;
	private int type;
	private String xmlStart;
	private String shapeName;
	private int level;
	private RegionPiece parent;
	private ArrayList<RegionPiece> children;
	private double arcStart,arcEnd;
	
	public static final int RECT_OR_OVAL=1;
	public static final int POLY_PT=2;
	public static final int END_POLYGON=3;
	public static final int SHAPE_3D=4;
	public static final int COMMAND_PIECE=5;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	// num is type, xmlData is string up to initial attribute
	public RegionPiece(int num,String xmlData,String shape,CmdViewer cmdDoc)
	{	doc = cmdDoc;
		type = num;
		xmlStart = xmlData;
		shapeName = shape;
		level = 0;
		parent = null;
		children = new ArrayList<RegionPiece>(10);
		arcStart = -1.;
		arcEnd = 0.;
	}
	
	public String getXmlStart() { return xmlStart; }
	public void appendXmlStart(String moreXml)
	{	xmlStart = xmlStart + moreXml;
	}
	
	// convert to XML string
	public String xmlString(String indent)
	{	// extra indent
		String levelIndent = "";
		for(int i=0;i<level;i++) levelIndent = levelIndent+"  ";
		
		// build xml buffer
		StringBuffer xmlText;
		if(type==POLY_PT)
		{	xmlText = new StringBuffer(indent + levelIndent + "  <Polygon>\n" + indent + levelIndent + xmlStart);
		
			// add children
			for(int i=0;i<children.size();i++)
			{	RegionPiece obj = children.get(i);
				xmlText.append(obj.xmlString(indent));
			}
		
			// finish up
			xmlText.append(indent + levelIndent + "  </Polygon>\n");
		}
		
		else
		{	xmlText = new StringBuffer(levelIndent+xmlStart);
		
			// and empty element
			if(children.size()==0 && arcStart<0.)
			{	xmlText.append("/>\n");
			}
			else
			{	// finish start
				xmlText.append(">\n");
				
				// arc angles
				if(arcStart>=0.)
				{	xmlText.append(indent+levelIndent+"    ");
					xmlText.append("<arc start='"+doc.formatDble(arcStart)+"' end='"+doc.formatDble(arcEnd)+"'/>\n");
				}
				
				// add children
				for(int i=0;i<children.size();i++)
				{	RegionPiece obj = children.get(i);
					xmlText.append(obj.xmlString(indent));
				}
			
				// terminate
				xmlText.append(indent+levelIndent+"  </"+shapeName+">\n");
			}
		}
		
		// return as string
		return xmlText.toString();
	}
	
	public RegionPiece getParent() { return parent; }
	public void setParent(RegionPiece obj) { parent = obj; }
	
	public void addChild(RegionPiece child) { children.add(child); }
	
	public int getType() { return type; }
	public int getLevel() { return level; }
	public void setLevel(int newLevel) { level = newLevel; }
	
	public void setArcAngles(double astart,double aend)
	{	arcStart = astart;
		arcEnd = aend;
	}

}
