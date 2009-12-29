/*******************************************************************
	MeshPlotView.java
	NairnFEAMPMViz

	Created by John Nairn on 2/27/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.geom.*;
import java.awt.font.*;
import javax.swing.*;
import java.awt.image.*;

public class MeshPlotView extends JPanel
{
	static final long serialVersionUID=9L;
	
	//----------------------------------------------------------------------------
	// variable and constants
	//----------------------------------------------------------------------------
	
	static final int BORDER=5;
	static final int JUST_DRAW=0;
	static final int CENTER_LABEL=0x01;
	static final int FRAME_LABEL=0x02;
	static final int FILL_FRAME=0x04;				// currently only if centered
	static final int MOVE_LABEL=CENTER_LABEL;		// add options that move before drawing
	static final int KEY_HEIGHT=25;
	
	public static final int MPMPARTICLE_PLOTS=0;
	public static final int MPMMESH_PLOTS=1;
	public static final int FEAMESH_PLOTS=2;
	
	static Color backColor=new Color((float)0.329,(float)0.349,(float)0.427);
	static Color meshLineColor=new Color((float)0.75,(float)0.75,(float)0.75);
	static Color meshNodesColor=new Color((float)0.83,(float)0.83,(float)0.83);
	static Color textColor=new Color((float)1.,(float)1.,(float)1.);

	private double xpt,ypt;
	public Rectangle2D.Double xyBounds;
	private double scale;
	private boolean showMatPts,showSquarePts,showCrackPlanes,showCrackSurfaces;
	private boolean showMesh,showMeshBCs,showNodeNums,showElemNums,showMatPtNums;
	private boolean showNodes,showDisplaced,transformPts;
	private double mpDiam=50.;
	public boolean repainting=false;
	private ResultsDocument resDoc;
	private Graphics2D g2Loc;
	public double dataMin,dataMax;
	public boolean dataLimitsSet=false;
	private int plotComponent,plotType;
	private boolean firstLoad=false;
	
	//----------------------------------------------------------------------------
	// initialize
	//----------------------------------------------------------------------------
	
	MeshPlotView(ResultsDocument gResDoc)
	{
		resDoc=gResDoc;
	}
	
	//----------------------------------------------------------------------------
	// plotting
	//----------------------------------------------------------------------------

	// convert component graphics to a buffered image
	public BufferedImage frameImage()
	{
		Dimension d=getSize();
		BufferedImage image=new BufferedImage(d.width,d.height,BufferedImage.TYPE_INT_RGB);
		Graphics2D g2d=image.createGraphics();
		paintComponent(g2d);
		g2d.dispose();
		return image;
	}
	
	// draw frame
	protected void paintComponent(Graphics g)
	{   repainting=true;
		int i;
		
		g2Loc=(Graphics2D)g;
		g2Loc.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
		getRangeAndScaleView();
		
		// fill background (and erase prior view)
		g2Loc.setColor(backColor);
		Dimension d=getSize();
		g2Loc.fill(new Rectangle2D.Double(0.,0.,(double)d.width,(double)d.height));
		
		if(!firstLoad) return;
		
		// fill elements
		if(plotType!=MPMPARTICLE_PLOTS && plotComponent!=PlotQuantity.MESHONLY)
		{	setLineWidth(3.f);			// for interfaces
			for(i=0;i<resDoc.elements.size();i++)
			{	(resDoc.elements.get(i)).fill(this);
			}
		}
		
		// elements
		g2Loc.setColor(meshLineColor);
		setLineWidth(ElementBase.lineWidth);
		if(showMesh)
		{	for(i=0;i<resDoc.elements.size();i++)
			{	(resDoc.elements.get(i)).stroke(this,false);
			}
			
			// mesh BCs
			if(showMeshBCs)
			{	for(i=0;i<resDoc.gridBCs.size();i++)
				{	(resDoc.gridBCs.get(i)).stroke(this,resDoc);
				}
			}
		}
		
		// displaced mesh - if requested and if available
		if(showDisplaced)
		{	if(resDoc.feaArchFormat[ReadArchive.ARCH_FEADisplacements]=='Y')
			for(i=0;i<resDoc.elements.size();i++)
			{	(resDoc.elements.get(i)).stroke(this,true);
			}
		}
		
		// nodes
		if(showNodes)
		{	g2Loc.setColor(meshNodesColor);
			for(i=0;i<resDoc.nodes.size();i++)
			{	(resDoc.nodes.get(i)).stroke(this,resDoc,showDisplaced);
			}
		}
		
		// node numbers
		Font font=new Font("sanserif",Font.PLAIN,9);
		g2Loc.setFont(font);
		g2Loc.setColor(textColor);
		if(showNodeNums)
		{	for(i=0;i<resDoc.nodes.size();i++)
			{	(resDoc.nodes.get(i)).number(this,resDoc,showDisplaced);
			}
		}
		
		// element numbers
		if(showElemNums)
		{	setLineWidth(ElementBase.lineWidth/(float)2.);
			for(i=0;i<resDoc.elements.size();i++)
			{	(resDoc.elements.get(i)).number(this,resDoc,showDisplaced);
			}
		}
		setLineWidth((float)1.);
		
		// material point
		if(showMatPts)
		{	for(i=0;i<resDoc.mpmPoints.size();i++)
			{	resDoc.mpmPoints.get(i).stroke(this,resDoc);
			}
		}
		
		// particle BCs
		g2Loc.setColor(meshLineColor);
		if(showMeshBCs && plotType==MPMPARTICLE_PLOTS)
		{	for(i=0;i<resDoc.particleBCs.size();i++)
			{	(resDoc.particleBCs.get(i)).stroke(this,resDoc);
			}
		}
		
		// material point numbers
		g2Loc.setColor(textColor);
		if(showMatPtNums)
		{	for(i=0;i<resDoc.mpmPoints.size();i++)
			{	resDoc.mpmPoints.get(i).number(this,resDoc);
			}
		}
		
		// crack surfaces
		if(showCrackPlanes || showCrackSurfaces)
		{	for(i=0;i<resDoc.mpmCracks.size();i++)
			{	resDoc.mpmCracks.get(i).stroke(this,resDoc,showCrackPlanes,showCrackSurfaces);
			}
		}
		
		// draw the key
		FontRenderContext frc=g2Loc.getFontRenderContext();
		String plotMin=MoviePlotWindow.formatDouble(dataMin);
		Rectangle2D bounds=g2Loc.getFont().getStringBounds(plotMin,frc);
		g2Loc.setColor(textColor);
		g2Loc.drawString(plotMin,BORDER,d.height-BORDER);
		double begin=BORDER+bounds.getWidth()+5.;
		
		String plotMax=MoviePlotWindow.formatDouble(dataMax,false);
		bounds=g2Loc.getFont().getStringBounds(plotMax,frc);
		double end=(double)d.width-BORDER-bounds.getWidth()-5.;
		g2Loc.drawString(plotMax,(float)(end+5.),(float)(d.height-BORDER));
		
		int num=20;
		double segment=(end-begin)/(double)num;
		for(i=0;i<num;i++)
		{	g2Loc.setColor(ColorPicker.PickRainbow((double)i/(double)(num-1)));
			g2Loc.fill(new Rectangle2D.Double(begin+i*segment,(double)d.height-KEY_HEIGHT+BORDER,segment,(double)(KEY_HEIGHT-2*BORDER)));
		}
		
		repainting=false;
	}
	
	// find scaling to pixels
	public void getRangeAndScaleView()
	{
		Rectangle2D.Double meshBounds=resDoc.getMeshBounds(showDisplaced);
		
		// add border for boundary conditions
		double bcExtra=Math.max(meshBounds.getWidth(),meshBounds.getHeight())*BoundaryCondition.BC_SIZE*1.15;
		xyBounds=new Rectangle2D.Double(meshBounds.getX()-bcExtra,meshBounds.getY()-bcExtra,
							meshBounds.getWidth()+2.*bcExtra,meshBounds.getHeight()+2.*bcExtra);
		
		// adjust aspect ratio to match frame aspect ratio and get scalling (pixels/mm)
		Dimension d=getSize();
		double xyAspect=xyBounds.getHeight()/xyBounds.getWidth();
		double pixAspect=(double)(d.height-KEY_HEIGHT-2*BORDER)/(double)(d.width-2*BORDER);
		if(pixAspect>xyAspect)
			scale=((double)(d.width-2*BORDER))/xyBounds.getWidth();
		else
			scale=((double)(d.height-2*BORDER-KEY_HEIGHT))/xyBounds.getHeight();
    
		// calculate border rectangle in mesh coordinates
		double xBorder=((double)d.width/scale - xyBounds.getWidth())/2.;
		double yBorder=((double)(d.height-KEY_HEIGHT)/scale - xyBounds.getHeight())/2.;
		xyBounds=new Rectangle2D.Double(xyBounds.getX()-xBorder,xyBounds.getY()-yBorder,
							(double)d.width/scale,(double)d.height/scale);
	}
	
	// move to a point
	public void moveTo(double x,double y)
	{	xpt=scale*(x-xyBounds.getX());
		ypt=scale*(xyBounds.getY()+xyBounds.getHeight()-y);
		ypt-=KEY_HEIGHT;
	}
	
	// point in pixels to location
	public Point2D.Double getCoords(Point pixs)
	{	return new Point2D.Double((double)pixs.x/scale+xyBounds.getX(),
						xyBounds.getY()+xyBounds.getHeight()-(double)(pixs.y+KEY_HEIGHT)/scale);
	}
	
	// nudge position by pixels
	public void nudge(double dx,double dy)
	{	xpt+=dx;
		ypt-=dy;
	}
	
	// line from previous point to new point
	public void lineTo(double x,double y)
	{	double xpt2=scale*(x-xyBounds.getX());
		double ypt2=scale*(xyBounds.getY()+xyBounds.getHeight()-y);
		ypt2-=KEY_HEIGHT;
		g2Loc.draw(new Line2D.Double(xpt,ypt,xpt2,ypt2));
		xpt=xpt2;
		ypt=ypt2;
	}
	
	// draw circle (because only 1 argument)
	public void drawOval(double radius)
	{	double diam=2.*scale*radius;
		g2Loc.draw(new Ellipse2D.Double(xpt-diam/2.,ypt-diam/2.,diam,diam));
	}
	
	// draw circle (because only 1 argument)
	public void fillOval(double radius)
	{	double diam=2.*scale*radius;
		g2Loc.fill(new Ellipse2D.Double(xpt-diam/2.,ypt-diam/2.,diam,diam));
	}
	
	// draw material point and fill with plot color
	public void drawMaterialPoint(Color theColor,double[] eps,double[] eplast,double erot)
	{	double diam=0.01*mpDiam*resDoc.cellMinSide*scale;
		double radiix=resDoc.xscale*diam/2.;
		double radiiy=resDoc.yscale*diam/2.;
		g2Loc.setColor(theColor);
		if(showSquarePts)
		{	if(transformPts)
			{	double dgrad00=0.01*(eps[MaterialPoint.XXID]+eplast[MaterialPoint.XXID]);
				double wrot=Math.PI*erot/180.;
				double dgrad01=0.005*(eps[MaterialPoint.XYID]+eplast[MaterialPoint.XYID]-wrot);
				double dgrad10=0.005*(eps[MaterialPoint.XYID]+eplast[MaterialPoint.XYID]+wrot);
				double dgrad11=0.01*(eps[MaterialPoint.YYID]+eplast[MaterialPoint.YYID]);
				
				GeneralPath quad=new GeneralPath();
				Point2D.Float pathPt0=new Point2D.Float((float)(xpt+radiix*(1.+dgrad00)+radiiy*dgrad01),
											(float)(ypt+radiix*dgrad10+radiiy*(1.+dgrad11)));
				quad.moveTo(pathPt0.x,pathPt0.y);
				quad.lineTo((float)(xpt+radiix*(1.+dgrad00)-radiiy*dgrad01),
											(float)(ypt+radiix*dgrad10-radiiy*(1.+dgrad11)));
				quad.lineTo((float)(xpt-radiix*(1.+dgrad00)-radiiy*dgrad01),
											(float)(ypt-radiix*dgrad10-radiiy*(1.+dgrad11)));
				quad.lineTo((float)(xpt-radiix*(1.+dgrad00)+radiiy*dgrad01),
											(float)(ypt-radiix*dgrad10+radiiy*(1.+dgrad11)));
				quad.lineTo(pathPt0.x,pathPt0.y);
				g2Loc.fill(quad);
			}
			else
				g2Loc.fill(new Rectangle2D.Double(xpt-radiix,ypt-radiiy,2.*radiix,2.*radiiy));
		}
		else
			g2Loc.fill(new Ellipse2D.Double(xpt-radiix,ypt-radiiy,2.*radiix,2.*radiiy));
	}
	
	// transform shape to coordinates and fill in current color
	public void fillShape(GeneralPath shape)
	{	AffineTransform tx=new AffineTransform(scale,0.,0.,-scale,
				-scale*xyBounds.getX(),scale*(xyBounds.getY()+xyBounds.getHeight())-KEY_HEIGHT);
		Shape transformed=shape.createTransformedShape(tx);
		g2Loc.draw(transformed);
		g2Loc.fill(transformed);
	}
	
	// transform shape to coordinates and fill in current color
	public void strokeShape(GeneralPath shape)
	{	AffineTransform tx=new AffineTransform(scale,0.,0.,-scale,
				-scale*xyBounds.getX(),scale*(xyBounds.getY()+xyBounds.getHeight())-KEY_HEIGHT);
		Shape transformed=shape.createTransformedShape(tx);
		g2Loc.draw(transformed);
	}
	
	public void drawColor(Color theColor) { g2Loc.setColor(theColor); }
	
	// draw string at current point
	public void drawString(String dString,int option)
	{	if((option&MOVE_LABEL)==0)
			g2Loc.drawString(dString,(float)xpt,(float)ypt);
		if(option==JUST_DRAW) return;
		
		FontRenderContext frc=g2Loc.getFontRenderContext();
		Rectangle2D bounds=g2Loc.getFont().getStringBounds(dString,frc);
		double width=bounds.getWidth();
		double height=bounds.getHeight();
				
		if((option&CENTER_LABEL)!=0)
		{	xpt-=width/2.;
			ypt+=height/2.;
			if((option&FILL_FRAME)!=0)
			{	g2Loc.setColor(backColor);
				g2Loc.fill(new Rectangle2D.Double(xpt-1.,ypt-height,width+2.,height+1.));
				g2Loc.setColor(textColor);
			}
			g2Loc.drawString(dString,(float)xpt,(float)ypt);
		}
				
		if((option&FRAME_LABEL)!=0)
		{	g2Loc.draw(new Rectangle2D.Double(xpt-1.,ypt-height,width+2.,height+1.));
		}
	}
	
	// set basic stroke of certain width
	public void setLineWidth(float lineWidth) { g2Loc.setStroke(new BasicStroke(lineWidth)); }
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// set plot options
	public void setOptions(boolean [] options,int theComponent,int theType)
	{	showMesh=options[PlotOptions.SHOW_MESH];
		showMeshBCs=options[PlotOptions.SHOW_MESHBCS];
		showNodeNums=options[PlotOptions.SHOW_NODENUMS];
		showElemNums=options[PlotOptions.SHOW_ELEMNUMS];
		showMatPts=options[PlotOptions.SHOW_MATPTS];
		showMatPtNums=options[PlotOptions.SHOW_MATPTNUMS];
		showSquarePts=options[PlotOptions.SHOW_SQUAREPTS];
		showCrackPlanes=options[PlotOptions.SHOW_CRACKPLANES];
		showCrackSurfaces=options[PlotOptions.SHOW_CRACKSURFACES];
		showNodes=options[PlotOptions.SHOW_NODES];
		showDisplaced=options[PlotOptions.SHOW_DISPLACEDMESH];
		transformPts=options[PlotOptions.TRANSFORM_PTS];
		
		plotComponent=theComponent;
		plotType=theType;
		
		if(plotType!=MPMPARTICLE_PLOTS)
		{	showMatPts=false;
			showMatPtNums=false;
			showCrackPlanes=false;
			showCrackSurfaces=false;
		}
		if(plotType!=FEAMESH_PLOTS)
		{	showDisplaced=false;
		}
	}
	
	public int getPlotType() { return plotType; }
	public int getPlotComponent() { return plotComponent; }
	public boolean inDisplaced() { return showDisplaced; }
	public boolean getFirstLoad() { return firstLoad; }
	public void setFirstLoad(boolean fload) { firstLoad=true; }
}
