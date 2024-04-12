/*******************************************************************
	ControlPanel.java
	NairnFEAMPMViz

	Created by John Nairn on Fri Mar 05 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;

import java.awt.geom.*;
import java.util.*;

public class ControlPanel extends JPanel
{
	static final long serialVersionUID=2L;
	
	//----------------------------------------------------------------------------
	// variables and constants
	//----------------------------------------------------------------------------
	static final int TOP_MARGIN=16;
	static final int LEFT_MARGIN=16;
	static final int ROW_SPACING=24;
	static final int COL_SPACING=16;
	static final int WIDTH=340;
	static final int VERTICAL_LINE=0;
	static final int HORIZONTAL_LINE=1;
	
	private LoadArchive load;
	private TimeSelector selectTime;
	public PlotQuantity quantity;
	private PlotOptions options;
	private TimePlotOptions timeoptions;
	private CrackSelector thecrack;
	private LimitsSelector thelimits;
	protected PlotLaunch launch;
	private DocViewer docCtrl;
	private int vStart;
	
	//----------------------------------------------------------------------------
	// initialize
	//----------------------------------------------------------------------------

	public ControlPanel(DocViewer dc)
	{   super();
		setLayout(null);
		//setPreferredSize(new Dimension(800, 100));
		setBackground(Color.gray);
		docCtrl=dc;
		
		int vloc=TOP_MARGIN;
		int hloc=LEFT_MARGIN;
		
		// load the archive
		load=new LoadArchive(docCtrl);
		load.setLocation(hloc,vloc);
		add(load);
		vloc += load.getHeight()+ROW_SPACING;
		
		// quantity selector
		quantity=new PlotQuantity(docCtrl);
		quantity.setLocation(hloc,vloc);
		add(quantity);
		vloc += quantity.getHeight()+ROW_SPACING;
		vStart = vloc;
		
		// time selector
		selectTime=new TimeSelector(null,docCtrl);
		selectTime.setLocation(hloc,vloc);
		add(selectTime);
		vloc += selectTime.getHeight()+ROW_SPACING;
		
		// movie and time options at same vertical position
		// they are never both active
		options=new PlotOptions(docCtrl);
		options.setLocation(hloc,vloc);
		add(options);
		
		// time options selector
		timeoptions=new TimePlotOptions(docCtrl);
		timeoptions.setLocation(hloc,vloc);
		add(timeoptions);
		vloc += Math.max(options.getHeight(),timeoptions.getHeight())+ROW_SPACING;
		
		// crack selector and limts at same vertical location
		// They are never used at the same time
		thecrack=new CrackSelector(docCtrl);
		thecrack.setLocation(hloc,vloc);
		add(thecrack);
		
		// limits selector
		thelimits=new LimitsSelector(docCtrl);
		thelimits.setLocation(hloc,vloc);
		add(thelimits);
		vloc += Math.max(thecrack.getHeight(),thelimits.getHeight())+ROW_SPACING;
		
		// plot button selector
		launch=new PlotLaunch(docCtrl);
		launch.setLocation(hloc,vloc);
		add(launch);
		vloc += launch.getHeight()+20;
		
		setPreferredSize(new Dimension(hloc+WIDTH+LEFT_MARGIN, vloc));
	}

	//----------------------------------------------------------------------------
	// methods
	//----------------------------------------------------------------------------

	// draw frame
	protected void paintComponent(Graphics g)
	{	
		Graphics2D g2=(Graphics2D)g;
		
		// fill background (and erase prior view)
		g2.setColor(new Color(0.617f,0.685f,0.791f));
		Dimension d=getSize();
		g2.fill(new Rectangle2D.Double(0.,0.,(double)d.width,(double)d.height));
		
		int selected=load.getSelected();
		int hloc = LEFT_MARGIN;
		int vloc = vStart;
		
		g2.setColor(new Color(1.0f,0.0f,0.0f));
		g2.setStroke(new BasicStroke((float)3.));
		PlotControl aboveLimits = options;
		switch(selected)
		{	case LoadArchive.PARTICLE_PLOT:
			case LoadArchive.MESH_PLOT:
				connectControls(load,quantity,g2);
				if(docCtrl.resDoc.isMPMAnalysis())
				{	connectControls(quantity,selectTime,g2);
					vloc += selectTime.getHeight()+ROW_SPACING;
					options.setLocation(hloc,vloc);
					connectControls(selectTime,options,g2);
				}
				else
				{	options.setLocation(hloc,vloc);
					connectControls(quantity,options,g2);
				}
				vloc += options.getHeight()+ROW_SPACING;
				thelimits.setLocation(hloc,vloc);
				connectControls(aboveLimits,thelimits,g2);
				vloc += thelimits.getHeight()+ROW_SPACING;
				launch.setLocation(hloc,vloc);
				connectControls(thelimits,launch,g2);
				break;
			
			case LoadArchive.TIME_PLOT:
				connectControls(load,quantity,g2);
				if(timeoptions.isVisible())
				{	timeoptions.setLocation(hloc,vloc);
					connectControls(quantity,timeoptions,g2);
					vloc += timeoptions.getHeight()+ROW_SPACING;
					launch.setLocation(hloc,vloc);
					connectControls(timeoptions,launch,g2);
				}
				else if(thecrack.isVisible())
				{	// time plot of crack data
					thecrack.setLocation(hloc,vloc);
					connectControls(quantity,thecrack,g2);
					vloc += thecrack.getHeight()+ROW_SPACING;
					launch.setLocation(hloc,vloc);
					connectControls(thecrack,launch,g2);
				}
				else
				{	// global results
					launch.setLocation(hloc,vloc);
					connectControls(quantity,launch,g2);					
				}
				break;
				
			case LoadArchive.MESH2D_PLOT:
				connectControls(load,quantity,g2);
				if(docCtrl.resDoc.isMPMAnalysis())
				{	if(timeoptions.isVisible())
					{	selectTime.setVisible(true);
						connectControls(quantity,selectTime,g2);
						vloc += selectTime.getHeight()+ROW_SPACING;
						timeoptions.setLocation(hloc,vloc);
						connectControls(selectTime,timeoptions,g2);
						vloc += timeoptions.getHeight()+ROW_SPACING;
						launch.setLocation(hloc,vloc);
						connectControls(timeoptions,launch,g2);
					}
					else if(thecrack.isVisible())
					{	// time plot of crack data
						selectTime.setVisible(true);
						connectControls(quantity,selectTime,g2);
						vloc += selectTime.getHeight()+ROW_SPACING;
						thecrack.setLocation(hloc,vloc);
						connectControls(selectTime,thecrack,g2);
						vloc += thecrack.getHeight()+ROW_SPACING;
						launch.setLocation(hloc,vloc);
						connectControls(thecrack,launch,g2);
					}
					else
					{	// must be import
						selectTime.setVisible(false);
						launch.setLocation(hloc,vloc);
						connectControls(quantity,launch,g2);
					}
				}
				else
				{	timeoptions.setLocation(hloc,vloc);
					connectControls(quantity,timeoptions,g2);
					vloc += timeoptions.getHeight()+ROW_SPACING;
					launch.setLocation(hloc,vloc);
					connectControls(timeoptions,launch,g2);
				}
				break;
			
			default:
				break;
		}
	}
	
	// connect two plot controls with a line
	protected void connectControls(PlotControl source,PlotControl dest,Graphics2D g2)
	{	Rectangle loc1=source.getControlRect();
		Rectangle loc2=dest.getControlRect();
		
		// same column (assume directly beneath and centered)
		if(loc1.x==loc2.x)
		{	connectingLine(loc1.x+(loc1.width/2),loc1.y+loc1.height,
							loc2.x+(loc2.width/2),loc2.y,g2,VERTICAL_LINE);
		}
		
		// control 2 in higher column
		else if(loc1.x<loc2.x)
		{	// horizontal to the right
			if(loc2.y<loc1.y+loc1.height)
			{	connectingLine(loc1.x+loc1.width,loc1.y+(loc1.height/2),
								loc2.x,loc2.y+(loc2.height/2),g2,HORIZONTAL_LINE);
			}
			
			// vertical to the right
			else
			{	connectingLine(loc1.x+(loc1.width/2),loc1.y+loc1.height,
								loc2.x+(loc2.width/2),loc2.y,g2,VERTICAL_LINE);
			}
		}
		
		// control 2 in lower column
		else
		{	// horizontal to the left
			if(loc2.y<loc1.y+loc1.height)
			{	connectingLine(loc1.x,loc1.y+(loc1.height/2),
								loc2.x+(loc2.width),loc2.y+(loc2.height/2),g2,HORIZONTAL_LINE);
			}
			
			// vertical to the right
			else
			{	connectingLine(loc1.x+(loc1.width/2),loc1.y+loc1.height,
									loc2.x+(loc2.width/2),loc2.y,g2,VERTICAL_LINE);
			}
		}
			
	}
	
	protected void connectingLine(int x1,int y1,int x2,int y2,Graphics2D g2,int direction)
	{
		if(x1==x2 || y1==y2)
			g2.draw(new Line2D.Double((double)x1,(double)y1,(double)x2,(double)y2));
		
		else if(direction==HORIZONTAL_LINE)
		{	int xmid= x2>x1 ? x2-COL_SPACING/2 : x2+COL_SPACING/2;
			g2.draw(new Line2D.Double((double)x1,(double)y1,(double)xmid,(double)y1));
			g2.draw(new Line2D.Double((double)xmid,(double)y1,(double)xmid,(double)y2));
			g2.draw(new Line2D.Double((double)xmid,(double)y2,(double)x2,(double)y2));
		}
		
		else if(direction==VERTICAL_LINE)
		{	int ymid= y2>y1 ? y2-ROW_SPACING/2 : y2+ROW_SPACING/2;
			g2.draw(new Line2D.Double((double)x1,(double)y1,(double)x1,(double)ymid));
			g2.draw(new Line2D.Double((double)x1,(double)ymid,x2,(double)ymid));
			g2.draw(new Line2D.Double((double)x2,(double)ymid,(double)x2,(double)y2));
		}
	}
		
	// called when new file is loaded
	public void fileHasLoaded()
	{	load.setEnabled();
		hiliteControls();
		selectTime.select.setValue(0);
	}
	
	// call when select new plot option and may need to adjust controls
	public void hiliteControls()
	{
		// type of plot 0,1,2,3 (-1 is no plot)
		int selected=load.getSelected();
		
		// quantity and time
		quantity.setEnabled(selected);
		int plotComponent=getPlotComponent(selected,false,null);
		selectTime.setEnabled(selected,plotComponent);
		
		// options or time options
		options.setEnabled(selected);
		timeoptions.setEnabled(selected,plotComponent,docCtrl.resDoc);
		
		// limits and crack
		thelimits.setEnabled(selected,plotComponent);
		thecrack.setEnabled(selected,plotComponent);
		
		plotOpened();
		repaint();
	}
	
	// call when plot window opened or close or plot type changed
	public void plotOpened() { launch.updateTitle(load.getSelected()); }
	
	// call when select new plot component
	public void changeComponent()
	{	int selected=load.getSelected();
		int plotComponent=getPlotComponent(selected,false,null);
		timeoptions.setEnabled(selected,plotComponent,docCtrl.resDoc);
		thecrack.setEnabled(selected,plotComponent);
		
		/*
		// include this to replot on each menu change in control panel, but user might
		// prefer to change both before replotting and it activated plot window too
		// note that getPlotType() no longer valid
		if(docCtrl.movieFrame!=null)
		{	if(getPlotType()==docCtrl.movieFrame.plotType)
				JNNotificationCenter.getInstance().postNotification("PlotQuantityChanged",docCtrl,null);
		}
		*/
	}
	
	//----------------------------------------------------------------------------
	// progress bar
	//----------------------------------------------------------------------------
	
	public void enableProgress(int numSteps)
	{	launch.progress.setMaximum(numSteps);
		launch.progress.setValue(0);
		launch.progress.setEnabled(true);
	}
	
	public void disableProgress()
	{	launch.progress.setValue(0);
		launch.progress.setEnabled(false);
	}
	
	public void setProgress(int thisSteps)
	{	launch.progress.setValue(thisSteps);
	}
	
	public boolean isPlotting()
	{	return launch.progress.isEnabled();
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// get plot wanted
	public int getPlotComponent(int selected,boolean getExpression,ISDictType settings)
	{	return quantity.getPlotComponent(selected,getExpression,settings);
	}
	public JComboBox<PlotMenuItem> getQuantityMenu() { return quantity.quant; }
	public JComboBox<String> getComponentMenu() { return quantity.cmpnt; }
	public String getExpression() { return quantity.getExpression(); }
	public double evaluateExpressionMP(MaterialPoint mptr,double angle,ResultsDocument doc)
	{	return quantity.evaluateExpressionMP(mptr,angle,doc);
	}
	public double evaluateExpressionNode(ElementBase eptr,int ndi,double angle,ResultsDocument doc)
	{	return quantity.evaluateExpressionNode(eptr,ndi,angle,doc);
	}
	
	// set to item for checking the mesh
	public void setCheckMeshItem() { quantity.setCheckMeshItem(); }

	// get list of options for the plot
	public boolean [] getOptions() { return options.getOptions(); }
	public double getParticleSize() { return (double)options.particleSize; }
	public JSlider getParticleSizeSlider() { return options.mpmParticleSize; }
	
	// get currently selected archive
	public int getArchiveIndex() { return selectTime.getArchiveIndex(); }
	public boolean setArchiveIndex(int newIndex) { return selectTime.setArchiveIndex(newIndex); }
	public void updateTimeDisplay() { selectTime.updateLabel(); }
	public boolean incrementArchiveIndex() { return selectTime.incrementArchiveIndex(); }
	public boolean decrementArchiveIndex() { return selectTime.decrementArchiveIndex(); }
	
	// get current plotting type
	public int getPlotType() { return load.getSelected(); }
	public void changePlotType(int newType) { load.changeSelected(newType); }
	
	// get particle number option for time plots
	public int getParticleNumber(ISDictType settings) throws Exception
	{	if(settings==null)
			return timeoptions.getParticleNumber();
	
		// default if xy plot
		String ptype = (String)settings.gcis_objectForKey("plottype");
		if(ptype.contentEquals("xyplot")) return 1;
	
		// get material option, then wif needed get material or point number
		int mopt = settings.gcis_integerForKey("materialoption");
		if(mopt==1 || mopt==3)
			return 0;
		else if(mopt==2 || mopt==4)
			return -settings.gcis_integerForKey("materialnumber");
		else
			return settings.gcis_integerForKey("initialpoint");
	}
	
	// get particle number option for time plots
	public boolean getAveraging(ISDictType settings) throws Exception
	{ 	if(settings==null)
			return timeoptions.getAveraged();
	
		// default if xy plot
		String ptype = (String)settings.gcis_objectForKey("plottype");
		if(ptype.contentEquals("xyplot")) return true;

		// default to true if needed
		int mopt = settings.gcis_integerForKey("materialoption");
		if(mopt==1 || mopt==2) return true;
		return false;
	}
	
	// get particle number option for time plots
	public int getCrackNumber(ISDictType settings) throws Exception
	{ 	if(settings==null)
			return thecrack.getCrackNumber();
		return settings.gcis_integerForKey("cracknumber");
	}
	
	// get particle number option for time plots
	public int getCrackTip(ISDictType settings) throws Exception
	{	if(settings==null)
			return thecrack.getCrackTip();
		return settings.gcis_integerForKey("tipnumber");
	}
	
	// adjust limits if desired
	public Point2D.Double adjustLimits(double dmin,double dmax)
	{	MoviePlotWindow movieFrame=docCtrl.getMovieFrame();
		return thelimits.adjustLimits(dmin,dmax,movieFrame.plotView.dataMin,
					movieFrame.plotView.dataMax,movieFrame.plotView.dataLimitsSet);
	}
	
	// for time plots, adjust for total options
	// skip if averaging has been set
	public int adjustComponent(int theComponent,ISDictType settings) throws Exception
	{	if(getParticleNumber(settings)>0 || getAveraging(settings)) return theComponent;
		switch(theComponent)
		{	case PlotQuantity.MPMSTRENERGY:
				theComponent=PlotQuantity.MPMTOTSTRENERGY;
				break;
			case PlotQuantity.MPMKINENERGY:
				theComponent=PlotQuantity.MPMTOTKINENERGY;
				break;
			case PlotQuantity.MPMENERGY:
				theComponent=PlotQuantity.MPMTOTENERGY;
				break;
			case PlotQuantity.MPMWORKENERGY:
				theComponent=PlotQuantity.MPMTOTWORKENERGY;
				break;
			case PlotQuantity.MPMPLASTICENERGY:
				theComponent=PlotQuantity.MPMTOTPLASTICENERGY;
				break;
			case PlotQuantity.MPMHEATENERGY:
				theComponent=PlotQuantity.MPMTOTHEATENERGY;
				break;
			case PlotQuantity.MPMELEMENTCROSSINGS:
				theComponent=PlotQuantity.MPMTOTELEMENTCROSSINGS;
				break;
			default:
				break;
		}
		return theComponent;
	}
	
	// get selected item for contour menu options
	// if settings, integer 0 to 6 (return -1 error on not set)
	public int getContour(ISDictType settings)
	{	if(settings==null)
			return timeoptions.getContour();
		try
		{	return settings.gcis_integerForKey("variable");
		}
		catch(Exception e) { }
		return -1;
	}
	
	// get contour expression
	// if settings get string (return "" error if not set)
	public String getContourFunction(ISDictType settings)
	{	if(settings==null)
			return timeoptions.getContourFunction();
		String ct = (String)settings.gcis_objectForKey("contour");
		if(ct==null) return "";
		return ct;
	}
	
	// get contour +/- range (or zero if empty)
	// if settings return 0 if none provided
	public double getPlusMinus(ISDictType settings) throws Exception
	{	if(settings==null)
			return timeoptions.getPlusMinus();
		try
		{	return settings.gcis_doubleForKey("contourrange");
		}
		catch(Exception e) { }
		return 0.;
	}
	
	//----------------------------------------------------------------------------
	// static methods
	//----------------------------------------------------------------------------
	
	// read a field which must be an integer
	public static int readInteger(JTextField theField,String fieldName) throws Exception
	{
		Scanner validate=new Scanner(theField.getText());		// will use user's locale
		if(validate.hasNextInt())
		{	int theInt = validate.nextInt();
			validate.close();
			return theInt;
		}
		validate.close();
		throw new Exception("The "+fieldName+" must be an integer");
	}
	
	// read a field which must be an integer
	public static double readDouble(JTextField theField,String fieldName) throws Exception
	{
		Scanner validate=new Scanner(theField.getText());		// will use user's locale
		if(validate.hasNextDouble())
		{	double theDble = validate.nextDouble();
			validate.close();
			return theDble;
		}
		validate.close();
		throw new Exception("The "+fieldName+" must be an floating point number");
	}
}
