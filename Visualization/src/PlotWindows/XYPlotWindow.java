/*******************************************************************
	XYPlotWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 1/22/08.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import geditcom.JNFramework.*;
import java.util.*;

import javax.swing.SwingUtilities;

import java.awt.Toolkit;
import java.awt.geom.*;
//import com.primalworld.math.*;

public class XYPlotWindow extends TwoDPlotWindow implements Runnable
{
	static final long serialVersionUID=25L;
	
	// constants
	static final int XorRContour=0;
	static final int YorZContour=1;
	static final int DContour=2;
	static final int TContour=3;
	static final int PPtAndDirectionLine=4;
	static final int P0PtAndDirectionLine=5;
	static final int CRACK_CONTOUR=10;
	static final int contourPoints=200;
	static final int AVGQUAD=10;
	private static double[] xq={-0.97390652,-0.8650633666,-0.6794095682,
		-0.4333953941,-0.1488743389,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.97390652};
	private static double[] wq={0.06667134,0.1494513491,0.2190863625,
		0.2692667193,0.2955242247,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.06667134};

	// variables
	private String var1=null;
	private String var2=null;
	private JNExpression contourExpr = null;
	HashMap<String,String> variables;
	private double contourRange;
	private int pmpts;
	private int functionType;
	private JNVector x0 = new JNVector(0.,0.,0.);
	private JNVector dir = new JNVector(0.,0.,0.);
	
	// for plot thread calculations
	private int component;
	private ControlPanel controls;
	private ISDictType settings;

	// initialize
	public XYPlotWindow(DocViewer parent)
	{	super(parent);
		setChildType("xyPlotFrame");
	}
	
	// add another plot
	public void addPlot(ControlPanel plotControls,ISDictType plotSettings) throws Exception
	{
		// the component to plot
		settings=plotSettings;
		if(settings!=null)
			settings.gcis_setObjectforKey("xyplot","plottype");
		controls=plotControls;
		component=controls.getPlotComponent(LoadArchive.MESH2D_PLOT,true,settings);
		if(component<0)
			throw new Exception("The result to be plotted was not recognized");
		component=controls.adjustComponent(component,settings);
		
		// plot file of global results
		if(component==PlotQuantity.IMPORTANDPLOTFILE)
		{	importAndPlotFile();			
			return;
		}
		
		// detach thread to gather plot information
		Thread plot2DThread=new Thread(this);
		plot2DThread.start();
	}
	
	public boolean isPlotting()
	{	return controls.isPlotting();
	}
	
	//----------------------------------------------------------------------------
	// detachable thread for loading time plot data
	//----------------------------------------------------------------------------

	public void run()
	{
		controls.enableProgress(contourPoints);
		
		String contourFxn=null;
		
		try
		{	
			switch(component)
			{	case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
				case PlotQuantity.MPMCRACKPROFILE:
				case PlotQuantity.MPMOPENINGFRACTION:
				case PlotQuantity.MPMSHEARFRACTION:
				case PlotQuantity.MPMCZMGI:
				case PlotQuantity.MPMCZMGII:
				case PlotQuantity.MPMTRACTION1:
				case PlotQuantity.MPMTRACTION2:
				case PlotQuantity.MPMTRACTION3:
				case PlotQuantity.MPMTRACTION4:
				case PlotQuantity.MPMTRACTION5:
				case PlotQuantity.MPMTRACTION6:
				case PlotQuantity.MPMTRACTION7:
				case PlotQuantity.MPMTRACTION8:
				case PlotQuantity.MPMTRACTION9:
				case PlotQuantity.MPMTRACTION10:
					functionType=CRACK_CONTOUR;
					break;
				
				default:
					functionType=controls.getContour(settings);
					contourFxn=controls.getContourFunction(settings);
					
					if(contourFxn.length()==0 || functionType<0)
						throw new Exception("The contour function cannot be empty");

					if(functionType==PPtAndDirectionLine || functionType==P0PtAndDirectionLine)
					{	// 2D needs 1, 2, or 4 comma separated values
						// 3D needs 1, 2, 3, or 6
						ResultsDocument resDoc=((DocViewer)document).resDoc;
						double[] args = new double[6];
						int numArgs = 0;
						int maxArgs = resDoc.is3D() ? 6 : 4 ;
						Scanner argScan=new Scanner(contourFxn);
						argScan.useDelimiter(",");
						while(argScan.hasNextDouble())
						{	// get full arguments and trim whitespace
							if(numArgs>=maxArgs)
							{	argScan.close();
								throw new Exception("More than "+maxArgs+"command-separated values not allowed");
							}
							args[numArgs] = argScan.nextDouble();
							numArgs++;
						}
						argScan.close();
						
						if(resDoc.is3D())
						{	if(numArgs==4)
							{	// plot along direction in item 4
								x0.Set(args[0],args[1],args[2]);
								dir.Set(0.,0.,0.);
								int axis = (int)(args[3]+.1)-1;
								if(axis<0 || axis>2)
									throw new Exception("Axis argument #4 must be 1, 2, or 3");
								dir.v[axis]=1.;
								if(axis==0)
									x0.v[0]=resDoc.xmin;
								else if(axis==1)
									x0.v[1]=resDoc.ymin;
								else
									x0.v[2]=resDoc.zmin;
							}
							else if(numArgs==6)
							{	// point and direction
								x0.Set(args[0],args[1],args[2]);;
								dir.Set(args[3],args[4],args[5]);
								double dirLength = dir.getLength();
								dir.Scale(1./dirLength);
							}
							else
								throw new Exception("3D particle plots need 4 or 6 comma-separated values");
						}
						else if(numArgs==1)
						{	// 2D plot for x a constant
							x0.Set(args[0],resDoc.ymin,0.);
							dir.Set(0.,1.,0.);
						}
						else if(numArgs==2)
						{	// 2D plot for y a constant
							x0.Set(resDoc.xmin,args[1],0.);
							dir.Set(1.,0.,0.);
						}
						else if(numArgs==4)
						{	// point and direction
							x0.Set(args[0],args[1],0.);
							dir.Set(args[2],args[3],0.);
							double dirLength = dir.getLength();
							dir.Scale(1./dirLength);
						}
						else
							throw new Exception("2D particle plots need 1, 2, or 4 comma-separated values");
							
						// get extra range
						contourRange=controls.getPlusMinus(settings);
					}
					
					else
					{	switch(functionType)
						{	case XorRContour:
								var1="y";
								var2="z";
								break;
							case YorZContour:
								var1="x";
								var2="r";
								break;
							case DContour:
								var1="T";
								break;
							case TContour:
								var1="D";
								break;
							default:
								break;
						}
						
						// test it
						contourExpr = new JNExpression(contourFxn,null);
						variables = new HashMap<String,String>(2);
						try
						{	evaluate(2.d);
						}
						catch(Exception e)
						{	throw new Exception("The contour function does not evaluate to a number");
						}
							
						// get range
						contourRange=controls.getPlusMinus(settings);
						if(ElementBase.DbleEqual(contourRange,0.))
							pmpts=1;
						else
							pmpts=AVGQUAD;
					}
					break;
			}
			
			// load the archive
			ResultsDocument resDoc=((DocViewer)document).resDoc;
			if(resDoc.isMPMAnalysis())
			{	int selArchive;
				if(settings!=null)
					selArchive = settings.gcis_integerForKey("archive")-1;
				else
					selArchive = controls.getArchiveIndex();
				resDoc.readSelectedArchive(selArchive);
			}
			
			// load element plot values
			if(functionType==PPtAndDirectionLine || functionType==P0PtAndDirectionLine)
				plotParticleResults();
			else if(functionType!=CRACK_CONTOUR)
			{	ElementBase.load2DPlotData(component,resDoc);
				plotXYResults();
			}
			else
				plotCrackResults();
		}
		catch(Exception pe)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(this,"X-Y plot error: "+pe.getMessage());
			if(plot2DView.getNumberOfPlots()==0) dispose();
			controls.disableProgress();
			return;
		}
		
		// finish up on EDT
		SwingUtilities.invokeLater(new Runnable()
		{	public void run()
			{	controls.disableProgress();
				setVisible(true);
				toFront();
			}
		});
	}
	
	// evaluate expression at tval
	public double evaluate(double tval) throws Exception
	{	variables.put(var1,JNUtilities.formatDouble(tval));
		if(var2!=null) variables.put(var2,JNUtilities.formatDouble(tval));
		contourExpr.evaluateWith(variables);
		return contourExpr.getNumericValue();
	}
	
	// entry point to get non-crack, x-y plot points
	public void plotXYResults() throws Exception
	{
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(contourPoints);
		ArrayList<Double> y=new ArrayList<Double>(contourPoints);
		
		// contour plots
		int i,iq;
		double pmin,pstep,separation,plotValue,rValue,pmdelx,pmdely;
		double mptX,mptY,prevX,prevY,thePtX,thePtY;
		double dval,aval;
		boolean hasValue;
		ElementBase theElem=null,prevElem;
		NodalPoint[] eNodes=new NodalPoint[9];
		
		//---------------------------------------------------------
		// begin plot array
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		Rectangle2D.Double mpmMeshBounds=resDoc.getMeshBounds(false);
		switch(functionType)
		{	case XorRContour:
				pmin=mpmMeshBounds.getY();
				pstep=mpmMeshBounds.getHeight()/(double)(contourPoints-1);
				break;
			case YorZContour:
				pmin=mpmMeshBounds.getX();
				pstep=mpmMeshBounds.getWidth()/(double)(contourPoints-1);
				break;
			case DContour:
			case TContour:
			default:
				throw new Exception("Countour type selected is not yet implemented");
		}
		double distance=pmin;
		
		// initialize a previous point for finding normals
		prevX=pmin-pstep;
		try
		{	prevY=evaluate(prevX);
		}
		catch(Exception e)
		{	throw new Exception("The contour function does not evaluate to a number at "+prevX);
		}
		if(functionType==XorRContour || functionType==DContour)
		{	double temp=prevX;
			prevX=prevY;
			prevY=temp;
		}
		if(functionType==DContour || functionType==TContour)
		{	// x is d and y is angle
			dval=prevX;
			aval=prevY;
			prevX=dval*Math.cos(aval);
			prevY=dval*Math.sin(aval);
		}
		
		//---------------------------------------------------------
		// loop plot points
		for(i=0;i<contourPoints;i++)
		{	// progress
			controls.setProgress(i);
			
			// plot
			mptX=pmin+((double)i)*pstep;
			try
			{	mptY = evaluate(mptX);
			}
			catch(Exception e)
			{	throw new Exception("The contour function does not evaluate to a number at "+mptX);
			}
			if(functionType==XorRContour || functionType==DContour)
			{	double temp=mptX;
				mptX=mptY;
				mptY=temp;
			}
			if(functionType==DContour || functionType==TContour)
			{	// x is d and y is angle
				dval=mptX;
				aval=mptY;
				mptX=dval*Math.cos(aval);
				mptY=dval*Math.sin(aval);
			}
			
			// possibly average over a range
			separation=Math.sqrt((mptX-prevX)*(mptX-prevX)+(mptY-prevY)*(mptY-prevY));
			plotValue=0.;
			rValue=0.;
			hasValue=false;
			for(iq=0;iq<pmpts;iq++)
			{	if(pmpts>1)
				{	pmdelx=(mptY-prevY)/separation;
					pmdely=(prevX-mptX)/separation;
					thePtX=mptX + xq[iq]*pmdelx*contourRange;
					thePtY=mptY + xq[iq]*pmdely*contourRange;
				}
				else
				{	thePtX=mptX;
					thePtY=mptY;
				}
			
				// find the element
				prevElem=theElem;
				theElem=ElementBase.findPtElement(thePtX,thePtY,theElem,resDoc);
			
				// find plot value if found an element
				if(theElem!=null)
				{	if(theElem!=prevElem) theElem.getNodalPoints(eNodes,resDoc.nodes);
					if(pmpts>1)
					{	if(resDoc.isAxisymmetric())
						{	plotValue+=wq[iq]*thePtX*theElem.findValueAtRealPt(thePtX,thePtY,eNodes);
							rValue+=wq[iq]*thePtX;
						}
						else
							plotValue+=wq[iq]*theElem.findValueAtRealPt(thePtX,thePtY,eNodes)/2.;
					}
					else
						plotValue=theElem.findValueAtRealPt(thePtX,thePtY,eNodes);
					hasValue=true;
				}
			}
			
			// add it now if found one
			if(hasValue)
			{	x.add(distance);
				if(pmpts>1 && resDoc.isAxisymmetric() && rValue>0.)
					y.add(plotValue/rValue);
				else
					y.add(plotValue);
			}
			
			// increment x axis
			if(i>0) distance+=separation;
			prevX=mptX;
			prevY=mptY;
		}

		//---------------------------------------------------------
		// Finish plot
		if(x.size()==0)
		{	throw new Exception("The selected plot contour did not pass through the mesh.");
		}
		else
		{	Hashtable<String,String> props = new Hashtable<String,String>();
			props.put("object.color",plot2DView.selectPlotColor());
			props.put("array.name",PlotQuantity.plotName(component,resDoc));
			plot2DView.plotData(x,y,props);
				
			// axis labels
			if(plot2DView.getNumberOfPlots()<2)
			{	plot2DView.setXTitle("Position ("+resDoc.units.lengthUnits()+")");
				plot2DView.setYTitle(PlotQuantity.plotLabel(component,resDoc));
			}
		}
	}
	
	// entry point to plot particle values near a line
	public void plotParticleResults() throws Exception
	{
		ResultsDocument resDoc=((DocViewer)document).resDoc;

		// go through points
		// get the plot
		ArrayList<PointData> pdata=new ArrayList<PointData>(100);
	    MaterialPoint mpm;
		for(int i=0;i<resDoc.mpmPoints.size();i++)
		{	mpm=resDoc.mpmPoints.get(i);

			// particle plots should not include rigid or reservoirparticles
			MaterialBase matl=resDoc.materials.get(mpm.materialIndex());
			if(matl.isRigid() || mpm.inReservoir()) continue;
			
			// check point
	        JNVector pt;
	        if(functionType==PPtAndDirectionLine)
	            pt = mpm.getPosition();
	        else
	            pt = mpm.getOrigPosition();
			Vector3 radius = mpm.getParticleRadius(resDoc);
			
			// find distance to the line

			// vector from x0 to point
			JNVector pmx = pt.Subtract(x0);
			
			// projection onto line (when dir normalized it is distance to intersection point
			double proj = pmx.Dot(dir);
			
			// intersection point
			JNVector cross = new JNVector(x0.v[0] + proj*dir.v[0],x0.v[1]+proj*dir.v[1],x0.v[2]+proj*dir.v[2]);
			
			// distance point to the line
			double dist = pt.getDistanceTo(cross);
			
			// is it close enough?
			if(dist<=Math.max(radius.x,contourRange))
			{	// add this point
				double plotValue = mpm.getForPlot(component,0.,resDoc);
				pdata.add(new PointData(proj,plotValue));
			}
		}
		
		if(pdata.size()==0)
		{	throw new Exception("The selected particle plot contour did not pass near any material points.");
		}
		
		// sort by location
		Collections.sort(pdata);

		//---------------------------------------------------------
		// Finish plot
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(contourPoints);
		ArrayList<Double> y=new ArrayList<Double>(contourPoints);
		for(int i=0;i<pdata.size();i++)
		{	PointData mpt = pdata.get(i);
			x.add(mpt.x);
			y.add(mpt.y);
		}
		
		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		String ptitle = PlotQuantity.plotName(component,resDoc) + " (";
		if(functionType==PPtAndDirectionLine)
			ptitle = ptitle + "p=";
		else
			ptitle = ptitle + "p0=";
		if(resDoc.is3D())
			ptitle = ptitle + x0.toString() + "+d"+ dir.toString();
		else
		{	JNVector x02D = new JNVector(x0.v[0],x0.v[1]);
			JNVector dir2D = new JNVector(dir.v[0],dir.v[1]);
			ptitle = ptitle + x02D.toString() + "+d"+ dir2D.toString();
		}
		if(contourRange>0.)
			ptitle = ptitle + ", +/- " + contourRange + ")";
		else
			ptitle = ptitle + ")";
		props.put("array.name",ptitle);
		plot2DView.plotData(x,y,props);
				
		// axis labels
		if(plot2DView.getNumberOfPlots()<2)
		{	plot2DView.setXTitle("Position ("+resDoc.units.lengthUnits()+")");
			plot2DView.setYTitle(PlotQuantity.plotLabel(component,resDoc));
		}
	}
	
	// entry point to get non-crack, x-y plot points
	public void plotCrackResults() throws Exception
	{
		// does crack exist?
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		int crackNum=controls.getCrackNumber(settings);
		if(crackNum>resDoc.mpmCracks.size())
			throw new Exception("The selected crack number does not exist.");
		CrackHeader header=resDoc.mpmCracks.get(crackNum-1);
		
		// read crack surface points
		// array for results
		int numseg=header.segments.size();
		ArrayList<Double> xa=new ArrayList<Double>(numseg);
		ArrayList<Double> ya=new ArrayList<Double>(numseg);
		header.getSurface(CrackSegment.ABOVE_POS,xa,ya);
		ArrayList<Double> xb=new ArrayList<Double>(numseg);
		ArrayList<Double> yb=new ArrayList<Double>(numseg);
		header.getSurface(CrackSegment.BELOW_POS,xb,yb);
		
		// traction history
		ArrayList<Double> hpts = null;
		if(component>=PlotQuantity.MPMTRACTION1 && component<=PlotQuantity.MPMTRACTION10)
		{	hpts=new ArrayList<Double>(numseg);
			header.getCrackTraction(component-PlotQuantity.MPMTRACTION1+1,hpts);
		}
	    else if(component==PlotQuantity.MPMCZMGI || component==PlotQuantity.MPMCZMGII)
		{	hpts=new ArrayList<Double>(numseg);
			header.getCrackZoneERR(component,hpts);
		}
		
		controls.setProgress((int)contourPoints/2);
		
		// crack profile
		if(component==PlotQuantity.MPMCRACKPROFILE)
		{	// top surface
			Hashtable<String,String> props = new Hashtable<String,String>();
			props.put("object.color",plot2DView.selectPlotColor());
			props.put("array.name","Crack "+crackNum+" top surface");
			plot2DView.plotData(xa,ya,props);
			
			// bottom surface
			props.put("object.color",plot2DView.selectPlotColor());
			props.put("array.name","Crack "+crackNum+" bottom surface");
			plot2DView.plotData(xb,yb,props);
			
			if(plot2DView.getNumberOfPlots()<3)
			{	plot2DView.setXTitle("Position ("+resDoc.units.lengthUnits()+")");
				plot2DView.setYTitle("Position ("+resDoc.units.lengthUnits()+")");
			}
		}
		
		// normal and shear COD or opening and sliding fraction or traction history
		else
		{	double distance=0.;
			double mag;
			
			// array for results
			ArrayList<Double> x=new ArrayList<Double>(contourPoints);
			ArrayList<Double> y=new ArrayList<Double>(contourPoints);
			
			// calculational points
			Point2D.Double endA=new Point2D.Double(xa.get(0),ya.get(0));
			Point2D.Double endB=new Point2D.Double(xb.get(0),yb.get(0));
			Point2D.Double planePt0=new Point2D.Double((endA.x+endB.x)/2.,(endA.y+endB.y)/2.);
			Point2D.Double planePt1=new Point2D.Double(0.,0.);
			Point2D.Double ctod0=new Point2D.Double(endA.x-endB.x,endA.y-endB.y);
			Point2D.Double ctod1=new Point2D.Double(0.,0.);
			Point2D.Double tangential=new Point2D.Double(0.,0.);
			Point2D.Double normal=new Point2D.Double(0.,0.);
			Point2D.Double modes=new Point2D.Double(0.,0.);
			
			// look at each segment
			int i;
			for(i=1;i<xa.size();i++)
			{	// end pt coordinates
				endA=new Point2D.Double(xa.get(i),ya.get(i));
				endB=new Point2D.Double(xb.get(i),yb.get(i));
				planePt1.setLocation((endA.x+endB.x)/2.,(endA.y+endB.y)/2.);
				ctod1.setLocation(endA.x-endB.x,endA.y-endB.y);
				
				// tangential vector is (dx,dy)
				tangential.setLocation(planePt1.x-planePt0.x,planePt1.y-planePt0.y);
				mag=Math.sqrt(tangential.x*tangential.x+tangential.y*tangential.y);
				tangential.setLocation(tangential.x/mag,tangential.y/mag);
				
				// normal is (-dy,dx)
				normal.setLocation(-tangential.y,tangential.x);
				
				// add point at each end of segment
				switch(component)
				{	case PlotQuantity.MPMNORMALCTOD:
						x.add(distance);
						y.add(ctod0.x*normal.x+ctod0.y*normal.y);
						distance+=PtSeparation2D(planePt1,planePt0);
						x.add(distance);
						y.add(ctod1.x*normal.x+ctod1.y*normal.y);
						break;
					
					case PlotQuantity.MPMSHEARCTOD:
						x.add(distance);
						y.add(ctod0.x*tangential.x+ctod0.y*tangential.y);
						distance+=PtSeparation2D(planePt1,planePt0);
						x.add(distance);
						y.add(ctod1.x*tangential.x+ctod1.y*tangential.y);
						break;
					
					case PlotQuantity.MPMOPENINGFRACTION:
					case PlotQuantity.MPMSHEARFRACTION:
						modes.setLocation(ctod0.x*normal.x+ctod0.y*normal.y,
												ctod0.x*tangential.x+ctod0.y*tangential.y);
						mag=Math.sqrt(modes.x*modes.x+modes.y*modes.y);
						modes.setLocation(modes.x/mag,modes.y/mag);
						x.add(distance);
						if(component==PlotQuantity.MPMOPENINGFRACTION)
							y.add(modes.x);
						else
							y.add(modes.y);
						distance+=PtSeparation2D(planePt1,planePt0);
						modes.setLocation(ctod1.x*normal.x+ctod1.y*normal.y,
												ctod1.x*tangential.x+ctod1.y*tangential.y);
						mag=Math.sqrt(modes.x*modes.x+modes.y*modes.y);
						modes.setLocation(modes.x/mag,modes.y/mag);
						x.add(distance);
						if(component==PlotQuantity.MPMOPENINGFRACTION)
							y.add(modes.x);
						else
							y.add(modes.y);
						break;
						
					case PlotQuantity.MPMCZMGI:
					case PlotQuantity.MPMCZMGII:
						// traction history
						x.add(distance);
						y.add(hpts.get(i));
						distance+=PtSeparation2D(planePt1,planePt0);
						break;
							
					default:
						// traction history
						x.add(distance);
						if(hpts!=null)
							y.add(hpts.get(i));
						else
							y.add(0.);
						distance+=PtSeparation2D(planePt1,planePt0);
						break;
				}
				
				// save end point values
				planePt0.setLocation(planePt1.x,planePt1.y);
				ctod0.setLocation(ctod1.x,ctod1.y);
			}
			
			// the plot
			Hashtable<String,String> props = new Hashtable<String,String>();
			props.put("object.color",plot2DView.selectPlotColor());
			props.put("array.name",PlotQuantity.plotName(component,resDoc)+" crack "+crackNum);
			plot2DView.plotData(x,y,props);
			
			if(plot2DView.getNumberOfPlots()<2)
			{	plot2DView.setXTitle("Distance ("+resDoc.units.lengthUnits()+")");
				if(component==PlotQuantity.MPMNORMALCTOD || component==PlotQuantity.MPMSHEARCTOD)
					plot2DView.setYTitle("Crack Opening ("+resDoc.units.lengthUnits()+")");
				else if(component>=PlotQuantity.MPMTRACTION1 && component<=PlotQuantity.MPMTRACTION10)
					plot2DView.setYTitle("Traction History");
				else if(component==PlotQuantity.MPMCZMGI || component==PlotQuantity.MPMCZMGI)
					plot2DView.setYTitle("Energy Released");
				else
					plot2DView.setYTitle("Mode Fraction");
			}
		}

		controls.setProgress(contourPoints);
	}
	
	// distance between two points
	public static double PtSeparation2D(Point2D.Double pt1,Point2D.Double pt2)
	{	double dx=pt2.x-pt1.x;
		double dy=pt2.y-pt1.y;
		return Math.sqrt(dx*dx+dy*dy);
	}
	
	// plot table of data
	public void plotPPTable(String plotTable) throws Exception
	{
		try
		{	plot2DView.readTable(plotTable);
			SwingUtilities.invokeLater(new Runnable()
			{	public void run()
				{	setVisible(true);
					toFront();
				}
			});
		}
		catch (Exception e)
		{	throw new Exception("Could not plot table of data:\n   " + e.getMessage());
		}
	}
}

// Store point location and value
// Used to sort by location only
class PointData implements Comparable<PointData>
{
	public double x;
	public double y;
	
	PointData(double h,double v)
	{	x = h;
		y = v;
	}
	
	public int compareTo(PointData pt)
	{
		if(x > pt.x)
			return 1;
		else if(x<pt.x)
			return -1;
		return 0;
	}
}