/*******************************************************************
	XYPlotWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 1/22/08.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import javax.swing.*;
import java.util.*;
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
	private MathEvaluator contourExpr=null;
	private double contourRange;
	private int pmpts;
	private int functionType;
	
	// for plot thread calculations
	private int component;
	private ControlPanel controls;

	// initialize
	public XYPlotWindow(ResultsDocument gResDoc)
	{	super(gResDoc);
	}
	
	// add another plot
	public void addPlot(ControlPanel plotControls) throws Exception
	{
		// the component to plot
		controls=plotControls;
		component=controls.adjustComponent(controls.getPlotComponent());
		
		// detach thread to gather plot information
		Thread plot2DThread=new Thread(this);
		plot2DThread.start();
	}
	
	//----------------------------------------------------------------------------
	// detachable thread for loading time plot data
	//----------------------------------------------------------------------------

	public void run()
	{
		controls.enableProgress(contourPoints);
		
		String contourFxn=null;
		
		try
		{	switch(component)
			{	case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
				case PlotQuantity.MPMCRACKPROFILE:
				case PlotQuantity.MPMOPENINGFRACTION:
				case PlotQuantity.MPMSHEARFRACTION:
					functionType=CRACK_CONTOUR;
					break;
				
				default:
					functionType=controls.getContour();
					contourFxn=controls.getContourFunction();
					
					if(contourFxn.length()==0)
						throw new Exception("The contour function cannot be empty");
					
					switch(functionType)
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
					contourExpr=new MathEvaluator(contourFxn);
					if(evaluate(2.d)==null)
						throw new Exception("The contour function does not evaluate to a number");
						
					// get range
					contourRange=controls.getPlusMinus();
					if(ElementBase.DbleEqual(contourRange,0.))
						pmpts=1;
					else
						pmpts=AVGQUAD;
					break;
			}
			
			// load the archive
			if(resDoc.isMPMAnalysis())
			{	resDoc.readSelectedArchive(controls.getArchiveIndex());
			}
			
			// load element plot values
			if(functionType!=CRACK_CONTOUR)
			{	ElementBase.load2DPlotData(component,resDoc);
				plotXYResults();
			}
			else
				plotCrackResults();
		}
		catch(Exception pe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(this,"X-Y plot error: "+pe.getMessage());
			if(plot2DView.getNumberOfPlots()==0) dispose();
			controls.disableProgress();
			return;
		}
		
		// finish up
		controls.disableProgress();
		setVisible(true);
		toFront();
	}
	
	// evaluate expression at tval
	public Double evaluate(double tval)
	{	contourExpr.reset();
		contourExpr.addVariable(var1,tval);
		if(var2!=null) contourExpr.addVariable(var2,tval);
		return contourExpr.getValue();
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
		Double nextVal=evaluate(prevX);
		if(nextVal==null)
			throw new Exception("The contour function does not evaluate to a number at "+prevX);
		prevY=nextVal.doubleValue();
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
			nextVal=evaluate(mptX);
			if(nextVal==null)
				throw new Exception("The contour function does not evaluate to a number at "+mptX);
			mptY=nextVal.doubleValue();
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
			props.put("array.name",PlotQuantity.plotName(component));
			plot2DView.plotData(x,y,props);
				
			// axis labels
			if(plot2DView.getNumberOfPlots()<2)
			{	plot2DView.setXTitle("Position (mm)");
				plot2DView.setYTitle(PlotQuantity.plotLabel(component,"mm","ms"));
			}
		}
	}
	
	// entry point to get non-crack, x-y plot points
	public void plotCrackResults() throws Exception
	{
		// does crack exist?
		int crackNum=controls.getCrackNumber();
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
			{	plot2DView.setXTitle("Position (mm)");
				plot2DView.setYTitle("Position (mm)");
			}
		}
		
		// normal and shear COD or opening and sliding fraction
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
					
					default:
						break;
				}
				
				// save end point values
				planePt0.setLocation(planePt1.x,planePt1.y);
				ctod0.setLocation(ctod1.x,ctod1.y);
			}
			
			// the plot
			Hashtable<String,String> props = new Hashtable<String,String>();
			props.put("object.color",plot2DView.selectPlotColor());
			props.put("array.name",PlotQuantity.plotName(component)+" "+crackNum);
			plot2DView.plotData(x,y,props);
			
			if(plot2DView.getNumberOfPlots()<2)
			{	plot2DView.setXTitle("Distance (mm)");
				if(component==PlotQuantity.MPMNORMALCTOD || component==PlotQuantity.MPMSHEARCTOD)
					plot2DView.setYTitle("Crack Opening (mm)");
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
}

