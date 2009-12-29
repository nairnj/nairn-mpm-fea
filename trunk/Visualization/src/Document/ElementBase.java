/*******************************************************************
	ElementBase.java
	NairnFEAMPMViz

	Created by John Nairn on 2/27/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.geom.*;
import java.util.*;

public class ElementBase
{
	// constants
	public static final int CST=1;
	public static final int FOUR_NODE_ISO=2;
	public static final int EIGHT_NODE_ISO=3;
	public static final int ISO_TRIANGLE=4;
	public static final int LINEAR_INTERFACE=5;
	public static final int QUAD_INTERFACE=6;
	public static final int EIGHT_NODE_ISO_BRICK=7;
	public static final int MAXITER=100;
	
	// variables
	public int num;
	public int[] node;
	protected double[] plotValues;
	protected double angle,thickness;
	protected double[] forcex;
	protected double[] forcey;
	protected double[] sigxx;
	protected double[] sigyy;
	protected double[] sigxy;
	protected double[] sigzz;
	protected int material;
	public double energy;		// in FEA
	
	protected ArrayList<GeneralPath> subPaths;
	protected ArrayList<Double> subValues;
	protected boolean subPathsDisplaced=false;
	protected ArrayList<Color> subColors;
	protected GeneralPath elemPath;
	protected GeneralPath displacedPath;
	
	public static int density=4;
	public static float lineWidth=(float)0.5;
	
	//--------------------------------------------------------------
	// Create element
	//--------------------------------------------------------------

	// initialize
	ElementBase(int elemNum,int [] nds)
	{
		num=elemNum;
		node=new int[getNumberNodes()+1];
		int i;
		for(i=0;i<getNumberNodes();i++)
			node[i]=nds[i]-1;
		node[getNumberNodes()]=nds[0]-1;			// to help some algorithms
		
		// for mesh plots
		plotValues=new double[getNumberNodes()];
		subPaths=new ArrayList<GeneralPath>(36);
		subValues=new ArrayList<Double>(36);
		subColors=new ArrayList<Color>(36);
		forcex=new double[getNumberNodes()];
		forcey=new double[getNumberNodes()];
		sigxx=new double[getNumberNodes()];
		sigyy=new double[getNumberNodes()];
		sigxy=new double[getNumberNodes()];
		sigzz=new double[getNumberNodes()];
	}
	
	//--------------------------------------------------------------
	// Plot element - general methods
	//--------------------------------------------------------------

	// draw the element (assumes elements trace boundary in order, override if not)
	public void stroke(MeshPlotView pv,boolean displaced)
	{	if(displaced)
			pv.strokeShape(displacedPath);
		else
			pv.strokeShape(elemPath);
	}
	
	// make element path
	public void setElemPath(ResultsDocument doc,boolean displaced)
	{	int ns=getNumberSides();
		GeneralPath newPath=new GeneralPath();
		
		// move to first point
		NodalPoint nd=doc.nodes.get(node[0]);
		if(displaced)
			newPath.moveTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
		else
			newPath.moveTo((float)nd.x,(float)nd.y);
			
		// if needed, line to mid side node
		if(hasMidSideNodes())
		{	nd=doc.nodes.get(node[ns]);
			if(displaced)
				newPath.lineTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
			else
				newPath.lineTo((float)nd.x,(float)nd.y);
		}
		
		// rest of the nodes
		for(int i=1;i<getNumberSides();i++)
		{	// line to next corner
			nd=doc.nodes.get(node[i]);
			if(displaced)
				newPath.lineTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
			else
				newPath.lineTo((float)nd.x,(float)nd.y);
				
			// if needed, line to mid side node
			if(hasMidSideNodes())
			{	nd=doc.nodes.get(node[i+ns]);
				if(displaced)
					newPath.lineTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
				else
					newPath.lineTo((float)nd.x,(float)nd.y);
			}
		}
		
		// back to origin
		nd=doc.nodes.get(node[0]);
		if(displaced)
			newPath.lineTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
		else
			newPath.lineTo((float)nd.x,(float)nd.y);
			
		// set the path
		if(displaced)
			displacedPath=newPath;
		else
			elemPath=newPath;
	}
	
	// fill all subelements
	public void fill(MeshPlotView pv)
	{	if(subPaths.size()==0) return;
		
		// draw filled paths at subelements
		int i;
		for(i=0;i<subPaths.size();i++)
		{	pv.drawColor(subColors.get(i));
			pv.fillShape(subPaths.get(i));
		}
	}
	
	// draw the number
	public void number(MeshPlotView pv,ResultsDocument doc,boolean displaced)
	{	int i;
		double xavg=0.,yavg=0.;
		for(i=0;i<getNumberNodes();i++)
		{	NodalPoint nd=doc.nodes.get(node[i]);
			xavg+=nd.x;
			yavg+=nd.y;
			if(displaced)
			{	xavg+=nd.dispx;
				yavg+=nd.dispy;
			}
		}
		xavg/=(double)getNumberNodes();
		yavg/=(double)getNumberNodes();
		pv.moveTo(xavg,yavg);
		pv.drawString(String.format("%d",num),MeshPlotView.CENTER_LABEL+MeshPlotView.FRAME_LABEL);
	}
	
	//--------------------------------------------------------------
	// Base element methods
	//--------------------------------------------------------------

	// find nearest node number to a point in this element
	public int NearestNode(Point2D.Double pt,ResultsDocument doc,boolean displaced)
	{	int i,nearest,numnds=getNumberNodes();
		double distance,approach;
		
		NodalPoint nd=doc.nodes.get(node[0]);
		if(displaced)
			distance=pt.distanceSq(nd.x+nd.dispx,nd.y+nd.dispy);
		else
			distance=pt.distanceSq(nd.x,nd.y);
		nearest=node[0];
		for(i=1;i<numnds;i++)
		{	nd=doc.nodes.get(node[i]);
			if(displaced)
				approach=pt.distanceSq(nd.x+nd.dispx,nd.y+nd.dispy);
			else
				approach=pt.distanceSq(nd.x,nd.y);
			if(approach<distance)
			{   distance=approach;
				nearest=node[i];
			}
		}
		return nearest+1;
	}
	
	// find minumum side this element (cellMinSide<0 to get for this cell or >=0 for global compare)
	// actually finds minimum extent in either direction - works all element types
	public double getMinSide(double cellMinSide,ResultsDocument doc)
	{	int i;
		NodalPoint nd=doc.nodes.get(node[0]);
		double xmin=nd.x,xmax=nd.x,ymin=nd.y,ymax=nd.y;
		for(i=1;i<getNumberNodes();i++)
		{	nd=doc.nodes.get(node[i]);
			xmin=Math.min(nd.x,xmin);
			xmax=Math.max(nd.x,xmax);
			ymin=Math.min(nd.y,ymin);
			ymax=Math.max(nd.y,ymax);
		}
		if(cellMinSide<0)
			cellMinSide=xmax-xmin;
		else
			cellMinSide=Math.min(cellMinSide,xmax-xmin);
		return Math.min(cellMinSide,ymax-ymin);
	}
	
	//--------------------------------------------------------------
	// Shape Function and (xi,eta) Methods needing overrides
	//--------------------------------------------------------------

	// shape functions only
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] ndpt) throws Exception
	{	throw new Exception("Element without shape function code was accessed");
	}

	// general shape function BMatrix evaluation
	public void getShapeBMatrix(Point2D.Double xiEta,double[] xiDeriv,double[] etaDeriv,double[] asbe,NodalPoint[] eNodes,boolean derivs) throws Exception
	{	throw new Exception("Element without shape function BMatrix code was accessed");
	}
	
	// Get dimensionless coordinate of node (0 based) - must override in isoparameteric elements
	public Point2D.Double getNodeXiEta(int i,NodalPoint[] eNodes)
	{	return new Point2D.Double(eNodes[i].x,eNodes[i].y);
	}
	
	//--------------------------------------------------------------
	// Shape Function usually not needing overrides
	//--------------------------------------------------------------

	// find dimensionless location of a point numerically
	public void getXiEta(Point2D.Double xiEta,Point2D.Double pt,NodalPoint[] ndpt) throws Exception
	{
		double xt,yt,dxxi,dxeta,dyxi,dyeta;
		double deter,dxi,deta,dist;
		double[] gfn=new double[9];
		double[] gdfnxi=new double[9];
		double[] gdfnet=new double[9];
		int i,j;
		
		/* solve for xi and eta using Newton-Rapheson (see FEA Notes)
			using shape functions and their derivatives */
		for(j=1;j<=MAXITER;j++)
		{	getShapeFunction(gfn,xiEta,ndpt);
			getShapeBMatrix(xiEta,gdfnxi,gdfnet,null,ndpt,true);
			xt=-pt.x;
			yt=-pt.y;
			dxxi=0.;
			dxeta=0.;
			dyxi=0.;
			dyeta=0.;
			for(i=0;i<getNumberNodes();i++)
			{	xt+=ndpt[i].x*gfn[i];
				yt+=ndpt[i].y*gfn[i];
				dxxi+=ndpt[i].x*gdfnxi[i];
				dxeta+=ndpt[i].x*gdfnet[i];
				dyxi+=ndpt[i].y*gdfnxi[i];
				dyeta+=ndpt[i].y*gdfnet[i];
			}
			deter=dxxi*dyeta-dxeta*dyxi;
			dxi=(-xt*dyeta+yt*dxeta)/deter;
			deta=(xt*dyxi-yt*dxxi)/deter;
			xiEta.x+=dxi;
			xiEta.y+=deta;
			dist=Math.sqrt(dxi*dxi+deta*deta);
			if(dist<.001) break;
		}
		
		// final answer is in xiEta
	}
	
	// find coordinates for given dimensionless position
	public Point2D.Double findCoords(Point2D.Double xiEta,NodalPoint[] ndpt,boolean displaced) throws Exception
	{	int i;
		double[] shape=new double[9];
		
		getShapeFunction(shape,xiEta,ndpt);
		Point2D.Double pt=new Point2D.Double(0.,0.);
		if(displaced)
		{	for(i=0;i<getNumberNodes();i++)
			{   pt.x+=shape[i]*(ndpt[i].x+ndpt[i].dispx);
				pt.y+=shape[i]*(ndpt[i].y+ndpt[i].dispy);
			}
		}
		else
		{	for(i=0;i<getNumberNodes();i++)
			{   pt.x+=shape[i]*ndpt[i].x;
				pt.y+=shape[i]*ndpt[i].y;
			}
		}
		return pt;
	}

	// find value for given dimensionless position
	public double findValueAt(Point2D.Double xiEta) throws Exception
	{
		int i;
		double result=0.;
		double [] shape=new double[9];
		
		getShapeFunction(shape,xiEta,null);			// assumes this element does need nodal points
		for(i=0;i<getNumberNodes();i++)
			result+=shape[i]*plotValues[i];
		
		return result;
	}
	
	// find value at an actual position in pt (nodal coordinates in ndpt)
	public double findValueAtRealPt(double x,double y,NodalPoint[] ndpt) throws Exception
	{
		Point2D.Double xiEta=getCentroid();
		getXiEta(xiEta,new Point2D.Double(x,y),ndpt);
		return findValueAt(xiEta);
	}
	
	//--------------------------------------------------------------
	// Plotting Support
	//--------------------------------------------------------------

	// copy nodal values to element values
	public void getPlotValuesFromNodes(ArrayList <NodalPoint>nodes)
	{	int i;
		for(i=0;i<getNumberNodes();i++)
			plotValues[i]=(nodes.get(node[i])).getPlotValue();
	}
	
	// get plot values for all nodes in this element
	public void getPlotValues(int comp,double stressAngle,ResultsDocument doc) throws Exception
	{	int i;
		for(i=0;i<getNumberNodes();i++)
			plotValues[i]=getNodeValue(i,comp,stressAngle,doc);
	}

	// get plot values from element data at one node. The first properties are calculated
	//   from element data. The later properties ar calculated from nodal
	//	 properties - these are all for plotting FEA results
	public double getNodeValue(int ndi,int comp,double stressAngle,ResultsDocument doc) throws Exception
	{	switch(comp)
		{	// --------------------------------------------
			//  Element Properties
			// --------------------------------------------
			case PlotQuantity.MESHMATERIAL:
				return (double)material;
			
			case PlotQuantity.MESHANGLE:
				return (double)angle;
				
			case PlotQuantity.MESHSTRAINX:
			case PlotQuantity.MESHSTRAINY:
			case PlotQuantity.MESHSTRAINXY:
			case PlotQuantity.MESHSTRAINZ:
			case PlotQuantity.MESHDUDY:
			case PlotQuantity.MESHDVDX:
				// Load nodal points
				NodalPoint[] eNodes=new NodalPoint[9];
				getNodalPoints(eNodes,doc.nodes);
				
				// find nodal value
				double[] dfnxi=new double[9];
				double[] dfnet=new double[9];
				double[] asbe=new double[9];
				double nodeValue=0.;
				getShapeBMatrix(getNodeXiEta(ndi,eNodes),dfnxi,dfnet,asbe,eNodes,false);
				if(DbleEqual(stressAngle,0.) || comp==PlotQuantity.MESHDUDY || comp==PlotQuantity.MESHDVDX || comp==PlotQuantity.MESHSTRAINZ)
				{	for(int k=0;k<getNumberNodes();k++)
					{	switch(comp)
						{	case PlotQuantity.MESHSTRAINX:
								nodeValue+=100.*dfnxi[k]*eNodes[k].dispx;
								break;
							case PlotQuantity.MESHSTRAINY:
								nodeValue+=100.*dfnet[k]*eNodes[k].dispy;
								break;
							case PlotQuantity.MESHSTRAINXY:
								nodeValue+=100.*(dfnet[k]*eNodes[k].dispx + dfnxi[k]*eNodes[k].dispy);
								break;
							case PlotQuantity.MESHSTRAINZ:
								// STRAINZ is only used for axisym gives STRAIN TT
								if(doc.isAxisymmetric())
								{	if(eNodes[ndi].x!=0.)
										nodeValue+=100.*asbe[k]*eNodes[k].dispx;
									else
										nodeValue+=100.*dfnxi[k]*eNodes[k].dispx;
								}
								else
									nodeValue=0.;
								break;
							case PlotQuantity.MESHDUDY:
								nodeValue+=100.*dfnet[k]*eNodes[k].dispx;
								break;
							case PlotQuantity.MESHDVDX:
								nodeValue+=100.*dfnxi[k]*eNodes[k].dispy;
								break;
							default:
								break;
						}
					}
				}
				else
				{   double radAngle=Math.PI*stressAngle/180.;
					double c=Math.cos(radAngle);
					double s=Math.sin(radAngle);
					double strx=0.;
					double stry=0.;
					double strxy=0.;
					for(int k=0;k<getNumberNodes();k++)
					{	strx+=100.*dfnxi[k]*eNodes[k].dispx;
						stry+=100.*dfnet[k]*eNodes[k].dispy;
						strxy+=100.*(dfnet[k]*eNodes[k].dispx + dfnxi[k]*eNodes[k].dispy);
					}
					switch(comp)
					{	case PlotQuantity.MESHSTRAINX:
							nodeValue=c*c*strx + s*s*stry - c*s*strxy;
							break;
						case PlotQuantity.MESHSTRAINY:
							nodeValue=s*s*strx + c*c*stry + c*s*strxy;
							break;
						case PlotQuantity.MESHSTRAINXY:
							nodeValue=2*c*s*(strx-stry) + (c*c-s*s)*strxy;
							break;
						default:
							nodeValue=0.;
							break;
					}
				}
				return nodeValue;
				
			case PlotQuantity.MESHELEMSIGMAX:
			case PlotQuantity.MESHELEMSIGMAY:
			case PlotQuantity.MESHELEMSIGMAXY:
				if(DbleEqual(stressAngle,0.))
				{   switch(comp)
					{   case PlotQuantity.MESHELEMSIGMAX:
							return sigxx[ndi];
						case PlotQuantity.MESHELEMSIGMAY:
							return sigyy[ndi];
						case PlotQuantity.MESHELEMSIGMAXY:
							return sigxy[ndi];
					}
				}
				else
				{   double radAngle=Math.PI*stressAngle/180.;
					double c=Math.cos(radAngle);
					double s=Math.sin(radAngle);
					switch(comp)
					{   case PlotQuantity.MESHELEMSIGMAX:
							return c*c*sigxx[ndi] + s*s*sigyy[ndi] - 2*c*s*sigxy[ndi];
						case PlotQuantity.MESHELEMSIGMAY:
							return s*s*sigxx[ndi] + c*c*sigyy[ndi] + 2*c*s*sigxy[ndi];
						case PlotQuantity.MESHELEMSIGMAXY:
							return c*s*(sigxx[ndi]-sigyy[ndi]) + (c*c-s*s)*sigxy[ndi];
					}
				}
				break;
				
			case PlotQuantity.MESHELEMSIGMAZ:
				return sigzz[ndi];
				
			case PlotQuantity.MESHFORCEX:
			case PlotQuantity.MESHFORCEY:
				if(DbleEqual(stressAngle,0.))
				{   switch(comp)
					{   case PlotQuantity.MESHFORCEX:
							return forcex[ndi];
						case PlotQuantity.MESHFORCEY:
							return forcey[ndi];
					}
				}
				else
				{   double radAngle=Math.PI*stressAngle/180.;
					double c=Math.cos(radAngle);
					double s=Math.sin(radAngle);
					switch(comp)
					{   case PlotQuantity.MESHFORCEX:
							return c*forcex[ndi] - s*forcey[ndi];
						case PlotQuantity.MESHFORCEY:
							return s*forcex[ndi] + c*forcey[ndi];
					}
				}
				break;
				
			case PlotQuantity.MESHSTRAINENERGY:
				return 1.e9*energy/(getArea(doc.nodes)*thickness);
			
			// --------------------------------------------
			//  Nodal Properties
			// --------------------------------------------
			case PlotQuantity.MESHSIGMAX:
				return (doc.nodes.get(node[ndi])).sigxx;
			case PlotQuantity.MESHSIGMAY:
				return (doc.nodes.get(node[ndi])).sigyy;
			case PlotQuantity.MESHSIGMAXY:
				return (doc.nodes.get(node[ndi])).sigxy;
			case PlotQuantity.MESHSIGMAZ:
				return (doc.nodes.get(node[ndi])).sigzz;
				
			case PlotQuantity.MESHDISPX:
				return (doc.nodes.get(node[ndi])).dispx;
			case PlotQuantity.MESHDISPY:
				return (doc.nodes.get(node[ndi])).dispy;
			case PlotQuantity.MESHNODEX:
				return (doc.nodes.get(node[ndi])).x;
			case PlotQuantity.MESHNODEY:
				return (doc.nodes.get(node[ndi])).y;
			
			case PlotQuantity.MESHNODEDISTANCE:
			case PlotQuantity.MESHNODEANGLE:
				NodalPoint eNode=doc.nodes.get(node[ndi]);
				double xpt=eNode.x;
				double ypt=eNode.y;
				double distance=Math.sqrt(xpt*xpt+ypt*ypt);
				if(comp==PlotQuantity.MESHNODEDISTANCE)
					return distance;
				else
				{	if(distance==0.)
						return 0.;
					else if(xpt>=0.)
						return Math.asin(ypt/distance);
					else
						return Math.PI-Math.asin(ypt/distance);
				}
			
			default:
				break;
		}
		
		return 0.;
	}

	// This is called when preparing a mesh plot. It should store sub elements
	// if needed for plot.
	public void allocateSubElements(ArrayList<NodalPoint>nodes,int density,boolean displaced) throws Exception
	{
		// exit if already have correct subelements
		if(subValues.size()==density*density && displaced==subPathsDisplaced) return;
		subPathsDisplaced=displaced;

		// get parameters
		NodalPoint[] ndpt=new NodalPoint[9];
		getNodalPoints(ndpt,nodes);
		double delta=2./(double)density;
		
		// clear all sub arrays
		int i,j,ii,size=density*density;
		subPaths.clear();
		subValues.clear();
		subColors.clear();
		for(i=0;i<size;i++)
		{	subValues.add(new Double(0.));
			subColors.add(Color.red);
		}
		
		// loop over sub elements
		Point2D.Double xiEta=new Point2D.Double();
		Point2D.Double[] pts=new Point2D.Double[4];
		for(i=0;i<density;i++)
		{   // find left edge
			xiEta.x=-1.;
			xiEta.y=-1.+i*delta;
			pts[0]=findCoords(xiEta,ndpt,displaced);
			xiEta.y+=delta;
			pts[3]=findCoords(xiEta,ndpt,displaced);
			xiEta.y-=delta;
			for(j=0;j<density;j++)
			{   xiEta.x+=delta;
				pts[1]=findCoords(xiEta,ndpt,displaced);
				xiEta.y+=delta;
				pts[2]=findCoords(xiEta,ndpt,displaced);
				xiEta.y-=delta;
				
				// make the path
				GeneralPath newPath=new GeneralPath();
				newPath.moveTo((float)pts[0].x,(float)pts[0].y);
				for(ii=1;ii<=3;ii++)
					newPath.lineTo((float)pts[ii].x,(float)pts[ii].y);
				newPath.lineTo((float)pts[0].x,(float)pts[0].y);
				
				// add to array
				subPaths.add(newPath);
				
				// reuse right side points
				pts[3]=pts[2];
				pts[0]=pts[1];
			}
		}
	}
	
	// create sub elements paths
	public void setPlotValues(int density) throws Exception
	{
		if(subValues.size()==0) return;			// can skip if not plotting

		// get parameters
		double delta=2./(double)density;
		
		// loop over sub elements
		int i,j,vindex=0;
		Point2D.Double xiEta=new Point2D.Double();
		for(i=0;i<density;i++)
		{   // find left edge
			xiEta.x=-1.+0.5*delta;
			xiEta.y=-1.+(i+0.5)*delta;
			for(j=0;j<density;j++)
			{   // get value at (xi,eta)
				subValues.set(vindex,new Double(findValueAt(xiEta)));
				
				// next point
				xiEta.x+=delta;
				vindex++;
			}
		}
	}
	
	// find range of plot values (non-plotting element return dlim.y < dlim.x)
	public Point2D.Double plotLimits()
	{	Point2D.Double dlim=new Point2D.Double(plotValues[0],plotValues[0]);
		for(int i=1;i<getNumberNodes();i++)
		{   dlim.x=Math.min(dlim.x,plotValues[i]);
			dlim.y=Math.max(dlim.y,plotValues[i]);
		}
		return dlim;
	}
	
	// scale should be 1/(dmax-dmin)
	public void setColors(double dmin,double scale)
	{
		if(subColors.size()==0) return;
		
		int i;
		double plotValue;
		for(i=0;i<subValues.size();i++)
		{	plotValue=subValues.get(i).doubleValue();
			subColors.set(i,ColorPicker.PickRainbow(scale*(plotValue-dmin)));
		}
	}
	
	// get value for point in this element
	// used in MeshPlotData to display value are current mouse point
	public double getValueAt(Point2D.Double pt)
	{	double zValue=0.;
		if(subPaths.size()==0) return zValue;
		
		int i;
		for(i=0;i<subPaths.size();i++)
		{	if((subPaths.get(i)).contains(pt))
			{	zValue=(subValues.get(i)).doubleValue();
				break;
			}
		}
		return zValue;
	}
	
	//--------------------------------------------------------------
	// Accessors usually not needing overrides
	//--------------------------------------------------------------

	// overide in element type
	public boolean PtInElement(Point2D.Double pt) { return elemPath.contains(pt); }
	public boolean PtInDispElement(Point2D.Double pt) { return displacedPath.contains(pt); }
	
	// load nodal points in zero-based ndpt[]
	public void getNodalPoints(NodalPoint[] ndpt,ArrayList <NodalPoint>nodes)
	{	int i;
		for(i=0;i<getNumberNodes();i++)
		{	ndpt[i]=nodes.get(node[i]);
		}
	}

	// Get area of the element (general for polygons and polygons with mid side nodes)
	public double getArea(ArrayList <NodalPoint>nodes)
	{
		int i,ns=getNumberSides();
		double area=0.;
		NodalPoint pt1,pt2,pt3;
		
		if(!hasMidSideNodes())
		{	for(i=0;i<ns;i++)
			{	pt1=nodes.get(node[i]);
				pt2=nodes.get(node[i+1]);
				area+=pt1.x*pt2.y-pt1.y*pt2.x;
			}
		}
		else
		{	for(i=0;i<ns-1;i++)
			{	pt1=nodes.get(node[i]);
				pt2=nodes.get(node[i+1]);
				pt3=nodes.get(node[i+ns]);
				area+=pt1.x*pt3.y-pt1.y*pt3.x
							+pt3.x*pt2.y-pt3.y*pt2.x;
			}
			ns--;
			int nn=getNumberNodes()-1;
			pt1=nodes.get(node[ns]);
			pt2=nodes.get(node[nn]);
			pt3=nodes.get(node[0]);
			area+=pt1.x*pt2.y-pt1.y*pt2.x
						+pt2.x*pt3.y-pt2.y*pt3.x;
		}
		return area/2.;
	}
	
	// coordinates of part of the face
	public void getFaceInfo(int face,int half,Point2D.Double origin,Point2D.Double slope,ArrayList<NodalPoint> nodes)
	{
		// got three points
		NodalPoint startPt=nodes.get(node[face]);
		int endNum = face==getNumberNodes()-1 ? 0 : face+1;
		NodalPoint endPt=nodes.get(node[endNum]);
		Point2D.Double midPt=new Point2D.Double(0.,0.);
		if(hasMidSideNodes())
		{	NodalPoint mid=nodes.get(node[face+getNumberSides()]);
			midPt.x=mid.x;
			midPt.y=mid.y;
		}
		else
		{	midPt.x=(startPt.x+endPt.x)/2.;
			midPt.y=(startPt.y+endPt.y)/2.;
		}
		
		// return coordinate for first or second half of the element
		if(half==1)
		{	origin.x=startPt.x;
			origin.y=startPt.y;
			slope.x=midPt.x-startPt.x;
			slope.y=midPt.y-startPt.y;
		}
		else
		{	origin.x=midPt.x;
			origin.y=midPt.y;
			slope.x=endPt.x-midPt.x;
			slope.y=endPt.y-midPt.y;
		}
	}
	
	// FEA elements need material, angle, and thickness
	public void setFEAProperties(int mat,double ang,double thick)
	{	material=mat;
		angle=ang;
		thickness=thick;
	}
	
	// FEA forces
	public void setForces(double fx,double fy,int i)
	{	forcex[i]=fx;
		forcey[i]=fy;
	}
	
	// FEA Element stresses
	public void setXYStresses(double sxx,double syy,double sxy,int i)
	{	sigxx[i]=sxx;
		sigyy[i]=syy;
		sigxy[i]=sxy;
	}
	
	// FEA Element stresses
	public void set3DStresses(double szz,double sxz,double syz,int i)
	{	sigzz[i]=szz;
	}
	
	//--------------------------------------------------------------
	// Accessors probably needing overrides
	//--------------------------------------------------------------

	// nodes, sides, and mid side nodes (override as needed)
	public int getNumberNodes() { return 4; }
	public int getNumberSides() { return 4; }
	public boolean hasMidSideNodes() { return false; }
	
	// centroid for this element and coordinates of the nodes
	public Point2D.Double getCentroid() { return new Point2D.Double(0.,0.); }

	//--------------------------------------------------------------
	// Static Class Methods
	//--------------------------------------------------------------

	// find minimum cell side all elements
	public static double getCellMinSide(ResultsDocument doc)
	{	int i;
		double cellMinSide=-1.;
		for(i=0;i<doc.elements.size();i++)
		{	cellMinSide=(doc.elements.get(i)).getMinSide(cellMinSide,doc);
		}
		return cellMinSide;
	}

	// load all material points with values for new plot component
	public static void loadPlotData(int component,ResultsDocument doc,int plotType,boolean inDisplaced) throws Exception
	{
		int i;
		double dmin=0.,dmax=1.;
		
		// exit if for FEA mesh only - no need to get data
		if(component==PlotQuantity.MESHONLY) return;
		
		// load plotValue of nodal points for MPM plots and transfer to element plotValues
		if(plotType==MeshPlotView.MPMMESH_PLOTS)
		{	// extrapolate to the grid
			mpmExtrapolation(component,0.0d,doc);
				
			// find range of the grid values
			NodalPoint nptr=doc.nodes.get(0);
			dmin=nptr.getPlotValue();
			dmax=dmin;
			for(i=1;i<doc.nodes.size();i++)
			{	nptr=doc.nodes.get(i);
				dmin=Math.min(dmin,nptr.getPlotValue());
				dmax=Math.max(dmax,nptr.getPlotValue());
			}
			
			// position/material will expand range slightly
			if(component==PlotQuantity.MPMPOS)
			{	dmin-=0.5;
				dmax+=0.5;
			}
			
			// transfer nodal values to elements
			for(i=0;i<doc.elements.size();i++)
				(doc.elements.get(i)).getPlotValuesFromNodes(doc.nodes);
		}
		
		// load plot values into element plotValues for FEA plots
		else
		{	// loop over elements
			for(i=0;i<doc.elements.size();i++)
				(doc.elements.get(i)).getPlotValues(component,0.0,doc);
				
			// find range from element values
			i=0;
			Point2D.Double dextent;
			
			// find first valid element
			while(i<doc.elements.size())
			{	dextent=(doc.elements.get(i)).plotLimits();
				i++;
				if(dextent.x>dextent.y) continue;
				dmin=dextent.x;
				dmax=dextent.y;
				break;
			}
			
			// search rest of elements
			while(i<doc.elements.size())
			{   dextent=(doc.elements.get(i)).plotLimits();
				i++;
				if(dextent.x>dextent.y) continue;
				dmin=Math.min(dmin,dextent.x);
				dmax=Math.max(dmax,dextent.y);
			}
			
			// material will expand range slightly
			if(component==PlotQuantity.MESHMATERIAL)
			{	dmin-=0.5;
				dmax+=0.5;
			}
		}
			
		// allocate subelements in elements for plotting
		for(i=0;i<doc.elements.size();i++)
			(doc.elements.get(i)).allocateSubElements(doc.nodes,density,inDisplaced);
		
		// calculate element plot values for the subelements based on element plot values
		for(i=0;i<doc.elements.size();i++)
			(doc.elements.get(i)).setPlotValues(density);
		
		// changed to fixed limits
		Point2D.Double limits = doc.docCtrl.controls.adjustLimits(dmin, dmax);
		dmin=limits.x;
		dmax=limits.y;
		
		// save global range and set the spectum
		ElementBase.setSpectrum(dmin,dmax,doc);
		
		// tell plot view its range
		doc.docCtrl.movieFrame.plotView.dataMin=dmin;
		doc.docCtrl.movieFrame.plotView.dataMax=dmax;
		doc.docCtrl.movieFrame.plotView.dataLimitsSet=true;
	}
	
	// get mesh data in preparatoin for 2D mesh plots
	public static void load2DPlotData(int component,ResultsDocument doc) throws Exception
	{
		int i;
		
		// load plotValue of nodal points for MPM plots and transfer to element plotValues
		if(doc.isMPMAnalysis())
		{	// extrapolate to the grid
			mpmExtrapolation(component,0.0d,doc);
			
			// transfer nodal values to elements
			for(i=0;i<doc.elements.size();i++)
				(doc.elements.get(i)).getPlotValuesFromNodes(doc.nodes);
		}
		
		// load plot values into element plotValues for FEA plots
		else
		{	// loop over elements
			for(i=0;i<doc.elements.size();i++)
				(doc.elements.get(i)).getPlotValues(component,0.0d,doc);
		}
	}
	
	// extrapolate MPM data to the grid plotValue
	public static void mpmExtrapolation(int component,double angle,ResultsDocument resDoc) throws Exception
	{
		int i,j;
		
		// zero nodal properties
		for(i=0;i<resDoc.nodes.size();i++)
		{	(resDoc.nodes.get(i)).zeroPlotValue();
		}
		
		// loop over particles
		MaterialPoint mptr;
		ElementBase eptr;
		NodalPoint[] ndpt=new NodalPoint[9];
		double plotValue,wt=1.;
		double[] fn=new double[9];
		for(i=0;i<resDoc.mpmPoints.size();i++)
		{	mptr=resDoc.mpmPoints.get(i);
		
			// mesh extrapolation should not include rigid particles
			MaterialBase matl=resDoc.materials.get(mptr.materialIndex());
			if(matl.isRigid()) continue;
			
			// get value to extrapolate
			plotValue=mptr.getForPlot(component,angle,resDoc);
			if(component!=PlotQuantity.MPMMASS) wt=mptr.mass;
			
			// Load element object
			eptr=resDoc.elements.get(mptr.inElem);
			
			// find dimensionsless location and shape functions
			eptr.getNodalPoints(ndpt,resDoc.nodes);
			Point2D.Double xiEta=eptr.getCentroid();
			eptr.getXiEta(xiEta,mptr.getPosition(),ndpt);
			eptr.getShapeFunction(fn,xiEta,ndpt);
			
			// Add particle property to each node in the element
			for(j=0;j<eptr.getNumberNodes();j++)
				ndpt[j].addPlotValue(plotValue,wt*fn[j]);
		}
		
		// normalize nodal properties
		if(component!=PlotQuantity.MPMMASS)
		{	for(i=0;i<resDoc.nodes.size();i++)
				(resDoc.nodes.get(i)).normPlotValue();
		}
	}
	
	// set color spectrum current range by changing colors
	//  of all element subelements.
	public static void setSpectrum(double dmin,double dmax,ResultsDocument doc)
	{
		// adjust for no range
		if(Math.abs(dmax-dmin)<1.e-15)
		{	dmax+=1.;
			dmin-=1.;
		}
		
		// Elements
		int i;
		double scale=1./(dmax-dmin);
		for(i=0;i<doc.elements.size();i++)
			(doc.elements.get(i)).setColors(dmin,scale);
	}
	
	// setting for change<small num can have an effect of some calculations
	public static boolean DbleEqual(double db1,double db2)
	{
		double ab1=Math.abs(db1);
		double ab2=Math.abs(db2);
		double change;
		
		if(db1==db2)
			return true;
		else if(ab1>ab2)
			change=Math.abs(db1-db2)/ab1;
		else
			change=Math.abs(db1-db2)/ab2;
				
		// Equal if different by less than 1 ppm
		if(change<1.e-6)
			return true;
		else
		{	if(ab1<1e-12 && ab2<1e-12)
				return true;
			else
				return false;
		}
	}
	
	// class method to get number of nodes before element object is created
	public static int NodesFromType(int elemID)
	{	switch(elemID)
		{	case CST:
				return 3;
			case FOUR_NODE_ISO:
			case LINEAR_INTERFACE:
				return 4;
			case EIGHT_NODE_ISO:
			case EIGHT_NODE_ISO_BRICK:
				return 8;
			case ISO_TRIANGLE:
			case QUAD_INTERFACE:
				return 6;
			default:
				break;
		}
		return 4;
	}
	
	// find element containing a point, check previous one first, then search all, return null if not found
	public static ElementBase findPtElement(double ptX,double ptY,ElementBase original,ResultsDocument doc)
	{
		Point2D.Double pt=new Point2D.Double(ptX,ptY);
		
		// is it in the same element
		if(original!=null)
		{	if(original.PtInElement(pt))
				return original;
		}
		
		// search for the element
		int i;
		ElementBase newElem;
		for(i=0;i<doc.elements.size();i++)
		{	newElem=doc.elements.get(i);
			if(newElem.PtInElement(pt))
				return newElem;
		}
		
		// none
		return null;
	}

}
