/*******************************************************************
	Interface.java
	NairnFEAMPMViz

	Created by John Nairn on 10/4/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.geom.*;
import java.util.*;

public class Interface extends ElementBase
{
	private boolean plottingTraction=false;
	
	// initialize
	Interface(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	// does not draw the mesh
	public void stroke(MeshPlotView pv,boolean displaced) {}
	
	// fill all subelements
	public void fill(MeshPlotView pv)
	{	if(subPaths.size()==0) return;
		
		// draw filled paths at subelements
		int i;
		for(i=0;i<subPaths.size();i++)
		{	pv.drawColor(subColors.get(i));
			pv.strokeShape(subPaths.get(i));
			pv.fillShape(subPaths.get(i));
		}
	}

	// no element path and not draw with mesh
	public void setElemPath(ResultsDocument doc,boolean displaced)
	{	if(!displaced) return;			// only make displaced path
		GeneralPath newPath=new GeneralPath();
		
		// first side
		NodalPoint nd=doc.nodes.get(node[0]);
		newPath.moveTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
		
		// interface element have nodes in order (and 0th node tacked on end)
		int i;
		for(i=1;i<=getNumberNodes();i++)
		{	nd=doc.nodes.get(node[i]);
			newPath.lineTo((float)(nd.x+nd.dispx),(float)(nd.y+nd.dispy));
		}
		
		displacedPath=newPath;
	}
	
	// copy nodal values to element values
	public void getPlotValuesFromNodes(ArrayList <NodalPoint>nodes)
	{	plottingTraction=false;
	}
	
	// get plot values for all nodes in this element
	public void getPlotValues(int comp,double stressAngle,ResultsDocument doc) throws Exception
	{	switch(comp)
		{	case PlotQuantity.INTERFACETRACTION_N:
			case PlotQuantity.INTERFACETRACTION_T:
				plottingTraction=true;
				super.getPlotValues(comp,stressAngle,doc);
				break;
			default:
				plottingTraction=false;
		}
	}
	
	// get plot values relevant for an interface element
	public double getNodeValue(int ndi,int comp,double stressAngle,ResultsDocument doc) throws Exception
	{	switch(comp)
		{	// --------------------------------------------
			//  Element Properties
			// --------------------------------------------
			case PlotQuantity.INTERFACETRACTION_N:
			case PlotQuantity.INTERFACETRACTION_T:
				if(ElementBase.DbleEqual(stressAngle,0.))
				{   switch(comp)
					{   case PlotQuantity.INTERFACETRACTION_N:
							return sigxx[ndi];
						case PlotQuantity.INTERFACETRACTION_T:
							return sigxy[ndi];
					}
				}
				else
				{   double radAngle=Math.PI*stressAngle/180.;
					double c=Math.cos(radAngle);
					double s=Math.sin(radAngle);
					switch(comp)
					{   case PlotQuantity.INTERFACETRACTION_N:
							return c*sigxx[ndi] - s*sigxy[ndi];
						case PlotQuantity.INTERFACETRACTION_T:
							return s*sigxx[ndi] + c*sigxy[ndi];
					}
				}
				break;
				
			default:
				break;
		}
		
		return 0.;
	}

	// find range of plot values
	public Point2D.Double plotLimits()
	{	if(!plottingTraction)
			return new Point2D.Double(1.,0.);
		else
			return super.plotLimits();
	}
	
	// This is called when preparing a mesh plot. It should store sub elements
	// if needed for plot.
	public void allocateSubElements(ArrayList<NodalPoint>nodes,int density,boolean displaced) throws Exception
	{
		// turn off plotting unless for tractions
		if(!plottingTraction)
		{	subPaths.clear();
			subValues.clear();
			subColors.clear();
			return;			// do not create new ones
		}

		// see if have correct number of subelements in right mesh
		if(subValues.size()==density && displaced==subPathsDisplaced) return;
		subPathsDisplaced=displaced;
				
		// get parameters
		NodalPoint[] ndpt=new NodalPoint[6];
		getNodalPoints(ndpt,nodes);
		double delta=2./(double)density;
		
		// clear all sub arrays
		int i,ii;
		subPaths.clear();
		subValues.clear();
		subColors.clear();
		for(i=0;i<density;i++)
		{	subValues.add(new Double(0.));
			subColors.add(Color.red);
		}
		
		Point2D.Double xiEta=new Point2D.Double(-1.,-1.);
		Point2D.Double[] pts=new Point2D.Double[4];
		for(i=0;i<density;i++)
		{	// first edge
			xiEta.y=-1.;
			pts[0]=findCoords(xiEta,ndpt,displaced);
			xiEta.x+=delta;
			pts[1]=findCoords(xiEta,ndpt,displaced);
			
			// second edge
			xiEta.y=1;
			pts[2]=findCoords(xiEta,ndpt,displaced);
			xiEta.x-=delta;
			pts[3]=findCoords(xiEta,ndpt,displaced);
				
			// make the path
			GeneralPath newPath=new GeneralPath();
			newPath.moveTo((float)pts[0].x,(float)pts[0].y);
			for(ii=1;ii<=3;ii++)
				newPath.lineTo((float)pts[ii].x,(float)pts[ii].y);
			newPath.lineTo((float)pts[0].x,(float)pts[0].y);
				
			// add to array
			subPaths.add(newPath);
				
			// on to next point
			xiEta.x+=delta;
		}
	}
	
	// create sub elements paths
	public void setPlotValues(int density) throws Exception
	{
		if(subValues.size()==0) return;			// can skip if not plotting

		// get parameters
		double delta=2./(double)density;

		// loop over sub elements
		int i;
		Point2D.Double xiEta=new Point2D.Double(-1.+0.5*delta,0.);
		for(i=0;i<density;i++)
		{   // get value at (xi)
			subValues.set(i,new Double(findValueAt(xiEta)));
			
			// next point
			xiEta.x+=delta;
		}
	}

	// find value for given dimensionless position
	// not isoparametric - thus use only half the nodes - only used to find traction
	// find value for given dimensionless position
	public double findValueAt(Point2D.Double xiEta) throws Exception
	{
		int i;
		double result=0.;
		double [] shape=new double[7];
		
		getShapeFunction(shape,xiEta,null);			//  this element does need nodal points
		result=0.;
		for(i=0;i<getNumberNodes()/2;i++)
			result+=shape[i]*plotValues[i];
		
		return result;
	}
	
	// xiEta.y<0 means first edge, xiEta.y>0 means second edge
	public Point2D.Double findCoords(Point2D.Double xiEta,NodalPoint[] ndpt,boolean displaced) throws Exception
	{	int i,i1=0,i2;
		double sign=1.;
		double[] shape=new double[7];
		
		getShapeFunction(shape,xiEta,ndpt);
		Point2D.Double pt=new Point2D.Double(0.,0.);
		if(xiEta.y<0)
			i2=getNumberNodes()/2;
		else
		{	i1=getNumberNodes()/2;
			i2=getNumberNodes();
			sign=-1.;
		}
		if(displaced)
		{	for(i=i1;i<i2;i++)
			{   pt.x+=sign*shape[i]*(ndpt[i].x+ndpt[i].dispx);
				pt.y+=sign*shape[i]*(ndpt[i].y+ndpt[i].dispy);
			}
		}
		else
		{	for(i=i1;i<i2;i++)
			{   pt.x+=sign*shape[i]*ndpt[i].x;
				pt.y+=sign*shape[i]*ndpt[i].y;
			}
		}
		return pt;
	}

	// accessors needing override
	public boolean PtInElement(Point2D.Double pt) { return false; }
	public int getNumberSides() { return 2; }
	public double getArea(ArrayList <NodalPoint>nodes) { return 0.; }
}
