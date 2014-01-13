/*******************************************************************
	MaterialPoint.java
	NairnFEAMPMViz

	Created by John Nairn on 2/28/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.nio.*;
import java.awt.*;
import java.awt.geom.*;


public class MaterialPoint
{
	//---------------------------------------------------------------------
	// variables and constants
	//---------------------------------------------------------------------
	
	public static final int XXID=0;
	public static final int YYID=1;
	public static final int ZZID=2;
	public static final int XYID=3;
	
	public int num;
	public int inElem;			// zero based element number
	public double mass,angleZ,thickness;
	public int material;
	public double x,y,origx,origy,velx,vely;
	public double deltaTemp,plastEnergy,dudy,dvdx,strainEnergy;
	public double history1,history2,history3,history4,concentration,dcdx,dcdy;
	public double workEnergy,heatEnergy;
	public int elementCrossings;
	public double[] sigma;
	public double[] eps;
	public double[] eplast;
	public double erot;
	
	private double plotValue;
	public Color plotColor;
	
	//---------------------------------------------------------------------
	// initialize
	//---------------------------------------------------------------------
	
	MaterialPoint(int ptNum)
	{	num=ptNum;
		plotValue=0.;
		sigma=new double[4];
		eps=new double[4];
		eplast=new double[4];
	}
	
	//---------------------------------------------------------------------
	// draw
	//---------------------------------------------------------------------
	
	// draw the material point
	public void stroke(MeshPlotView pv,ResultsDocument doc)
	{	pv.moveTo(x,y);
		pv.drawMaterialPoint(plotColor,this);
	}
	
	// draw the material point
	public void addToClip(MeshPlotView pv,ResultsDocument doc,GeneralPath theClip)
	{	pv.moveTo(x,y);
		pv.clipMaterialPoint(this,theClip);
	}
	
	// return shape for the particle at current location
	// (xpt,ypt) are pixel locations and (radiix,radiiy) are pixel sizes
	public Shape particleShape(ResultsDocument resDoc,double xpt,double ypt,double radiix,double radiiy,
								boolean showSquarePts,boolean transformPts)
	{	// depends on setting
		if(showSquarePts)
		{	if(transformPts)
			{	double dgrad00,dgrad01,dgrad10,dgrad11;
				if(resDoc.materials.get(materialIndex()).hasPlasticStrainForGradient())
				{	dgrad00=0.01*(eps[MaterialPoint.XXID]+eplast[MaterialPoint.XXID]);
					double wrot=Math.PI*(angleZ-erot)/180.;
					dgrad01=0.005*(eps[MaterialPoint.XYID]+eplast[MaterialPoint.XYID])+wrot;
					dgrad10=0.005*(eps[MaterialPoint.XYID]+eplast[MaterialPoint.XYID])-wrot;
					dgrad11=0.01*(eps[MaterialPoint.YYID]+eplast[MaterialPoint.YYID]);
				}
				else
				{	dgrad00=0.01*eps[MaterialPoint.XXID];
					double wrot=Math.PI*(angleZ-erot)/180.;
					dgrad01=0.005*eps[MaterialPoint.XYID]+wrot;
					dgrad10=0.005*eps[MaterialPoint.XYID]-wrot;
					dgrad11=0.01*eps[MaterialPoint.YYID];
				}
				
				// This works in Java 1.5
				GeneralPath quad=new GeneralPath();
				Point2D.Float pathPt0=new Point2D.Float((float)(xpt+radiix*(1.+dgrad00)-radiiy*dgrad10),
											(float)(ypt+radiix*dgrad01-radiiy*(1.+dgrad11)));
				quad.moveTo(pathPt0.x,pathPt0.y);
				quad.lineTo((float)(xpt+radiix*(1.+dgrad00)+radiiy*dgrad10),
											(float)(ypt+radiix*dgrad01+radiiy*(1.+dgrad11)));
				quad.lineTo((float)(xpt-radiix*(1.+dgrad00)+radiiy*dgrad10),
											(float)(ypt-radiix*dgrad01+radiiy*(1.+dgrad11)));
				quad.lineTo((float)(xpt-radiix*(1.+dgrad00)-radiiy*dgrad10),
											(float)(ypt-radiix*dgrad01-radiiy*(1.+dgrad11)));
				quad.lineTo(pathPt0.x,pathPt0.y);
				/*
				// This requires Java 1.6
				Path2D.Double quad=new Path2D.Double();
				Point2D.Double pathPt0=new Point2D.Double((double)(xpt+radiix*(1.+dgrad00)-radiiy*dgrad10),
											(double)(ypt+radiix*dgrad01-radiiy*(1.+dgrad11)));
				quad.moveTo(pathPt0.x,pathPt0.y);
				quad.lineTo((double)(xpt+radiix*(1.+dgrad00)+radiiy*dgrad10),
											(double)(ypt+radiix*dgrad01+radiiy*(1.+dgrad11)));
				quad.lineTo((double)(xpt-radiix*(1.+dgrad00)+radiiy*dgrad10),
											(double)(ypt-radiix*dgrad01+radiiy*(1.+dgrad11)));
				quad.lineTo((double)(xpt-radiix*(1.+dgrad00)-radiiy*dgrad10),
											(double)(ypt-radiix*dgrad01-radiiy*(1.+dgrad11)));
				quad.lineTo(pathPt0.x,pathPt0.y);
				*/
				return quad;
			}
			else
				return new Rectangle2D.Double(xpt-radiix,ypt-radiiy,2.*radiix,2.*radiiy);
		}
		else
			return new Ellipse2D.Double(xpt-radiix,ypt-radiiy,2.*radiix,2.*radiiy);

	}
	
	// draw the number
	public void number(MeshPlotView pv,ResultsDocument doc)
	{
		pv.moveTo(x,y);
		pv.drawString(String.format("%d",num),
				MeshPlotView.CENTER_LABEL+MeshPlotView.FRAME_LABEL+MeshPlotView.FILL_FRAME);
	}
	
	//---------------------------------------------------------------------
	// read record from and archive file into this material point
	//---------------------------------------------------------------------
	
	public void readRecord(ByteBuffer bb,char[] mpmOrder,double lengthScale,double timeScale,boolean has3D)
	{	// required elements
		inElem=bb.getInt()-1;					// in element number (zero based)
		mass=bb.getDouble();					// mass in g
		material=bb.getShort();					// material number
		bb.getShort();							// skip 2 bytes
		
		// 3D will have 3 angles, but not read in this app yet
		if(has3D)
		{	angleZ=bb.getDouble();					// angle in degrees about Z axis
			bb.getDouble();
			bb.getDouble();
			thickness=1.;
		}
		else
		{	angleZ=bb.getDouble();					// angle in degrees about Z axis
			thickness=bb.getDouble()*lengthScale;	// thickness in length units
		}
		
		x=bb.getDouble()*lengthScale;			// x position
		y=bb.getDouble()*lengthScale;			// y position
		// 3D will have third position\
		if(has3D) bb.getDouble();
		
		// new default puts original position (in length units) here
		if(mpmOrder[ReadArchive.ARCH_Defaults]=='Y')
		{	origx=bb.getDouble()*lengthScale;
			origy=bb.getDouble()*lengthScale;
			// 3D will have third position
			if(has3D) bb.getDouble();
		}

		// velocity (length units/time units), but archive has mm/sec so need to
		// divide by 1000 to get default units here which are ms
		if(mpmOrder[ReadArchive.ARCH_Velocity]=='Y')
		{	double rescale=lengthScale/(1000.*timeScale);
			velx=bb.getDouble()*rescale;
			vely=bb.getDouble()*rescale;
			// 3D will have third velocity
			if(has3D) bb.getDouble();
		}
		else
		{	velx=0.;
			vely=0.;
			// 3D will have third velocity
		}

		// stress (in MPa) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_Stress]=='Y')
		{	sigma[XXID]=1.e-6*bb.getDouble();
			sigma[YYID]=1.e-6*bb.getDouble();
			sigma[ZZID]=1.e-6*bb.getDouble();
			sigma[XYID]=1.e-6*bb.getDouble();
			// 3D will have two more components
			if(has3D)
			{	bb.getDouble();
				bb.getDouble();
			}
		}
		else
		{	sigma[XXID]=0.;
			sigma[YYID]=0.;
			sigma[ZZID]=0.;
			sigma[XYID]=0.;
			// 3D will have two more components
		}
		
		// strain (in %)  (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_Strain]=='Y')
		{	eps[XXID]=100.*bb.getDouble();
			eps[YYID]=100.*bb.getDouble();
			eps[ZZID]=100.*bb.getDouble();
			eps[XYID]=100.*bb.getDouble();
			// 3D will have two more components
			if(has3D)
			{	bb.getDouble();
				bb.getDouble();
			}
		}
		else
		{	eps[XXID]=0.;
			eps[YYID]=0.;
			eps[ZZID]=0.;
			eps[XYID]=0.;
			// 3D will have two more components
		}
		
		// plastic strain (in %) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_PlasticStrain]=='Y')
		{	eplast[XXID]=100.*bb.getDouble();
			eplast[YYID]=100.*bb.getDouble();
			eplast[ZZID]=100.*bb.getDouble();
			eplast[XYID]=100.*bb.getDouble();
			// 3D will have two more components
			if(has3D)
			{	bb.getDouble();
				bb.getDouble();
			}
		}
		else
		{	eplast[XXID]=0.;
			eplast[YYID]=0.;
			eplast[ZZID]=0.;
			eplast[XYID]=0.;
			// 3D will have two more components
		}
				
		// old method for original positions (in length units)
		if(mpmOrder[ReadArchive.ARCH_OldOrigPosition]=='Y')
		{	origx=bb.getDouble()*lengthScale;
			origy=bb.getDouble()*lengthScale;
			// 3D will have one more position
			if(has3D) bb.getDouble();
		}
		
		// external work (cumulative) in J (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_WorkEnergy]=='Y')
			workEnergy=bb.getDouble();
		else
			workEnergy=0.;
		
		// temperature (C) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_DeltaTemp]=='Y')
			deltaTemp=bb.getDouble();
		else
			deltaTemp=0.;
				
		// total plastic energy (Volume*energy) in J (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_PlasticEnergy]=='Y')
			plastEnergy=bb.getDouble();
		else
			plastEnergy=0.;
		
		// shear components (dimensionless) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_ShearComponents]=='Y')
		{	dudy=bb.getDouble();
			dvdx=bb.getDouble();
		}
		else
		{	dudy=0.;
			dvdx=0.;
		}
			
		// total strain energy (Volume*energy) in J (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_StrainEnergy]=='Y')
			strainEnergy=bb.getDouble();
		else
			strainEnergy=0.;
			
		// particle history (variable units) (not converted for units)
		history1=0.;
		history2=0.;
		history3=0.;
		history4=0.;
		if(mpmOrder[ReadArchive.ARCH_History]=='Y')
			history1=bb.getDouble();
		else if(mpmOrder[ReadArchive.ARCH_History]!='N')
		{	int history=(int)mpmOrder[ReadArchive.ARCH_History];
			if((history & 0x01) !=0) history1=bb.getDouble();
			if((history & 0x02) !=0) history2=bb.getDouble();
			if((history & 0x04) !=0) history3=bb.getDouble();
			if((history & 0x08) !=0) history4=bb.getDouble();
		}
			
		// concentration and gradients
		if(mpmOrder[ReadArchive.ARCH_Concentration]=='Y')
		{	concentration=bb.getDouble();
			dcdx=bb.getDouble()/lengthScale;
			dcdy=bb.getDouble()/lengthScale;
			// 3D will have one more
			if(has3D) bb.getDouble();
		}
		else
		{	concentration=0.;
			dcdx=0.;
			dcdy=0.;
			// 3D will have one more
		}
		
		// thermal energy (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_HeatEnergy]=='Y')
			heatEnergy=bb.getDouble();
		else
			heatEnergy=0.;
			
		// element crossings (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_ElementCrossings]=='Y')
			elementCrossings=bb.getInt();
		else
			elementCrossings=0;
			
		// rot strain (in degrees) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_RotStrain]=='Y')
			erot=bb.getDouble();
		else
			erot=angleZ;	// cause rotational strain to be zero since not known
	}
	
	//---------------------------------------------------------------------
	// To plot material point property, this method must set the plotValue
	// variable according to the plot component. The angle variable
	// will rotate the component if it makes sense
	//---------------------------------------------------------------------
	
	public void loadForPlot(int component,double angle,ResultsDocument doc)
	{	plotValue=getForPlot(component,angle,doc);
	}
	
	public double getForPlot(int component,double angle,ResultsDocument doc)
	{
		double radAngle,sigx,sigy,sigxy,c,s;
		//unsigned hist;
		double theValue=0.;
		
		switch(component)
		{	// Stresses
			case PlotQuantity.MPMSIGMAX:
			case PlotQuantity.MPMSIGMAY:
			case PlotQuantity.MPMSIGMAXY:
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMSIGMAX:
							theValue=sigma[XXID];
							break;
						case PlotQuantity.MPMSIGMAY:
							theValue=sigma[YYID];
							break;
						case PlotQuantity.MPMSIGMAXY:
							theValue=sigma[XYID];
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=sigma[XXID];
					sigy=sigma[YYID];
					sigxy=sigma[XYID];
					switch(component)
					{   case PlotQuantity.MPMSIGMAX:
							theValue=c*c*sigx + s*s*sigy - 2*c*s*sigxy;
							break;
						case PlotQuantity.MPMSIGMAY:
							theValue=s*s*sigx + c*c*sigy + 2*c*s*sigxy;
							break;
						case PlotQuantity.MPMSIGMAXY:
							theValue=c*s*(sigx-sigy) + (c*c-s*s)*sigxy;
							break;
					}
				}
				break;
				
			case PlotQuantity.MPMSIGMAZ:
				theValue=sigma[ZZID];
				break;
		
			// Strains
			case PlotQuantity.MPMEPSX:
			case PlotQuantity.MPMEPSY:
			case PlotQuantity.MPMEPSXY:
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMEPSX:
							theValue=eps[XXID];
							break;
						case PlotQuantity.MPMEPSY:
							theValue=eps[YYID];
							break;
						case PlotQuantity.MPMEPSXY:
							theValue=eps[XYID];
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=eps[XXID];
					sigy=eps[YYID];
					sigxy=eps[XYID];
					switch(component)
					{   case PlotQuantity.MPMEPSX:
							theValue=c*c*sigx + s*s*sigy - c*s*sigxy;
							break;
						case PlotQuantity.MPMEPSY:
							theValue=s*s*sigx + c*c*sigy + c*s*sigxy;
							break;
						case PlotQuantity.MPMEPSXY:
							theValue=2*c*s*(sigx-sigy) + (c*c-s*s)*sigxy;
							break;
					}
				}
				break;
				
			case PlotQuantity.MPMEPSZ:
				theValue=eps[ZZID];
				break;
		
			// Plastic Strains
			case PlotQuantity.MPMPLEPSX:
			case PlotQuantity.MPMPLEPSY:
			case PlotQuantity.MPMPLEPSXY:
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMPLEPSX:
							theValue=eplast[XXID];
							break;
						case PlotQuantity.MPMPLEPSY:
							theValue=eplast[YYID];
							break;
						case PlotQuantity.MPMPLEPSXY:
							theValue=eplast[XYID];
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=eplast[XXID];
					sigy=eplast[YYID];
					sigxy=eplast[XYID];
					switch(component)
					{   case PlotQuantity.MPMPLEPSX:
							theValue=c*c*sigx + s*s*sigy - c*s*sigxy;
							break;
						case PlotQuantity.MPMPLEPSY:
							theValue=s*s*sigx + c*c*sigy + c*s*sigxy;
							break;
						case PlotQuantity.MPMPLEPSXY:
							theValue=2*c*s*(sigx-sigy) + (c*c-s*s)*sigxy;
							break;
					}
				}
				break;
				
			case PlotQuantity.MPMPLEPSZ:
				theValue=eplast[ZZID];
				break;
			
				// Strains
			case PlotQuantity.MPMEPSTOTX:
			case PlotQuantity.MPMEPSTOTY:
			case PlotQuantity.MPMEPSTOTXY:
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMEPSTOTX:
							theValue=eps[XXID]+eplast[XXID];
							break;
						case PlotQuantity.MPMEPSTOTY:
							theValue=eps[YYID]+eplast[YYID];
							break;
						case PlotQuantity.MPMEPSTOTXY:
							theValue=eps[XYID]+eplast[XYID];
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=eps[XXID]+eplast[XXID];
					sigy=eps[YYID]+eplast[YYID];
					sigxy=eps[XYID]+eplast[XYID];
					switch(component)
					{   case PlotQuantity.MPMEPSTOTX:
							theValue=c*c*sigx + s*s*sigy - c*s*sigxy;
							break;
						case PlotQuantity.MPMEPSTOTY:
							theValue=s*s*sigx + c*c*sigy + c*s*sigxy;
							break;
						case PlotQuantity.MPMEPSTOTXY:
							theValue=2*c*s*(sigx-sigy) + (c*c-s*s)*sigxy;
							break;
					}
				}
				break;
				
			case PlotQuantity.MPMEPSTOTZ:
				theValue=eps[ZZID]+eplast[ZZID];
				break;
		
			// Energy (totals are getting this point only)
			case PlotQuantity.MPMTOTENERGY:
			case PlotQuantity.MPMENERGY:
			case PlotQuantity.MPMTOTSTRENERGY:
			case PlotQuantity.MPMSTRENERGY:
				// stress in MPa, strain in %, mass in g, rho in g/cm^3 -> Joules
				theValue=strainEnergy;
				if(component==PlotQuantity.MPMSTRENERGY || component==PlotQuantity.MPMTOTSTRENERGY) break;
			case PlotQuantity.MPMTOTKINENERGY:
			case PlotQuantity.MPMKINENERGY:
				if(component==PlotQuantity.MPMKINENERGY || component==PlotQuantity.MPMTOTKINENERGY) theValue=0.;
				// mass in g, vel in mm/msec (once converted by current scales) -> Joules
				theValue+=0.5e-3*mass*(velx*velx+vely*vely)*doc.timeScale/doc.lengthScale;
				break;
			
			case PlotQuantity.MPMTOTWORKENERGY:
			case PlotQuantity.MPMWORKENERGY:
				theValue=workEnergy;
				break;
			
			case PlotQuantity.MPMTOTPLASTICENERGY:
			case PlotQuantity.MPMPLASTICENERGY:
				theValue=plastEnergy;
				break;
			
			case PlotQuantity.MPMTOTHEATENERGY:
			case PlotQuantity.MPMHEATENERGY:
				theValue=heatEnergy;
				break;
			
			case PlotQuantity.MPMTEMPERATURE:
				theValue=deltaTemp;
				break;
			
			// Velocity
			case PlotQuantity.MPMVELX:
				theValue=velx;
				break;
			case PlotQuantity.MPMVELY:
				theValue=vely;
				break;
			
			// Displacements
			case PlotQuantity.MPMDISPX:
				theValue=x-origx;
				break;
			case PlotQuantity.MPMDISPY:
				theValue=y-origy;
				break;
			
			// Position and angle
			case PlotQuantity.MPMPOS:
				theValue=(double)material;
				break;
			
			case PlotQuantity.MPMANGLEZ:
				theValue=angleZ;
				break;
			
			case PlotQuantity.MPMPOSX:
				theValue=x;
				break;
			
			case PlotQuantity.MPMPOSY:
				theValue=y;
				break;
			
			// shear components
			case PlotQuantity.MPMDUDY:
				theValue=dudy;
				break;
			case PlotQuantity.MPMDVDX:
				theValue=dvdx;
				break;
			
			// concentration
			case PlotQuantity.MPMCONCENTRATION:
				theValue=concentration;
				break;
			case PlotQuantity.MPMDCDY:
				theValue=dcdy;
				break;
			case PlotQuantity.MPMDCDX:
				theValue=dcdx;
				break;
			
			// histrory variables
			case PlotQuantity.MPMHISTORY1:
				theValue=history1;
				break;
			case PlotQuantity.MPMHISTORY2:
				theValue=history2;
				break;
			case PlotQuantity.MPMHISTORY3:
				theValue=history3;
				break;
			case PlotQuantity.MPMHISTORY4:
				theValue=history4;
				break;
			
			case PlotQuantity.MPMMASS:
				theValue=mass;
				break;
			
			case PlotQuantity.MPMTOTELEMENTCROSSINGS:
			case PlotQuantity.MPMELEMENTCROSSINGS:
				theValue=(double)elementCrossings;
				break;
				
			// Unknown
			default:
				break;
		}
		return theValue;
	}
	
	//---------------------------------------------------------------------
	// Class methods dealing with all material points
	//---------------------------------------------------------------------
	
	// load all material points with values for new plot component
	public static void loadPlotData(int component,ResultsDocument doc)
	{
		double dmin,dmax;
		double angle=0.;

		// load plotValue of each material point and get range of data
		int i;
		MaterialPoint mp;
		mp=doc.mpmPoints.get(0);
		mp.loadForPlot(component,angle,doc);
		dmin=mp.getPlotValue();
		dmax=mp.getPlotValue();
		for(i=1;i<doc.mpmPoints.size();i++)
		{	mp=doc.mpmPoints.get(i);
			mp.loadForPlot(component,angle,doc);
			dmin=Math.min(dmin,mp.getPlotValue());
			dmax=Math.max(dmax,mp.getPlotValue());
		}
		
		// position will expand range slightly
		if(component==PlotQuantity.MPMPOS)
		{	dmin=0.0;
			dmax=Math.max(1.0,dmax);
		}
		
		// changed to fixed limits
		Point2D.Double limits = doc.docCtrl.controls.adjustLimits(dmin, dmax);
		dmin=limits.x;
		dmax=limits.y;
		
		// set the spectum
		setSpectrum(dmin,dmax,doc);
		
		// tell plot view its range
		MoviePlotWindow movieFrame=doc.docCtrl.getMovieFrame();
		movieFrame.plotView.dataMin=dmin;
		movieFrame.plotView.dataMax=dmax;
		movieFrame.plotView.dataLimitsSet=true;
	}
	
	// set color spectrum current range by changing colors
	//  of all material points or element subelements. No need
	//  to reload all the data.
	public static void setSpectrum(double dmin,double dmax,ResultsDocument doc)
	{
		// adjust for no range
		if(Math.abs(dmax-dmin)<1.e-15)
		{	dmax+=1.;
			dmin-=1.;
		}
		double scale=1./(dmax-dmin);
		
		// particle plots
		int i;
		MaterialPoint mp;
		for(i=0;i<doc.mpmPoints.size();i++)
		{	mp=doc.mpmPoints.get(i);
			mp.setPlotColor(dmin,scale);
		}
	}
	
	//-----------------------------------------------------------------
	// Accessors
	//-----------------------------------------------------------------
	
	// set color from a rainbow
	public void setPlotColor(double dmin,double scale)
	{	plotColor=ColorPicker.PickRainbow(scale*(plotValue-dmin));
	}
	
	// return the plotValue
	public double getPlotValue() { return plotValue; }
	
	// get position as a point
	public Point2D.Double getPosition() { return new Point2D.Double(x,y); }
	
	// index to material type (zero based)
	public int materialIndex() { return material-1; }

}
