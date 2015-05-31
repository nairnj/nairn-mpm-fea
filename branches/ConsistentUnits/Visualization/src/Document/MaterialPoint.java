/*******************************************************************
	MaterialPoint.java
	NairnFEAMPMViz

	Created by John Nairn on 2/28/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.nio.*;
import java.awt.*;
import java.awt.geom.*;
import javax.vecmath.Vector3d;

public class MaterialPoint
{
	//---------------------------------------------------------------------
	// variables and constants
	//---------------------------------------------------------------------
	
	public static final int XXID=0;
	public static final int YYID=1;
	public static final int ZZID=2;
	public static final int XYID=3;
	public static final int XZID=3;
	public static final int YZID=3;
	
	public int num;
	public int inElem;			// zero based element number
	public double mass,angleZ,thickness,angleY,angleX;
	public int material;
	public double x,y,z,origx,origy,origz,velx,vely,velz;
	public double deltaTemp,plastEnergy,strainEnergy;
	public double history1,history2,history3,history4,concentration,dcdx,dcdy,dcdz;
	public double workEnergy,heatEnergy;
	public int elementCrossings;
	public double[] sigma;
	public double[] eps;
	public double[] eplast;
	public double[] erot;
	public double sizeX,sizeY,sizeZ;
	
	private double plotValue;
	public Color plotColor;
	
	//---------------------------------------------------------------------
	// initialize
	//---------------------------------------------------------------------
	
	MaterialPoint(int ptNum)
	{	num=ptNum;
		plotValue=0.;
		sigma=new double[6];
		eps=new double[6];
		eplast=new double[6];
		erot=new double[3];
		sizeX = -1.;
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
								boolean showSquarePts,boolean transformPts,double mpDiam)
	{	// depends on setting
		
		// adjust size
		if(sizeX>0.)
		{	radiix = radiix*sizeX/(0.01*mpDiam);
			radiiy = radiiy*sizeY/(0.01*mpDiam);
		}
	    
		if(showSquarePts)
		{	if(transformPts)
			{	// This works in Java 1.5
				// transformed (radiix,0)
				Matrix3 F = getDeformationGradient(resDoc);
				double r1x = radiix*F.get(0,0);
				double r1y = -radiix*F.get(1,0);
				// transformed (0,radiiy)
				double r2x = radiiy*F.get(0,1);
				double r2y = -radiiy*F.get(1,1);
				GeneralPath quad=new GeneralPath();
				Point2D.Float pathPt0=new Point2D.Float((float)(xpt+r1x+r2x),(float)(ypt+r1y+r2y));
				quad.moveTo(pathPt0.x,pathPt0.y);
				quad.lineTo((float)(xpt-r1x+r2x),(float)(ypt-r1y+r2y));
				quad.lineTo((float)(xpt-r1x-r2x),(float)(ypt-r1y-r2y));
				quad.lineTo((float)(xpt+r1x-r2x),(float)(ypt+r1y-r2y));
				quad.lineTo(pathPt0.x,pathPt0.y);
				/*
				// This requires Java 1.6
				Path2D.Double quad=new Path2D.Double();
				Point2D.Double pathPt0=new Point2D.Double((double)(xpt+r1x+r2x),
											(double)(ypt+r1y+r2y));
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
	
	public void readRecord(ByteBuffer bb,char[] mpmOrder,JNUnits units,boolean has3D)
	{	// required elements
		inElem=bb.getInt()-1;					// in element number (zero based)
		mass=bb.getDouble()*units.massScale();	// mass in g
		material=bb.getShort();					// material number
		bb.getShort();							// skip 2 bytes
		
		// 3D will have 3 angles, but not read in this app yet
		if(has3D)
		{	angleZ=bb.getDouble();					// angle in degrees about Z axis
			angleY=bb.getDouble();
			angleX=bb.getDouble();
			thickness=1.;
		}
		else
		{	angleZ=bb.getDouble();					// angle in degrees about Z axis
			thickness=bb.getDouble()*units.lengthScale();	// thickness in length units
		}
		
		x=bb.getDouble()*units.lengthScale();			// x position
		y=bb.getDouble()*units.lengthScale();			// y position
		// 3D will have third position
		if(has3D) z=bb.getDouble()*units.lengthScale();
		
		// new default puts original position (in length units) here
		if(mpmOrder[ReadArchive.ARCH_Defaults]=='Y')
		{	origx=bb.getDouble()*units.lengthScale();
			origy=bb.getDouble()*units.lengthScale();
			// 3D will have third position
			if(has3D) origz=bb.getDouble()*units.lengthScale();
		}

		// velocity (length units/time units), but archive has mm/sec so need to
		// divide by 1000 to get default units here which are ms
		if(mpmOrder[ReadArchive.ARCH_Velocity]=='Y')
		{	velx=bb.getDouble()*units.velocityScale();
			vely=bb.getDouble()*units.velocityScale();
			// 3D will have third velocity
			if(has3D) velz=bb.getDouble()*units.velocityScale();
		}
		else
		{	velx=0.;
			vely=0.;
			velz=0.;
		}

		// stress
		if(mpmOrder[ReadArchive.ARCH_Stress]=='Y')
		{	sigma[XXID]=bb.getDouble()*units.mpmStressScale();
			sigma[YYID]=bb.getDouble()*units.mpmStressScale();
			sigma[ZZID]=bb.getDouble()*units.mpmStressScale();
			sigma[XYID]=bb.getDouble()*units.mpmStressScale();
			// 3D will have two more components
			if(has3D)
			{	sigma[XZID]=bb.getDouble()*units.mpmStressScale();
				sigma[YZID]=bb.getDouble()*units.mpmStressScale();
			}
		}
		else
		{	sigma[XXID]=0.;
			sigma[YYID]=0.;
			sigma[ZZID]=0.;
			sigma[XYID]=0.;
			sigma[XZID]=0.;
			sigma[YZID]=0.;
		}
		
		// strain (in %)
		if(mpmOrder[ReadArchive.ARCH_Strain]=='Y')
		{	eps[XXID]=bb.getDouble()*units.strainScale();
			eps[YYID]=bb.getDouble()*units.strainScale();
			eps[ZZID]=bb.getDouble()*units.strainScale();
			eps[XYID]=bb.getDouble()*units.strainScale();
			// 3D will have two more components
			if(has3D)
			{	eps[XZID]=bb.getDouble()*units.strainScale();
				eps[YZID]=bb.getDouble()*units.strainScale();
			}
		}
		else
		{	eps[XXID]=0.;
			eps[YYID]=0.;
			eps[ZZID]=0.;
			eps[XYID]=0.;
			eps[XZID]=0.;
			eps[YZID]=0.;
		}
		
		// plastic strain (in %) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_PlasticStrain]=='Y')
		{	eplast[XXID]=bb.getDouble()*units.strainScale();
			eplast[YYID]=bb.getDouble()*units.strainScale();
			eplast[ZZID]=bb.getDouble()*units.strainScale();
			eplast[XYID]=bb.getDouble()*units.strainScale();
			// 3D will have two more components
			if(has3D)
			{	eplast[XZID]=bb.getDouble()*units.strainScale();
				eplast[YZID]=bb.getDouble()*units.strainScale();
			}
		}
		else
		{	eplast[XXID]=0.;
			eplast[YYID]=0.;
			eplast[ZZID]=0.;
			eplast[XYID]=0.;
			eplast[XZID]=0.;
			eplast[YZID]=0.;
		}
				
		// old method for original positions (in length units)
		if(mpmOrder[ReadArchive.ARCH_OldOrigPosition]=='Y')
		{	origx=bb.getDouble()*units.lengthScale();
			origy=bb.getDouble()*units.lengthScale();
			// 3D will have one more position
			if(has3D) origz=bb.getDouble()*units.lengthScale();
		}
		
		// external work (cumulative) in J (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_WorkEnergy]=='Y')
			workEnergy=bb.getDouble()*units.energyScale();
		else
			workEnergy=0.;
		
		// temperature (C) (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_DeltaTemp]=='Y')
			deltaTemp=bb.getDouble();
		else
			deltaTemp=0.;
				
		// total plastic energy (Volume*energy) in J (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_PlasticEnergy]=='Y')
			plastEnergy=bb.getDouble()*units.energyScale();
		else
			plastEnergy=0.;
		
		// shear components (dimensionless) no ignores
		if(mpmOrder[ReadArchive.ARCH_ShearComponents]=='Y')
		{	bb.getDouble();
			bb.getDouble();
		}
			
		// total strain energy (Volume*energy) in J (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_StrainEnergy]=='Y')
			strainEnergy=bb.getDouble()*units.energyScale();
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
			dcdx=bb.getDouble()/units.lengthScale();
			dcdy=bb.getDouble()/units.lengthScale();
			// 3D will have one more
			if(has3D) dcdz=bb.getDouble()/units.lengthScale();
		}
		else
		{	concentration=0.;
			dcdx=0.;
			dcdy=0.;
			dcdz=0.;
		}
		
		// thermal energy (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_HeatEnergy]=='Y')
			heatEnergy=bb.getDouble()*units.energyScale();
		else
			heatEnergy=0.;
			
		// element crossings (not converted for units)
		if(mpmOrder[ReadArchive.ARCH_ElementCrossings]=='Y')
			elementCrossings=bb.getInt();
		else
			elementCrossings=0;
			
		// rot strain (in degrees) (not converted for units)
		// actually initial angle to subtraction of material angle to
		// get rotational strain
		if(mpmOrder[ReadArchive.ARCH_RotStrain]=='Y')
		{	erot[0]=bb.getDouble();
			if(has3D)
			{	erot[1]=bb.getDouble();
				erot[2]=bb.getDouble();
			}
		}
		else
		{	erot[0]=0.;		// assume initial angles were all zero
			erot[1]=0.;
			erot[2]=0.;
		}
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
		
			case PlotQuantity.MPMSIGMAXZ:
				theValue=sigma[XZID];
				break;
		
			case PlotQuantity.MPMSIGMAYZ:
				theValue=sigma[YZID];
				break;
		
			// Strains
			case PlotQuantity.MPMEPSX:
			case PlotQuantity.MPMEPSY:
			case PlotQuantity.MPMEPSXY:
			{	Matrix3 biot = getElasticStrain(doc);
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMEPSX:
							theValue=biot.get(0,0);
							break;
						case PlotQuantity.MPMEPSY:
							theValue=biot.get(1,1);
							break;
						case PlotQuantity.MPMEPSXY:
							theValue=2.*biot.get(0,1);
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=biot.get(0,0);
					sigy=biot.get(1,1);
					sigxy=2.*biot.get(0,1);
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
				theValue *= doc.units.strainScale();
				break;
			}
				
			case PlotQuantity.MPMEPSZ:
			{	Matrix3 biot = getElasticStrain(doc);
				theValue=biot.get(1,1)*doc.units.strainScale();
				break;
			}
		
			case PlotQuantity.MPMEPSXZ:
			{	Matrix3 biot = getElasticStrain(doc);
				theValue=2.*biot.get(0,2)*doc.units.strainScale();
				break;
			}
		
			case PlotQuantity.MPMEPSYZ:
			{	Matrix3 biot = getElasticStrain(doc);
				theValue=2.*biot.get(1,2)*doc.units.strainScale();
				break;
			}
		
			// Plastic Strains
			case PlotQuantity.MPMPLEPSX:
			case PlotQuantity.MPMPLEPSY:
			case PlotQuantity.MPMPLEPSXY:
			{	Matrix3 biot = getPlasticStrain(doc);
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMPLEPSX:
							theValue=biot.get(0,0);
							break;
						case PlotQuantity.MPMPLEPSY:
							theValue=biot.get(1,1);
							break;
						case PlotQuantity.MPMPLEPSXY:
							theValue=2.*biot.get(0,1);
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=biot.get(0,0);
					sigy=biot.get(1,1);
					sigxy=2.*biot.get(0,1);
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
				theValue *= doc.units.strainScale();
				break;
			}
				
			case PlotQuantity.MPMPLEPSZ:
			{	Matrix3 biot = getPlasticStrain(doc);
				theValue=biot.get(1,1)*doc.units.strainScale();
				break;
			}
			
			case PlotQuantity.MPMPLEPSXZ:
			{	Matrix3 biot = getPlasticStrain(doc);
				theValue=2.*biot.get(0,2)*doc.units.strainScale();
				break;
			}
			
			case PlotQuantity.MPMPLEPSYZ:
			{	Matrix3 biot = getPlasticStrain(doc);
				theValue=2.*biot.get(1,2)*doc.units.strainScale();
				break;
			}
			
			// Strains
			case PlotQuantity.MPMEPSTOTX:
			case PlotQuantity.MPMEPSTOTY:
			case PlotQuantity.MPMEPSTOTXY:
			{	Matrix3 biot = getBiotStrain(doc);
				if(angle==0.)
				{   switch(component)
					{   case PlotQuantity.MPMEPSTOTX:
							theValue=biot.get(0,0);
							break;
						case PlotQuantity.MPMEPSTOTY:
							theValue=biot.get(1,1);
							break;
						case PlotQuantity.MPMEPSTOTXY:
							theValue=2.*biot.get(0,1);
							break;
						default:
							break;
					}
				}
				else
				{   radAngle=Math.PI*angle/180.;
					c=Math.cos(radAngle);
					s=Math.sin(radAngle);
					sigx=biot.get(0,0);
					sigy=biot.get(1,1);
					sigxy=2.*biot.get(0,1);
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
				theValue *= doc.units.strainScale();
				break;
			}
				
			case PlotQuantity.MPMEPSTOTZ:
				if(doc.is3D())
				{	Matrix3 biot = getBiotStrain(doc);
					theValue = biot.get(0,0)*doc.units.strainScale();
				}
				else
				{	// 2D only
					if(doc.materials.get(materialIndex()).hasPlasticStrainForGradient(doc))
						theValue=eps[ZZID]+eplast[ZZID];
					else
						theValue=eps[ZZID];
				}
				break;
		
			case PlotQuantity.MPMEPSTOTXZ:
			{	Matrix3 biot = getBiotStrain(doc);
				theValue = 2.*biot.get(0,2)*doc.units.strainScale();
				break;
			}
				
			case PlotQuantity.MPMEPSTOTYZ:
			{	Matrix3 biot = getBiotStrain(doc);
				theValue = 2.*biot.get(1,2)*doc.units.strainScale();
				break;
			}
				
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
				theValue+=0.5*mass*(velx*velx+vely*vely)*doc.units.calcVelocityScale();
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
			case PlotQuantity.MPMVELZ:
				theValue=velz;
				break;
			
			// Displacements
			case PlotQuantity.MPMDISPX:
				theValue=x-origx;
				break;
			case PlotQuantity.MPMDISPY:
				theValue=y-origy;
				break;
			case PlotQuantity.MPMDISPZ:
				theValue=y-origz;
				break;
			
			// Position and angle
			case PlotQuantity.MPMPOS:
				theValue=(double)material;
				break;
			
			case PlotQuantity.MPMANGLEZ:
				theValue=angleZ;
				break;
			case PlotQuantity.MPMANGLEY:
				theValue=angleY;
				break;
			case PlotQuantity.MPMANGLEX:
				theValue=angleX;
				break;
			
			case PlotQuantity.MPMPOSX:
				theValue=x;
				break;
			
			case PlotQuantity.MPMPOSY:
				theValue=y;
				break;
			
			case PlotQuantity.MPMPOSZ:
				theValue=z;
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
			case PlotQuantity.MPMDCDZ:
				theValue=dcdz;
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
	
	// set dimensionless size for membrane particles
	public void setDimensionlessSize(double sx,double sy,double sz)
	{	sizeX = sx;
		sizeY = sy;
		sizeZ = sz;
	}
	
	// get biot strain from total strain in biot matrix
	public Matrix3 getBiotStrain(ResultsDocument doc)
	{	Matrix3 F = getDeformationGradient(doc);
		Matrix3 biot = F.LeftDecompose(null);
		biot.set(0,0,biot.get(0,0)-1.);
		biot.set(1,1,biot.get(1,1)-1.);
		biot.set(2,2,biot.get(2,2)-1.);
		return biot;
	}
	
	// get biot elastic strain tensor from total strain and maybe needs plastic strain too
	public Matrix3 getElasticStrain(ResultsDocument doc)
	{	Matrix3 biot;
		if(doc.materials.get(materialIndex()).AltStrainContains()==MaterialBase.ENG_BIOT_PLASTIC_STRAIN)
		{	// elastic strain = total strain minus plastic strain
			biot = getBiotStrain(doc);
			double scale = 1./doc.units.strainScale();
			biot.set(0,0,biot.get(0,0)-scale*eplast[XXID]);
			biot.set(1,1,biot.get(1,1)-scale*eplast[YYID]);
			biot.set(2,2,biot.get(2,2)-scale*eplast[ZZID]);
			biot.set(0,1,biot.get(0,1)-0.5*scale*eplast[XYID]);
			biot.set(1,0,biot.get(0,1));
			if(doc.is3D())
			{	biot.set(0,2,biot.get(0,2)-0.5*scale*eplast[XZID]);
				biot.set(2,0,biot.get(0,2));
				biot.set(1,2,biot.get(1,2)-0.5*scale*eplast[YZID]);
				biot.set(2,1,biot.get(1,2));
			}
		}
		else if(doc.materials.get(materialIndex()).AltStrainContains()==MaterialBase.LEFT_CAUCHY_ELASTIC_B_STRAIN)
		{	// get elastic strain from elastic B
			biot = getElasticStrainFromB(doc);
		}
		else
		{	// elastic strain = total strain
			biot = getBiotStrain(doc);
		}
		return biot;
	}
	
	// get biot plastic strain tensor from plastic strain or from elastic and plastis
	// zero if not archived or not plastic strain
	public Matrix3 getPlasticStrain(ResultsDocument doc)
	{	Matrix3 biot;
		if(doc.materials.get(materialIndex()).AltStrainContains()==MaterialBase.ENG_BIOT_PLASTIC_STRAIN)
		{	// plastic strain all ready
			double scale = 1./doc.units.strainScale();
			biot = new Matrix3();
			biot.set(0,0,scale*eplast[XXID]);
			biot.set(1,1,scale*eplast[YYID]);
			biot.set(0,0,scale*eplast[XYID]);
			biot.set(0,1,0.5*scale*eplast[XYID]);
			biot.set(1,0,biot.get(0,1));
			if(doc.is3D())
			{	biot.set(0,2,0.5*scale*eplast[XZID]);
				biot.set(2,0,biot.get(0,2));
				biot.set(1,2,0.5*scale*eplast[YZID]);
				biot.set(2,1,biot.get(1,2));
			}
		}
		else if(doc.materials.get(materialIndex()).AltStrainContains()==MaterialBase.LEFT_CAUCHY_ELASTIC_B_STRAIN)
		{	// plastic strain = total strain minus elastic strain (from B) (needs both)
			biot = getBiotStrain(doc);
			Matrix3 biotElast = getElasticStrainFromB(doc);
			biot.subtract(biotElast);
		}
		else
		{	// others have zero plastic strain
			biot = new Matrix3();
		}
		return biot;
	}
	
	// get elastic biot strain tensor from elastic B tensor
	public Matrix3 getElasticStrainFromB(ResultsDocument doc)
	{
		// Get sqrt(B)
		double scale = 1./doc.units.strainScale();
		Matrix3 B = null;
		Vector3d lam = null;
		if(doc.is3D())
		{	B = new Matrix3(scale*eplast[XXID],scale*eplast[XYID],scale*eplast[XZID],
				scale*eplast[XYID],scale*eplast[YYID],scale*eplast[XZID],
				scale*eplast[XZID],scale*eplast[XZID],scale*eplast[ZZID]);
			lam = new Vector3d();		// 3D does not need precaculated eigenvalues
		}
		else
		{	B = new Matrix3(scale*eplast[XXID],scale*eplast[XYID],scale*eplast[XYID],
					scale*eplast[YYID],scale*eplast[ZZID]);
			lam = B.Eigenvalues();		// 2D needs eigenvalue inputs
		}
		Matrix3 Q = B.Eigenvectors(lam);
		
		// put U = sqrt(B) in output matrix
		Matrix3 biot = new Matrix3(Math.sqrt(lam.x),0.,0.,Math.sqrt(lam.y),Math.sqrt(lam.z));
		biot.RMRT(Q);
		
		// Get rotation matrix
		Matrix3 F = getDeformationGradient(doc);
		Matrix3 R = new Matrix3();
		F.LeftDecompose(R);
		
		// Rotate U to V and subtract I
		biot.RMRT(R);
		biot.set(0,0,biot.get(0,0)-1.);
		biot.set(1,1,biot.get(1,1)-1.);
		biot.set(2,2,biot.get(2,2)-1.);
		
		// release and done
		return biot;
	}

	
	// get deformation gradient
	public Matrix3 getDeformationGradient(ResultsDocument resDoc)
	{	Matrix3 F = new Matrix3();
		double scale = 1./resDoc.units.strainScale();
		
		if(resDoc.is3D())
		{	// tensorial rotational strain archives are found from rotation (zero if not archived)
			double wxy,wxz,wyz;
			wxy = Math.PI*(erot[0]-angleZ)/180.;			// 0.5(dv/dx-du/dy)
			wxz = -Math.PI*(erot[1]-angleY)/180.;		// 0.5(dw/dx-du/dz)
			wyz = Math.PI*(erot[2]-angleX)/180.;			// 0.5(dw/dy-dv/dz)
			
			if(resDoc.materials.get(materialIndex()).hasPlasticStrainForGradient(resDoc))
			{	F.set(0,0,1.+scale*(eps[XXID]+eplast[XXID]));
				F.set(0,1,0.5*scale*(eps[XYID]+eplast[XYID])-wxy);
				F.set(0,2,0.5*scale*(eps[XZID]+eplast[XZID])-wxz);
				F.set(1,0,0.5*scale*(eps[XYID]+eplast[XYID])+wxy);
				F.set(1,1,1.+scale*(eps[YYID]+eplast[YYID]));
				F.set(1,2,0.5*scale*(eps[YZID]+eplast[YZID])-wyz);
				F.set(2,0,0.5*scale*(eps[XZID]+eplast[XZID])+wxz);
				F.set(2,1,0.5*scale*(eps[YZID]+eplast[YZID])+wyz);
				F.set(2,2,1.+scale*(eps[ZZID]+eplast[ZZID]));
			}
			else
			{	F.set(0,0,1.+scale*eps[XXID]);
				F.set(0,1,0.5*scale*eps[XYID]-wxy);
				F.set(0,2,0.5*scale*eps[XZID]-wxz);
				F.set(1,0,0.5*scale*eps[XYID]+wxy);
				F.set(1,1,1.+scale*eps[YYID]);
				F.set(1,2,0.5*scale*eps[YZID]-wyz);
				F.set(2,0,0.5*scale*eps[XZID]+wxz);
				F.set(2,1,0.5*scale*eps[YZID]+wyz);
				F.set(2,2,1.+scale*eps[ZZID]);
			}
			F.setIs2D(false);

		}
		else
		{	double wrot=Math.PI*(erot[0]-angleZ)/180.;
			if(resDoc.materials.get(materialIndex()).hasPlasticStrainForGradient(resDoc))
			{	F.set(0,0,1.+scale*(eps[MaterialPoint.XXID]+eplast[MaterialPoint.XXID]));
				F.set(0,1,0.5*scale*(eps[MaterialPoint.XYID]+eplast[MaterialPoint.XYID]) - wrot);
				F.set(1,0,0.5*scale*(eps[MaterialPoint.XYID]+eplast[MaterialPoint.XYID]) + wrot);
				F.set(1,1,1.+scale*(eps[MaterialPoint.YYID]+eplast[MaterialPoint.YYID]));
				F.set(2,2,1.+scale*(eps[MaterialPoint.ZZID]+eplast[MaterialPoint.ZZID]));
			}
			else
			{	F.set(0,0,1.+scale*eps[MaterialPoint.XXID]);
				F.set(0,1,0.5*scale*eps[MaterialPoint.XYID] - wrot);
				F.set(1,0,0.5*scale*eps[MaterialPoint.XYID] + wrot);
				F.set(1,1,1.+scale*eps[MaterialPoint.YYID]);
				F.set(2,2,1.+scale*eps[MaterialPoint.ZZID]);
			}
		}
		return F;
	}


}
