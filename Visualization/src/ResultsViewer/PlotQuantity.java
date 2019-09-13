/*******************************************************************
	PlotQuantity.java
	NairnFEAMPMViz

	Created by John Nairn on Wed Mar 10 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class PlotQuantity extends PlotControl
{
	static final long serialVersionUID=18L;
	
	// MPM plotting options
	static final int MPMSIGMAX=1;
	static final int MPMSIGMAY=2;
	static final int MPMSIGMAXY=3;
	static final int MPMSIGMAZ=4;
	static final int MPMSIGMAXZ=5;
	static final int MPMSIGMAYZ=6;
	
	static final int MPMEPSX=7;
	static final int MPMEPSY=8;
	static final int MPMEPSXY=9;
	static final int MPMEPSZ=10;
	static final int MPMEPSXZ=11;
	static final int MPMEPSYZ=12;
	
	static final int MPMPLEPSX=13;
	static final int MPMPLEPSY=14;
	static final int MPMPLEPSXY=15;
	static final int MPMPLEPSZ=16;
	static final int MPMPLEPSXZ=17;
	static final int MPMPLEPSYZ=18;
	
	static final int MPMVELX=19;
	static final int MPMVELY=20;
	static final int MPMVELZ=21;
	
	static final int MPMVELVEC=22;
	static final int MPMSTRENERGY=23;
	static final int MPMKINENERGY=24;
	static final int MPMENERGY=25;
	static final int MPMTOTSTRENERGY=26;
	static final int MPMTOTKINENERGY=27;
	static final int MPMTOTENERGY=28;
	static final int MPMPOS=29;
	
	static final int MPMDISPX=30;
	static final int MPMDISPY=31;
	static final int MPMDISPZ=32;
	
	static final int MPMPOSX=33;
	static final int MPMPOSY=34;
	static final int MPMPOSZ=35;
	
	static final int MPMEXPRESSION=36;
	static final int MPMWORKENERGY=37;
	static final int MPMPLASTICENERGY=38;
	static final int MPMTEMPERATURE=39;
	static final int MPMDVDX=40;
	static final int MPMDUDY=41;
	static final int MPMHISTORYOLD=42;
	static final int MPMTOTWORKENERGY=43;
	static final int MPMTOTPOTENERGY=44;				// no longer used
	static final int MPMTOTPLASTICENERGY=45;
    static final int MPMJ1=46;
	static final int MPMJ2=47;
	static final int MPMKI=48;
	static final int MPMKII=49;
	static final int MPMLENGTH=50;
	static final int MPMMASS=51;
	static final int MPMARCHIVETIME=52;
	static final int MPMGLOBALRESULTS=53;
	static final int MPMCONCENTRATION=54;
	static final int MPMDCDX=55;
	static final int MPMDCDY=56;
	static final int MPMDCDZ=57;
	static final int MPMCRACKRELEASE=58;
	static final int MPMCRACKABSORB=59;
	static final int MPMNORMALCTOD=60;
	static final int MPMSHEARCTOD=61;
	static final int MPMHEATENERGY=62;
	static final int MPMTOTHEATENERGY=63;
	static final int MPMCRACKPROFILE=64;
	static final int MPMOPENINGFRACTION=65;
	static final int MPMSHEARFRACTION=66;
	static final int MPMANGLEZ=67;
	static final int MPMANGLEY=68;
	static final int MPMANGLEX=69;
	static final int MPMELEMENTCROSSINGS=70;
	static final int MPMTOTELEMENTCROSSINGS=71;
	static final int MPMDEBONDLENGTH=72;
	static final int MPMDEBONDNCTOD=73;
	static final int MPMDEBONDSCTOD=74;
	static final int MPMHISTORY1=75;
	static final int MPMHISTORY2=76;
	static final int MPMHISTORY3=77;
	static final int MPMHISTORY4=78;
	
	static final int MPMEPSTOTX=90;
	static final int MPMEPSTOTY=91;
	static final int MPMEPSTOTXY=92;
	static final int MPMEPSTOTZ=93;
	static final int MPMEPSTOTXZ=94;
	static final int MPMEPSTOTYZ=95;
	
	static final int MPMEQUIVSTRESS=105;
	static final int MPMEQUIVSTRAIN=106;
	static final int MPMPRESSURE=107;

	static final int MPMSPINMOMENTUMX=110;
	static final int MPMSPINMOMENTUMY=111;
	static final int MPMSPINMOMENTUMZ=112;
	static final int MPMSPINVELOCITYX=113;
	static final int MPMSPINVELOCITYY=114;
	static final int MPMSPINVELOCITYZ=115;
	
	static final int MESHONLY=1001;
	static final int MESHSIGMAX=1002;
	static final int MESHSIGMAY=1003;
	static final int MESHSIGMAXY=1004;
	static final int MESHSIGMAZ=1005;
	static final int MESHDISPX=1006;
	static final int MESHDISPY=1007;
	static final int MESHMATERIAL=1008;
	static final int MESHSTRAINX=1009;
	static final int MESHSTRAINY=1010;
	static final int MESHSTRAINXY=1011;
	static final int MESHSTRAINZ=1012;
	static final int MESHDVDX=1013;
	static final int MESHDUDY=1014;
	static final int MESHELEMSIGMAX=1015;
	static final int MESHELEMSIGMAY=1016;
	static final int MESHELEMSIGMAXY=1017;
	static final int MESHELEMSIGMAZ=1018;
	static final int MESHFORCEX=1019;
	static final int MESHFORCEY=1020;
	static final int MESHSTRAINENERGY=1021;
	static final int FEAEXPRESSION=1022;
	static final int MESHNODEX=1023;		// in expressions
	static final int MESHNODEY=1024;		// in expressions
	static final int INTERFACETRACTION_N=1025;
	static final int INTERFACETRACTION_T=1026;
	static final int MESHANGLE=1027;
	static final int MESHNODEDISTANCE=1028;		// in expressions - x^2+y^2 
	static final int MESHNODEANGLE=1029;		// in expressions - ccw angle in radians from x axis
	static final int MESHPRESSURE=1032;
	
	public static int mpmMovieQuant=-1;
	public static int mpmMovieComp=-1;
	public static int mpmTimeQuant=-1;
	public static int mpmTimeComp=-1;
	public static int mpmMesh2DQuant=-1;
	public static int mpmMesh2DComp=-1;
	public static int feaMeshQuant=-1;
	public static int feaMeshComp=-1;
	
	static final int IMPORTANDPLOTFILE=5000;
	
	// pop-up menus
	public JComboBox<PlotMenuItem> quant=new JComboBox<PlotMenuItem>();
	public JComboBox<String> cmpnt=new JComboBox<String>();
	private int checkMeshItem;
	private boolean suspendStoreSelection = false;
	
	// axes
	private String xchar="x";
	private String ychar="y";
	private String zchar="z";
	
	// initialize
	PlotQuantity(DocViewer dc)
	{   super(ControlPanel.WIDTH,64,dc);
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		//setLayout(new GridLayout(2,2));
		setLayout(gridbag);

		// label and quantity menu
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 0.0;
		JLabel label=new JLabel("Plot:",JLabel.RIGHT);
		gridbag.setConstraints(label, c);
		add(label);
		
		// quantity combo box
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(quant, c);
		quant.setToolTipText("Selected the archived quantity to be plotter");
		add(quant);
		
		// when quantity changes, update component menu
		quant.addItemListener(new ItemListener()
		{	public void itemStateChanged(ItemEvent e)
			{	setComponentMenu();
				docCtrl.controls.changeComponent();
			}
		});
		
		// label and quantity menu
		c.gridwidth = 1;
		//c.fill = GridBagConstraints.BOTH;
		c.weightx = 0.0;
		label=new JLabel("Component:",JLabel.RIGHT);
		gridbag.setConstraints(label, c);
		add(label);
		
		// quantity combo box
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(cmpnt, c);
		cmpnt.addItem("--");
		cmpnt.setToolTipText("Select component of currently selected quantity, which will depend on the type of quantity");
		add(cmpnt);
		
		setEnabled(LoadArchive.NO_PLOT);
	}
	
	// called when new file loaded
	public void setEnabled(int selected)
	{
		if(selected!=LoadArchive.NO_PLOT)
		{   quant.setEnabled(true);
			cmpnt.setEnabled(true);
		
			// update the items
			quant.removeAllItems();
			
			if(docCtrl.resDoc.isAxisymmetric())
			{	xchar="r";
				ychar="z";
				zchar="t";
			}
			
			int preselect = -1;
			int preComp = -1;
			suspendStoreSelection = true;
			
			if(docCtrl.resDoc.isMPMAnalysis())
			{	byte [] arch=docCtrl.resDoc.archFormat.getBytes();
				
				if(arch[ReadArchive.ARCH_Stress]=='Y')
				{	quant.addItem(new PlotMenuItem("Stress",MPMSIGMAX));
					quant.addItem(new PlotMenuItem("Pressure",MPMPRESSURE));
					quant.addItem(new PlotMenuItem("Equiv. Stress",MPMEQUIVSTRESS));
				}
				if(arch[ReadArchive.ARCH_Strain]=='Y')
				{	if(arch[ReadArchive.ARCH_PlasticStrain]=='Y')
					{	// need both to support all materials
						quant.addItem(new PlotMenuItem("Total Strain",MPMEPSTOTX));
						quant.addItem(new PlotMenuItem("Elastic Strain",MPMEPSX));
						quant.addItem(new PlotMenuItem("Plastic Strain",MPMPLEPSX));
						quant.addItem(new PlotMenuItem("Equiv. Strain",MPMEQUIVSTRAIN));
					}
					else if(docCtrl.resDoc.units.getNewUnitsVersion())
					{	// new version only needs strain to get total strain
						quant.addItem(new PlotMenuItem("Total Strain",MPMEPSTOTX));
						quant.addItem(new PlotMenuItem("Equiv. Strain",MPMEQUIVSTRAIN));
					}
				}
					
				if(arch[ReadArchive.ARCH_StrainEnergy]=='Y')
				{   quant.addItem(new PlotMenuItem("Strain Energy",MPMSTRENERGY));
					if(arch[ReadArchive.ARCH_Velocity]=='Y')
					{	quant.addItem(new PlotMenuItem("Kinetic Energy",MPMKINENERGY));
						quant.addItem(new PlotMenuItem("Energy",MPMENERGY));
					}
				}
				else if(arch[ReadArchive.ARCH_Velocity]=='Y')
					quant.addItem(new PlotMenuItem("Kinetic Energy",MPMKINENERGY));
					
				if(arch[ReadArchive.ARCH_PlasticEnergy]=='Y')
					quant.addItem(new PlotMenuItem("Plastic Energy",MPMPLASTICENERGY));
				if(arch[ReadArchive.ARCH_WorkEnergy]=='Y')
					quant.addItem(new PlotMenuItem("Work Energy",MPMWORKENERGY));
				if(arch[ReadArchive.ARCH_HeatEnergy]=='Y')
					quant.addItem(new PlotMenuItem("Heat Energy",MPMHEATENERGY));
					
				if(arch[ReadArchive.ARCH_Velocity]=='Y')
					quant.addItem(new PlotMenuItem("Velocity",MPMVELX));
				if(arch[ReadArchive.ARCH_SpinVelocity]=='Y')
					quant.addItem(new PlotMenuItem("Angular Velocity",MPMSPINVELOCITYX));
				if(arch[ReadArchive.ARCH_SpinMomentum]=='Y')
					quant.addItem(new PlotMenuItem("Angular Momentum",MPMSPINMOMENTUMX));
				quant.addItem(new PlotMenuItem("Displacement",MPMDISPX));
				checkMeshItem=quant.getItemCount();
				quant.addItem(new PlotMenuItem("Material",MPMPOS));
				quant.addItem(new PlotMenuItem("Material Angle",MPMANGLEZ));
				quant.addItem(new PlotMenuItem("Density",MPMMASS));
				if(arch[ReadArchive.ARCH_DeltaTemp]=='Y')
					quant.addItem(new PlotMenuItem("Temperature",MPMTEMPERATURE));
					
				if(arch[ReadArchive.ARCH_Concentration]=='Y')
				{   if(docCtrl.resDoc.hasPorePressure)
					{	quant.addItem(new PlotMenuItem("Pore Pressure",MPMCONCENTRATION));
						quant.addItem(new PlotMenuItem("Pore Press Gradient",MPMDCDX));
					}
					else
					{	quant.addItem(new PlotMenuItem("Concentration",MPMCONCENTRATION));
						quant.addItem(new PlotMenuItem("Conc Gradient",MPMDCDX));
					}
				}
				
				if(arch[ReadArchive.ARCH_History]=='Y')
					quant.addItem(new PlotMenuItem("History 1",MPMHISTORY1));
				else if(arch[ReadArchive.ARCH_History]!='N')
				{	int history=(int)arch[ReadArchive.ARCH_History];
					if((history & 0x01) !=0)
						quant.addItem(new PlotMenuItem("History 1",MPMHISTORY1));
					if((history & 0x02) !=0)
						quant.addItem(new PlotMenuItem("History 2",MPMHISTORY2));
					if((history & 0x04) !=0)
						quant.addItem(new PlotMenuItem("History 3",MPMHISTORY3));
					if((history & 0x08) !=0)
						quant.addItem(new PlotMenuItem("History 4",MPMHISTORY4));
				}
				
				if(arch[ReadArchive.ARCH_ElementCrossings]=='Y')
					quant.addItem(new PlotMenuItem("Element Crossings",MPMELEMENTCROSSINGS));
				
				// additional time plot options
				if(selected==LoadArchive.TIME_PLOT)
				{	// global results
					if(docCtrl.resDoc.globalArchive!=null)
						quant.addItem(new PlotMenuItem("Global Results",MPMGLOBALRESULTS));
					
					if(!docCtrl.resDoc.is3D())
					{	byte [] carch=docCtrl.resDoc.crackFormat.getBytes();
					
						if(carch[ReadArchive.ARCH_JIntegral]=='Y')
						{	quant.addItem(new PlotMenuItem("J1",MPMJ1));
							quant.addItem(new PlotMenuItem("J2",MPMJ2));
						}
					
						if(carch[ReadArchive.ARCH_StressIntensity]=='Y')
						{	quant.addItem(new PlotMenuItem("KI",MPMKI));
							quant.addItem(new PlotMenuItem("KII",MPMKII));
						}
					
						if(carch[ReadArchive.ARCH_BalanceResults]=='Y')
						{	quant.addItem(new PlotMenuItem("Global Released",MPMCRACKRELEASE));
							quant.addItem(new PlotMenuItem("Global Dissipated",MPMCRACKABSORB));
						}
					
						quant.addItem(new PlotMenuItem("Crack Length",MPMLENGTH));
						quant.addItem(new PlotMenuItem("Debonded Crack Length",MPMDEBONDLENGTH));
						quant.addItem(new PlotMenuItem("Normal CTOD",MPMNORMALCTOD));
						quant.addItem(new PlotMenuItem("Shear CTOD",MPMSHEARCTOD));
						quant.addItem(new PlotMenuItem("Debond Tip Normal COD",MPMDEBONDNCTOD));
						quant.addItem(new PlotMenuItem("Debond Tip Shear COD",MPMDEBONDSCTOD));
					}
					
					preselect = mpmTimeQuant;
					preComp = mpmTimeComp;
				}
				
				// additional time plot options
				else if(selected==LoadArchive.MESH2D_PLOT)
				{	// x-y crack results
					quant.addItem(new PlotMenuItem("Crack Profile",MPMCRACKPROFILE));
					quant.addItem(new PlotMenuItem("Crack Normal CTOD",MPMNORMALCTOD));
					quant.addItem(new PlotMenuItem("Crack Tangential CTOD",MPMSHEARCTOD));
					quant.addItem(new PlotMenuItem("Crack Opening Fraction",MPMOPENINGFRACTION));
					quant.addItem(new PlotMenuItem("Crack Sliding Fraction",MPMSHEARFRACTION));
					
					preselect = mpmMesh2DQuant;
					preComp = mpmMesh2DComp;
				}
				
				else
				{	preselect = mpmMovieQuant;
					preComp = mpmMovieComp;
				}
				
				// import and plot
				quant.addItem(new PlotMenuItem("Import...",IMPORTANDPLOTFILE));
			}
			
			// FEA plot quantities
			else
			{	char [] arch=docCtrl.resDoc.feaArchFormat;
			
				checkMeshItem=quant.getItemCount();
				quant.addItem(new PlotMenuItem("Mesh Only",MESHONLY));
				if(arch[ReadArchive.ARCH_FEAAvgStress]=='Y')
				{	quant.addItem(new PlotMenuItem("Stress",MESHSIGMAX));
					quant.addItem(new PlotMenuItem("Pressure",MESHPRESSURE));
				}
				if(arch[ReadArchive.ARCH_FEAElemStress]=='Y')
					quant.addItem(new PlotMenuItem("Element Stress",MESHELEMSIGMAX));
				if(arch[ReadArchive.ARCH_FEADisplacements]=='Y')
					quant.addItem(new PlotMenuItem("Strain",MESHSTRAINX));
				if(arch[ReadArchive.ARCH_FEAElemEnergy]=='Y')
					quant.addItem(new PlotMenuItem("Energy",MESHSTRAINENERGY));
				if(arch[ReadArchive.ARCH_FEAElemForce]=='Y')
					quant.addItem(new PlotMenuItem("Element Force",MESHFORCEX));
				if(arch[ReadArchive.ARCH_Interfaces]=='Y' && arch[ReadArchive.ARCH_FEAElemStress]=='Y')
					quant.addItem(new PlotMenuItem("Interface Traction",INTERFACETRACTION_N));
				if(arch[ReadArchive.ARCH_FEADisplacements]=='Y')
				{	quant.addItem(new PlotMenuItem("Displacement",MESHDISPX));
					quant.addItem(new PlotMenuItem("Shear Component",MESHDVDX));
				}
				quant.addItem(new PlotMenuItem("Material",MESHMATERIAL));
				quant.addItem(new PlotMenuItem("Material Angle",MESHANGLE));
				
				// import and plot
				quant.addItem(new PlotMenuItem("Import...",IMPORTANDPLOTFILE));

				preselect = feaMeshQuant;
				preComp = feaMeshComp;
			}
			
			// preselect and item if has one
			if(preselect>=0)
			{	for(int i=0;i<quant.getItemCount();i++)
				{	PlotMenuItem pm=(PlotMenuItem)quant.getItemAt(i);
					if(pm.getTag()==preselect)
					{	quant.setSelectedIndex(i);
						if(preComp>=0 && preComp<cmpnt.getItemCount())
							cmpnt.setSelectedIndex(preComp);
						break;
					}
				}
			}
			suspendStoreSelection = false;
		}
	    else
		{   quant.setEnabled(false);
			cmpnt.setEnabled(false);
		}
	}
	
	// called when state changed in quantity menu
	public void setComponentMenu()
	{
		PlotMenuItem pm=(PlotMenuItem)quant.getSelectedItem();
		if(pm==null) return;
		
		int numItems = cmpnt.getItemCount();
		switch(pm.getTag())
		{   case MPMSIGMAX:
			case MPMEPSX:
			case MPMPLEPSX:
			case MPMEPSTOTX:
				if(docCtrl.resDoc.is3D())
				{	if(numItems!=6 || !cmpnt.getItemAt(0).equals(xchar+xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar+xchar);
						cmpnt.addItem(ychar+ychar);
						cmpnt.addItem(xchar+ychar);
						cmpnt.addItem(zchar+zchar);
						cmpnt.addItem(xchar+zchar);
						cmpnt.addItem(ychar+zchar);
					}
				}
				else
				{	if(numItems!=4 || !cmpnt.getItemAt(0).equals(xchar+xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar+xchar);
						cmpnt.addItem(ychar+ychar);
						cmpnt.addItem(xchar+ychar);
						cmpnt.addItem(zchar+zchar);
					}
				}
				cmpnt.setEnabled(true);
				break;
				
			case MESHSIGMAX:
			case MESHSTRAINX:
			case MESHELEMSIGMAX:
				if(numItems!=4 || !cmpnt.getItemAt(0).equals(xchar+xchar))
				{	cmpnt.removeAllItems();
					cmpnt.addItem(xchar+xchar);
					cmpnt.addItem(ychar+ychar);
					cmpnt.addItem(xchar+ychar);
					cmpnt.addItem(zchar+zchar);
				}
				cmpnt.setEnabled(true);
				break;
			
			case MPMVELX:
			case MPMDISPX:
				if(docCtrl.resDoc.is3D())
				{	if(numItems!=3 || !cmpnt.getItemAt(0).equals(xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar);
						cmpnt.addItem(ychar);
						cmpnt.addItem(zchar);
					}
				}
				else
				{	if(numItems!=2 || !cmpnt.getItemAt(0).equals(xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar);
						cmpnt.addItem(ychar);
					}
				}
				cmpnt.setEnabled(true);
				break;
				
			case MPMSPINVELOCITYX:
			case MPMSPINMOMENTUMX:
				if(docCtrl.resDoc.is3D())
				{	if(numItems!=3 || !cmpnt.getItemAt(0).equals(xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar);
						cmpnt.addItem(ychar);
						cmpnt.addItem(zchar);
					}
				}
				else
				{	if(numItems!=1 || !cmpnt.getItemAt(0).equals(zchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(zchar);
					}
				}
				cmpnt.setEnabled(true);
				break;
				
			case MESHDISPX:
			case MESHFORCEX:
				if(numItems!=2 || !cmpnt.getItemAt(0).equals(xchar))
				{	cmpnt.removeAllItems();
					cmpnt.addItem(xchar);
					cmpnt.addItem(ychar);
				}
				cmpnt.setEnabled(true);
				break;
			
				
			case MPMDCDX:
				if(docCtrl.resDoc.is3D())
				{	if(docCtrl.resDoc.hasPorePressure)
					{	if(numItems!=3 || !cmpnt.getItemAt(0).equals("dp/d"+xchar))
						{	cmpnt.removeAllItems();
							cmpnt.addItem("dp/d"+xchar);
							cmpnt.addItem("dp/d"+ychar);
							cmpnt.addItem("dp/d"+zchar);
						}
					}
					else
					{	if(numItems!=3 || !cmpnt.getItemAt(0).equals("dc/d"+xchar))
						{	cmpnt.removeAllItems();
							cmpnt.addItem("dc/d"+xchar);
							cmpnt.addItem("dc/d"+ychar);
							cmpnt.addItem("dc/d"+zchar);
						}
					}
				}
				else
				{	if(docCtrl.resDoc.hasPorePressure)
					{	if(numItems!=2 || !cmpnt.getItemAt(0).equals("dp/d"+xchar))
						{	cmpnt.removeAllItems();
							cmpnt.addItem("dp/d"+xchar);
							cmpnt.addItem("dp/d"+ychar);
						}
					}
					else
					{	if(numItems!=2 || !cmpnt.getItemAt(0).equals("dc/d"+xchar))
						{	cmpnt.removeAllItems();
							cmpnt.addItem("dc/d"+xchar);
							cmpnt.addItem("dc/d"+ychar);
						}
					}
				}
				cmpnt.setEnabled(true);
				break;
				
			case MESHDVDX:
				if(numItems!=2 || !cmpnt.getItemAt(0).equals("dv/d"+xchar))
				{	cmpnt.removeAllItems();
					cmpnt.addItem("dv/d"+xchar);
					cmpnt.addItem("du/d"+ychar);
				}
				cmpnt.setEnabled(true);
				break;
			
			case INTERFACETRACTION_N:
				if(numItems!=2 || !cmpnt.getItemAt(0).equals("normal"))
				{	cmpnt.removeAllItems();
					cmpnt.addItem("normal");
					cmpnt.addItem("tangential");
				}
				cmpnt.setEnabled(true);
				break;
				
			default:
				cmpnt.setEnabled(false);
				break;
		}
	}
	
	// get plot component from current selection
	public int getPlotComponent(int selected)
	{	PlotMenuItem pm=(PlotMenuItem)quant.getSelectedItem();
		if(pm==null) return -1;
		int plotComponent=pm.getTag();
		
		// adjust component menus - add component selected in component menu
		int extra = 0;
		switch(plotComponent)
		{   case MPMSIGMAX:
			case MPMEPSX:
			case MPMPLEPSX:
			case MPMEPSTOTX:
			case MPMVELX:
			case MPMDCDX:
			case MPMDISPX:
			case MESHSIGMAX:
			case MESHDISPX:
			case MESHSTRAINX:
			case MESHDVDX:
			case MESHELEMSIGMAX:
			case MESHFORCEX:
			case INTERFACETRACTION_N:
				extra=cmpnt.getSelectedIndex();
				break;
			case PlotQuantity.MPMSPINVELOCITYX:
			case PlotQuantity.MPMSPINMOMENTUMX:
				if(docCtrl.resDoc.is3D())
					extra=cmpnt.getSelectedIndex();
				else
					extra = 2;	// only z in 2D
				break;
			default:
				break;
		}
		
		// save selection
		if(plotComponent>=0 && !suspendStoreSelection)
		{	if(docCtrl.resDoc.isMPMAnalysis())
			{	if(selected==LoadArchive.TIME_PLOT)
				{	mpmTimeQuant=plotComponent;
					mpmTimeComp=extra;
				}
				else if(selected==LoadArchive.MESH2D_PLOT)
				{	mpmMesh2DQuant=plotComponent;
					mpmMesh2DComp=extra;
				}
				else 
				{	mpmMovieQuant=plotComponent;
					mpmMovieComp=extra;
				}
			}
			else
			{	feaMeshQuant=plotComponent;
				feaMeshComp=extra;
			}
		}
		
		return plotComponent+extra;
	}
	
	// Label for plot axes (i.e., just generic name and plot units)
	public static String plotLabel(int component,ResultsDocument resDoc)
	{
		JNUnits units = resDoc.units;
		switch(component)
		{	// Stresses
			case MPMSIGMAX:
			case MPMSIGMAY:
			case MPMSIGMAXY:
			case MPMSIGMAZ:
			case MPMSIGMAXZ:
			case MPMSIGMAYZ:
			case MESHSIGMAX:
			case MESHSIGMAY:
			case MESHSIGMAXY:
			case MESHSIGMAZ:
			case MESHELEMSIGMAX:
			case MESHELEMSIGMAY:
			case MESHELEMSIGMAXY:
			case MESHELEMSIGMAZ:
			case MPMPRESSURE:
			case MPMEQUIVSTRESS:
			case MESHPRESSURE:
				return "Stress ("+units.stressUnits()+")";
		
			// Strains
			case MPMEPSX:
			case MPMEPSY:
			case MPMEPSXY:
			case MPMEPSZ:
			case MPMEPSXZ:
			case MPMEPSYZ:
			case MPMPLEPSX:
			case MPMPLEPSY:
			case MPMPLEPSXY:
			case MPMPLEPSZ:
			case MPMPLEPSXZ:
			case MPMPLEPSYZ:
			case MPMEPSTOTX:
			case MPMEPSTOTY:
			case MPMEPSTOTXY:
			case MPMEPSTOTZ:
			case MPMEPSTOTXZ:
			case MPMEPSTOTYZ:
			case MESHSTRAINX:
			case MESHSTRAINY:
			case MESHSTRAINXY:
			case MESHSTRAINZ:
			case MPMEQUIVSTRAIN:
				return "Strain ("+units.strainUnits()+")";
			
			// Energy Density
			case MPMENERGY:
			case MPMSTRENERGY:
			case MPMKINENERGY:
			case MPMWORKENERGY:
			case MPMPLASTICENERGY:
			case MPMHEATENERGY:
				return "Energy Density ("+units.energyUnits()+"/"+units.lengthUnits()+"^3)";

			// Energy
			case MPMTOTSTRENERGY:
			case MPMTOTKINENERGY:
			case MPMTOTENERGY:
			case MPMTOTWORKENERGY:
			case MPMTOTPLASTICENERGY:
			case MPMTOTHEATENERGY:
			case MESHSTRAINENERGY:
				return "Energy ("+units.energyUnits()+")";
			
			case MPMTEMPERATURE:
				return "Temperature (K)";
			
			// Velocity
			case MPMVELVEC:
			case MPMVELX:
			case MPMVELY:
			case MPMVELZ:
				return "Velocity ("+units.velocityUnits()+")";
				
			// Spin velocity
			case MPMSPINVELOCITYX:
			case MPMSPINVELOCITYY:
			case MPMSPINVELOCITYZ:
				return "Ang. Velocity (1/"+units.timeUnits()+")";
			
			// Spin momentum
			case MPMSPINMOMENTUMX:
			case MPMSPINMOMENTUMY:
			case MPMSPINMOMENTUMZ:
				return "Ang. Momentum ("+units.energyUnits()+"-"+units.timeUnits()+")";
			
			// Displacements
			case MPMDISPX:
			case MPMDISPY:
			case MPMDISPZ:
			case MPMNORMALCTOD:
			case MPMSHEARCTOD:
			case MPMDEBONDNCTOD:
			case MPMDEBONDSCTOD:
			case MESHDISPX:
			case MESHDISPY:
				return "Displacement ("+units.lengthUnits()+")";
			
			// Position or material
			case MESHANGLE:
			case MPMANGLEZ:
			case MPMANGLEY:
			case MPMANGLEX:
				return "Material Angle";
				
			// Position or material
			//case MESHMATERIAL:
			case MPMPOS:
				return "Material";
					 
			// shear components
			case MESHDVDX:
			case MESHDUDY:
				return "Shear Component ("+units.strainUnits()+")";
			
			// histrory variables
			case MPMHISTORY1:
			case MPMHISTORY2:
			case MPMHISTORY3:
			case MPMHISTORY4:
				return "Material History";
			
			// concentration
			case MPMCONCENTRATION:
				if(resDoc.hasPorePressure)
					return "Pore Pressure ("+units.stressUnits()+")";
				else
					return "Concentration";
			
			// concentration gradient
			case MPMDCDX:
			case MPMDCDY:
			case MPMDCDZ:
				if(resDoc.hasPorePressure)
					return "Pore Press Grad ("+units.stressUnits()+"/"+units.lengthUnits()+")";
				else
					return "Conc Grad (1/"+units.lengthUnits()+")";
			
			// an Expression
			case MPMEXPRESSION:
			case FEAEXPRESSION:
				return "?";
			
			case MPMJ1:
			case MPMJ2:
				return "J Integral (J/m^2)";
			
			case MPMKI:
			case MPMKII:
				return "Stress Intensity (MPa m^0.5)";
			
			case MPMLENGTH:
			case MPMDEBONDLENGTH:
				return "Length ("+units.lengthUnits()+")";
			
			case MPMCRACKRELEASE:
			case MPMCRACKABSORB:
				return "Energy Rate (J/m^2)";
			
			case MESHFORCEX:
			case MESHFORCEY:
				return "Force ("+units.forceUnits()+")";
			
			case INTERFACETRACTION_N:
			case INTERFACETRACTION_T:
				return "Interface Traction ("+units.stressUnits()+")";
				
			case MPMPOSY:
			case MPMPOSX:
			case MPMPOSZ:
			case MPMCRACKPROFILE:
				return "Position ("+units.lengthUnits()+")";
			
			case MPMOPENINGFRACTION:
			case MPMSHEARFRACTION:
				return "Fraction";
			
			case MPMMASS:
				return "Density ("+units.massUnits()+"/"+units.lengthUnits()+"^3)";
			
			case MPMARCHIVETIME:
				return "Time ("+units.timeUnits()+")";
			
			case MPMGLOBALRESULTS:
				return "Global Results";
			
			case IMPORTANDPLOTFILE:
				return "Plot Data";
			
			case MPMELEMENTCROSSINGS:
			case MPMTOTELEMENTCROSSINGS:
				return "Element Crossings";
			
			// Unknown
			default:
				break;
		}
		
		return "Y Axis";
	}
	
	// Units for plot quantity
	public static String plotUnits(int component,ResultsDocument resDoc)
	{
		JNUnits units = resDoc.units;
		switch(component)
		{	// Stresses
			case MPMSIGMAX:
			case MPMSIGMAY:
			case MPMSIGMAXY:
			case MPMSIGMAZ:
			case MPMSIGMAXZ:
			case MPMSIGMAYZ:
			case MESHSIGMAX:
			case MESHSIGMAY:
			case MESHSIGMAXY:
			case MESHSIGMAZ:
			case MESHELEMSIGMAX:
			case MESHELEMSIGMAY:
			case MESHELEMSIGMAXY:
			case MESHELEMSIGMAZ:
			case INTERFACETRACTION_N:
			case INTERFACETRACTION_T:
				return units.stressUnits();
		
			// Strains
			case MPMEPSX:
			case MPMEPSY:
			case MPMEPSXY:
			case MPMEPSZ:
			case MPMEPSXZ:
			case MPMEPSYZ:
			case MPMPLEPSX:
			case MPMPLEPSY:
			case MPMPLEPSXY:
			case MPMPLEPSZ:
			case MPMPLEPSXZ:
			case MPMPLEPSYZ:
			case MPMEPSTOTX:
			case MPMEPSTOTY:
			case MPMEPSTOTXY:
			case MPMEPSTOTZ:
			case MPMEPSTOTXZ:
			case MPMEPSTOTYZ:
			case MESHSTRAINX:
			case MESHSTRAINY:
			case MESHSTRAINXY:
			case MESHSTRAINZ:
			case MESHDVDX:
			case MESHDUDY:
				return units.strainUnits();
			
			// Energy Density
			case MPMENERGY:
			case MPMSTRENERGY:
			case MPMKINENERGY:
			case MPMWORKENERGY:
			case MPMPLASTICENERGY:
			case MPMHEATENERGY:
				return units.energyUnits()+"/"+units.lengthUnits()+"^3";
				
			// Energy
			case MPMTOTSTRENERGY:
			case MPMTOTKINENERGY:
			case MPMTOTENERGY:
			case MPMTOTWORKENERGY:
			case MPMTOTPLASTICENERGY:
			case MPMTOTHEATENERGY:
			case MESHSTRAINENERGY:
				return units.energyUnits();
			
			case MPMTEMPERATURE:
				return "K";
			
			// Velocity
			case MPMVELVEC:
			case MPMVELX:
			case MPMVELY:
			case MPMVELZ:
				return units.velocityUnits();
			
			// Spin velocity
			case MPMSPINVELOCITYX:
			case MPMSPINVELOCITYY:
			case MPMSPINVELOCITYZ:
				return "1/"+units.timeUnits();
			
			// Spin momentum
			case MPMSPINMOMENTUMX:
			case MPMSPINMOMENTUMY:
			case MPMSPINMOMENTUMZ:
				return units.energyUnits()+"-"+units.timeUnits();
				
			// Displacements
			case MPMDISPX:
			case MPMDISPY:
			case MPMDISPZ:
			case MPMNORMALCTOD:
			case MPMSHEARCTOD:
			case MPMDEBONDNCTOD:
			case MPMDEBONDSCTOD:
			case MESHDISPX:
			case MESHDISPY:
			case MPMLENGTH:
			case MPMDEBONDLENGTH:
			case MPMPOSY:
			case MPMPOSX:
			case MPMPOSZ:
			case MPMCRACKPROFILE:
				return units.lengthUnits();
			
			case MPMJ1:
			case MPMJ2:
			case MPMCRACKRELEASE:
			case MPMCRACKABSORB:
				return "J/m^2";
			
			case MPMKI:
			case MPMKII:
				return "MPa m^0.5";
			
			// concentration gradient
			case MPMDCDX:
			case MPMDCDY:
			case MPMDCDZ:
				if(resDoc.hasPorePressure)
					return units.stressUnits()+"/"+units.lengthUnits();
				else
					return "1/"+units.lengthUnits();
			
			// concentration and pore pressure
			case MPMCONCENTRATION:
				if(resDoc.hasPorePressure)
					return units.stressUnits();
				else
					return "";
				
			case MESHFORCEX:
			case MESHFORCEY:
				return units.forceUnits();
			
			case MPMMASS:
				return units.massUnits()+"/"+units.lengthUnits()+"^3";
			
			case MPMARCHIVETIME:
				return units.timeUnits();
			
			// rest of no (or unknown) units
			default:
				break;
		}
		
		return "";
	}

	
	// select item for checking the mesh
	public void setCheckMeshItem() { quant.setSelectedIndex(checkMeshItem); }
	
	// name for plot component
	public static String plotName(int component,ResultsDocument resDoc)
	{
		char xc='x', yc='y', zc='z';
		if(resDoc!=null)
		{	if(resDoc.isAxisymmetric())
			{	xc = 'r';
				yc = 'z';
				zc = 't';
			}
		}
		switch(component)
		{	case MPMSIGMAX:
			case MESHSIGMAX:
			case MESHELEMSIGMAX:
				return "Stress "+xc+xc;
				
			case MPMSIGMAY:
			case MESHSIGMAY:
			case MESHELEMSIGMAY:
				return "Stress "+yc+yc;
				
			case MPMSIGMAXY:
			case MESHSIGMAXY:
			case MESHELEMSIGMAXY:
				return "Stress "+xc+yc;
				
			case MPMSIGMAZ:
			case MESHSIGMAZ:
			case MESHELEMSIGMAZ:
				return "Stress "+zc+zc;
				
			case MPMSIGMAXZ:
				return "Stress "+xc+zc;
				
			case MPMSIGMAYZ:
				return "Stress "+yc+zc;
			
			case MPMPRESSURE:
			case MESHPRESSURE:
				return "Pressure";
			
			case MPMEQUIVSTRESS:
				return "Equiv Stress";
				
			case MPMEPSX:
			case MPMEPSTOTX:
			case MESHSTRAINX:
				return "Strain "+xc+xc;
				
			case MPMEPSY:
			case MPMEPSTOTY:
			case MESHSTRAINY:
				return "Strain "+yc+yc;
				
			case MPMEPSXY:
			case MPMEPSTOTXY:
			case MESHSTRAINXY:
				return "Strain "+xc+yc;
				
			case MPMEPSZ:
			case MPMEPSTOTZ:
			case MESHSTRAINZ:
				return "Strain "+zc+zc;
				
			case MPMEPSXZ:
			case MPMEPSTOTXZ:
				return "Strain "+xc+zc;

			case MPMEPSYZ:
			case MPMEPSTOTYZ:
				return "Strain "+yc+zc;

			case MPMEQUIVSTRAIN:
				return "Equiv Strain";
				
			case MPMPLEPSX:
				return "Plastic Strain "+xc+xc;
				
			case MPMPLEPSY:
				return "Plastic Strain "+yc+yc;
				
			case MPMPLEPSXY:
				return "Plastic Strain "+xc+yc;
				
			case MPMPLEPSZ:
				return "Plastic Strain "+zc+zc;

			case MPMPLEPSXZ:
				return "Plastic Strain "+xc+zc;

			case MPMPLEPSYZ:
				return "Plastic Strain "+yc+zc;

			case MPMENERGY:
				return "Energy";

			case MPMSTRENERGY:
				return "Strain Energy Density";
				
			case MESHSTRAINENERGY:
				return "Strain Energy";

			case MPMKINENERGY:
				return "Kinetic Energy Density";

			case MPMWORKENERGY:
				return "Work Energy Density";

			case MPMPLASTICENERGY:
				return "Plastic Energy Density";

			case MPMHEATENERGY:
				return "Heat Energy Density";

			case MPMTOTSTRENERGY:
				return "Total Strain Energy";
				
			case MPMTOTKINENERGY:
				return "Total Kinetic Energy";
				
			case MPMTOTENERGY:
				return "Total Energy";
				
			case MPMTOTWORKENERGY:
				return "Total Work Work";
				
			case MPMTOTPLASTICENERGY:
				return "Total Plastic Energy";
						
			case MPMTOTHEATENERGY:
				return "Total Thermal Energy";
						
			case MPMTEMPERATURE:
				return "Temperature";
			
			case MPMVELX:
				return "Velocity "+xc;
			
			case MPMVELY:
				return "Velocity "+yc;
			
			case MPMVELZ:
				return "Velocity "+zc;
			
			case MPMSPINVELOCITYX:
				return "Ang. Velocity "+xc;
			
			case MPMSPINVELOCITYY:
				return "Ang. Velocity "+yc;
			
			case MPMSPINVELOCITYZ:
				return "Ang. Velocity "+zc;
			
			case MPMSPINMOMENTUMX:
				return "Ang. Momentum "+xc;
			
			case MPMSPINMOMENTUMY:
				return "Ang. Momentum "+yc;
			
			case MPMSPINMOMENTUMZ:
				return "Ang. Momentum "+zc;
			
			case MPMVELVEC:
				return "Velocity";
			
			case MPMDISPX:
			case MESHDISPX:
				return "Displacment "+xc;
			
			case MPMDISPY:
			case MESHDISPY:
				return "Displacment "+yc;
			
			case MPMDISPZ:
				return "Displacment "+zc;
			
			case MPMPOS:
				return "Position";
			
			case MPMPOSX:
				return "Position "+xc;
			
			case MPMPOSY:
				return "Position "+yc;
			
			case MPMPOSZ:
				return "Position "+zc;
			
			case MESHDUDY:
				return "Strain du/d"+yc;
			
			case MESHDVDX:
				return "Strain dv/d"+xc;
			
			case MPMHISTORY1:
			case MPMHISTORY2:
			case MPMHISTORY3:
			case MPMHISTORY4:
				return "Material History";
			
			case MPMDCDY:
				if(resDoc.hasPorePressure)
					return "Pore Press Grad dp/d"+yc;
				else
					return "Conc Grad dc/d"+yc;
			
			case MPMDCDX:
				if(resDoc.hasPorePressure)
					return "Pore Press Grad dp/d"+xc;
				else
					return "Conc Grad dc/d"+xc;
			
			case MPMDCDZ:
				if(resDoc.hasPorePressure)
					return "Pore Press Grad dp/d"+zc;
				else
					return "Conc Grad dc/d"+zc;
			
			case MPMCONCENTRATION:
				if(resDoc.hasPorePressure)
					return "Pore Pressure";
				else
					return "Concentration";
			
			case MPMEXPRESSION:
			case FEAEXPRESSION:
				return "?";
			
			case MPMJ1:
				return "J1 Integral";
			
			case MPMJ2:
				return "J2 Integral";
			
			case MPMKI:
				return "KI";
			
			case MPMKII:
				return "KII";
			
			case MPMLENGTH:
				return "Crack Length";
			
			case MPMDEBONDLENGTH:
				return "Debonded Crack Length";
			
			case MPMNORMALCTOD:
				return "Normal CTOD";
			
			case MPMSHEARCTOD:
				return "Shear CTOD";
			
			case MPMDEBONDNCTOD:
				return "Debond Tip Normal COD";
			
			case MPMDEBONDSCTOD:
				return "Debond Tip Shear COD";
			
			case MPMCRACKRELEASE:
				return "Energy Release Rate";
			
			case MPMCRACKABSORB:
				return "Energy Dissipation Rate";
			
			case MESHMATERIAL:
				return "Material";
			
			case MESHANGLE:
			case MPMANGLEZ:
			case MPMANGLEY:
			case MPMANGLEX:
				return "Material Angle";
			
			case MESHFORCEX:
				return "Force "+xc;
				
			case MESHFORCEY:
				return "Force "+yc;
			
			case INTERFACETRACTION_N:
				return "Interface Normal Traction";
				
			case INTERFACETRACTION_T:
				return "Interface Tangential Traction";
				
			case MPMMASS:
				return "Density";
			
			case MPMARCHIVETIME:
				return "Time";
			
			case MPMCRACKPROFILE:
				return "Crack Profile";
			
			case MPMOPENINGFRACTION:
				return "Opening Mode Fraction";
				
			case MPMSHEARFRACTION:
				return "Shear Mode Fraction";
				
			case MPMELEMENTCROSSINGS:
				return "Element Crossings";
				
			case MPMTOTELEMENTCROSSINGS:
				return "Total Element Crossings";
				
			default:
				break;
		}
		return "Unknown";
	}
}
