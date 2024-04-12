/*******************************************************************
	PlotQuantity.java
	NairnFEAMPMViz

	Created by John Nairn on Wed Mar 10 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.*;

import geditcom.JNFramework.JNDialog;
import geditcom.JNFramework.JNExpression;
import geditcom.JNFramework.JNUtilities;

public class PlotQuantity extends PlotControl
{
	static final long serialVersionUID=18L;
	
	// MPM plotting options
	static final int SHIFT4MAGNITUDE=10000;
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
	static final int MPMVELS=MPMVELX+SHIFT4MAGNITUDE;
	static final int MPMVELVEC=22;				// extra in Java app
	
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
	static final int MPMDISPS=MPMDISPX+SHIFT4MAGNITUDE;
	
	static final int MPMPOSX=33;
	static final int MPMPOSY=34;
	static final int MPMPOSZ=35;
	
	static final int PLOTEXPRESSION=36;
	//static final int MPMEXPRESSION=36;
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
	static final int MPMDENSITY=51;
	static final int MPMARCHIVETIME=52;
	static final int MPMGLOBALRESULTS=53;
	static final int MPMCONCENTRATION=54;
	static final int MPMDCDX=55;
	static final int MPMDCDY=56;
	static final int MPMDCDZ=57;
	static final int MPMMODEIFB=58;
	static final int MPMMODEIIFB=59;
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
	
	static final int MPMPOSDISTANCE=72;				// expression only
	static final int MPMPOSANGLE=73;
	static final int MPMORIGPOSX=74;
	static final int MPMORIGPOSY=75;
	static final int MPMORIGPOSZ=76;
	
	static final int MPMDEBONDLENGTH=77;
	static final int MPMDEBONDNCTOD=78;
	static final int MPMDEBONDSCTOD=79;
	static final int MPMHISTORY1=80;
	static final int MPMHISTORY2=81;
	static final int MPMHISTORY3=82;
	static final int MPMHISTORY4=83;
	
	static final int MPMEPSTOTX=87;
	static final int MPMEPSTOTY=88;
	static final int MPMEPSTOTXY=89;
	static final int MPMEPSTOTZ=90;
	static final int MPMEPSTOTXZ=91;
	static final int MPMEPSTOTYZ=92;
	
	static final int MPMWXY=93;
	static final int MPMWXZ=94;
	static final int MPMWYZ=95;
	
	static final int MPMMAXSTRESS=102;
	static final int MPMMINSTRESS=103;
	static final int MPMSTRESSDIR=104;
	
	static final int MPMEQUIVSTRESS=105;
	static final int MPMEQUIVSTRAIN=106;
	static final int MPMPRESSURE=107;

	static final int MPMSPINMOMENTUMX=110;
	static final int MPMSPINMOMENTUMY=111;
	static final int MPMSPINMOMENTUMZ=112;
	static final int MPMSPINVELOCITYX=113;
	static final int MPMSPINVELOCITYY=114;
	static final int MPMSPINVELOCITYZ=115;
	
	static final int MPMHISTORY5=116;
	static final int MPMHISTORY6=117;
	static final int MPMHISTORY7=118;
	static final int MPMHISTORY8=119;
	static final int MPMHISTORY9=120;
	static final int MPMHISTORY10=121;
	static final int MPMHISTORY11=122;
	static final int MPMHISTORY12=123;
	static final int MPMHISTORY13=124;
	static final int MPMHISTORY14=125;
	static final int MPMHISTORY15=126;
	static final int MPMHISTORY16=127;
	static final int MPMHISTORY17=128;
	static final int MPMHISTORY18=129;
	static final int MPMHISTORY19=130;
	static final int MPMTRACTION1=131;
	static final int MPMTRACTION2=132;
	static final int MPMTRACTION3=133;
	static final int MPMTRACTION4=134;
	static final int MPMTRACTION5=135;
	static final int MPMTRACTION6=136;
	static final int MPMTRACTION7=137;
	static final int MPMTRACTION8=138;
	static final int MPMTRACTION9=139;
	static final int MPMTRACTION10=140;
	static final int MPMCZMGI=141;
	static final int MPMCZMGII=142;
	static final int MPMCZLENGTH=143;

	static final int MPMMASS=144;
	static final int MPMSIZE=145;
	static final int MPMMAXCODE=146;

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
	//static final int FEAEXPRESSION=1022;
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
	private String totalchar="magnitude";
	
	// plot expression
	private String expression = "0";
	private static String expressionTitle = null;
	private JNExpression expr = new JNExpression("0");
	private HashMap<String,String> exprVars = new HashMap<String,String>(1);
	private ArrayList<String> keys = new ArrayList<String>(10);
	private int[] keyCode = new int[MPMMAXCODE];
	
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
		quant.setFocusable(false);
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
		cmpnt.setFocusable(false);
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
					quant.addItem(new PlotMenuItem("Max Principal Stress",MPMMAXSTRESS));
					quant.addItem(new PlotMenuItem("Min Principal Stress",MPMMINSTRESS));
					quant.addItem(new PlotMenuItem("Max Stress Angle",MPMSTRESSDIR));
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
				
				quant.addItem(new PlotMenuItem("Expression...",PLOTEXPRESSION));
					
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
					
				if(arch[ReadArchive.ARCH_Velocity]=='Y') {
					quant.addItem(new PlotMenuItem("Velocity",MPMVELX));
				}
					
				if(arch[ReadArchive.ARCH_SpinVelocity]=='Y')
					quant.addItem(new PlotMenuItem("Angular Velocity",MPMSPINVELOCITYX));
				if(arch[ReadArchive.ARCH_SpinMomentum]=='Y')
					quant.addItem(new PlotMenuItem("Angular Momentum",MPMSPINMOMENTUMX));
				quant.addItem(new PlotMenuItem("Displacement",MPMDISPX));
				checkMeshItem=quant.getItemCount();
				quant.addItem(new PlotMenuItem("Material",MPMPOS));
				quant.addItem(new PlotMenuItem("Rotational Strain",MPMANGLEZ));
				quant.addItem(new PlotMenuItem("Density",MPMDENSITY));
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
				
				if(arch[ReadArchive.ARCH_History59]=='Y')
					quant.addItem(new PlotMenuItem("History 5",MPMHISTORY5));
				else if(arch[ReadArchive.ARCH_History59]!='N')
				{	int history=(int)arch[ReadArchive.ARCH_History59];
					if((history & 0x01) !=0)
						quant.addItem(new PlotMenuItem("History 5",MPMHISTORY5));
					if((history & 0x02) !=0)
						quant.addItem(new PlotMenuItem("History 6",MPMHISTORY6));
					if((history & 0x04) !=0)
						quant.addItem(new PlotMenuItem("History 7",MPMHISTORY7));
					if((history & 0x08) !=0)
						quant.addItem(new PlotMenuItem("History 8",MPMHISTORY8));
					if((history & 0x10) !=0)
						quant.addItem(new PlotMenuItem("History 9",MPMHISTORY9));
				}
				
				if(arch[ReadArchive.ARCH_History1014]=='Y')
					quant.addItem(new PlotMenuItem("History 10",MPMHISTORY10));
				else if(arch[ReadArchive.ARCH_History1014]!='N')
				{	int history=(int)arch[ReadArchive.ARCH_History1014];
					if((history & 0x01) !=0)
						quant.addItem(new PlotMenuItem("History 10",MPMHISTORY10));
					if((history & 0x02) !=0)
						quant.addItem(new PlotMenuItem("History 11",MPMHISTORY11));
					if((history & 0x04) !=0)
						quant.addItem(new PlotMenuItem("History 12",MPMHISTORY12));
					if((history & 0x08) !=0)
						quant.addItem(new PlotMenuItem("History 13",MPMHISTORY13));
					if((history & 0x10) !=0)
						quant.addItem(new PlotMenuItem("History 14",MPMHISTORY14));
				}
				
				if(arch[ReadArchive.ARCH_History1519]=='Y')
					quant.addItem(new PlotMenuItem("History 15",MPMHISTORY15));
				else if(arch[ReadArchive.ARCH_History1519]!='N')
				{	int history=(int)arch[ReadArchive.ARCH_History1519];
					if((history & 0x01) !=0)
						quant.addItem(new PlotMenuItem("History 15",MPMHISTORY15));
					if((history & 0x02) !=0)
						quant.addItem(new PlotMenuItem("History 16",MPMHISTORY16));
					if((history & 0x04) !=0)
						quant.addItem(new PlotMenuItem("History 17",MPMHISTORY17));
					if((history & 0x08) !=0)
						quant.addItem(new PlotMenuItem("History 18",MPMHISTORY18));
					if((history & 0x10) !=0)
						quant.addItem(new PlotMenuItem("History 19",MPMHISTORY19));
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
					
						quant.addItem(new PlotMenuItem("Crack Length",MPMLENGTH));
						quant.addItem(new PlotMenuItem("Debonded Crack Length",MPMDEBONDLENGTH));
						if(carch[ReadArchive.ARCH_CZMDeltaG]=='Y')
						{	quant.addItem(new PlotMenuItem("Cohesive Damage Length",MPMCZLENGTH));
						}
						
						quant.addItem(new PlotMenuItem("Normal CTOD",MPMNORMALCTOD));
						quant.addItem(new PlotMenuItem("Shear CTOD",MPMSHEARCTOD));
						quant.addItem(new PlotMenuItem("Debond Tip Normal COD",MPMDEBONDNCTOD));
						quant.addItem(new PlotMenuItem("Debond Tip Shear COD",MPMDEBONDSCTOD));
						
						if(carch[ReadArchive.ARCH_CZMDeltaG]=='Y')
						{	quant.addItem(new PlotMenuItem("CZM Mode I Force",MPMMODEIFB));
							quant.addItem(new PlotMenuItem("CZM Mode II Force",MPMMODEIIFB));
						}
					}
					
					preselect = mpmTimeQuant;
					preComp = mpmTimeComp;
				}
				
				// additional time plot options
				else if(selected==LoadArchive.MESH2D_PLOT)
				{	byte [] carch=docCtrl.resDoc.crackFormat.getBytes();
				
					// x-y crack results
					quant.addItem(new PlotMenuItem("Crack Profile",MPMCRACKPROFILE));
					quant.addItem(new PlotMenuItem("Crack Normal CTOD",MPMNORMALCTOD));
					quant.addItem(new PlotMenuItem("Crack Tangential CTOD",MPMSHEARCTOD));
					quant.addItem(new PlotMenuItem("Crack Opening Fraction",MPMOPENINGFRACTION));
					quant.addItem(new PlotMenuItem("Crack Sliding Fraction",MPMSHEARFRACTION));
					if(carch[ReadArchive.ARCH_CZMDeltaG]=='Y')
					{	quant.addItem(new PlotMenuItem("CZM GI",MPMCZMGI));
						quant.addItem(new PlotMenuItem("CZM GII",MPMCZMGII));
					}
					
					char histChar = (char)carch[ReadArchive.ARCH_Traction15];
					if(histChar=='Y')
						quant.addItem(new PlotMenuItem("Traction 1",MPMTRACTION1));
					else if(histChar!='N')
					{	int history=(int)histChar;
						if((history & 0x01) !=0)
							quant.addItem(new PlotMenuItem("Traction 1",MPMTRACTION1));
						if((history & 0x02) !=0)
							quant.addItem(new PlotMenuItem("Traction 2",MPMTRACTION2));
						if((history & 0x04) !=0)
							quant.addItem(new PlotMenuItem("Traction 3",MPMTRACTION3));
						if((history & 0x08) !=0)
							quant.addItem(new PlotMenuItem("Traction 4",MPMTRACTION4));
						if((history & 0x10) !=0)
							quant.addItem(new PlotMenuItem("Traction 5",MPMTRACTION5));
					}
					histChar = (char)carch[ReadArchive.ARCH_Traction610];
					if(histChar=='Y')
						quant.addItem(new PlotMenuItem("Traction 6",MPMTRACTION6));
					else if(histChar!='N')
					{	int history=(int)histChar;
						if((history & 0x01) !=0)
							quant.addItem(new PlotMenuItem("Traction 6",MPMTRACTION6));
						if((history & 0x02) !=0)
							quant.addItem(new PlotMenuItem("Traction 7",MPMTRACTION7));
						if((history & 0x04) !=0)
							quant.addItem(new PlotMenuItem("Traction 8",MPMTRACTION8));
						if((history & 0x08) !=0)
							quant.addItem(new PlotMenuItem("Traction 9",MPMTRACTION9));
						if((history & 0x10) !=0)
							quant.addItem(new PlotMenuItem("Traction 10",MPMTRACTION10));
					}
					
					preselect = mpmMesh2DQuant;
					preComp = mpmMesh2DComp;
				}
				
				else
				{	preselect = mpmMovieQuant;
					preComp = mpmMovieComp;
				}
				
				// import and plot
				if(selected==LoadArchive.TIME_PLOT || selected==LoadArchive.MESH2D_PLOT)
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
				
				// an expression
				quant.addItem(new PlotMenuItem("Expression...",PLOTEXPRESSION));
				
				// import and plot
				if(selected==LoadArchive.MESH2D_PLOT)
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
	
	// convert menutext to plot quantity ID
	// but not "Import"
	public int menuTextToPlotTag(String menutext)
	{
		int cmpnt = -1;
		
		if(docCtrl.resDoc.isMPMAnalysis())
		{	// MPM plots
			if(menutext.equals("Stress"))
				cmpnt = MPMSIGMAX;
			else if(menutext.equals("Pressure"))
				cmpnt = MPMPRESSURE;
			else if(menutext.equals("Equiv. Stress") || menutext.equals("Equivalent Stress"))
				cmpnt = MPMEQUIVSTRESS;
			else if(menutext.equals("Max Principal Stress"))
				cmpnt = MPMMAXSTRESS;
			else if(menutext.equals("Min Principal Stress"))
				cmpnt = MPMMINSTRESS;
			else if(menutext.equals("Max Stress Angle"))
				cmpnt = MPMSTRESSDIR;
			else if(menutext.equals("Total Strain"))
				cmpnt = MPMEPSTOTX;
			else if(menutext.equals("Elastic Strain"))
				cmpnt = MPMEPSX;
			else if(menutext.equals("Plastic Strain"))
				cmpnt = MPMPLEPSX;
			else if(menutext.equals("Equiv. Strain") || menutext.equals("Equivalent Strain"))
				cmpnt = MPMEQUIVSTRAIN;
			else if(menutext.equals("Expression...") || menutext.equals("Expression"))
				cmpnt = PLOTEXPRESSION;
			else if(menutext.equals("Strain Energy"))
				cmpnt = MPMSTRENERGY;
			else if(menutext.equals("Kinetic Energy"))
				cmpnt = MPMKINENERGY;
			else if(menutext.equals("Energy"))
				cmpnt = MPMENERGY;
			else if(menutext.equals("Plastic Energy"))
				cmpnt = MPMPLASTICENERGY;
			else if(menutext.equals("Work Energy"))
				cmpnt = MPMWORKENERGY;
			else if(menutext.equals("Heat Energy"))
				cmpnt = MPMHEATENERGY;
			else if(menutext.equals("Velocity"))
				cmpnt = MPMVELX;
			else if(menutext.equals("Angular Velocity"))
				cmpnt = MPMSPINVELOCITYX;
			else if(menutext.equals("Angular Momentum"))
				cmpnt = MPMSPINMOMENTUMX;
			else if(menutext.equals("Displacement"))
				cmpnt = MPMDISPX;
			else if(menutext.equals("Material"))
				cmpnt = MPMPOS;
			else if(menutext.equals("Rotational Strain"))
				cmpnt = MPMANGLEZ;
			else if(menutext.equals("Density"))
				cmpnt = MPMDENSITY;
			else if(menutext.equals("Temperature"))
				cmpnt = MPMTEMPERATURE;
			else if(menutext.equals("Pore Pressure"))
				cmpnt = MPMCONCENTRATION;
			else if(menutext.equals("Pore Press Gradient"))
				cmpnt = MPMDCDX;
			else if(menutext.equals("Concentration"))
				cmpnt = MPMCONCENTRATION;
			else if(menutext.equals("Conc Gradient") || menutext.equals("Conc. Gradient"))
				cmpnt = MPMDCDX;
			else if(menutext.equals("History") || menutext.equals("History Data"))
				cmpnt = MPMHISTORY1;
			else if(menutext.equals("Element Crossings"))
				cmpnt = MPMELEMENTCROSSINGS;
			else if(menutext.equals("Global Results"))
				cmpnt = MPMGLOBALRESULTS;
			else if(menutext.equals("J1"))
				cmpnt = MPMJ1;
			else if(menutext.equals("J2"))
				cmpnt = MPMJ2;
			else if(menutext.equals("KI"))
				cmpnt = MPMKI;
			else if(menutext.equals("KII"))
				cmpnt = MPMKII;
			else if(menutext.equals("CZM Mode I Force"))
				cmpnt = MPMMODEIFB;
			else if(menutext.equals("CZM Mode II Force"))
				cmpnt = MPMMODEIIFB;
			else if(menutext.equals("Cohesive Damage Length"))
				cmpnt = MPMCZLENGTH;
			else if(menutext.equals("Crack Length"))
				cmpnt = MPMLENGTH;
			else if(menutext.equals("Debonded Crack Length"))
				cmpnt = MPMDEBONDLENGTH;
			else if(menutext.equals("Normal CTOD"))
				cmpnt = MPMNORMALCTOD;
			else if(menutext.equals("Shear CTOD"))
				cmpnt = MPMSHEARCTOD;
			else if(menutext.equals("Debond Tip Normal CTOD"))
				cmpnt = MPMDEBONDNCTOD;
			else if(menutext.equals("Debond Shear Normal CTOD"))
				cmpnt = MPMDEBONDSCTOD;
			else if(menutext.equals("Crack Profile"))
				cmpnt = MPMCRACKPROFILE;
			else if(menutext.equals("Crack Normal CTOD"))
				cmpnt = MPMNORMALCTOD;
			else if(menutext.equals("Crack Tangential CTOD"))
				cmpnt = MPMSHEARCTOD;
			else if(menutext.equals("Crack Opening Fraction"))
				cmpnt = MPMOPENINGFRACTION;
			else if(menutext.equals("Crack Sliding Fraction"))
				cmpnt = MPMSHEARFRACTION;
			else if(menutext.equals("CZM GI"))
				cmpnt = MPMCZMGI;
			else if(menutext.equals("CZM GII"))
				cmpnt = MPMCZMGII;
			else if(menutext.equals("Traction") || menutext.equals("Crack Traction Data"))
				cmpnt = MPMTRACTION1;
		}
		else
		{	// FEA plots
			if(menutext.equals("Stress"))
				cmpnt = MESHSIGMAX;
			else if(menutext.equals("Pressure"))
				cmpnt = MESHPRESSURE;
			else if(menutext.equals("Element Stress"))
				cmpnt = MESHELEMSIGMAX;
			else if(menutext.equals("Strain"))
				cmpnt = MESHSTRAINX;
			else if(menutext.equals("Energy"))
				cmpnt = MESHSTRAINENERGY;
			else if(menutext.equals("Element Force"))
				cmpnt = MESHFORCEX;
			else if(menutext.equals("Interface Traction"))
				cmpnt = INTERFACETRACTION_N;
			else if(menutext.equals("Displacement"))
				cmpnt = MESHDISPX;
			else if(menutext.equals("Shear Component"))
				cmpnt = MESHDVDX;
			else if(menutext.equals("Material"))
				cmpnt = MESHMATERIAL;
			else if(menutext.equals("Material Angle"))
				cmpnt = MESHANGLE;
			else if(menutext.equals("Expression...") || menutext.equals("Expression"))
				cmpnt = PLOTEXPRESSION;
		}
		
		return cmpnt;
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
				{	if(numItems!=4 || !cmpnt.getItemAt(0).equals(xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar);
						cmpnt.addItem(ychar);
						cmpnt.addItem(zchar);
						cmpnt.addItem(totalchar);
					}
				}
				else
				{	if(numItems!=3 || !cmpnt.getItemAt(0).equals(xchar))
					{	cmpnt.removeAllItems();
						cmpnt.addItem(xchar);
						cmpnt.addItem(ychar);
						cmpnt.addItem(totalchar);
					}
				}
				cmpnt.setEnabled(true);
				break;
			
			case MPMANGLEZ:
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
	public int getPlotComponent(int selected,boolean getExpression,ISDictType settings)
	{	int vindex=0,tindex=0,plotComponent=-1;
		if(settings==null)
		{	PlotMenuItem pm=(PlotMenuItem)quant.getSelectedItem();
			if(pm==null) return -1;
			plotComponent=pm.getTag();
			vindex=cmpnt.getSelectedIndex();
			tindex=vindex;
		}
		else
		{	String mt = (String)settings.gcis_objectForKey("menutext");
			if(mt==null) return -1;
			plotComponent = menuTextToPlotTag(mt);
			try
			{	vindex = settings.gcis_integerForKey("vectorindex");
			}
			catch(Exception e) {}
			try
			{	tindex = settings.gcis_integerForKey("tensorindex");
			}
			catch(Exception e) {}
		}
		if(plotComponent<0 || vindex<0 || tindex<0) return -1;
		
		// adjust component menus - add component selected in component menu
		int extra = 0;
		switch(plotComponent)
		{   
			case MPMVELX:
			{
				extra=vindex;
				int nextVal = docCtrl.resDoc.is3D()?3:2;
				if(extra == nextVal) extra=SHIFT4MAGNITUDE;
				break;
			}
			case MPMDISPX:
			{
				extra=vindex;
				int nextVal = docCtrl.resDoc.is3D()?3:2;
				if(extra == nextVal) extra=SHIFT4MAGNITUDE;
				break;
			}
			case MPMSIGMAX:
			case MPMEPSX:
			case MPMPLEPSX:
			case MPMEPSTOTX:
			case MESHSIGMAX:
			case MESHSTRAINX:
			case MESHELEMSIGMAX:
				extra=tindex;
				break;
				
			case MPMDCDX:
			case MESHDVDX:
			case MESHDISPX:
			case MESHFORCEX:
			case INTERFACETRACTION_N:
				extra=vindex;
				break;
	
			case MPMSPINVELOCITYX:
			case MPMSPINMOMENTUMX:
				if(docCtrl.resDoc.is3D())
					extra=vindex;
				else
					extra = 2;	// only z in 2D
				break;
			
			case PLOTEXPRESSION:
				if(getExpression) setupForExpression(settings);
				break;
				
			case MPMHISTORY1:
				if(settings!=null)
				{	if(tindex<4)
					{	// tindex 0 to 3 select history 1 to 4
						extra = tindex;
					}
					else
					{	// tindex 4 on selects history variable 5 on
						extra = MPMHISTORY5-MPMHISTORY1+tindex-4;
					}
				}
				break;
				
			case MPMTRACTION1:
				if(settings!=null)
					extra=tindex;
				break;
				
			default:
				break;
		}
		
		// save selection
		if(plotComponent>=0 && !suspendStoreSelection && settings==null)
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
			case MPMMAXSTRESS:
			case MPMMINSTRESS:
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
			case MPMVELS:
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
			case MPMDISPS:
			case MPMNORMALCTOD:
			case MPMSHEARCTOD:
			case MPMDEBONDNCTOD:
			case MPMDEBONDSCTOD:
			case MESHDISPX:
			case MESHDISPY:
				return "Displacement ("+units.lengthUnits()+")";
			
			// Position or material
			case MESHANGLE:
				if(resDoc.isFEAAnalysis()) return "Material Angle";
			case MPMANGLEZ:
			case MPMANGLEY:
			case MPMANGLEX:
				return "Rotational Strain";
				
			// Position or material
			//case MESHMATERIAL:
			case MPMPOS:
				return "Material";
					 
			// shear components
			case MESHDVDX:
			case MESHDUDY:
				return "Shear Component ("+units.strainUnits()+")";
			
			case MPMWXY:
			case MPMWXZ:
			case MPMWYZ:
				return "Rotational Strain ("+units.strainUnits()+")";
				
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
			
			case PLOTEXPRESSION:
				if(expressionTitle!=null) return expressionTitle;
				return "Expression";
			
			case MPMJ1:
			case MPMJ2:
				return "J Integral (J/m^2)";
			
			case MPMKI:
			case MPMKII:
				return "Stress Intensity (MPa m^{0.5})";
			
			case MPMLENGTH:
			case MPMDEBONDLENGTH:
			case MPMCZLENGTH:
				return "Length ("+units.lengthUnits()+")";
			
	        case MPMCZMGI:
	        case MPMCZMGII:
				return "Energy Release Rate ("+units.forceUnits()+"/"+units.lengthUnits()+")";
	            
			case MPMMODEIFB:
			case MPMMODEIIFB:
				return "Energy Release Force ("+units.forceUnits()+")";
			
			case MESHFORCEX:
			case MESHFORCEY:
				return "Force ("+units.forceUnits()+")";
			
			case INTERFACETRACTION_N:
			case INTERFACETRACTION_T:
				return "Interface Traction ("+units.stressUnits()+")";
				
			case MPMPOSY:
			case MPMPOSX:
			case MPMPOSZ:
			case MPMORIGPOSY:
			case MPMORIGPOSX:
			case MPMORIGPOSZ:
			case MPMCRACKPROFILE:
				return "Position ("+units.lengthUnits()+")";
			
			case MPMOPENINGFRACTION:
			case MPMSHEARFRACTION:
				return "Fraction";
			
			case MPMDENSITY:
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
				return "Traction History";
						
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
			case MPMMAXSTRESS:
			case MPMMINSTRESS:
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
			case MPMWXY:
			case MPMWXZ:
			case MPMWYZ:
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
			case MPMVELS:
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
			case MPMDISPS:
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
			case MPMORIGPOSY:
			case MPMORIGPOSX:
			case MPMORIGPOSZ:
			case MPMCRACKPROFILE:
				return units.lengthUnits();
			
			case MPMJ1:
			case MPMJ2:
				return "J/m^2";
				
			case MPMMODEIFB:
			case MPMMODEIIFB:
			case MPMCZLENGTH:
				return units.forceUnits();
			
	        case MPMCZMGI:
	        case MPMCZMGII:
				return units.forceUnits()+"/"+units.lengthUnits();
				
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
			
			case MPMDENSITY:
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
	
	// get expression
	public String getExpression() { return expression; }
	
	// set up for using current epxresison
	public void setupForExpression(ISDictType settings)
	{
		// loop until enter a valid expression
		while(true)
		{	String newExpr;
			expressionTitle = null;
			if(settings==null)
			{	//String newExpr = (String)JOptionPane.showInputDialog(null,msg,"Plot Expression",
				//	JOptionPane.PLAIN_MESSAGE,null,null,expression);
				GetExpression exprDialog = new GetExpression(docCtrl,expression);
				if(exprDialog.getClickedButton()==JNDialog.CANCEL_BUTTON) break;
				
				// get entered expression
				newExpr = exprDialog.getEnteredExpression();
			}
			else
			{	newExpr = (String)settings.gcis_objectForKey("expression");
				if(newExpr==null) newExpr = "0";
				expressionTitle = (String)settings.gcis_objectForKey("expressionTitle");
			}
			
			
			// read variables
			try
			{	JNExpression testExpr = new JNExpression(newExpr,null);
				HashMap<String,String> testVars = testExpr.usedVariables();
				// if no exception, we are done
				expression = newExpr;
				expr = testExpr;
				exprVars = testVars;
				break;
			}
			catch(Exception err)
			{	JNUtilities.showMessage(docCtrl,"Invalid expression: "+err.getLocalizedMessage());
			}
		}
		
		// Get list of keys
		exprVars.put("#PI","3.141592653589793");
		
		// FEA or MPM variables
		keys.clear();
		if(docCtrl.resDoc.isMPMAnalysis())
		{
		    setKeyCode("#wpx",MPMSPINVELOCITYX);
		    setKeyCode("#wpy",MPMSPINVELOCITYY);
		    setKeyCode("#wpz",MPMSPINVELOCITYZ);
		    
		    setKeyCode("#c",MPMCONCENTRATION);
		    setKeyCode("#dcdx",MPMDCDX);
		    setKeyCode("#dcdy",MPMDCDY);
		    setKeyCode("#dcdz",MPMDCDZ);
		    
		    setKeyCode("#rho",MPMDENSITY);
		    
		    setKeyCode("#dispx",MPMDISPX);
		    setKeyCode("#dispR",MPMDISPX);
		    setKeyCode("#dispy",MPMDISPY);
		    setKeyCode("#dispZ",MPMDISPY);
		    setKeyCode("#dispz",MPMDISPZ);

		    setKeyCode("#exxe",MPMEPSX);
		    setKeyCode("#eRRe",MPMEPSX);
		    setKeyCode("#eyye",MPMEPSY);
		    setKeyCode("#eZZe",MPMEPSY);
		    setKeyCode("#ezze",MPMEPSZ);
		    setKeyCode("#eTTe",MPMEPSZ);
		    setKeyCode("#exye",MPMEPSXY);
		    setKeyCode("#eRZe",MPMEPSXY);
		    setKeyCode("#exze",MPMEPSXZ);
		    setKeyCode("#eyze",MPMEPSYZ);
		    
		    setKeyCode("#xing",MPMELEMENTCROSSINGS);

		    setKeyCode("#ener",MPMENERGY);
		    
		    setKeyCode("#h1",MPMHISTORY1);
		    setKeyCode("#h2",MPMHISTORY2);
		    setKeyCode("#h3",MPMHISTORY3);
		    setKeyCode("#h4",MPMHISTORY4);
		    setKeyCode("#h5",MPMHISTORY5);
		    setKeyCode("#h6",MPMHISTORY6);
		    setKeyCode("#h7",MPMHISTORY7);
		    setKeyCode("#h8",MPMHISTORY8);
		    setKeyCode("#h9",MPMHISTORY9);
		    setKeyCode("#h10",MPMHISTORY10);
		    setKeyCode("#h11",MPMHISTORY11);
		    setKeyCode("#h12",MPMHISTORY12);
		    setKeyCode("#h13",MPMHISTORY13);
		    setKeyCode("#h14",MPMHISTORY14);
		    setKeyCode("#h15",MPMHISTORY15);
		    setKeyCode("#h16",MPMHISTORY16);
		    setKeyCode("#h17",MPMHISTORY17);
		    setKeyCode("#h18",MPMHISTORY18);
		    setKeyCode("#h19",MPMHISTORY19);
		    
		    setKeyCode("#kine",MPMKINENERGY);

		    setKeyCode("#m",MPMMASS);

		    setKeyCode("#mat",MPMPOS);

		    setKeyCode("#x0",MPMORIGPOSX);
		    setKeyCode("#y0",MPMORIGPOSY);
		    setKeyCode("#z0",MPMORIGPOSZ);
		    
		    setKeyCode("#plaste",MPMPLASTICENERGY);
		    
		    setKeyCode("#exxp",MPMPLEPSX);
		    setKeyCode("#eRRp",MPMPLEPSX);
		    setKeyCode("#eyyp",MPMPLEPSY);
		    setKeyCode("#eZZp",MPMPLEPSY);
		    setKeyCode("#ezzp",MPMPLEPSZ);
		    setKeyCode("#eTTp",MPMPLEPSZ);
		    setKeyCode("#exyp",MPMPLEPSXY);
		    setKeyCode("#eRZp",MPMPLEPSXY);
		    setKeyCode("#exzp",MPMPLEPSXZ);
		    setKeyCode("#eyzp",MPMPLEPSYZ);

		    setKeyCode("#pp",MPMCONCENTRATION);
		    setKeyCode("#dpdx",MPMDCDX);
		    setKeyCode("#dpdy",MPMDCDY);
		    setKeyCode("#dpdz",MPMDCDZ);

		    setKeyCode("#x",MPMPOSX);
		    setKeyCode("#R",MPMPOSX);
		    setKeyCode("#y",MPMPOSY);
		    setKeyCode("#Z",MPMPOSY);
		    setKeyCode("#z",MPMPOSZ);

		    setKeyCode("#D",MPMPOSDISTANCE);

		    setKeyCode("#wxy",MPMWXY);
		    setKeyCode("#wxz",MPMWXZ);
		    setKeyCode("#wyz",MPMWYZ);

		    setKeyCode("#exx",MPMEPSTOTX);
		    setKeyCode("#eRR",MPMEPSTOTX);
		    setKeyCode("#eyy",MPMEPSTOTY);
		    setKeyCode("#eZZ",MPMEPSTOTY);
		    setKeyCode("#ezz",MPMEPSTOTZ);
		    setKeyCode("#eTT",MPMEPSTOTZ);
		    setKeyCode("#exy",MPMEPSTOTXY);
		    setKeyCode("#eRZ",MPMEPSTOTXY);
		    setKeyCode("#exz",MPMEPSTOTXZ);
		    setKeyCode("#eyz",MPMEPSTOTYZ);

		    setKeyCode("#stre",MPMSTRENERGY);

		    setKeyCode("#sxx",MPMSIGMAX);
			setKeyCode("#sRR",MPMSIGMAX);
		    setKeyCode("#syy",MPMSIGMAY);
			setKeyCode("#sZZ",MPMSIGMAY);
		    setKeyCode("#szz",MPMSIGMAZ);
		    setKeyCode("#sTT",MPMSIGMAZ);
		    setKeyCode("#sxy",MPMSIGMAXY);
		    setKeyCode("#sRZ",MPMSIGMAXY);
		    setKeyCode("#sxz",MPMSIGMAXZ);
		    setKeyCode("#syz",MPMSIGMAYZ);
		    
		    setKeyCode("#T",MPMPOSANGLE);

		    setKeyCode("#temp",MPMTEMPERATURE);

		    setKeyCode("#heate",MPMHEATENERGY);

		    setKeyCode("#t",MPMARCHIVETIME);

		    setKeyCode("#velx",MPMVELX);
		    setKeyCode("#velR",MPMVELX);
		    setKeyCode("#vely",MPMVELY);
		    setKeyCode("#velZ",MPMVELY);
		    setKeyCode("#velz",MPMVELZ);

		    setKeyCode("#wrke",MPMWORKENERGY);
			
		}
		else
		{
		    setKeyCode("#dispx",MESHDISPX);
		    setKeyCode("#dispR",MESHDISPX);
		    setKeyCode("#dispy",MESHDISPY);
		    setKeyCode("#dispZ",MESHDISPY);

		    setKeyCode("#fx",MESHFORCEX);
		    setKeyCode("#fy",MESHFORCEY);

		    setKeyCode("#esxx",MESHELEMSIGMAX);
		    setKeyCode("#esRR",MESHELEMSIGMAX);
		    setKeyCode("#esyy",MESHELEMSIGMAY);
		    setKeyCode("#esZZ",MESHELEMSIGMAY);
		    setKeyCode("#eszz",MESHELEMSIGMAZ);
		    setKeyCode("#esTT",MESHELEMSIGMAZ);
		    setKeyCode("#esxy",MESHELEMSIGMAXY);
		    setKeyCode("#esRZ",MESHELEMSIGMAXY);

		    setKeyCode("#ener",MESHSTRAINENERGY);
		    
		    setKeyCode("#x",MESHNODEX);
		    setKeyCode("#R",MESHNODEX);
		    setKeyCode("#y",MESHNODEY);
		    setKeyCode("#Z",MESHNODEY);
		    
		    setKeyCode("#D",MESHNODEDISTANCE);

		    setKeyCode("#exx",MESHSTRAINX);
		    setKeyCode("#eRR",MESHSTRAINX);
		    setKeyCode("#eyy",MESHSTRAINY);
		    setKeyCode("#eZZ",MESHSTRAINY);
		    setKeyCode("#ezz",MESHSTRAINZ);
		    setKeyCode("#eTT",MESHSTRAINZ);
		    setKeyCode("#exy",MESHSTRAINXY);
		    setKeyCode("#eRZ",MESHSTRAINXY);
		    
		    // shear components not mentioned in help files
		    setKeyCode("#dudy",MESHDUDY);
		    setKeyCode("#dvdx",MESHDVDX);

		    setKeyCode("#sxx",MESHSIGMAX);
		    setKeyCode("#sRR",MESHSIGMAX);
		    setKeyCode("#syy",MESHSIGMAY);
		    setKeyCode("#sZZ",MESHSIGMAY);
		    setKeyCode("#szz",MESHSIGMAZ);
		    setKeyCode("#sTT",MESHSIGMAZ);
		    setKeyCode("#sxy",MESHSIGMAXY);
		    setKeyCode("#sRZ",MESHSIGMAXY);

		    setKeyCode("#T",MESHNODEANGLE);
			
		}
	}
	
	// if variable is used, add to list of used variables
	public void setKeyCode(String key,int code)
	{	if(exprVars.get(key)!=null)
	    {	keys.add(key);
	        keyCode[keys.size()-1]=code;
	    }
	}
	
	// evaluate expression for a material point
	public double evaluateExpressionMP(MaterialPoint mptr,double angle,ResultsDocument doc)
	{
		for(int i=0;i<keys.size();i++)
		{	double keyValue = mptr.getForPlot(keyCode[i],angle,doc);
			exprVars.put(keys.get(i),JNUtilities.formatDouble(keyValue));
		}
		String errMsg = expr.evaluateWith(exprVars);
		double theValue = 1.e-20;
		if(errMsg==null)
		{	try
			{	theValue = expr.getNumericValue();
			}
			catch(Exception e)
			{	
			}
		}
		return theValue;
	}
	
	// evaluate expression for a material point
	public double evaluateExpressionNode(ElementBase eptr,int ndi,double angle,ResultsDocument doc)
	{
		for(int i=0;i<keys.size();i++)
		{	double keyValue;
			try
			{	keyValue = eptr.getNodeValue(ndi,keyCode[i],angle,doc);
			}
			catch(Exception e)
			{	keyValue = 0.;
			}
			exprVars.put(keys.get(i),JNUtilities.formatDouble(keyValue));
		}
		String errMsg = expr.evaluateWith(exprVars);
		double theValue = 1.e-20;
		if(errMsg==null)
		{	try
			{	theValue = expr.getNumericValue();
			}
			catch(Exception e)
			{	
			}
		}
		return theValue;
	}
	
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
				
			case MPMMAXSTRESS:
				return "Max Principal Stress";
				
			case MPMMINSTRESS:
				return "Min Principal Stress";
				
			case MPMSTRESSDIR:
				return "Stress Direction";
						
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
			
			case MPMVELS:
				return "||Velocity|| ";
			
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
				
			case MPMDISPS:
				return "||Displacment||";
			
			case MPMPOS:
				return "Position";
			
			case MPMPOSX:
			case MPMORIGPOSX:
				return "Position "+xc;
			case MPMPOSY:
			case MPMORIGPOSY:
				return "Position "+yc;
			case MPMPOSZ:
			case MPMORIGPOSZ:
				return "Position "+zc;
			
			case MESHDUDY:
				return "Strain du/d"+yc;
			
			case MESHDVDX:
				return "Strain dv/d"+xc;
			
			case MPMHISTORY1:
			case MPMHISTORY2:
			case MPMHISTORY3:
			case MPMHISTORY4:
				return "Material History "+(component-MPMHISTORY1+1);
			
			case MPMHISTORY5:
			case MPMHISTORY6:
			case MPMHISTORY7:
			case MPMHISTORY8:
			case MPMHISTORY9:
			case MPMHISTORY10:
			case MPMHISTORY11:
			case MPMHISTORY12:
			case MPMHISTORY13:
			case MPMHISTORY14:
			case MPMHISTORY15:
			case MPMHISTORY16:
			case MPMHISTORY17:
			case MPMHISTORY18:
			case MPMHISTORY19:
				return "Material History "+(component-MPMHISTORY5+5);
			
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
			
			case PLOTEXPRESSION:
				if(expressionTitle!=null) return expressionTitle;
				return resDoc.docCtrl.controls.getExpression();
			
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
			
			case MPMMODEIFB:
				return "Mode I Force";
			
			case MPMMODEIIFB:
				return "Mode II Force";
			
			case MPMCZLENGTH:
				return "Cohesive Damage Length";
			
	        case MPMCZMGI:
	            return "Mode I Energy Released";
	        
	        case MPMCZMGII:
	            return "Mode II Energy Released";
	            
			case MESHMATERIAL:
				return "Material";
			
			case MESHANGLE:
				if(resDoc.isFEAAnalysis()) return "Material Angle";
			case MPMANGLEZ:
			case MPMWXY:
				return "Rotational Strain "+xc+yc;
			case MPMANGLEY:
			case MPMWXZ:
				return "Rotational Strain "+xc+zc;
			case MPMANGLEX:
			case MPMWYZ:
				return "Rotational Strain "+yc+zc;
			
			
			case MESHFORCEX:
				return "Force "+xc;
				
			case MESHFORCEY:
				return "Force "+yc;
			
			case INTERFACETRACTION_N:
				return "Interface Normal Traction";
				
			case INTERFACETRACTION_T:
				return "Interface Tangential Traction";
				
			case MPMDENSITY:
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
				
			case MPMTRACTION1:
			case MPMTRACTION2:
			case MPMTRACTION3:
			case MPMTRACTION4:
			case MPMTRACTION5:
			case MPMTRACTION6:
			case MPMTRACTION7:
			case MPMTRACTION8:
			case MPMTRACTION9:
			case MPMTRACTION10:
				return "Traction history "+(component-MPMTRACTION1+1);
			
			default:
				break;
		}
		return "Unknown";
	}
}
