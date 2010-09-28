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
	static final int MPMEPSX=5;
	static final int MPMEPSY=6;
	static final int MPMEPSXY=7;
	static final int MPMEPSZ=8;
	static final int MPMPLEPSX=9;
	static final int MPMPLEPSY=10;
	static final int MPMPLEPSXY=11;
	static final int MPMPLEPSZ=12;
	static final int MPMVELX=13;
	static final int MPMVELY=14;
	static final int MPMVELVEC=15;
	static final int MPMSTRENERGY=16;
	static final int MPMKINENERGY=17;
	static final int MPMENERGY=18;
	static final int MPMTOTSTRENERGY=19;
	static final int MPMTOTKINENERGY=20;
	static final int MPMTOTENERGY=21;
	static final int MPMPOS=22;
	static final int MPMDISPX=23;
	static final int MPMDISPY=24;
	static final int MPMPOSX=25;
	static final int MPMPOSY=26;
	static final int MPMEXPRESSION=27;
	static final int MPMEXTWORK=28;
	static final int MPMPLASTICENERGY=29;
	static final int MPMTEMPERATURE=30;
	static final int MPMDVDX=31;
	static final int MPMDUDY=32;
	static final int MPMHISTORYOLD=33;
	static final int MPMTOTEXTWORK=34;
	static final int MPMTOTPOTENERGY=35;
	static final int MPMTOTPLASTICENERGY=36;
    static final int MPMJ1=37;
	static final int MPMJ2=38;
	static final int MPMKI=39;
	static final int MPMKII=40;
	static final int MPMLENGTH=41;
	static final int MPMMASS=42;
	static final int MPMARCHIVETIME=43;
	static final int MPMGLOBALRESULTS=44;
	static final int MPMCONCENTRATION=45;
	static final int MPMDCDX=46;
	static final int MPMDCDY=47;
	static final int MPMCRACKRELEASE=48;
	static final int MPMCRACKABSORB=49;
	static final int MPMNORMALCTOD=50;
	static final int MPMSHEARCTOD=51;
	static final int MPMTHERMALENERGY=52;
	static final int MPMTOTTHERMALENERGY=53;
	static final int MPMCRACKPROFILE=54;
	static final int MPMOPENINGFRACTION=55;
	static final int MPMSHEARFRACTION=56;
	static final int MPMANGLEZ=57;
	static final int MPMELEMENTCROSSINGS=58;
	static final int MPMTOTELEMENTCROSSINGS=59;
	static final int MPMDEBONDLENGTH=60;
	static final int MPMDEBONDNCTOD=61;
	static final int MPMDEBONDSCTOD=62;
	static final int MPMHISTORY1=63;
	static final int MPMHISTORY2=64;
	static final int MPMHISTORY3=65;
	static final int MPMHISTORY4=66;
	
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
	
	// pop-up menus
	JComboBox quant=new JComboBox();
	JComboBox cmpnt=new JComboBox();
	private int checkMeshItem;
	
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
			
			if(docCtrl.resDoc.isMPMAnalysis())
			{	byte [] arch=docCtrl.resDoc.archFormat.getBytes();
				
				if(arch[ReadArchive.ARCH_Stress]=='Y')
					quant.addItem(new PlotMenuItem("Stress",MPMSIGMAX));
				if(arch[ReadArchive.ARCH_Strain]=='Y')
					quant.addItem(new PlotMenuItem("Strain",MPMEPSX));
				if(arch[ReadArchive.ARCH_PlasticStrain]=='Y')
					quant.addItem(new PlotMenuItem("Plastic Strain",MPMPLEPSX));
					
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
				if(arch[ReadArchive.ARCH_ExtWork]=='Y')
					quant.addItem(new PlotMenuItem("External Work",MPMEXTWORK));
				if(arch[ReadArchive.ARCH_ThermalEnergy]=='Y')
					quant.addItem(new PlotMenuItem("Thermal Energy",MPMTHERMALENERGY));
					
				if(arch[ReadArchive.ARCH_Velocity]=='Y')
					quant.addItem(new PlotMenuItem("Velocity",MPMVELX));
				quant.addItem(new PlotMenuItem("Displacement",MPMDISPX));
				checkMeshItem=quant.getItemCount();
				quant.addItem(new PlotMenuItem("Material",MPMPOS));
				quant.addItem(new PlotMenuItem("Material Angle",MPMANGLEZ));
				quant.addItem(new PlotMenuItem("Mass",MPMMASS));
				if(arch[ReadArchive.ARCH_DeltaTemp]=='Y')
					quant.addItem(new PlotMenuItem("Temperature",MPMTEMPERATURE));
					
				if(arch[ReadArchive.ARCH_Concentration]=='Y')
				{   quant.addItem(new PlotMenuItem("Concentration",MPMCONCENTRATION));
					quant.addItem(new PlotMenuItem("Conc Gradient",MPMDCDX));
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
				
				if(arch[ReadArchive.ARCH_ShearComponents]=='Y')
					quant.addItem(new PlotMenuItem("Shear Component",MPMDVDX));
					
				if(arch[ReadArchive.ARCH_ElementCrossings]=='Y')
					quant.addItem(new PlotMenuItem("Element Crossings",MPMELEMENTCROSSINGS));
				
				// additional time plot options
				if(selected==LoadArchive.TIME_PLOT)
				{	// global results
					if(docCtrl.resDoc.globalArchive!=null)
						quant.addItem(new PlotMenuItem("Global Results",MPMGLOBALRESULTS));
						
					byte [] carch=docCtrl.resDoc.crackFormat.getBytes();
					
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
				
				// additional time plot options
				else if(selected==LoadArchive.MESH2D_PLOT)
				{	// x-y crack results
					quant.addItem(new PlotMenuItem("Crack Profile",MPMCRACKPROFILE));
					quant.addItem(new PlotMenuItem("Crack Normal CTOD",MPMNORMALCTOD));
					quant.addItem(new PlotMenuItem("Crack Tangential CTOD",MPMSHEARCTOD));
					quant.addItem(new PlotMenuItem("Crack Opening Fraction",MPMOPENINGFRACTION));
					quant.addItem(new PlotMenuItem("Crack Sliding Fraction",MPMSHEARFRACTION));
				}
			}
			
			// FEA plot quantities
			else
			{	char [] arch=docCtrl.resDoc.feaArchFormat;
			
				checkMeshItem=quant.getItemCount();
				quant.addItem(new PlotMenuItem("Mesh Only",MESHONLY));
				if(arch[ReadArchive.ARCH_FEAAvgStress]=='Y')
					quant.addItem(new PlotMenuItem("Stress",MESHSIGMAX));
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
			}
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
		switch(pm.getTag())
		{   case MPMSIGMAX:
			case MPMEPSX:
			case MPMPLEPSX:
			case MESHSIGMAX:
			case MESHSTRAINX:
			case MESHELEMSIGMAX:
				if(!cmpnt.getItemAt(0).equals(xchar+xchar))
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
			case MESHDISPX:
			case MESHFORCEX:
				if(!cmpnt.getItemAt(0).equals(xchar))
				{	cmpnt.removeAllItems();
					cmpnt.addItem(xchar);
					cmpnt.addItem(ychar);
				}
				cmpnt.setEnabled(true);
				break;
			
				
			case MPMDCDX:
				if(!cmpnt.getItemAt(0).equals("dc/d"+xchar))
				{	cmpnt.removeAllItems();
					cmpnt.addItem("dc/d"+xchar);
					cmpnt.addItem("dc/d"+ychar);
				}
				cmpnt.setEnabled(true);
				break;
				
			case MPMDVDX:
			case MESHDVDX:
				if(!cmpnt.getItemAt(0).equals("dv/d"+xchar))
				{	cmpnt.removeAllItems();
					cmpnt.addItem("dv/d"+xchar);
					cmpnt.addItem("du/d"+ychar);
				}
				cmpnt.setEnabled(true);
				break;
			
			case INTERFACETRACTION_N:
				if(!cmpnt.getItemAt(0).equals("normal"))
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
	public int getPlotComponent()
	{	PlotMenuItem pm=(PlotMenuItem)quant.getSelectedItem();
		if(pm==null) return -1;
		int plotComponent=pm.getTag();
		
		// adjust component menus - add component selected in component menu
		int extra;
		switch(plotComponent)
		{   case MPMSIGMAX:
			case MPMEPSX:
			case MPMPLEPSX:
			case MPMVELX:
			case MPMDCDX:
			case MPMDISPX:
			case MPMDVDX:
			case MESHSIGMAX:
			case MESHDISPX:
			case MESHSTRAINX:
			case MESHDVDX:
			case MESHELEMSIGMAX:
			case MESHFORCEX:
			case INTERFACETRACTION_N:
				extra=cmpnt.getSelectedIndex();
				if(extra>=0) plotComponent+=extra;
				break;
			default:
				break;
		}
		return plotComponent;
	}
	
	// Label for plot axes (i.e., just generic name and plot units)
	public static String plotLabel(int component,String distU,String timeU)
	{
		switch(component)
		{	// Stresses
			case MPMSIGMAX:
			case MPMSIGMAY:
			case MPMSIGMAXY:
			case MPMSIGMAZ:
			case MESHSIGMAX:
			case MESHSIGMAY:
			case MESHSIGMAXY:
			case MESHSIGMAZ:
			case MESHELEMSIGMAX:
			case MESHELEMSIGMAY:
			case MESHELEMSIGMAXY:
			case MESHELEMSIGMAZ:
				return "Stress (MPa)";
		
			// Strains
			case MPMEPSX:
			case MPMEPSY:
			case MPMEPSXY:
			case MPMEPSZ:
			case MPMPLEPSX:
			case MPMPLEPSY:
			case MPMPLEPSXY:
			case MPMPLEPSZ:
			case MESHSTRAINX:
			case MESHSTRAINY:
			case MESHSTRAINXY:
			case MESHSTRAINZ:
				return "Strain (%)";
			
			// Energy
			case MPMENERGY:
			case MPMSTRENERGY:
			case MPMKINENERGY:
			case MPMEXTWORK:
			case MPMPLASTICENERGY:
			case MPMTHERMALENERGY:
			case MPMTOTSTRENERGY:
			case MPMTOTKINENERGY:
			case MPMTOTENERGY:
			case MPMTOTEXTWORK:
			case MPMTOTPOTENERGY:
			case MPMTOTPLASTICENERGY:
			case MPMTOTTHERMALENERGY:
			case MESHSTRAINENERGY:
				return "Energy (J)";
			
			case MPMTEMPERATURE:
				return "Temperature (K)";
			
			// Velocity
			case MPMVELVEC:
			case MPMVELX:
			case MPMVELY:
				return "Velocity ("+distU+"/"+timeU+")";
			
			// Displacements
			case MPMDISPX:
			case MPMDISPY:
			case MPMNORMALCTOD:
			case MPMSHEARCTOD:
			case MPMDEBONDNCTOD:
			case MPMDEBONDSCTOD:
			case MESHDISPX:
			case MESHDISPY:
				return "Displacement ("+distU+")";
			
			// Position or material
			case MESHANGLE:
			case MPMANGLEZ:
				return "Material Angle";
				
			// Position or material
			//case MESHMATERIAL:
			case MPMPOS:
				return "Material";
					 
			// shear components
			case MPMDUDY:
			case MPMDVDX:
			case MESHDVDX:
			case MESHDUDY:
				return "Shear Component (%)";
			
			// histrory variables
			case MPMHISTORY1:
			case MPMHISTORY2:
			case MPMHISTORY3:
			case MPMHISTORY4:
				return "Material History";
			
			// concentration
			case MPMCONCENTRATION:
				return "Concentration";
			
			// concentration gradient
			case MPMDCDX:
			case MPMDCDY:
				return "Displacement (1/"+distU+")";
			
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
				return "Length ("+distU+")";
			
			case MPMCRACKRELEASE:
			case MPMCRACKABSORB:
				return "Energy Rate (J/m^2)";
			
			case MESHFORCEX:
			case MESHFORCEY:
				return "Force (N)";
			
			case INTERFACETRACTION_N:
			case INTERFACETRACTION_T:
				return "Interface Traction (MPa)";
				
			case MPMPOSY:
			case MPMPOSX:
			case MPMCRACKPROFILE:
				return "Position ("+distU+")";
			
			case MPMOPENINGFRACTION:
			case MPMSHEARFRACTION:
				return "Fraction";
			
			case MPMMASS:
				return "Mass (g)";
			
			case MPMARCHIVETIME:
				return "Time ("+timeU+")";
			
			case MPMGLOBALRESULTS:
				return "Global Results";
			
			case MPMELEMENTCROSSINGS:
			case MPMTOTELEMENTCROSSINGS:
				return "Element Crossings";
			
			// Unknown
			default:
				break;
		}
		
		return "Y Axis";
	}
	
	// select item for checking the mesh
	public void setCheckMeshItem() { quant.setSelectedIndex(checkMeshItem); }
	
	// name for plot component
	public static String plotName(int component)
	{
		switch(component)
		{	case MPMSIGMAX:
			case MESHSIGMAX:
			case MESHELEMSIGMAX:
				return "Stress xx";
				
			case MPMSIGMAY:
			case MESHSIGMAY:
			case MESHELEMSIGMAY:
				return "Stress yy";
				
			case MPMSIGMAXY:
			case MESHSIGMAXY:
			case MESHELEMSIGMAXY:
				return "Stress xy";
				
			case MPMSIGMAZ:
			case MESHSIGMAZ:
			case MESHELEMSIGMAZ:
				return "Stress zz";
				
			case MPMEPSX:
			case MESHSTRAINX:
				return "Strain xx";
				
			case MPMEPSY:
			case MESHSTRAINY:
				return "Strain yy";
				
			case MPMEPSXY:
			case MESHSTRAINXY:
				return "Strain xy";
				
			case MPMEPSZ:
			case MESHSTRAINZ:
				return "Strain zz";
				
			case MPMPLEPSX:
				return "Plastic Strain xx";
				
			case MPMPLEPSY:
				return "Plastic Strain yy";
				
			case MPMPLEPSXY:
				return "Plastic Strain xy";
				
			case MPMPLEPSZ:
				return "Plastic Strain zz";

			case MPMENERGY:
				return "Energy";

			case MPMSTRENERGY:
			case MESHSTRAINENERGY:
				return "Strain Energy";

			case MPMKINENERGY:
				return "KineticEnergy";

			case MPMEXTWORK:
				return "External Work";

			case MPMPLASTICENERGY:
				return "Plastic Energy";

			case MPMTHERMALENERGY:
				return "Thermal Energy";

			case MPMTOTSTRENERGY:
				return "Total Strain Energy";
				
			case MPMTOTKINENERGY:
				return "Total Kinetic Energy";
				
			case MPMTOTENERGY:
				return "Total Energy";
				
			case MPMTOTEXTWORK:
				return "Total External Work";
				
			case MPMTOTPOTENERGY:
				return "Total Potential Energy";
				
			case MPMTOTPLASTICENERGY:
				return "Total Plastic Energy";
						
			case MPMTOTTHERMALENERGY:
				return "Total Thermal Energy";
						
			case MPMTEMPERATURE:
				return "Temperature";
			
			case MPMVELX:
				return "Velocity x";
			
			case MPMVELY:
				return "Velocity y";
			
			case MPMVELVEC:
				return "Velocity";
			
			case MPMDISPX:
			case MESHDISPX:
				return "Displacment x";
			
			case MPMDISPY:
			case MESHDISPY:
				return "Displacment y";
			
			case MPMPOS:
				return "Position";
			
			case MPMPOSX:
				return "Position x";
			
			case MPMPOSY:
				return "Position y";
			
			case MPMDUDY:
			case MESHDUDY:
				return "Strain du/dy";
			
			case MPMDVDX:
			case MESHDVDX:
				return "Strain dv/dx";
			
			case MPMHISTORY1:
			case MPMHISTORY2:
			case MPMHISTORY3:
			case MPMHISTORY4:
				return "Material History";
			
			case MPMDCDY:
				return "Conc Gradient dc/dy";
			
			case MPMDCDX:
				return "Conc Gradient dc/dx";
			
			case MPMCONCENTRATION:
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
				return "Material Angle";
			
			case MESHFORCEX:
				return "Force x";
				
			case MESHFORCEY:
				return "Force y";
			
			case INTERFACETRACTION_N:
				return "Interface Normal Traction";
				
			case INTERFACETRACTION_T:
				return "Interface Tangential Traction";
				
			case MPMMASS:
				return "Mass";
			
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
