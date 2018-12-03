/*******************************************************************
	ResultsDocument.java
	NairnFEAMPMViz

	Created by John Nairn on Sat Mar 06 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import geditcom.JNFramework.JNUtilities;

import java.awt.geom.*;
import java.io.*;
import java.nio.*;

import javax.swing.table.*;

import java.util.*;

public class ResultsDocument extends AbstractTableModel
{
	static final long serialVersionUID=19L;
	
	//---------------------------------------------------------
	// variables and constants
	//---------------------------------------------------------
	public static final int PLANE_STRAIN=0;
	public static final int PLANE_STRESS=1;
	public static final int AXI_SYM=2;
	public static final int THREE_DIM=3;
	public static final int BEGIN_MPM_TYPES=9;
	public static final int PLANE_STRAIN_MPM=10;
	public static final int PLANE_STRESS_MPM=11;
	public static final int THREED_MPM=12;
	public static final int AXI_SYM_MPM=13;
	
	public String name;
	public String fullPath;
	public File path;
	public ArrayList<String> sectionText;
	public ArrayList<String> sectionTitle;
	public ArrayList<NodalPoint> nodes;
	public ArrayList<GridDispBC> gridBCs;
	public ArrayList<ParticleBC> particleBCs;
	public ArrayList<ElementBase> elements;
	public ArrayList<File> archives;
	public File globalArchive;
	public ArrayList<Double> archiveTimes;
	public ArrayList<MaterialPoint> mpmPoints;
	public ArrayList<CrackHeader> mpmCracks;
	public ArrayList<MaterialBase> materials;
	public ArrayList<Vector3> pointDims;
	public String archFormat,crackFormat;
	public char[] feaArchFormat = {'N','N','N','N','N','N' };
	public double xmin,xmax,ymin,ymax,zmin,zmax;			// mesh bounds
	public double dxmin,dxmax,dymin,dymax;					// displaced mesh bounds
	// cell dimension in any direction will be cellMinSide*iscale
	public double cellMinSide,xscale=-1.,yscale=-1.,zscale=-1;
	public int np;
	public DocViewer docCtrl;
	public int recSize;
	public double totalEnergy;
	public boolean hasPorePressure;
	
	// units
	public JNUnits units = new JNUnits();
	
	private int currentArchive;
	
	//---------------------------------------------------------
	// Initialize
	//---------------------------------------------------------
	
	public ResultsDocument()
	{
		sectionTitle=new ArrayList<String>(20);
		sectionText=new ArrayList<String>(20);
		nodes=new ArrayList<NodalPoint>(20000);
		elements=new ArrayList<ElementBase>(20000);
		gridBCs=new ArrayList<GridDispBC>(200);
		particleBCs=new ArrayList<ParticleBC>(200);
		mpmPoints=new ArrayList<MaterialPoint>(20000);
		mpmCracks=new ArrayList<CrackHeader>(10);
		archives=new ArrayList<File>(100);
		archiveTimes=new ArrayList<Double>(100);
		materials=new ArrayList<MaterialBase>(10);
		pointDims=null;
		currentArchive=-1;
	}
	
	//-----------------------------------------------------------------
	// After file is read into section (in TextDisplay), this method
	// is called to decode all sections and extract data needed
	// for plotting. Throw exception if an errors in the data
	//-----------------------------------------------------------------
	
	public void DecodeFileSections(File file) throws Exception
	{   int endIndex,lineStart;
		String line,archDir;
		char strChar;
		String [] words;
		Scanner s=null,sline=null;
		int numMps=1;
		
		//----------------------------------------------------------
		// Header
		String header=section("TITLE");
		String fileUnits = "_none_";
		try
		{	s=new Scanner(header);
			s.useDelimiter("\\r\\n|\\n|\\r");
			while(s.hasNext())
			{	line = s.next();
				if(line.startsWith("Units:"))
				{	fileUnits = line.substring(7);
					break;
				}
			}
			s.close();
			s=null;
		}
		catch(NoSuchElementException e)
		{	s.close();
			throw new Exception("Could not read the header for this file");
		}
		units.setOutputUnits(fileUnits);
		
		//----------------------------------------------------------
		// Expected number of nodes and elements
		String summary=section("NODES AND ELEMENTS");
		int nnodes,nelems;
		Scanner nline=null;
		try
		{	// get nodes and elements
			s=new Scanner(summary);
			s.useDelimiter("\\r\\n|\\n|\\r");
			s.useLocale(Locale.US);
			s.next();
			s.next();
			
			nline=new Scanner(s.next());
			nline.useLocale(Locale.US);
			nline.next();
			nnodes=nline.nextInt();
			nline.next();
			nelems=nline.nextInt();
			nline.close();
		
			// Options are 2D Plane Strain Analysis, 2D Plane Stress Analysis, Axisymmetric Analysis,
			//	2D Plane Strain MPM Analysis, 2D Plane Stress MPM Analysis, 3D MPM Analysis
			nline=new Scanner(s.next());
			nline.useLocale(Locale.US);
			nline.next();
			nline.next();
			nline.next();
			nline.next();
			String word=nline.next();
			if(word.equals("3D"))
				np=THREED_MPM;
			else if(word.equals("Axisymmetric"))
			{	if(nline.next().equals("MPM"))
					np=AXI_SYM_MPM;
				else
					np=AXI_SYM;
			}
			else
			{	nline.next();
				if(nline.next().equals("Strain"))
					np=PLANE_STRAIN;
				else
					np=PLANE_STRESS;
				if(nline.next().equals("MPM"))
					np = (np==PLANE_STRAIN) ? PLANE_STRAIN_MPM : PLANE_STRESS_MPM;
			}
			nline.close();
			s.close();
			s=null;
		}
		catch(NoSuchElementException e)
		{	if(s!=null) s.close();
			if(nline!=null) nline.close();
			throw new Exception("Could not decode analysis type for this file");
		}
		//if(is3D())
		//	throw new Exception("This tool cannot visualize 3D results; see help information for other options.");
		
		//----------------------------------------------------------
		// Nodal Point Coordinates
		String ndst=section("NODAL POINT COORDINATES");
		lineStart=findNextLine(ndst,"-----");
		if(lineStart<0)
			throw new Exception("Error: missing nodal point coordinates section");
		s=new Scanner(ndst.substring(lineStart,ndst.length()-1));
		s.useLocale(Locale.US);
		int prevNodeNum=0,nodeNum;
		double xpt,ypt,zpt;
		if(is3D())
		{	while(s.hasNextInt())
			{	nodeNum=s.nextInt();
				xpt=s.nextDouble()*units.lengthScale();
				ypt=s.nextDouble()*units.lengthScale();
				zpt=s.nextDouble()*units.lengthScale();
				if(nodeNum!=prevNodeNum+1)
				{	s.close();
					throw new Exception("Some node numbers are missing");
				}
				addNode(nodeNum,xpt,ypt,zpt);
				prevNodeNum=nodeNum;
			}
		}
		else
		{	while(s.hasNextInt())
			{	nodeNum=s.nextInt();
				xpt=s.nextDouble()*units.lengthScale();
				ypt=s.nextDouble()*units.lengthScale();
				if(nodeNum!=prevNodeNum+1)
				{	s.close();
					throw new Exception("Some node numbers are missing");
				}
				addNode(nodeNum,xpt,ypt);
				prevNodeNum=nodeNum;
			}
		}
		s.close();
		s=null;
		if(prevNodeNum!=nnodes)
			throw new Exception("Number of nodes does not match expected number of nodes.");
		
		// bounds now in xmin, xmax, etc., transferred to displaced bounds
		dxmin=xmin;
		dxmax=xmax;
		dymin=ymin;
		dymax=ymax;
		
		//----------------------------------------------------------
		// Element Definitions
		String elems=section("ELEMENT DEFINITIONS");
		lineStart=findNextLine(elems,"-----");
		if(lineStart<0)
			throw new Exception("Error decoding element definitions");
		s=new Scanner(elems.substring(lineStart,elems.length()-1));
		s.useLocale(Locale.US);
		int prevElemNum=0,elemNum,elemID,i,matID=0;
		double elemThickness=1.,elemAngle=0.;
		int[] nds;
		nds=new int[9];
		while(s.hasNextInt())
		{	elemNum=s.nextInt();				// element number
			if(elemNum!=prevElemNum+1)
			{	s.close();
				throw new Exception("Some element numbers are missing");
			}
			elemID=s.nextInt();					// element ID
			
			// material and thickness
			if(isFEAAnalysis())
			{	matID=s.nextInt();
				elemAngle=s.nextDouble();
				if(!isAxisymmetric())
					elemThickness=s.nextDouble();
			}
			
			// nodes
			switch(elemID)
			{	case ElementBase.FOUR_NODE_ISO:
				case ElementBase.EIGHT_NODE_ISO:
				case ElementBase.CST:
				case ElementBase.ISO_TRIANGLE:
				case ElementBase.LINEAR_INTERFACE:
				case ElementBase.QUAD_INTERFACE:
				case ElementBase.EIGHT_NODE_ISO_BRICK:
				case ElementBase.LAGRANGE_2D:
					// read all nodes
					for(i=0;i<ElementBase.NodesFromType(elemID);i++)
					{	nds[i]=s.nextInt();
					}
					break;
				
				default:
					throw new Exception("Element type found ("+elemID+") that is not yet supported in this tool.");
			}
			
			addElement(elemNum,elemID,nds,matID,elemAngle,elemThickness*units.lengthScale());
			prevElemNum=elemNum;
		}
		s.close();
		s=null;
		if(prevElemNum!=nelems)
			throw new Exception("Number of elements does not match expected number of elements.");
		
		//----------------------------------------------------------
		// defined materials
		String matls=section("DEFINED MATERIALS");
		s=new Scanner(matls);
		s.useDelimiter("\\r\\n|\\n|\\r");
		s.useLocale(Locale.US);
		while(true)
		{	// scan to start of next material and get its name
			String matLine,matName=null;
			while(matName==null && s.hasNext())
			{	matLine=s.next();
				if(matLine.length()<9) continue;
				if(matLine.substring(0,9).equals("Material "))
				{	int beginIndex=matLine.indexOf(":");
					if(beginIndex<0)
						matName="Unknown Material Name";
					else
						matName=matLine.substring(beginIndex+2,matLine.length());
				}
			}
			if(matName==null) break;
			
			// get material type
			sline=new Scanner(s.next());
			sline.useLocale(Locale.US);
			String word1=sline.next();
			String word2=sline.next();
			MaterialBase matl;
			if(word1.equals("Isotropic"))
			{	if(word2.equals("Dugdale"))
					matl=new MaterialBase(matName,MaterialBase.DUGDALE);
        		else if(word2.equals("Elastic-Plastic"))
        			matl=new MaterialBase(matName,MaterialBase.ISOPLASTICITY);
            	else if(word2.equals("Softening"))
            		matl=new MaterialBase(matName,MaterialBase.ISOSOFTENING);
            	else if(word2.equals("Plastic"))
            		matl=new MaterialBase(matName,MaterialBase.ISOPLASTICSOFTENING);
				else
					matl=new IsotropicMat(matName);
			}
			else if(word1.equals("Tranversely"))
			{	sline.next();
				sline.next();
				word1 = sline.next();
				// old style
				if(word1.equals("normal"))
					matl=new MaterialBase(matName,MaterialBase.TRANSISO1);
				else if(word1.equals("in"))
					matl=new MaterialBase(matName,MaterialBase.TRANSISO2);
				else
				{	// new style
					sline.next();
					word1=sline.next();
					word2=sline.next();
					if(word1.equals("z"))
						matl=new MaterialBase(matName,MaterialBase.TRANSISO1);
					else if(word1.equals("y"))
						matl=new MaterialBase(matName,MaterialBase.TRANSISO2);
					else if(word2.equals("z"))
						matl=new MaterialBase(matName,MaterialBase.TRANSISOSOFTENING1);
					else
						matl=new MaterialBase(matName,MaterialBase.TRANSISOSOFTENING2);
				}
			}
			else if(word1.equals("Orthotropic"))
				matl=new MaterialBase(matName,MaterialBase.ORTHO);
			else if(word1.equals("Bistable"))
				matl=new MaterialBase(matName,MaterialBase.BISTABLEISO);
			else if(word1.equals("Von"))
				matl=new MaterialBase(matName,MaterialBase.VONMISES);
			else if(word1.equals("Mooney-Rivlin"))
				matl=new MaterialBase(matName,MaterialBase.MOONEYRIVLIN);
			else if(word1.equals("Viscoelastic"))
				matl=new MaterialBase(matName,MaterialBase.VISCOELASTIC);
			else if(word1.equals("Interface"))
				matl=new MaterialBase(matName,MaterialBase.INTERFACEPARAMS);
			else if(word1.equals("Rigid"))
			{	if(word2.equals("Block"))
					matl=new RigidMaterial(matName,MaterialBase.RIGIDBLOCKMATERIAL);
				else if(word2.equals("Contact"))
					matl=new RigidMaterial(matName,MaterialBase.RIGIDCONTACTMATERIAL);
				else
					matl=new RigidMaterial(matName,MaterialBase.RIGIDBCMATERIAL);
			}
			else if(word1.equals("Triangular"))
				matl=new MaterialBase(matName,MaterialBase.COHESIVEZONEMATERIAL);
			else if(word1.equals("Linear"))
			{	if(word2.equals("Imperfect"))
					matl=new MaterialBase(matName,MaterialBase.LINEARIMPERFECT);
				else
					matl=new MaterialBase(matName,MaterialBase.LINEARTRACTIONMATERIAL);
			}
			else if(word1.equals("Cubic"))
				matl=new MaterialBase(matName,MaterialBase.CUBICTRACTIONMATERIAL);
			else if(word1.equals("Elastic-Plastic"))
				matl=new MaterialBase(matName,MaterialBase.HILLPLASTIC);
	        else if(word1.equals("M-G,"))
				matl=new MaterialBase(matName,MaterialBase.MGSCGLMATERIAL);
			else if(word1.equals("Steinberg-Lund"))
				matl=new MaterialBase(matName,MaterialBase.SLMATERIAL);
			else if(word1.equals("Johnson-Cook"))
				matl=new MaterialBase(matName,MaterialBase.JOHNSONCOOK);
			else if(word1.equals("Wood"))
				matl=new MaterialBase(matName,MaterialBase.WOODMATERIAL);
			else if(word1.equals("Ideal"))
				matl=new MaterialBase(matName,MaterialBase.IDEALGAS);
			else if(word1.equals("Hyperelastic"))
	        {   if(word2.equals("Isotropic"))
					matl=new MaterialBase(matName,MaterialBase.HEISOTROPIC);
				else if(word2.equals("MGEOS"))
					matl=new MaterialBase(matName,MaterialBase.HEMGEOSMATERIAL);
				else
					matl=new MaterialBase(matName,MaterialBase.HEANISOTROPIC);
	        }
			else if(word1.equals("Phase"))
				matl=new MaterialBase(matName,MaterialBase.PHASETRANSITION);
			else if(word1.equals("Contact"))
				matl=new MaterialBase(matName,MaterialBase.IGNORECONTACT);
			else if(word1.equals("Coulomb"))
				matl=new MaterialBase(matName,MaterialBase.COULOMBFRICTION);
			else if(word1.equals("Adhesion"))
				matl=new MaterialBase(matName,MaterialBase.ADHESIVEFRICTION);
			else if(word1.equals("Liquid/Wall"))
				matl=new MaterialBase(matName,MaterialBase.LIQUIDCONTACT);
			else
			{	// try to continue with unknown material type
				matl=new MaterialBase(matName,MaterialBase.UNKNOWNMATERIAL);
			}
			sline.close();
			sline=null;
			
			// decode material properteis
			matl.decodeData(s);
			materials.add(matl);
		}
		s.close();
		s=null;

		//----------------------------------------------------------
		// mesh boundary conditions
		String bcs=section("NODAL POINTS WITH FIXED DISPLACEMENTS");
		lineStart=findNextLine(bcs,"-----");
		if(lineStart>0 && lineStart<bcs.length())
		{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
			s.useLocale(Locale.US);
			int dof,bcID;
			double bcVal,bcArg,bcAngle=0.,bcAngle2=0.;
			while(s.hasNextInt())
			{	nodeNum=s.nextInt();
				dof=s.nextInt();
				
				if(isMPMAnalysis())
				{	bcID=s.nextInt();
					bcVal=s.nextDouble();
					bcArg=s.nextDouble();
					bcAngle=0.;
					
					// keep reading as long as next item is not integer on the next line
					// could be a problem if function value looks like an integer
					if(!s.hasNextInt())
					{	if(s.hasNextDouble())
							bcAngle=s.nextDouble();
					}
					if(!s.hasNextInt())
					{	if(s.hasNextDouble())
							bcAngle2=s.nextDouble();
					}
					if(bcID==BoundaryCondition.FUNCTION_VALUE && !s.hasNextInt())
						s.next();
					
					addGridBC(nodeNum,dof,bcID,bcVal*units.lengthScale(),bcArg*units.altTimeScale(),bcAngle,bcAngle2);
				}
				else
				{	bcVal=s.nextDouble();
					bcArg=0.;
					bcAngle=0.;
					if(dof==3)
						bcID=BoundaryCondition.FEA_ROTATEBC;
					else
					{	bcID=BoundaryCondition.FEA_DISPBC;
						// find out if this condition is rotated
						for(i=0;i<gridBCs.size();i++)
						{	GridDispBC obj=gridBCs.get(i);
							if(obj.rotates(nodeNum))
								bcAngle+=obj.getValue();
						}
					}
					addGridBC(nodeNum,dof,bcID,bcVal*units.lengthScale(),bcArg*units.altTimeScale(),bcAngle,bcAngle2);
				}
			}
			s.close();
			s=null;
		}
		
		//----------------------------------------------------------
		// particle boundary conditions and grid info
		if(isMPMAnalysis())
		{	bcs=section("MATERIAL POINTS WITH EXTERNAL FORCES");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int dof,bcID;
				double bcLoad,bcArg;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();		// particle number here
					dof=s.nextInt();
					bcID=s.nextInt();
					bcLoad=s.nextDouble();
					bcArg=s.nextDouble();
					if(bcID==BoundaryCondition.FUNCTION_VALUE)
						s.next();
					
					// value not scaled, but current not visualizes
					addParticleBC(nodeNum,dof,0,bcID,bcLoad,bcArg*units.altTimeScale());
				}
				s.close();
				s=null;
			}
			
			bcs=section("MATERIAL POINTS WITH TRACTIONS");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int dof,bcID,bcFace;
				double bcLoad,bcArg;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();		// particle number here
					dof=s.nextInt();
					bcFace=s.nextInt();
					bcID=s.nextInt();
					bcLoad=s.nextDouble();
					bcArg=s.nextDouble();
					if(bcID==BoundaryCondition.FUNCTION_VALUE)
						s.next();
					addParticleBC(nodeNum,dof,bcFace,bcID,bcLoad,bcArg*units.altTimeScale());
				}
				s.close();
				s=null;
			}
			
			bcs=section("NODAL POINTS WITH FIXED CONCENTRATIONS");
			lineStart=findNextLine(bcs,"--------");
			if(lineStart<=0)
			{	bcs=section("NODAL POINTS WITH FIXED PORE PRESSURE");
				lineStart=findNextLine(bcs,"--------");
			}
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int bcID;
				double bcVal,bcTime;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();
					bcID=s.nextInt();
					bcVal=s.nextDouble();
					bcTime=s.nextDouble();
					if(bcID==BoundaryCondition.FUNCTION_VALUE)
						s.next();
					addGridBC(nodeNum,BoundaryCondition.CONCENTRATION_DIR,
								bcID,bcVal,bcTime*units.altTimeScale(),0.,0.);
				}
				s.close();
				s=null;
			}
			
			bcs=section("MATERIAL POINTS WITH CONCENTRATION FLUX");
			lineStart=findNextLine(bcs,"--------");
			if(lineStart<=0)
			{	bcs=section("NODAL POINTS WITH PORE PRESSURE FLUX");
				lineStart=findNextLine(bcs,"--------");
			}
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int bcID,bcDir,bcFace;
				double bcVal,bcTime;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();
					bcDir=s.nextInt();
					bcFace=s.nextInt();
					bcID=s.nextInt();
					bcVal=s.nextDouble();
					bcTime=s.nextDouble();
					if(bcID==BoundaryCondition.FUNCTION_VALUE)
						s.next();
					addParticleBC(nodeNum,bcDir,bcFace,bcID,bcVal,bcTime*units.altTimeScale());
					particleBCs.get(particleBCs.size()-1).setBCType(BoundaryCondition.CONCENTRATION_DIR);
				}
				s.close();
				s=null;
			}
			
			bcs=section("NODAL POINTS WITH FIXED TEMPERATURES");
			lineStart=findNextLine(bcs,"--------");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int bcID;
				double bcVal,bcTime;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();
					bcID=s.nextInt();
					bcVal=s.nextDouble();
					bcTime=s.nextDouble();
					if(bcID==BoundaryCondition.FUNCTION_VALUE)
						s.next();
					addGridBC(nodeNum,BoundaryCondition.TEMPERATURE_DIR,
								bcID,bcVal,bcTime*units.altTimeScale(),0.,0.);
				}
				s.close();
				s=null;
			}

			bcs=section("MATERIAL POINTS WITH HEAT FLUX");
			lineStart=findNextLine(bcs,"--------");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int bcID,bcDir,bcFace;
				double bcVal,bcTime;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();
					bcDir=s.nextInt();
					bcFace=s.nextInt();
					bcID=s.nextInt();
					bcVal=s.nextDouble();
					bcTime=s.nextDouble();
					if(bcID==BoundaryCondition.FUNCTION_VALUE)
						s.next();
					addParticleBC(nodeNum,bcDir,bcFace,bcID,bcVal,bcTime*units.altTimeScale());
					particleBCs.get(particleBCs.size()-1).setBCType(BoundaryCondition.TEMPERATURE_DIR);
				}
				s.close();
				s=null;
			}

			bcs=section("FULL MASS MATRIX");
			s=new Scanner(bcs);
			s.useDelimiter("\\r\\n|\\n|\\r");
			s.useLocale(Locale.US);
			// scan to grid info
			String gridInfo=null;
			while(s.hasNext() && gridInfo==null)
			{	String gridLine=s.next();
				if(gridLine.length()<23) continue;
				if(gridLine.substring(0,21).equals("Structured orthogonal") || 
						gridLine.substring(0,23).equals("Unstructured orthogonal") ||
						gridLine.substring(0,10).equals("Orthogonal"))
				{	int beginIndex=gridLine.indexOf(":");
					if(beginIndex<0) break;
					gridInfo=gridLine.substring(beginIndex+1,gridLine.length());
				}
				else if(gridLine.substring(0,10).equals("Number of "))
				{	int beginIndex=gridLine.indexOf(":");
					if(beginIndex<0) break;
					String numStr = gridLine.substring(beginIndex+1,gridLine.length());
					numMps = Integer.parseInt(numStr.trim());
				}
			}
			s.close();
			s=null;
			
			if(gridInfo!=null)
			{	// read grid info and find relative cell sides
				// minimum side scale is 1 and other are relative to that side
				sline=new Scanner(gridInfo);
				sline.useDelimiter("[ :]");
				sline.useLocale(Locale.US);
				if(sline.hasNextDouble())
					xscale=sline.nextDouble()*units.lengthScale();
				else
					xscale=1.;
				if(sline.hasNext()) sline.next();
				if(sline.hasNext()) sline.next();
				if(sline.hasNextDouble())
					yscale=sline.nextDouble()*units.lengthScale();
				else
					yscale=xscale;
				if(is3D())
				{	if(sline.hasNext()) sline.next();
					if(sline.hasNextDouble())
						zscale=sline.nextDouble()*units.lengthScale();
					else
						zscale=xscale;
					if(zscale<xscale && zscale<yscale)
					{	xscale/=zscale;
						yscale/=zscale;
						cellMinSide = zscale;
						zscale=1;
					}
					else if(yscale<xscale)
					{	xscale/=yscale;
						zscale/=yscale;
						cellMinSide = yscale;
						yscale=1.;
					}
					else
					{	yscale/=xscale;
						zscale/=xscale;
						cellMinSide = xscale;
						xscale=1.;
					}
				}
				else
				{	if(xscale>yscale)
					{	xscale/=yscale;
						cellMinSide = yscale;
						yscale=1.;
					}
					else
					{	yscale/=xscale;
						cellMinSide = xscale;
						xscale=1.;
					}
				}
				
				sline.close();
				sline=null;
			}
			
			bcs=section("SCHEDULED CUSTOM TASKS");
			s=new Scanner(bcs);
			s.useDelimiter("\\r\\n|\\n|\\r");
			s.useLocale(Locale.US);
			// scan THROUGH CUSTOM TASKS
			hasPorePressure = false;
			while(s.hasNext())
			{	String gridLine=s.next();
				if(gridLine.length()>12)
				{	if(gridLine.substring(0,12).equals("Coupled pore"))
					{	hasPorePressure = true;
						break;
					}
				}
			}
		}
		
		//----------------------------------------------------------
		// nodal load conditions
		else
		{	bcs=section("NODAL POINTS WITH APPLIED LOADS");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int dof;
				double bcLoad,bcAngle;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();		// node number here
					if(nodeNum>nnodes || nodeNum<1)
					{	s.close();
						throw new Exception("Found nodal bondary condition with unexpected node number.");
					}
					dof=s.nextInt();
					bcLoad=s.nextDouble();
					bcAngle=0.;
					
					// looks for angle, but it shouldnot be read
					//for(i=0;i<gridBCs.size();i++)
					//{	GridDispBC obj=gridBCs.get(i);
					//	if(obj.rotates(nodeNum))
					//		bcAngle+=obj.getValue();
					//}
					
					// value not scaled, but not currently visualized
					addNodalLoadBC(nodeNum,dof,bcLoad,bcAngle);
				}
				s.close();
				s=null;
			}
		}
		
		//---------------------------------------------------------------
		// element face stresses (FEA)
		if(isFEAAnalysis())
		{	bcs=section("FACES WITH APPLIED STRESS");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int face;
				String orient;
				double str1,str2,str3;
				while(s.hasNextInt())
				{	elemNum=s.nextInt();		// element number here
					if(elemNum>nelems || elemNum<1)
					{	s.close();
						throw new Exception("Found element boundary condition with unexpected element number.");
					}
					face=s.nextInt();			// face number
					orient=s.next();			// "Normal" or "Shear"
					str1=s.nextDouble();		// 3stresses
					str2=s.nextDouble();
					if(elements.get(elemNum-1).hasMidSideNodes())
						str3=s.nextDouble();
					else
						str3=0.;
					
					// stress not scaled, but currently not visualized
					addElementBC(elemNum,face,orient,str1,str2,str3);
				}
				s.close();
				s=null;
			}
		}
		
		//----------------------------------------------------------
		// global archive file
		if(isMPMAnalysis())
		{	String globalResults=section("ARCHIVED GLOBAL RESULTS");
			if(globalResults.length()>0)
			{	Scanner ss=new Scanner(globalResults);
				ss.useDelimiter("\\r\\n|\\n|\\r");
				ss.useLocale(Locale.US);
				ss.next();
				ss.next();
				
				// global results file name
				Scanner ssline=new Scanner(ss.next());
				ssline.useDelimiter(": ");
				ssline.useLocale(Locale.US);
				ssline.next();
				line=ssline.next();
				setGlobalPath(file.getParentFile(),line);
				
				ss.close();
				ssline.close();
			}
			else
				globalArchive=null;
		}
		
		//----------------------------------------------------------
		// read archive list
		if(isMPMAnalysis())
		{	String archives=section("ARCHIVED ANALYSIS RESULTS");
			s=new Scanner(archives);
			s.useDelimiter("\\r\\n|\\n|\\r");
			s.useLocale(Locale.US);
			s.next();
			s.next();
			
			// root file name
			sline=new Scanner(s.next());
			sline.useDelimiter(": ");
			sline.useLocale(Locale.US);
			sline.next();
			line=sline.next();
			sline.close();
			sline=null;
			endIndex=line.lastIndexOf('/');
			if(endIndex>=0)
				archDir=line.substring(0,endIndex+1);
			else
				archDir="";
			setPath(file.getParentFile(),archDir);
			String ptsPath = line.substring(0, line.length()-1)+"_PtSizes.txt";
			String dimPath = line.substring(0, line.length()-1)+"_PtDims.txt";
			
			// archive format and check it
			sline=new Scanner(s.next());
			sline.useDelimiter(": ");
			sline.useLocale(Locale.US);
			sline.next();
			setArchFormat(sline.next());
			sline.close();
			sline=null;
			if(archFormat.length()>ReadArchive.ARCH_MAXMPMITEMS)
			{	for(int ii=ReadArchive.ARCH_MAXMPMITEMS;ii<archFormat.length();ii++)
				{	if(archFormat.charAt(ii)=='Y')
					{	s.close();
						throw new Exception("This archive includes data not supported by this version of NairnFEAMPMViz");
					}
				}
			}
			
			// crack archive format
			sline=new Scanner(s.next());
			sline.useDelimiter(": ");
			sline.useLocale(Locale.US);
			sline.next();
			setCrackFormat(sline.next());
			if(crackFormat.length()>ReadArchive.ARCH_MAXCRACKITEMS)
			{	for(int ii=ReadArchive.ARCH_MAXCRACKITEMS;ii<crackFormat.length();ii++)
				{	if(crackFormat.charAt(ii)=='Y')
					{	s.close();
						sline.close();
						throw new Exception("This archive includes data not supported by this version of NairnFEAMPMViz");
					}
				}
			}
			sline.close();
			sline=null;
			
			// archive files
			s.next();
			s.next();
			s.next();
			while(s.hasNext())
			{   line=s.next();
				if(line.length()==0) break;
				strChar=line.charAt(0);
				if(strChar=='#') continue;
				
				// get time and file name (split any number of spaces)
				words=line.trim().split(" +");
				if(words.length<3) continue;
				
				// errors might be non-comment text
				try
				{	Double atime=new Double(words[1]);
					addArchiveFile(atime.doubleValue()*units.altTimeScale(),words[2]);
				}
				catch (Exception e)
				{	System.out.println("Invalid line with archived files:\n   "+line);
				}
			}
			s.close();
			s=null;
			
			// error in no files were found
			if(archiveTimes.size()<1)
				throw new Exception("None of the archived results files could be found.");
			
			// look for point dimensions (new code)
			pointDims = new ArrayList<Vector3>(numMps);
			File dimFile = new File(file.getParentFile(),dimPath);
			if(dimFile.exists())
			{	// read the file
				FileReader fr=new FileReader(dimFile);
				char [] buffer=new char [(int)dimFile.length()];
				fr.read(buffer);
				String ptsSection=new String(buffer);
				fr.close();
				
				if(ptsSection.length()>10)
				{	s=new Scanner(ptsSection );
					s.useDelimiter("\\r\\n|\\n|\\r");
					s.useLocale(Locale.US);
					s.next();
					s.next();
					s.next();
					
					double lscale = 0.5*units.lengthScale();		// 0.5 to get radius
					double sx,sy,sz;
					
					while(s.hasNext())
					{	Scanner pline=new Scanner(s.next());
						pline.useLocale(Locale.US);
						pline.nextInt();
						// store point number, and x-y-z sizes
						sx = lscale*pline.nextDouble();
						sy = lscale*pline.nextDouble();
						sz = lscale*pline.nextDouble();
						pointDims.add(new Vector3(sx,sy,sz));
						pline.close();
					}
					s.close();
					s=null;
				}
			}
			
			// read point sizes (old code membranes only)
			if(pointDims==null)
			{	// store particle sizes (assumes a regular grid was found
				double msx=(float)0.25*cellMinSide;			// radius of points when 2 per axis
				double msy=msx,msz=msx;
				if(xscale>0)
				{	msx *= xscale;
					msy *= yscale;
					msz *= zscale;
				}
				for(int p=0;p<numMps;p++)
					pointDims.add(new Vector3(msx,msy,msz));

				// look for file for membrane particles
				File ptsFile = new File(file.getParentFile(),ptsPath);
				if(ptsFile.exists())
				{	// read the file
					FileReader fr=new FileReader(ptsFile);
					char [] buffer=new char [(int)ptsFile.length()];
					fr.read(buffer);
					String ptsSection=new String(buffer);
					fr.close();
				
					if(ptsSection.length()>10)
					{	s=new Scanner(ptsSection );
						s.useDelimiter("\\r\\n|\\n|\\r");
						s.useLocale(Locale.US);
						s.next();
						s.next();
						s.next();
					
						double sx,sy,sz;
						while(s.hasNext())
						{	Scanner pline=new Scanner(s.next());
							pline.useLocale(Locale.US);
							// store point number, and x-y-z sizes
							int pnum=pline.nextInt();
							sx = 2.*msx*pline.nextDouble();
							sy = 2.*msy*pline.nextDouble();
							sz = 2.*msz*pline.nextDouble();
							pointDims.set(pnum-1, new Vector3(sx,sy,sz));
							pline.close();
						}
						s.close();
						s=null;
					}
				}
			}
		}
			
		//---------------------------------------------------------------
		// FEA specific reslts
		if(isFEAAnalysis())
		{	// FEA nodal displacements
			bcs=section("NODAL DISPLACEMENTS");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int numFound;
				NodalPoint anode;
				while(s.hasNextInt())
				{	numFound=s.nextInt();		// node number
					if(numFound>nnodes || numFound<1)
					{	s.close();
						throw new Exception("Found nodal displacement with unexpected node number.");
					}
					anode=nodes.get(numFound-1);
					anode.dispx=s.nextDouble()*units.lengthScale();
					anode.dispy=s.nextDouble()*units.lengthScale();
					dxmin=Math.min(dxmin,anode.x+anode.dispx);
					dxmax=Math.max(dxmax,anode.x+anode.dispx);
					dymin=Math.min(dymin,anode.y+anode.dispy);
					dymax=Math.max(dymax,anode.y+anode.dispy);
				}
				feaArchFormat[ReadArchive.ARCH_FEADisplacements]='Y';
				s.close();
				s=null;
			}
			
			for(i=0;i<elements.size();i++)
				elements.get(i).setElemPath(this,true);
				
			// FEA nodal stresses
			bcs=section("AVERAGE NODAL STRESSES");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				int numFound;
				NodalPoint anode;
				while(s.hasNextInt())
				{	numFound=s.nextInt();		// node number
					if(numFound>nnodes || numFound<1)
					{	s.close();
						throw new Exception("Found nodal displacement with unexpected node number.");
					}
					anode=nodes.get(numFound-1);
					anode.sigxx=s.nextDouble()*units.feaStressScale();
					anode.sigyy=s.nextDouble()*units.feaStressScale();
					anode.sigzz=s.nextDouble()*units.feaStressScale();
					anode.sigxy=s.nextDouble()*units.feaStressScale();
				}
				feaArchFormat[ReadArchive.ARCH_FEAAvgStress]='Y';
				s.close();
				s=null;
			}
			
			// FEA element energies
			bcs=section("STRAIN ENERGIES IN ELEMENTS");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				ElementBase aelem;
				totalEnergy = 0.;
				while(s.hasNextInt())
				{	elemNum=s.nextInt();		// element number
					if(elemNum>nelems || elemNum<1)
					{	s.close();
						throw new Exception("Found element energy with unexpected element number.");
					}
					aelem=elements.get(elemNum-1);
					aelem.energy=s.nextDouble()*units.energyScale();
					totalEnergy += aelem.energy;
				}
				feaArchFormat[ReadArchive.ARCH_FEAElemEnergy]='Y';
				s.close();
				s=null;
			}
			
			// NODAL FORCES AND ELEMENT STRESSES (FEA)
			bcs=section("NODAL FORCES");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				s.useLocale(Locale.US);
				ElementBase aelem=null;
				boolean hasForces=false,hasStresses=false;
				while(s.hasNext())
				{	// check first line for what is here
					String word1=s.next();
					if(word1.equals("Element"))
					{	aelem=null;
						s.next();
						word1=s.next();
					}
					else
					{	if(aelem==null)
						{	s.close();
							throw new Exception("Found element results outside element definition.");
						}
					}
					
					// skip to first number
					while(!s.hasNextDouble()) s.next();
					
					// word1 is "F_" or "sig_"
					if(word1.charAt(0)=='F')
					{	hasForces=true;
					
						// read element is needed
						if(aelem==null)
						{	elemNum=s.nextInt();
							if(elemNum<1 || elemNum>nelems)
							{	s.close();
								throw new Exception("Found element results with unexpected element number.");
							}
							aelem=elements.get(elemNum-1);
						}
						
						// read nodal forces
						for(i=0;i<aelem.getNumberNodes();i++)
						{	s.nextInt();			// skip node number
							aelem.setForces(s.nextDouble()*units.forceScale(),
									s.nextDouble()*units.forceScale(),i);
						}
					}
					
					else if(word1.charAt(0)=='s')
					{	hasStresses=true;
					
						// read element is needed
						if(aelem==null)
						{	elemNum=s.nextInt();
							if(elemNum<1 || elemNum>nelems)
							{	s.close();
								throw new Exception("Found element results with unexpected element number.");
							}
							aelem=elements.get(elemNum-1);
						}
						
						// read nodal stresses
						for(i=0;i<aelem.getNumberNodes();i++)
						{	s.nextInt();			// skip node number
							if(word1.charAt(3)=='x' || word1.charAt(3)=='r')
							{	aelem.setXYStresses(s.nextDouble()*units.feaStressScale(),
										s.nextDouble()*units.feaStressScale(),
										s.nextDouble()*units.feaStressScale(),i);
							}
							else
							{	aelem.set3DStresses(s.nextDouble()*units.feaStressScale(),
									s.nextDouble()*units.feaStressScale(),
									s.nextDouble()*units.feaStressScale(),i);
							}
						}
					}
					
					// skip underlying "----"
					if(s.hasNext()) s.next();
				}
				
				// save what was found (if even on selected ones were found)
				if(hasForces) feaArchFormat[ReadArchive.ARCH_FEAElemForce]='Y';
				if(hasStresses) feaArchFormat[ReadArchive.ARCH_FEAElemStress]='Y';
				
				s.close();
				s=null;
			}
		}
		
		// all decoded correctly
		setName(file);
		
		// other calculations needed for eventual plotting
		if(isMPMAnalysis()) recSize=ReadArchive.getRecordSize(this);
		cellMinSide=ElementBase.getCellMinSide(this);
	}
	
	// read the currently selected archive (unless it is already loaded)
	public void readSelectedArchive(int newArchive) throws Exception
	{	// if already read then done
		if(newArchive==currentArchive) return;
		
		// read archive file
		currentArchive=newArchive;
		ReadArchive.load(this,archives.get(currentArchive));
		
	}
	
	// open archive by index and return its byte buffer
	public ByteBuffer openSelectedArchive(int newArchve) throws Exception
	{	return ReadArchive.openBuffer(this,archives.get(newArchve));
	}
	
	//-----------------------------------------------------------------
	// Most parts of a document are held in a series of ArrayLists
	// These methods add items to a list or clear all lists
	//-----------------------------------------------------------------
	
	// clear all sections
	// if doAll is false don't clear sections (used on rescaling results)
	public void clear(boolean doAll)
	{	if(doAll)
		{   sectionTitle.clear();
			sectionText.clear();
		}
		nodes.clear();
		elements.clear();
		gridBCs.clear();
		particleBCs.clear();
		mpmPoints.clear();
		archives.clear();
		archiveTimes.clear();
		currentArchive=-1;
	}
	
	// add new text section
	public void add(String title,String section)
	{   sectionTitle.add(title);
		sectionText.add(section);
	}
	
	// replace named text section
	public void replace(String title,String section)
	{   int index;
		for(index=0;index<sectionTitle.size();index++)
		{   if(title.equals(sectionTitle.get(index)))
				break;
		}
		if(index>=sectionTitle.size())
			add(title,section);
		else
		{	sectionText.remove(index);
			sectionText.add(index,section);
		}
	}
	
	// add 2D nodal point
	public void addNode(int nodeNum,double xpt,double ypt)
	{	if(nodeNum==1)
		{	xmin=xpt;
			xmax=xpt;
			ymin=ypt;
			ymax=ypt;
		}
		else
		{	if(xpt<xmin) xmin=xpt;
			else if(xpt>xmax) xmax=xpt;
			if(ypt<ymin) ymin=ypt;
			else if(ypt>ymax) ymax=ypt;
		}
		nodes.add(new NodalPoint(nodeNum,xpt,ypt));
	}

	// add 3D nodal point
	public void addNode(int nodeNum,double xpt,double ypt,double zpt)
	{	if(nodeNum==1)
		{	xmin=xpt;
			xmax=xpt;
			ymin=ypt;
			ymax=ypt;
			zmin=zpt;
			zmax=zpt;
		}
		else
		{	if(xpt<xmin) xmin=xpt;
			else if(xpt>xmax) xmax=xpt;
			if(ypt<ymin) ymin=ypt;
			else if(ypt>ymax) ymax=ypt;
			if(zpt<zmin) zmin=ypt;
			else if(zpt>zmax) zmax=ypt;
		}
		nodes.add(new NodalPoint(nodeNum,xpt,ypt,zpt));
	}
	// add element
	public void addElement(int elemNum,int elemID,int [] nds,int matID,double angle,double thickness)
	{
		ElementBase newElem;
		
		switch(elemID)
		{	case ElementBase.FOUR_NODE_ISO:
				newElem=new FourNodeIsoparam(elemNum,nds);
				break;
			
			case ElementBase.EIGHT_NODE_ISO:
				newElem=new EightNodeIsoparam(elemNum,nds);
				break;
			
			case ElementBase.CST:
				newElem=new CSTriangle(elemNum,nds);
				break;
			
			case ElementBase.ISO_TRIANGLE:
				newElem=new IsoTriangle(elemNum,nds);
				break;
				
			case ElementBase.EIGHT_NODE_ISO_BRICK:
				newElem=new EightNodeBrick(elemNum,nds);
				break;				
			
			case ElementBase.LINEAR_INTERFACE:
				feaArchFormat[ReadArchive.ARCH_Interfaces]='Y';
				newElem=new LinearInt(elemNum,nds);
				break;
			
			case ElementBase.QUAD_INTERFACE:
				feaArchFormat[ReadArchive.ARCH_Interfaces]='Y';
				newElem=new QuadInt(elemNum,nds);
				break;
			
			case ElementBase.LAGRANGE_2D:
				newElem=new Lagrange2D(elemNum,nds);
				break;
				
			default:
				return;
		}
		
		elements.add(newElem);
		if(isFEAAnalysis()) newElem.setFEAProperties(matID,angle,thickness);
		newElem.setElemPath(this,false);
	}
	
	// add grid BC
	public void addGridBC(int nodeNum,int dof,int bcID,double bcVal,double bcArg,double bcAngle,double bcAngle2)
	{	
		gridBCs.add(new GridDispBC(nodeNum,dof,bcID,bcVal,bcArg,bcAngle,bcAngle2));
	}
	
	// add nodal load BC
	public void addNodalLoadBC(int nodeNum,int dof,double bcVal,double bcAngle)
	{	
		gridBCs.add(new NodalLoadBC(nodeNum,dof,bcVal,bcAngle));
	}
	
	// add element face BC
	public void addElementBC(int elemNum,int face,String orient,double str1,double str2,double str3)
	{
		gridBCs.add(new ElementStressBC(elemNum,face,orient,str1,str2,str3));
	}

	// add particle BC
	public void addParticleBC(int partNum,int dof,int bcFace,int bcID,double bcVal,double bcArg)
	{	
		particleBCs.add(new ParticleBC(partNum,dof,bcFace,bcID,bcVal,bcArg));
	}
	
	// add archive file if it exists
	public void addArchiveFile(double time,String name)
	{   // do not add unless the file exists
	    File archive=new File(path,name);
		if(!archive.exists()) return;
		archives.add(archive);
	    archiveTimes.add(new Double(time));
	}
	
	//-----------------------------------------------------------------
	// Utilities for decoding text when reading a file
	//-----------------------------------------------------------------
	
	// find index to start of line after line containing substring
	public int findNextLine(String data,String lookStr)
	{
		int endIndex=data.indexOf(lookStr);
		if(endIndex<0) return endIndex;
		
		// find end of the line
		char strChar='\n';
		while(endIndex<data.length())
		{   strChar=data.charAt(endIndex);
			if(strChar=='\n' || strChar=='\r') break;
			endIndex++;
		}
		
		// check for carriage return and line feed
		endIndex++;
		if(strChar=='\r' && endIndex<data.length())
		{   strChar=data.charAt(endIndex);
			if(strChar=='\n') endIndex++;
		}
		return endIndex;
	}
	
	//-----------------------------------------------------------------
	// Accessor Methods
	//-----------------------------------------------------------------
	
	// return string of a section by index of the section
	public String section(int index)
	{   if(index<0 || index>=sectionText.size()) return "";
		return sectionText.get(index);
	}
	
	// return string of section whose fill name starts with variable name
	public String section(String name)
	{   int index;
		for(index=0;index<sectionTitle.size();index++)
		{	if(sectionTitle.get(index).startsWith(name))
				return section(index);
		}
		return "";
	}
	
	// return full name of a section
	
	// set name of the mpm file
	public void setName(File mpmFile)
	{	name=new String(mpmFile.getName());
		fullPath=new String(mpmFile.getPath());
	}
	
	// set root path to archive files, parent is directory of mpm file
	public void setPath(File parent,String archDir) { path=new File(parent,archDir); }
	
	// set path to global file, but only if it exists
	public void setGlobalPath(File parent,String globalPath)
	{	globalArchive=new File(parent,globalPath);
		if(!globalArchive.exists())
		{	globalArchive=null;
			return;
		}
	}
	
	// set archive formats
	public void setArchFormat(String format) { archFormat=new String(format+"NNNNNNNNNNN"); }
	public void setCrackFormat(String format) { crackFormat=new String(format+"NNNNNNNNNNN"); }
	
	// mesh bounds
	public Rectangle2D.Double getMeshBounds(boolean displaced)
	{	if(displaced)
			return new Rectangle2D.Double(dxmin,dymin,dxmax-dxmin,dymax-dymin);
		else
			return new Rectangle2D.Double(xmin,ymin,xmax-xmin,ymax-ymin);
	}
	
	// return energy as string (for scripts)
	public String getEnergy()
	{	if(feaArchFormat[ReadArchive.ARCH_FEAElemEnergy]=='Y')
			return JNUtilities.formatDouble(totalEnergy);
		else
			return null;
	}
	
	public boolean isFEAAnalysis() { return np<BEGIN_MPM_TYPES; }
	public boolean isMPMAnalysis() { return np>=BEGIN_MPM_TYPES; }
	public boolean isAxisymmetric() { return (np==AXI_SYM) || (np==AXI_SYM_MPM) ; }
	public boolean is3D() { return np==THREED_MPM; }
	public double currentTime() { return currentArchive>=0 ? archiveTimes.get(currentArchive) : 0.; }
	
	// set controller
	public void setDocController(DocViewer dc) { docCtrl=dc; }
	public DocViewer getDocController() { return docCtrl; }
	
	// gather time plot data as a scripting object call
	// Expression starts @obj.timeplot.component.
	//    for crack plots add crackNum.tipNum, except length and debond length add only crackNum
	public String collectTimePlotData(String[] atoms,HashMap<String, Double> variables,
			HashMap<String, String> variablesStrs) throws Exception
	{
		// get component
		if(atoms.length<3) throw new Exception("timeplot property has too few parameters");
		
		// get plot component from quote string, string, or string variable
		atoms[2] = CmdViewer.atomString(atoms[2],variablesStrs);
		
		// trap global data plot
		if(atoms[2].equals("global"))
		{	if(atoms.length<4) throw new Exception("timeplot property missing global quantity");
			String quant = CmdViewer.atomString(atoms[3],variablesStrs);
			int dataLines=0;
			String setName=null;
			ArrayList<String> table=new ArrayList<String>(100);
			try
			{	FileReader fr=new FileReader(globalArchive);
				char [] buffer=new char [(int)globalArchive.length()];
				fr.read(buffer);
				fr.close();
				
				// scanner for lines
				Scanner s=new Scanner(new String(buffer));
				s.useDelimiter("[\\n\\r]");
				while(s.hasNext())
				{	String line=s.next();
					if(line.charAt(0)!='#')
					{	dataLines++;
						table.add(line);
					}
					else if(line.length()>7)
					{	String tcmd = line.substring(0,8);
						if(tcmd.equals("#setName")) setName=line;
					}
				}
				s.close();
			}
			catch (Exception e)
			{	throw new Exception("Could not load global results file:\n   " + e.getMessage());
			}
			
			// check table header
			if(setName==null || dataLines==0)
				throw new Exception("global file column labels not found or file has no data lines");
			
			// find column to read
			int column=-1,col=0;
			Scanner s=new Scanner(setName);
			s.useDelimiter("\t");
			while(s.hasNext())
			{	String colVal = s.next();
				// remove quotes
				int clen = colVal.length();
				if(clen>1)
				{	if(colVal.charAt(0)=='"' && colVal.charAt(clen-1)=='"')
						colVal = colVal.substring(1,clen-1);
				}
				if(quant.equals(colVal))
				{	column = col;
					break;
				}
				col++;
			}
			s.close();
			if(column<0)
				throw new Exception("The quantity '"+quant+"' not found in global results");
			
			// collect
			StringBuffer gdata = new StringBuffer();
			int row;
			for(row=0;row<table.size();row++)
			{	String dataRow=table.get(row);
				Scanner r=new Scanner(dataRow);
				r.useDelimiter("\t");
				
				// get x value
				if(!r.hasNext())
				{	r.close();
					continue;
				}
				double xg = r.nextDouble();
				
				// skip intervening columns
				int j=column;
				while(j>1)
				{	if(r.hasNext()) r.next();
					j--;
				}
				
				// get x value
				if(!r.hasNext())
				{	r.close();
					continue;
				}
				double yg = r.nextDouble();
				
				// append data
				gdata.append(xg+" "+yg+"\n");
				
				r.close();
				
			}

			return gdata.toString();
		}
		
		// other plot components
		int cmpnt = -1;
		if(atoms[2].equals("J1"))
			cmpnt = PlotQuantity.MPMJ1;
		else if(atoms[2].equals("J2"))
			cmpnt = PlotQuantity.MPMJ2;
		else if(atoms[2].equals("KI"))
			cmpnt = PlotQuantity.MPMKI;
		else if(atoms[2].equals("KII"))
			cmpnt = PlotQuantity.MPMKII;
		else if(atoms[2].equals("CrackRelease"))
			cmpnt = PlotQuantity.MPMCRACKRELEASE;
		else if(atoms[2].equals("CrackAbsorb"))
			cmpnt = PlotQuantity.MPMCRACKABSORB;
		else if(atoms[2].equals("Length"))
			cmpnt = PlotQuantity.MPMLENGTH;
		else if(atoms[2].equals("DebondLength"))
			cmpnt = PlotQuantity.MPMDEBONDLENGTH;
		else if(atoms[2].equals("TipNCOD"))
			cmpnt = PlotQuantity.MPMNORMALCTOD;
		else if(atoms[2].equals("TipSCOD"))
			cmpnt = PlotQuantity.MPMSHEARCTOD;
		else if(atoms[2].equals("DebondNCOD"))
			cmpnt = PlotQuantity.MPMDEBONDNCTOD;
		else if(atoms[2].equals("DebondSCOD"))
			cmpnt = PlotQuantity.MPMDEBONDSCTOD;
		
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(archives.size());
		ArrayList<Double> y=new ArrayList<Double>(archives.size());
		StringBuffer pdata = new StringBuffer();
		int i,crackNum,tipNum;
		switch(cmpnt)
		{	case PlotQuantity.MPMJ1:
			case PlotQuantity.MPMJ2:
			case PlotQuantity.MPMKI:
			case PlotQuantity.MPMKII:
			case PlotQuantity.MPMCRACKRELEASE:
			case PlotQuantity.MPMCRACKABSORB:
			case PlotQuantity.MPMNORMALCTOD:
			case PlotQuantity.MPMSHEARCTOD:
			case PlotQuantity.MPMDEBONDNCTOD:
			case PlotQuantity.MPMDEBONDSCTOD:
				// get crack numberand tip number
				if(atoms.length<5) throw new Exception("timeplot property has too few parameters");
				crackNum = CmdViewer.atomInt(atoms[3],variables);
				tipNum = CmdViewer.atomInt(atoms[4],variables);
				if(tipNum<0 || tipNum>1)
					throw new Exception("timeplot crack tip number must be 0 (start) or 1 (end)");
				getTimeCrackData(null,cmpnt,crackNum,tipNum,x,y);
				for(i=0;i<x.size();i++)
					pdata.append(x.get(i)+" "+y.get(i)+"\n");
				break;
				
			case PlotQuantity.MPMLENGTH:
			case PlotQuantity.MPMDEBONDLENGTH:
				// these crack values don't need tip number
				if(atoms.length<4) throw new Exception("timeplot property has too few parameters");
				crackNum = CmdViewer.atomInt(atoms[3],variables);
				getTimeCrackData(null,cmpnt,crackNum,0,x,y);
				for(i=0;i<x.size();i++)
					pdata.append(x.get(i)+" "+y.get(i)+"\n");
				break;
			
			default:
				throw new Exception(atoms[2]+" is an unsupported timeplot quantity");
		}
		
		return pdata.toString();
	}
	
	// Get crack data for a time plot
	// ctrls for progress bar (or null to skip), crkCmpnt (crack time plot quantity),
	// crack and tip number.
	// Plot returned in x and y
	public void getTimeCrackData(ControlPanel ctrls,int crkCmpnt,int crackNum,int tipNum,
				ArrayList<Double> x,ArrayList<Double> y) throws Exception
	{
		// settings
		int i,npts=archives.size();
		ArrayList<Integer> crackEnds=new ArrayList<Integer>(20);

		// variables while decoding
		byte[] version=new byte[4];
		ByteBuffer bb;
		CrackSegment seg=new CrackSegment();
		int lastoffset,offset,tipOffset,endOffset;
		boolean foundTip;
		short matnum;
		double yvalue;
		Point2D.Double pt,lastPt,cod;
		
		// format
		char[] crackOrder=new char[ReadArchive.ARCH_MAXCRACKITEMS];
		crackFormat.getChars(0,ReadArchive.ARCH_MAXCRACKITEMS,crackOrder,0);
		
		// loop over all archives
		for(i=0;i<npts;i++)
		{	// adjust progress bar
			if(ctrls!=null) ctrls.setProgress(i+1);
		
			// open file (try to continue on errors
			try
			{	bb=openSelectedArchive(i);
			}
			catch(Exception bbe)
			{	continue;
			}
			
			// check version
			bb.get(version);
			int headerLength=4;
			int vernum=version[3]-'0';
			if(vernum>=4)
				headerLength=64;
			else if(vernum!=3)
				throw new Exception("Archive file is too old for this tool");
			int nummpms=(int)((bb.remaining()+4-headerLength)/recSize);
			
			// find the first crack (matnum>0), remember ends
			offset=headerLength+(nummpms-1)*recSize+ReadArchive.sizeofInt+ReadArchive.sizeofDouble;
			lastoffset=offset;
			crackEnds.clear();
			crackEnds.add(new Integer(lastoffset));
			while(offset>0)
			{	bb.position(offset);
				matnum=bb.getShort();
				if(matnum>0) break;
				offset-=recSize;
				if(matnum==-1) crackEnds.add(new Integer(offset));
			}
			
			// find start or end of desired crack (if it exists)
			if(crackNum>crackEnds.size()-1) continue;
			if(tipNum==CrackSelector.CRACK_START || crkCmpnt==PlotQuantity.MPMLENGTH
						|| crkCmpnt==PlotQuantity.MPMDEBONDLENGTH)
			{	Integer offObj=crackEnds.get(crackEnds.size()-crackNum);
				offset=offObj.intValue()+recSize;
				offObj=crackEnds.get(crackEnds.size()-crackNum-1);
				endOffset=offObj.intValue();
			}
			else
			{	Integer offObj=crackEnds.get(crackEnds.size()-crackNum-1);
				offset=offObj.intValue();
				offObj=crackEnds.get(crackEnds.size()-crackNum);
				endOffset=offObj.intValue()+recSize;
			}
			offset-=(ReadArchive.sizeofInt+ReadArchive.sizeofDouble);
			endOffset-=(ReadArchive.sizeofInt+ReadArchive.sizeofDouble);
			
			// read segment
			bb.position(offset);
			seg.readRecord(bb,crackOrder,units);
			
			// crack tip properties
			switch(crkCmpnt)
			{   case PlotQuantity.MPMJ1:
					yvalue=seg.J1;
					break;
					
				case PlotQuantity.MPMJ2:
					yvalue=seg.J2;
					break;
					
				case PlotQuantity.MPMKI:
					yvalue=seg.KI;
					break;

				case PlotQuantity.MPMKII:
					yvalue=seg.KII;
					break;

				case PlotQuantity.MPMLENGTH:
				case PlotQuantity.MPMDEBONDLENGTH:
					yvalue=0.;
					double bonded=0.;
					lastPt=seg.getMedianPosition();
					while(true)
					{   offset+=recSize;
						if(offset>lastoffset) break;
						bb.position(offset);
						seg.readRecord(bb,crackOrder,units);
						if(seg.startFlag==-1) break;
						pt=seg.getMedianPosition();
						double segLength=Math.sqrt((pt.x-lastPt.x)*(pt.x-lastPt.x) +
										(pt.y-lastPt.y)*(pt.y-lastPt.y));
						yvalue+=segLength;
						if(seg.tractionMaterial>0) bonded+=segLength;
						lastPt=pt;
					}
					if(crkCmpnt==PlotQuantity.MPMDEBONDLENGTH) yvalue-=bonded;
					break;
				
				case PlotQuantity.MPMCRACKRELEASE:
					yvalue=seg.release;
					break;
				
				case PlotQuantity.MPMCRACKABSORB:
					yvalue=seg.absorb;
					break;
				
				case PlotQuantity.MPMDEBONDNCTOD:
				case PlotQuantity.MPMDEBONDSCTOD:
					tipOffset=offset;
					foundTip=true;
					pt=seg.getPt();
					int tlCount=0;
					lastPt=new Point2D.Double(0.,0.);
					cod=new Point2D.Double(0.,0.);
					// scan to end of debond zone from this tip
					while(seg.tractionMaterial>0)
					{	lastPt=pt;
						cod=seg.getCOD();
						pt=seg.getPt();
						tlCount++;
						if(tipNum==CrackSelector.CRACK_START)
						{	tipOffset+=recSize;
							if(tipOffset>endOffset)
							{	foundTip=false;
								break;
							}
						}
						else
						{	tipOffset-=recSize;
							if(tipOffset<endOffset)
							{	foundTip=false;
								break;
							}
						}
						bb.position(tipOffset);
						seg.readRecord(bb,crackOrder,units);
					}
					
					// not found in the crack (entire crack is traction law)
					if(!foundTip)
					{	yvalue=0.;
						break;
					}
					
					// if tlCount==0, then no traction at this crack tip, so just fall through and use regular crack tip cod
					//   otherwise do calculations
					//	Here pt and cod are at the debond tip. lastPt is at previous traction law or at tip if tlCount==1
					if(tlCount>0)
					{	pt=seg.getPt();
						double dx=lastPt.x-pt.x;
						double dy=lastPt.y-pt.y;
						double norm=Math.sqrt(dx*dx+dy*dy);
					
						if(crkCmpnt==PlotQuantity.MPMDEBONDNCTOD)
							yvalue=(-cod.x*dy + cod.y*dx)/norm;
						else
							yvalue=(cod.x*dx + cod.y*dy)/norm;
						break;
					}
				
				case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
					cod=seg.getCOD();
					pt=seg.getMedianPosition();

					// read previous segment
					offset = tipNum==CrackSelector.CRACK_START ? offset+recSize : offset-recSize;
					bb.position(offset);
					seg.readRecord(bb,crackOrder,units);
					lastPt=seg.getMedianPosition();
					double dx=pt.x-lastPt.x;
					double dy=pt.y-lastPt.y;
					double norm=Math.sqrt(dx*dx+dy*dy);
					
					if(crkCmpnt==PlotQuantity.MPMNORMALCTOD || crkCmpnt==PlotQuantity.MPMDEBONDNCTOD)
						yvalue=(-cod.x*dy + cod.y*dx)/norm;
					else
						yvalue=(cod.x*dx + cod.y*dy)/norm;
					break;
				
				default:
					yvalue=0.;
					break;
			}
			
			// add crack value to the  plot
			x.add(new Double(archiveTimes.get(i)));
			y.add(new Double(yvalue));
		}
		
		if(x.size()==0)
			throw new Exception("No data found for that plot quantity");
	}

	//-----------------------------------------------------------------
	// Standard methods to support JTable to display file sections
	//-----------------------------------------------------------------
	public int getRowCount() { return sectionTitle.size();}
	public int getColumnCount() { return 1; }
	public Object getValueAt(int row,int column) { return sectionTitle.get(row); }
	public String getColumnName(int column) { return new String("Sections"); }
	public Class<?> getColumnClass(int column) { return String.class; }
	
}

