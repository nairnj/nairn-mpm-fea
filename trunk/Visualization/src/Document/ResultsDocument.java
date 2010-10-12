/*******************************************************************
	ResultsDocument.java
	NairnFEAMPMViz

	Created by John Nairn on Sat Mar 06 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

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
	public String archFormat,crackFormat;
	public char[] feaArchFormat = {'N','N','N','N','N','N' };
	public double xmin,xmax,ymin,ymax;					// mesh bounds
	public double dxmin,dxmax,dymin,dymax;				// displaced mesh bounds
	public double cellMinSide,xscale=1.,yscale=1.;
	public int np;
	public DocViewer docCtrl;
	public int recSize;
	
	// units
	public double lengthScale=1.0;			// mm
	public String distU="mm";
	public double timeScale=1.0;			// ms
	public String timeU="ms";
	
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
		Scanner s,sline;
		
		//----------------------------------------------------------
		// Expected number of nodes and elements
		String summary=section("NODES AND ELEMENTS (Background Grid)");
		if(summary.length()==0)
			summary=section("NODES AND ELEMENTS");
		int nnodes,nelems;
		try
		{	// get nodes and elements
			s=new Scanner(summary).useDelimiter("[\\n\\r]");
			s.next();
			s.next();
			sline=new Scanner(s.next());
			sline.next();
			nnodes=sline.nextInt();
			sline.next();
			nelems=sline.nextInt();
		
			// Options are 2D Plane Strain Analysis, 2D Plane Stress Analysis, Axisymmetric Analysis,
			//	2D Plane Strain MPM Analysis, 2D Plane Stress MPM Analysis, 3D MPM Analysis
			sline=new Scanner(s.next());
			sline.next();
			sline.next();
			sline.next();
			sline.next();
			String word=sline.next();
			if(word.equals("3D"))
				np=THREED_MPM;
			else if(word.equals("Axisymmetric"))
				np=AXI_SYM;
			else
			{	sline.next();
				if(sline.next().equals("Strain"))
					np=PLANE_STRAIN;
				else
					np=PLANE_STRESS;
				if(sline.next().equals("MPM"))
					np = (np==PLANE_STRAIN) ? PLANE_STRAIN_MPM : PLANE_STRESS_MPM;
			}
		}
		catch(NoSuchElementException e)
		{	throw new Exception("Could not decode analysis type for this file");
		}
		if(is3D())
			throw new Exception("This tool cannot visualize 3D results; see help information for other options.");
		
		//----------------------------------------------------------
		// Nodal Point Coordinates
		String ndst=section("NODAL POINT COORDINATES (in mm)");
		lineStart=findNextLine(ndst,"-----");
		if(lineStart<0)
			throw new Exception("Error decoding nodal point coordinates");
		s=new Scanner(ndst.substring(lineStart,ndst.length()-1));
		int prevNodeNum=0,nodeNum;
		double xpt,ypt;
		while(s.hasNextInt())
		{	nodeNum=s.nextInt();
			xpt=s.nextDouble()*lengthScale;
			ypt=s.nextDouble()*lengthScale;
			if(nodeNum!=prevNodeNum+1)
				throw new Exception("Some node numbers are missing");
			addNode(nodeNum,xpt,ypt);
			prevNodeNum=nodeNum;
		}
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
		int prevElemNum=0,elemNum,elemID,i,matID=0;
		double elemThickness=1.,elemAngle=0.;
		int[] nds;
		nds=new int[8];
		while(s.hasNextInt())
		{	elemNum=s.nextInt();				// element number
			if(elemNum!=prevElemNum+1)
				throw new Exception("Some element numbers are missing");
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
					// read all nodes
					for(i=0;i<ElementBase.NodesFromType(elemID);i++)
						nds[i]=s.nextInt();
					break;
				
				default:
					throw new Exception("Element type found ("+elemID+") that is not yet supported in this tool.");
			}
			
			addElement(elemNum,elemID,nds,matID,elemAngle,elemThickness*lengthScale);
			prevElemNum=elemNum;
		}
		if(prevElemNum!=nelems)
			throw new Exception("Number of elements does not match expected number of elements.");
		
		//----------------------------------------------------------
		// defined materials
		String matls=section("DEFINED MATERIALS");
		s=new Scanner(matls).useDelimiter("[\\n\\r]");
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
			String word1=sline.next();
			String word2=sline.next();
			MaterialBase matl;
			if(word1.equals("Isotropic"))
			{	if(word2.equals("Dugdale"))
					matl=new MaterialBase(matName,MaterialBase.DUGDALE);
				else
					matl=new IsotropicMat(matName);
			}
			else if(word1.equals("Tranversely"))
			{	sline.next();
				sline.next();
				if(sline.next().equals("normal"))
					matl=new MaterialBase(matName,MaterialBase.TRANSISO1);
				else
					matl=new MaterialBase(matName,MaterialBase.TRANSISO2);
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
				matl=new RigidMaterial(matName);
			else if(word1.equals("Triangular"))
				matl=new MaterialBase(matName,MaterialBase.COHESIVEZONEMATERIAL);
			else if(word1.equals("Linear"))
				matl=new MaterialBase(matName,MaterialBase.LINEARTRACTIONMATERIAL);
			else if(word1.equals("Cubic"))
				matl=new MaterialBase(matName,MaterialBase.CUBICTRACTIONMATERIAL);
			else if(word1.equals("Elastic-Plastic"))
				matl=new MaterialBase(matName,MaterialBase.HILLPLASTIC);
			else
			{	// try to continue with unknown material type
				matl=new MaterialBase(matName,MaterialBase.UNKNOWNMATERIAL);
			}
			
			// decode material properteis
			matl.decodeData(s);
			materials.add(matl);
		}

		//----------------------------------------------------------
		// mesh boundary conditions
		String bcs=section("NODAL POINTS WITH FIXED DISPLACEMENTS");
		lineStart=findNextLine(bcs,"-----");
		if(lineStart>0 && lineStart<bcs.length())
		{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
			int dof,bcID;
			double bcVal,bcArg,bcAngle;
			while(s.hasNextInt())
			{	nodeNum=s.nextInt();
				dof=s.nextInt();
				
				if(isMPMAnalysis())
				{	bcID=s.nextInt();
					bcVal=s.nextDouble();
					bcArg=s.nextDouble();
					bcAngle=0.;
					if(!s.hasNextInt())
					{	if(s.hasNextDouble())
						{	bcAngle=s.nextDouble();
							if(bcID==BoundaryCondition.FUNCTION_VALUE)
								s.next();
						}
					}
					addGridBC(nodeNum,dof,bcID,bcVal*lengthScale,bcArg*timeScale,bcAngle);
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
					addGridBC(nodeNum,dof,bcID,bcVal*lengthScale,bcArg*timeScale,bcAngle);
				}
			}
		}
		
		//----------------------------------------------------------
		// particle boundary conditions and grid info
		if(isMPMAnalysis())
		{	bcs=section("MATERIAL POINTS WITH EXTERNAL FORCES");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
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
					addParticleBC(nodeNum,dof,bcID,bcLoad,bcArg);
				}
			}
			
			bcs=section("FULL MASS MATRIX");
			s=new Scanner(bcs).useDelimiter("[\\n\\r]");
			// scan to grid info
			String gridInfo=null;
			while(s.hasNext() && gridInfo==null)
			{	String gridLine=s.next();
				if(gridLine.length()<10) continue;
				if(gridLine.substring(0,10).equals("Orthogonal"))
				{	int beginIndex=gridLine.indexOf(":");
					if(beginIndex<0) break;
					gridInfo=gridLine.substring(beginIndex+1,gridLine.length());
				}
			}
			if(gridInfo!=null)
			{	sline=new Scanner(gridInfo).useDelimiter("[ :]");
				if(sline.hasNextDouble()) xscale=sline.nextDouble()*lengthScale;
				if(sline.hasNext()) sline.next();
				if(sline.hasNextDouble()) yscale=sline.nextDouble()*lengthScale;
				if(xscale>yscale)
				{	xscale/=yscale;
					yscale=1.;
				}
				else
				{	yscale/=xscale;
					xscale=1.;
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
				int dof;
				double bcLoad,bcAngle;
				while(s.hasNextInt())
				{	nodeNum=s.nextInt();		// node number here
					if(nodeNum>nnodes || nodeNum<1)
						throw new Exception("Found nodal bondary condition with unexpected node number.");
					dof=s.nextInt();
					bcLoad=s.nextDouble();
					bcAngle=0.;
					for(i=0;i<gridBCs.size();i++)
					{	GridDispBC obj=gridBCs.get(i);
						if(obj.rotates(nodeNum))
							bcAngle+=obj.getValue();
					}
					addNodalLoadBC(nodeNum,dof,bcLoad,bcAngle);
				}
			}
		}
		
		//---------------------------------------------------------------
		// element face stresses (FEA)
		if(isFEAAnalysis())
		{	bcs=section("FACES WITH APPLIED STRESS (MPa)");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				int face;
				String orient;
				double str1,str2,str3;
				while(s.hasNextInt())
				{	elemNum=s.nextInt();		// element number here
					if(elemNum>nelems || elemNum<1)
						throw new Exception("Found element boundary condition with unexpected element number.");
					face=s.nextInt();			// face number
					orient=s.next();			// "Normal" or "Shear"
					str1=s.nextDouble();		// 3stresses
					str2=s.nextDouble();
					if(elements.get(elemNum-1).hasMidSideNodes())
						str3=s.nextDouble();
					else
						str3=0.;
					addElementBC(elemNum,face,orient,str1,str2,str3);
				}
			}
		}
		
		//----------------------------------------------------------
		// global archive file
		if(isMPMAnalysis())
		{	String globalResults=section("ARCHIVED GLOBAL RESULTS");
			if(globalResults.length()>0)
			{	s=new Scanner(globalResults).useDelimiter("[\\n\\r]");
				s.next();
				s.next();
				
				// global results file name
				sline=new Scanner(s.next()).useDelimiter(": ");
				sline.next();
				line=sline.next();
				setGlobalPath(file.getParentFile(),line);
			}
			else
				globalArchive=null;
		}
		
		//----------------------------------------------------------
		// read archive list
		if(isMPMAnalysis())
		{	String archives=section("ARCHIVED ANALYSIS RESULTS");
			s=new Scanner(archives).useDelimiter("[\\n\\r]");
			s.next();
			s.next();
			
			// root file name
			sline=new Scanner(s.next()).useDelimiter(": ");
			sline.next();
			line=sline.next();
			endIndex=line.lastIndexOf('/');
			if(endIndex>=0)
				archDir=line.substring(0,endIndex+1);
			else
				archDir="";
			setPath(file.getParentFile(),archDir);
			System.out.println(file.getParentFile()+","+archDir);
			
			// archive format and check it
			sline=new Scanner(s.next()).useDelimiter(": ");
			sline.next();
			setArchFormat(sline.next());
			if(archFormat.length()>ReadArchive.ARCH_MAXMPMITEMS)
			{	for(int ii=ReadArchive.ARCH_MAXMPMITEMS;ii<archFormat.length();ii++)
				{	if(archFormat.charAt(ii)=='Y')
						throw new Exception("This archive includes data not supported by this version of NairnFEAMPMViz");
				}
			}
			
			// crack archive format
			sline=new Scanner(s.next()).useDelimiter(": ");
			sline.next();
			setCrackFormat(sline.next());
			if(crackFormat.length()>ReadArchive.ARCH_MAXCRACKITEMS)
			{	for(int ii=ReadArchive.ARCH_MAXCRACKITEMS;ii<crackFormat.length();ii++)
				{	if(crackFormat.charAt(ii)=='Y')
						throw new Exception("This archive includes data not supported by this version of NairnFEAMPMViz");
				}
			}
			
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
				Double atime=new Double(words[1]);
				addArchiveFile(atime.doubleValue()*timeScale,words[2]);
			}
			
			// error in no files were found
			if(archiveTimes.size()<1)
				throw new Exception("None of the archived results files could be found.");
		}
			
		//---------------------------------------------------------------
		// FEA specific reslts
		if(isFEAAnalysis())
		{	// FEA nodal displacements
			bcs=section("NODAL DISPLACEMENTS (in mm)");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				int numFound;
				NodalPoint anode;
				while(s.hasNextInt())
				{	numFound=s.nextInt();		// node number
					if(numFound>nnodes || numFound<1)
						throw new Exception("Found nodal displacement with unexpected node number.");
					anode=nodes.get(numFound-1);
					anode.dispx=s.nextDouble()*lengthScale;
					anode.dispy=s.nextDouble()*lengthScale;
					dxmin=Math.min(dxmin,anode.x+anode.dispx);
					dxmax=Math.max(dxmax,anode.x+anode.dispx);
					dymin=Math.min(dymin,anode.y+anode.dispy);
					dymax=Math.max(dymax,anode.y+anode.dispy);
				}
				feaArchFormat[ReadArchive.ARCH_FEADisplacements]='Y';
			}
			
			for(i=0;i<elements.size();i++)
				elements.get(i).setElemPath(this,true);
				
			// FEA nodal stresses
			bcs=section("AVERAGE NODAL STRESSES (in MPa)");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				int numFound;
				NodalPoint anode;
				while(s.hasNextInt())
				{	numFound=s.nextInt();		// node number
					if(numFound>nnodes || numFound<1)
						throw new Exception("Found nodal displacement with unexpected node number.");
					anode=nodes.get(numFound-1);
					anode.sigxx=s.nextDouble();
					anode.sigyy=s.nextDouble();
					anode.sigzz=s.nextDouble();
					anode.sigxy=s.nextDouble();
				}
				feaArchFormat[ReadArchive.ARCH_FEAAvgStress]='Y';
			}
			
			// FEA element energies
			bcs=section("STRAIN ENERGIES IN ELEMENTS (in J)");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
				ElementBase aelem;
				while(s.hasNextInt())
				{	elemNum=s.nextInt();		// element number
					if(elemNum>nelems || elemNum<1)
						throw new Exception("Found element energy with unexpected element number.");
					aelem=elements.get(elemNum-1);
					aelem.energy=s.nextDouble();
				}
				feaArchFormat[ReadArchive.ARCH_FEAElemEnergy]='Y';
			}
			
			// NODAL FORCES AND ELEMENT STRESSES (FEA)
			bcs=section("NODAL FORCES (in N) AND STRESSES (in MPa) IN EACH ELEMENT");
			lineStart=findNextLine(bcs,"-----");
			if(lineStart>0 && lineStart<bcs.length())
			{	s=new Scanner(bcs.substring(lineStart,bcs.length()-1));
			
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
							throw new Exception("Found element results outside element definition.");
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
								throw new Exception("Found element results with unexpected element number.");
							aelem=elements.get(elemNum-1);
						}
						
						// read nodal forces
						for(i=0;i<aelem.getNumberNodes();i++)
						{	s.nextInt();			// skip node number
							aelem.setForces(s.nextDouble(),s.nextDouble(),i);
						}
					}
					
					else if(word1.charAt(0)=='s')
					{	hasStresses=true;
					
						// read element is needed
						if(aelem==null)
						{	elemNum=s.nextInt();
							if(elemNum<1 || elemNum>nelems)
								throw new Exception("Found element results with unexpected element number.");
							aelem=elements.get(elemNum-1);
						}
						
						// read nodal forces
						for(i=0;i<aelem.getNumberNodes();i++)
						{	s.nextInt();			// skip node number
							if(word1.charAt(3)=='x' || word1.charAt(3)=='r')
								aelem.setXYStresses(s.nextDouble(),s.nextDouble(),s.nextDouble(),i);
							else
								aelem.set3DStresses(s.nextDouble(),s.nextDouble(),s.nextDouble(),i);
						}
					}
					
					// skip underlying "----"
					if(s.hasNext()) s.next();
				}
				
				// save what was found (if even on selected ones were found)
				if(hasForces) feaArchFormat[ReadArchive.ARCH_FEAElemForce]='Y';
				if(hasStresses) feaArchFormat[ReadArchive.ARCH_FEAElemStress]='Y';
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
	
	// add nodal point
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
			
			case ElementBase.LINEAR_INTERFACE:
				feaArchFormat[ReadArchive.ARCH_Interfaces]='Y';
				newElem=new LinearInt(elemNum,nds);
				break;
			
			case ElementBase.QUAD_INTERFACE:
				feaArchFormat[ReadArchive.ARCH_Interfaces]='Y';
				newElem=new QuadInt(elemNum,nds);
				break;
				
			default:
				return;
		}
		
		elements.add(newElem);
		if(isFEAAnalysis()) newElem.setFEAProperties(matID,angle,thickness);
		newElem.setElemPath(this,false);
	}
	
	// add grid BC
	public void addGridBC(int nodeNum,int dof,int bcID,double bcVal,double bcArg,double bcAngle)
	{	
		gridBCs.add(new GridDispBC(nodeNum,dof,bcID,bcVal,bcArg,bcAngle));
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
	public void addParticleBC(int partNum,int dof,int bcID,double bcVal,double bcArg)
	{	
		particleBCs.add(new ParticleBC(partNum,dof,bcID,bcVal,bcArg));
	}
	
	// add archive file if it exists
	public void addArchiveFile(double time,String name)
	{   // do not add unless the file exists
	    File archive=new File(path,name);
		if(!archive.exists()) return;
		archives.add(archive);
		System.out.println((new Double(time))+","+name);
	    archiveTimes.add(new Double(time));
	}
	
	//-----------------------------------------------------------------
	// Utilities for decoding text when reading a file
	//-----------------------------------------------------------------
	
	// find index to start of line after line containg substring
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
	
	// return string of section by name of the section
	public String section(String name)
	{   int index;
		for(index=0;index<sectionTitle.size();index++)
		{   if(name.equals(sectionTitle.get(index)))
				return section(index);
		}
		return "";
	}
	
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
	
	public boolean isFEAAnalysis() { return np<BEGIN_MPM_TYPES; }
	public boolean isMPMAnalysis() { return np>=BEGIN_MPM_TYPES; }
	public boolean isAxisymmetric() { return (np==AXI_SYM) || (np==AXI_SYM_MPM) ; }
	public boolean is3D() { return np==THREED_MPM; }
	public double currentTime() { return currentArchive>=0 ? archiveTimes.get(currentArchive) : 0.; }
	
	// set controller
	public void setDocController(DocViewer dc) { docCtrl=dc; }
	public DocViewer getDocController() { return docCtrl; }
	
	//-----------------------------------------------------------------
	// Standard methods to support JTable to display file sections
	//-----------------------------------------------------------------
	public int getRowCount() { return sectionTitle.size();}
	public int getColumnCount() { return 1; }
	public Object getValueAt(int row,int column) { return sectionTitle.get(row); }
	public String getColumnName(int column) { return new String("Sections"); }
	public Class<?> getColumnClass(int column) { return String.class; }
	
}

