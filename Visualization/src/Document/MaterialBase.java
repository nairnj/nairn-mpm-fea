/*******************************************************************
	MaterialBase.java
	NairnFEAMPMViz

	Created by John Nairn on 9/6/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.Color;
import java.util.*;

public class MaterialBase
{
	public static final int ISOTROPIC=1;
	public static final int TRANSISO1=2;
	public static final int TRANSISO2=3;
	public static final int ORTHO=4;
	public static final int INTERFACEPARAMS=5;
	public static final int DUGDALE=6;
	public static final int VISCOELASTIC=7;
	public static final int MOONEYRIVLIN=8;
	public static final int ISOPLASTICITY=9;
	public static final int VONMISES=9;
	public static final int BISTABLEISO=10;
	public static final int RIGIDBCMATERIAL=11;
	public static final int COHESIVEZONEMATERIAL=12;
	public static final int LINEARTRACTIONMATERIAL=13;
	public static final int CUBICTRACTIONMATERIAL=14;
	public static final int HILLPLASTIC=15;
	public static final int JOHNSONCOOK=16;
	public static final int MGSCGLMATERIAL=17;
	public static final int MGEOSMATERIAL=17;
	public static final int SLMATERIAL=18;
	public static final int WOODMATERIAL=19;
	public static final int TRILINEARTRACTIONLAW=20;
	public static final int HEANISOTROPIC=21;
	public static final int IDEALGAS=22;
	public static final int COUPLESAWTOOTHTRACTIONLAW=23;
	public static final int HEISOTROPIC=24;
	public static final int HEMGEOSMATERIAL=25;
	public static final int TAITLIQUID=27;
	public static final int NEOHOOKEAN=28;
	public static final int CLAMPEDNEOHOOKEAN=29;
	public static final int PHASETRANSITION=30;
	public static final int RIGIDCONTACTMATERIAL=35;
	public static final int RIGIDBLOCKMATERIAL=36;
	public static final int ISOSOFTENING=50;
	public static final int TRANSISOSOFTENING1=51;
	public static final int TRANSISOSOFTENING2=52;
	public static final int ISOPLASTICSOFTENING=53;
	public static final int IGNORECONTACT=60;
	public static final int COULOMBFRICTION=61;
	public static final int LINEARIMPERFECT=62;
	public static final int ADHESIVEFRICTION=63;
	public static final int LIQUIDCONTACT=64;
	public static final int UNKNOWNMATERIAL=65;
	
	public static final int NOTHING=0;
	public static final int ENG_BIOT_PLASTIC_STRAIN=1;
	public static final int LEFT_CAUCHY_ELASTIC_B_STRAIN=2;
	public static final int LEFT_CAUCHY_TOTAL_B_STRAIN=3;
	public static final int MEMBRANE_DEFORMATION=4;

	protected int type;
	protected String name;
	protected double rho;
	protected Color matClr;
	
	// initialize
	MaterialBase(String matName,int matType)
	{	name=matName;
		type=matType;
		matClr = null;
	}
	
	// decode data
	public void decodeData(Scanner s)
	{	// scan to end
		while(s.hasNext())
		{	String matLine = s.next();
			Scanner sline=new Scanner(matLine);
			sline.useLocale(Locale.US);
			if(!sline.hasNext())
			{	sline.close();
				break;
			}
			String prop = sline.next();
			if(prop.equals("rho="))
				rho=sline.nextDouble();
			else if(prop.equals("rho"))
			{	if(sline.hasNext())
				{	prop = sline.next();
					if(sline.hasNext())
						rho=sline.nextDouble();
				}
			}
			else if(prop.equals("color="))
			{	String [] cstr = matLine.substring(7).split(", ");
				float red=0.f,green=0.f,blue=0.f,alpha=1.f;
				if(cstr.length>0)
					red = Float.parseFloat(cstr[0]);
				if(cstr.length>1)
					green = Float.parseFloat(cstr[1]);
				if(cstr.length>2)
					blue = Float.parseFloat(cstr[2]);
				if(cstr.length>3)
					alpha = Float.parseFloat(cstr[3]);
				matClr = new Color(red,green,blue,alpha);
			}
			sline.close();
		}
	}
	
	// most are not rigid
	public boolean isRigid() { return false; }
	
	// true it this material archives plastic strain that
	// should be added to elastic strain when finding deformation
	// gradient. It is onlu low-strain plasticity materials
	public boolean hasPlasticStrainForGradient(ResultsDocument doc)
	{
		// new style always has full deformation gradietn
		if(doc.units.getNewUnitsVersion()) return false;
		
		// legacy files have it for some materials
		switch(type)
		{	case DUGDALE:
			case VONMISES:
			case HILLPLASTIC:
			case JOHNSONCOOK:
			case MGSCGLMATERIAL:
			case SLMATERIAL:
			case WOODMATERIAL:
				return true;
			default:
				break;
		}
		return false;
	}

	// return contents of plastic strain for this material
	public int AltStrainContains()
	{
		switch(type)
		{	case ISOTROPIC:
			case TRANSISO1:
			case TRANSISO2:
			case INTERFACEPARAMS:
			case VISCOELASTIC:
			case BISTABLEISO:
			case RIGIDBCMATERIAL:
			case RIGIDCONTACTMATERIAL:
			case RIGIDBLOCKMATERIAL:
				return NOTHING;
				
			case DUGDALE:
			case ISOPLASTICITY:
			case HILLPLASTIC:
			case JOHNSONCOOK:
			case MGEOSMATERIAL:
			case SLMATERIAL:
			case WOODMATERIAL:
			case ISOSOFTENING:
			case TRANSISOSOFTENING1:
			case TRANSISOSOFTENING2:
			case ISOPLASTICSOFTENING:
				return ENG_BIOT_PLASTIC_STRAIN;
			
			case MOONEYRIVLIN:
			case IDEALGAS:
			case TAITLIQUID:
			case NEOHOOKEAN:
			case CLAMPEDNEOHOOKEAN:
				return LEFT_CAUCHY_TOTAL_B_STRAIN;
			
			case HEISOTROPIC:
			case HEMGEOSMATERIAL:
			case PHASETRANSITION:
				return LEFT_CAUCHY_ELASTIC_B_STRAIN;
			
			case HEANISOTROPIC:
				return MEMBRANE_DEFORMATION;
				
			default:
				break;
		}
		
		// unknown materials default to plastic strain
		// would be better to look up in an editable resource
		return ENG_BIOT_PLASTIC_STRAIN;
		
	}

}
