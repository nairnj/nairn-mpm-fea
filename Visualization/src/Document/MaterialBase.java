/*******************************************************************
	MaterialBase.java
	NairnFEAMPMViz

	Created by John Nairn on 9/6/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

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
	public static final int RIGIDMATERIAL=11;
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
	public static final int ISOSOFTENING=50;
	public static final int IGNORECONTACT=60;
	public static final int COULOMBFRICTION=61;
	public static final int LINEARIMPERFECT=62;
	public static final int ADHESIVEFRICTION=63;
	public static final int UNKNOWNMATERIAL=64;
	
	public static final int NOTHING=0;
	public static final int ENG_BIOT_PLASTIC_STRAIN=1;
	public static final int LEFT_CAUCHY_ELASTIC_B_STRAIN=2;
	public static final int LEFT_CAUCHY_TOTAL_B_STRAIN=3;
	public static final int MEMBRANE_DEFORMATION=4;

	protected int type;
	protected String name;
	protected double rho;
	
	// initialize
	MaterialBase(String matName,int matType)
	{	name=matName;
		type=matType;
	}
	
	// decode data
	public void decodeData(Scanner s)
	{	// scan to end
		while(s.hasNext())
		{	Scanner sline=new Scanner(s.next());
			sline.useLocale(Locale.US);
			if(!sline.hasNext())
			{	sline.close();
				break;
			}
			if(sline.next().equals("rho="))
				rho=sline.nextDouble();
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
			case RIGIDMATERIAL:
				return NOTHING;
				
			case DUGDALE:
			case ISOPLASTICITY:
			case HILLPLASTIC:
			case JOHNSONCOOK:
			case MGEOSMATERIAL:
			case SLMATERIAL:
			case WOODMATERIAL:
			case ISOSOFTENING:
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
