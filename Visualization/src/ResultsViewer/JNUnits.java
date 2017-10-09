/*
 * JNUNits
 * NairnFEAMPMViz
 * 
 * Created 4/22/2015, John A. Nairn
 * Copyright 2015, RSAC Software, All Rights Reserved
 * 
 * Notes
 * 1. Default viz units are
 *		Length = mm, mass = g, time = ms
 *		Velocity = m/s (mm/ms), force = N, stress = MPa
 *		Energy = J
 * 2. The output...Scale constants convert the output results of the
 *		the current document to the default viz units
 * 3. The ...Scale constants convert the output units of the current
 *		document to the current viz units (i.e., they are output...Scale
 *		times conversion from default viz units to current viz units)
 * 4. To convert current data to Legacy use output...Scale/...Scale
 * 
 */

public class JNUnits
{
	private boolean newUnitsVersion;
	
	private String lengthUnits;
	private int lengthScaleIndex;
	private double lengthScale;
	private double ccLengthScale;
	private double outputLengthScale;
	
	private String massUnits;
	private int massScaleIndex;
	private double massScale;
	private double outputMassScale;
	
	private String timeUnits;
	private String globalTimeUnits;
	private int timeScaleIndex;
	private double timeScale;
	private double altTimeScale;
	private double outputTimeScale;
	private double outputAltTimeScale;
	
	private String velocityUnits;
	private int velocityScaleIndex;
	private double velocityScale;
	private double outputVelocityScale;
	
	private String forceUnits;
	private int forceScaleIndex;
	private double forceScale;
	private double ccForceScale;
	private double outputForceScale;
	
	private String stressUnits;
	private int stressScaleIndex;
	private double feaStressScale;
	private double mpmStressScale;
	private double outputFeaStressScale;
	private double outputMpmStressScale;
	
	private String strainUnits;
	private int strainScaleIndex;
	private double strainScale;
	private double outputStrainScale;
	
	private String energyUnits;
	private int energyScaleIndex;
	double energyScale;
	private double outputEnergyScale;
	
	// initialize
	public JNUnits()
	{	lengthScaleIndex = 2;			// mm
		massScaleIndex = 2;				// g
		timeScaleIndex = 1;				// ms
		velocityScaleIndex = 2;			// m/s = mm/msec
		forceScaleIndex = 3;			// N
		stressScaleIndex = 1;			// MPa
		energyScaleIndex = 3;			// J
		strainScaleIndex = 1;			// %
		setOutputUnits("Legacy");
	}

	// set scaling to convert output results in this document to the
	// default visualization units
	public void setOutputUnits(String style)
	{
		newUnitsVersion = true;
		String[] parts = style.split("-");
		
		if(parts.length<3)
		{	// Legacy units
			outputLengthScale = 1.;				// mm to mm
			outputTimeScale = 1000.;			// s to ms
			outputAltTimeScale = 1.;			// ms to s
			outputMassScale = 1.;				// g to g
			outputVelocityScale = 0.001;		// mm/s to m/s
			outputForceScale = 1.;				// N to N
			outputFeaStressScale = 1.;			// MPa to MPa
			outputMpmStressScale = 1.e-6;		// Pa to MPa
			outputEnergyScale = 1.;				// J to J
			
			// check for old version
			if(style.equals("_none_"))
				newUnitsVersion = false;
		}
		
		else
		{	// length (to m)
			int lexp = 0;
			if(parts[0].equals("L")) parts[0] = "m";
			if(parts[0].length()>1) lexp = readPrefix(parts[0].charAt(0));
			
			// mass (to g)
			int mexp = 0;
			if(parts[1].equals("M")) parts[1] = "kg";
			if(parts[1].length()>1) mexp = readPrefix(parts[1].charAt(0));

			// time (to s)
			int texp = 0;
			if(parts[2].equals("L")) parts[2] = "s";
			if(parts[2].length()>1) texp = readPrefix(parts[2].charAt(0));
			
			// get scale to convert file to defauilt viz units
			if(lexp<100)
				outputLengthScale = 1000.*Math.pow(10.,(double)lexp);		// length to mm
			else
				outputLengthScale = 1.;
			
			if(mexp<100)
				outputMassScale = Math.pow(10.,(double)mexp);				// mass to g
			else
				outputMassScale = 1.;
			
			if(texp<100)
			{	outputTimeScale = 1000.*Math.pow(10.,(double)texp);			// time to ms
				outputAltTimeScale = outputTimeScale;
			}
			else
			{	outputTimeScale = 1.;
				outputAltTimeScale = 1.;
			}
			
			// velocity L/T to m/s
			if(lexp<100 && texp<100)
				outputVelocityScale = Math.pow(10.,(double)(lexp-texp));
			else
				outputVelocityScale = 1.;
			
			// Force M-L/T^2 to N, and Stress M/(L-T^2) to MPa, and energy F-L to J
			if(lexp<100 && texp<100 && mexp<100)
			{	outputForceScale = Math.pow(10.,(double)(mexp-3+lexp-2*texp));
				outputFeaStressScale = Math.pow(10.,(double)(mexp-3-lexp-2*texp-6));
				outputMpmStressScale = outputFeaStressScale;
				outputEnergyScale = Math.pow(10.,(double)(mexp-3+2*lexp-2*texp));
			}
			else
			{	outputForceScale = 1.;
				outputFeaStressScale = 1.;
				outputMpmStressScale = 1.;
				outputEnergyScale = 1.;
			}
			
		}
		
		// strain absolute in files and Legacy viz in %
		outputStrainScale = 100.;
		
		// adjust viz units
		setLengthScaleIndex(lengthScaleIndex);
		setMassScaleIndex(massScaleIndex);
		setTimeScaleIndex(timeScaleIndex);
		setVelocityScaleIndex(velocityScaleIndex);
		setForceScaleIndex(forceScaleIndex);
		setStressScaleIndex(stressScaleIndex);
		setEnergyScaleIndex(energyScaleIndex);
		setStrainScaleIndex(strainScaleIndex);
		
		// get global time
		if(outputAltTimeScale>999.)
			setGlobalTimeUnits("s");
		else if(outputAltTimeScale<.1)
			setGlobalTimeUnits("\u03BCs");			// microseconds
		else
			setGlobalTimeUnits("ms");
	}
	
	// return flag to now version with units or not
	public boolean getNewUnitsVersion() { return newUnitsVersion; }

	// read prefixes for length, mass, and time
	public int readPrefix(char prefix)
	{	switch(prefix)
		{	case 'k':
				return 3;
			case 'd':
				return -1;
			case 'c':
				return -2;
			case 'm':
				return -3;
			case 'u':
				return -6;
			case 'n':
				return -9;
			default:
				break;
		}
		return 100;
	}

	// set new scale, but return YES or NO if changed
	// ccLengthScale is from viz scale to mm
	public boolean setLengthScaleIndex( int newScale)
	{	boolean change = (newScale != lengthScaleIndex);
		lengthScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setLengthUnits("nm");
				lengthScale = outputLengthScale*1.e6;
				ccLengthScale = 1.e-6;
				break;
			case 1:
				setLengthUnits("\u03BCm");					// micro meters
				lengthScale = outputLengthScale*1.e3;
				ccLengthScale = 1.e-3;
				break;
			case 2:
				setLengthUnits("mm");
				lengthScale = outputLengthScale;
				ccLengthScale = 1.;
				break;
			case 3:
				setLengthUnits("cm");
				lengthScale = outputLengthScale*0.1;
				ccLengthScale = 10.;
				break;
			case 4:
				setLengthUnits("m");
				lengthScale = outputLengthScale*0.001;
				ccLengthScale = 1000.;
				break;
			case 5:
				setLengthUnits("mil");
				lengthScale = outputLengthScale*1000./25.4;
				ccLengthScale = 25.4/1000.;
				break;
			case 6:
				setLengthUnits("in");
				lengthScale = outputLengthScale/25.4;
				ccLengthScale = 25.4;
				break;
			default:
				setLengthUnits("ft");
				lengthScale = outputLengthScale/(12.*25.4);
				ccLengthScale = 112.*25.4;
				break;
		}
		return change;
	}
	
	public int lengthScaleIndex() { return lengthScaleIndex; }

	// mass units
	public boolean setMassScaleIndex(int newScale)
	{	boolean change = (newScale != massScaleIndex);
		massScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setMassUnits("\u03BCg");						// micrograms
				massScale = outputMassScale*1.e6;
				break;
			case 1:
				setMassUnits("mg");
				massScale = outputMassScale*1000.;
				break;
			case 2:
				setMassUnits("g");
				massScale = outputMassScale;
				break;
			case 3:
				setMassUnits("kg");
				massScale = outputMassScale*0.001;
				break;
			case 4:
				setMassUnits("oz");
				massScale = outputMassScale*0.001*0.028349;
				break;
			default:
				setMassUnits("lb");
				massScale = outputMassScale*0.001*0.453592376254;
				break;
		}
		return change;
	}
	
	public int massScaleIndex() { return massScaleIndex; }

	// set time scale nd return YES or NO if really changed
	public boolean setTimeScaleIndex(int newScale)
	{	boolean change = (newScale != timeScaleIndex);
		timeScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setTimeUnits("\u03BCs");						// microseconds
				timeScale = outputTimeScale*1000.;
				altTimeScale = outputAltTimeScale*1000.;
				break;
			case 1:
				setTimeUnits("ms");
				timeScale = outputTimeScale;
				altTimeScale = outputAltTimeScale;
				break;
			default:
				setTimeUnits("s");
				timeScale = outputTimeScale*0.001;
				altTimeScale = outputAltTimeScale*0.001;
				break;
		}
		return change;
	}
	
	public int timeScaleIndex() { return timeScaleIndex; }

	// velocity units
	public boolean setVelocityScaleIndex(int newScale)
	{	boolean change = (newScale != velocityScaleIndex);
		velocityScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setVelocityUnits("mm/s");
				velocityScale = outputVelocityScale*1000.;
				break;
			case 1:
				setVelocityUnits("cm/s");
				velocityScale = outputVelocityScale*100.;
				break;
			case 2:
				setVelocityUnits("m/s");			// also mm/msec
				velocityScale = outputVelocityScale;
				break;
			case 3:
				setVelocityUnits("km/s");
				velocityScale = outputVelocityScale*0.001;
				break;
			case 4:
				setVelocityUnits("in/s");
				velocityScale = outputVelocityScale*1000/25.4;
				break;
			case 5:
			default:
				setVelocityUnits("ft/s");
				velocityScale = outputVelocityScale*1000/(12.*25.4);
				break;
		}
		return change;
	}
	
	public int velocityScaleIndex() { return velocityScaleIndex; }

	// force units
	public boolean setForceScaleIndex(int newScale)
	{	boolean change = (newScale != forceScaleIndex);
		forceScaleIndex = newScale;
		// ccForceScale is viz force to N
		switch(newScale)
		{	case 0:
				setForceUnits("nN");
				forceScale = outputForceScale*1.e9;
				ccForceScale = 1.e-9;
				break;
			case 1:
				setForceUnits("\u03BCN");				// micro newtons
				forceScale = outputForceScale*1.e6;
				ccForceScale = 1.e-6;
				break;
			case 2:
				setForceUnits("mN");
				forceScale = outputForceScale*1000.;
				ccForceScale = 1.e-3;
				break;
			case 3:
				setForceUnits("N");
				forceScale = outputForceScale;
				ccForceScale = 1.;
				break;
			case 4:
				setForceUnits("kN");
				forceScale = outputForceScale*0.001;
				ccForceScale = 1000.;
				break;
			case 5:
				setForceUnits("dyne");
				forceScale = outputForceScale*1.e5;
				ccForceScale = 1.e-5;
				break;
			case 6:
				setForceUnits("kgf");
				forceScale = outputForceScale*0.101972;
				ccForceScale = 1./0.101972;
				break;
			case 7:
			default:
				setForceUnits("lb");
				forceScale = outputForceScale*0.224809;
				ccForceScale = 1./0.224809;
				break;
		}
		return change;
	}
	
	public int forceScaleIndex() { return forceScaleIndex; }

	// stress units
	public boolean setStressScaleIndex(int newScale)
	{	boolean change = (newScale != stressScaleIndex);
		stressScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setStressUnits("Pa");
				feaStressScale = outputFeaStressScale*1.e6;
				mpmStressScale = outputMpmStressScale*1.e6;
				break;
			case 1:
				setStressUnits("MPa");
				feaStressScale = outputFeaStressScale;
				mpmStressScale = outputMpmStressScale;
				break;
			case 2:
				setStressUnits("GPa");
				feaStressScale = outputFeaStressScale*1.e-3;
				mpmStressScale = outputMpmStressScale*1.e-3;
				break;
			case 3:
				setStressUnits("Ba");
				feaStressScale = outputFeaStressScale*1.e7;
				mpmStressScale = outputMpmStressScale*1.e7;
				break;
			case 4:
				setStressUnits("psi");
				feaStressScale = outputFeaStressScale*145.0377;
				mpmStressScale = outputMpmStressScale*145.0377;
				break;
			case 5:
			default:
				setStressUnits("kpsi");
				feaStressScale = outputFeaStressScale*.1450377;
				mpmStressScale = outputMpmStressScale*.1450377;
				break;
		}
		return change;
	}
	
	public int stressScaleIndex() { return stressScaleIndex; }

	// energy units
	public boolean setEnergyScaleIndex(int newScale)
	{	boolean change = (newScale != energyScaleIndex);
		energyScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setEnergyUnits("nJ");
				energyScale = outputEnergyScale*1.e9;
				break;
			case 1:
				setEnergyUnits("\u03BCJ");					// micro joules
				energyScale = outputEnergyScale*1.e6;
				break;
			case 2:
				setEnergyUnits("mJ");
				energyScale = outputEnergyScale*1.e3;
				break;
			case 3:
				setEnergyUnits("J");
				energyScale = outputEnergyScale;
				break;
			case 4:
				setEnergyUnits("kJ");
				energyScale = outputEnergyScale*1.e-3;
				break;
			case 5:
				setEnergyUnits("ergs");
				energyScale = outputEnergyScale*1.e7;
				break;
			case 6:
			default:
				setEnergyUnits("ft-lbs");
				energyScale = outputEnergyScale*0.737563;
				break;
		}
		return change;
	}
	
	public int energyScaleIndex() { return energyScaleIndex; }

	// strain units
	public boolean setStrainScaleIndex(int newScale)
	{	boolean change = (newScale != strainScaleIndex);
		strainScaleIndex = newScale;
		switch(newScale)
		{	case 0:
				setStrainUnits("");
				strainScale = outputStrainScale*.01;
				break;
			case 1:
				setStrainUnits("%");
				strainScale = outputStrainScale;
				break;
			case 2:
			default:
				setStrainUnits("\u03BCstr");				// micro strain
				strainScale = outputStrainScale*1.e4;
				break;
		}
		return change;
	}
	public int strainScaleIndex() { return strainScaleIndex; }

	// find scaling of current velocity mass and velocity in viz units
	// and ten to current energy viz units
	public double calcVelocityScale()
	{	// convert mass to kg and velocity to m/sec
		double mkg = 0.001*outputMassScale/massScale;
		double velmpers = outputVelocityScale/velocityScale;
		
		// convert to Legacy J and then to current energy units
		return mkg*velmpers*velmpers*energyScale/outputEnergyScale;
	}

	public String lengthUnits() { return lengthUnits; }
	public void setLengthUnits(String argStr) { lengthUnits = new String(argStr); }
	public double lengthScale() { return lengthScale; }
	public double ccLengthScale() { return ccLengthScale; }

	public String massUnits() { return massUnits; }
	public void setMassUnits(String argStr) { massUnits = new String(argStr); }
	public double massScale() { return massScale; }
	public void setMassScale(double arg) { massScale = arg; }

	public String timeUnits() { return timeUnits; }
	public void setTimeUnits(String argStr) { timeUnits = new String(argStr); }
	public String globalTimeUnits() { return globalTimeUnits; }
	public void setGlobalTimeUnits(String argStr) { globalTimeUnits = new String(argStr); }
	public double timeScale() { return timeScale; }
	public double altTimeScale() { return altTimeScale; }

	public String velocityUnits() { return velocityUnits; }
	public void setVelocityUnits(String argStr) { velocityUnits = new String(argStr); }
	public double velocityScale() { return velocityScale; }

	public String forceUnits() { return forceUnits; }
	public void setForceUnits(String argStr) { forceUnits = new String(argStr); }
	public double forceScale() { return forceScale; }
	public double ccForceScale() { return ccForceScale; }

	public String stressUnits() { return stressUnits; }
	public void setStressUnits(String argStr) { stressUnits = new String(argStr); }
	public double feaStressScale() { return feaStressScale; }
	public double mpmStressScale() { return mpmStressScale; }

	public String strainUnits() { return strainUnits; }
	public void setStrainUnits(String argStr) { strainUnits = new String(argStr); }
	public double strainScale() { return strainScale; }

	public String energyUnits() { return energyUnits; }
	public void setEnergyUnits(String argStr) { energyUnits = new String(argStr); }
	public double energyScale() { return energyScale; }

}
