/********************************************************************************
	FailureSurface.cpp
	nairn-mpm-fea

	Created by John Nairn, June 26, 2015.
	Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/FailureSurface.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"

#include "NairnMPM_Class/NairnMPM.hpp"

#define SQRTTWOMINUSONE 0.4142135623730950

#pragma mark FailureSurface::Constructors and Destructors

// Constructors
FailureSurface::FailureSurface() {}

// Main Constructor
FailureSurface::FailureSurface(MaterialBase *pair)
{
	criticalNormal = 1.e40;
	criticalShear = 1.e40;
	parent = pair;
	failureSurface = 0;
    InitPDTerms();
}

#pragma mark FailureSurface::Initialize

// Read hardening law properties
char *FailureSurface::InputInitationProperty(char *xName,int &input,double &gScaling)
{
    // principle stress failure loads
    if(strcmp(xName,"sigmac")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalNormal,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"tauc")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalShear,gScaling,1.e6);
    }

    // rest are to handle pressure dependence
    return InputPDProperty(xName,input,gScaling);
}

// verify settings and some initial calculations
const char *FailureSurface::VerifyAndLoadProperties(int np)
{
	// reduced failure stresses
	double rho = parent->GetRho(NULL);
    sigmacRed = criticalNormal/rho;
    taucRed = criticalShear/rho;
	
    // pressure dependent terms
    return VerifyAndLoadPDTerms(rho);
}

// print just initiation properties to output window
void FailureSurface::PrintInitiationProperties(void) const
{
    cout << GetInitiationLawName() << endl;
    MaterialBase::PrintProperty("sigc",criticalNormal*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("tauc",criticalShear*UnitsController::Scaling(1.e-6),"");
    PrintPDTerms();
}

#pragma mark FailureSurface::Methods

// check stress state in initial axis system to see if failed
// If failed, find the normal vector in the same material axis system, but for
//     3D, find Euler angles in normal, change relStrength to
//          (relStrength)/(failure_stress)
//     for scaling when using smoothed stresses, and return failure mode
// If not failed, return 0
int FailureSurface::ShouldInitiateFailure(Tensor *str,Vector *normal,int np,double &relStrength) const
{
	// initialize
	int failureMode = NO_FAILURE;

// To eliminiate some 3D transition zones, uncomment
//#define NO3DTRANSITIONS
	
	// get principal stress
	if(np==THREED_MPM)
	{
		// get principle stress
		Matrix3 strmat = TensorToMatrix(str,true);
		Vector s = strmat.Eigenvalues();
		
		// Six options - convert to x>y>z
		if(s.y<s.z)
		{	// swap y and z
			double temp = s.y;
			s.y = s.z;
			s.z = temp;
		}
		// Now have x>y>z, y>x>z, or y>z>x
		if(s.x<s.z)
		{	// rotate y->x,z->y,x->z
			Vector t = s;
			s.x = t.y;
			s.y = t.z;
			s.z = t.x;
		}
		else if(s.x<s.y)
		{	// swap x and y
			double temp = s.y;
			s.y = s.x;
			s.x = temp;
		}
		
		// get pressure, and find region
		double P = -(s.x+s.y+s.z)/3.;
		double sigi = relStrength*sigmacRed;
		double taui = relStrength*taucRed;
#ifndef NO3DTRANSITIONS
		double nx = 1.,ny = 0.,nz = 0.;
#endif
		double dt=0.,ds=0.;

        // pressure dependence
        taui = GetPDStrength(taui,P,relStrength);

		if(failureSurface==OVOID_SURFACE && sigi>taui)
		{	// Handle this surface using Mohr's stress plot
			// always possible when sigi>taui
 			
			// Mohr's circle radius and center
			double r = 0.5*(s.x-s.z);
			double p = 0.5*(s.x+s.z);
			
			// Three cases depending on r and ellipse properties
			double tipRcrit = taui*taui/sigi;
			if(r<=tipRcrit)
			{	// Small radius curvature fails in tension when p>sigi-r or r>sigi-p
				if(r>sigi-p)
				{	// fails in tension, watch for multiple tension
					if(s.y>sigi)
					{	if(s.z<sigi)
						{	// only sig.y is greater than sigi}
							nx = s.x-sigi;
							ny = s.y-sigi;
							failureMode = TENSILE_TENSILE_FAILURE;
						}
						else
						{	// s.y and s.z are greater than sigi
							nx = s.x-sigi;
							ny = s.y-sigi;
							nz = s.z-sigi;
                            failureMode = TENSILE_TENSILE_FAILURE;
						}
					}
					else
					{	failureMode = TENSILE_FAILURE;
					}
				}
			}
			
			else if(r>=taui)
			{	// shear failure
				nz = -1.;
				failureMode = SHEAR_FAILURE;
			}
			
			else
			{	// ovoid failure
				double c = sqrt(sigi*sigi-taui*taui);
				double arg = sqrt(taui*taui-r*r);
				double pc = c*arg/taui;
				if(p>pc)
				{	double theta = 0.5*acos(taui*arg/(c*r));
					nx = cos(theta);
					nz = -sin(theta);
					failureMode = TENSILE_SHEAR_FAILURE;
				}
			}
			
		}
        
		else
		{	// Analyze projection on the plane normal to hydrostatic line
			// For details see JANOSU-019-14
            // This code only works for P_INDEPENDENT or P_LINEAR
			// special case for hydrostatic
			int region = 0;
			if(s.x+P==0.)
			{	// if s.x==P, then all must equal P
				if(s.x>sigi)
				{
#ifndef NO3DTRANSITIONS
					failureMode = TENSILE_TENSILE_FAILURE;
					ny = 1.;
					nz = 1.;
#else
					failureMode = TENSILE_FAILURE;
#endif
				}
			}
			else
			{	// pick a region
				double range = 1.5*(sigi+P);
				if(range>=2.*taui)
					region = 1.;
				else if(range>=taui)
					region = 2;
				else if(range>=0.)
					region = 3;
				else
					region = 4;
				
				// key points
				dt = (sigi+P)/(s.x+P);
				ds = 2.*taui/(s.x-s.z);
			}

			// start with regions 2, 3, 4
			if(region==2)
			{	if(ds<1. || dt<1.)
				{	if(ds<dt)
					{	// Region I failure, calculations done below
						region = 1;
					}
					else if((s.z<sigi-2.*taui) && CloserToShear(&s,sigi,P,taui))
					{	// in gap between tension and shear - but closer to shear
						region = 1;
					}
					else
					{	// Region II tensile failure or in gap but closer to tensile
						failureMode = TENSILE_FAILURE;
					}
				}
			}
			
#ifndef NO3DTRANSITIONS
			else if(region==3)
			{	if(dt<1)
				{	// failure, check for mode
					double dt2 = s.y+P > 0. ? (sigi+P)/(s.y+P) : 2.;
					if(dt2>=1.)
					{	// Region III tensile failure
						failureMode = TENSILE_FAILURE;
					}
					else if((s.z<sigi-2.*taui) && CloserToShear(&s,sigi,P,taui))
					{	// in gap between tension and shear - but closer to shear
						region = 1;
					}
					else
					{	// Region III tensile-tensile failure or in gap but closer to tensile
						// changed
						if(s.y>sigi)
						{	failureMode = TENSILE_TENSILE_FAILURE;
							nx = s.x-sigi;
							ny = s.y-sigi;
						}
						else
							failureMode = TENSILE_FAILURE;
					}
				}
			}
			
			else if(region==4)
			{	// region IV is always failed, but how?
				if(s.z>sigi)
				{	// then all are above sigi
					// Region IV tensile-tensile-tensile failure
					failureMode = TENSILE_TENSILE_FAILURE;
					// changed
					nx = s.x-sigi;
					ny = s.y-sigi;
					nz = s.z-sigi;
				}
				else if(s.y>sigi)
				{	// changed
					// Region IV tensile-tensile failure
					failureMode = TENSILE_TENSILE_FAILURE;
					if(s.y>sigi)
					{	nx = s.x-sigi;
						ny = s.y-sigi;
					}
				}
				else
				{	// Region IV tensile failure
					failureMode = TENSILE_FAILURE;
				}
			}
			
			// final region (or above converted to it)
			if(region==1)
			{	if(ds<1.)
				{	// check for shear-shear gap
					double edge = 1.5*(s.y+P);
					if(edge>taui)
					{	// Region I shear-shear gap on right edge
						failureMode = SHEAR_SHEAR_FAILURE;
						nx = s.x+P-2.*taui/3.;
						if(nx>0.)
						{	ny = s.y+P-2.*taui/3.;
							nz = s.z+P+4.*taui/3.;
						}
						else
						{	nx = 1.;
							nz = -1.;
						}
					}
					else if(edge<-taui)
					{	// Region I shear-shear gap on left edge
						failureMode = SHEAR_SHEAR_FAILURE;
						nx = s.x+P-4.*taui/3.;
						if(nx>0.)
						{	ny = s.y+P+2.*taui/3.;
							nz = s.z+P+2.*taui/3.;
						}
						else
						{	nx = 1.;
							nz = -1.;
						}
					}
					else
					{	// Region I shear failure
						failureMode = SHEAR_FAILURE;
						nz = -1.;
					}
				}
			}
#else
			else if(region==3)
			{	if(dt<1)
				{	// failure, check for mode
					double dt2 = (sigi+P)/(s.y+P);
					if(fabs(dt2)>=1.)
					{	// Region III tensile failure
						failureMode = TENSILE_FAILURE;
					}
					else if((s.z<sigi-2.*taus) && CloserToShear(&s,sigi,P,taus))
					{	// in gap between tension and shear - but closer to shear
						region = 1;
					}
					else
					{	// Region III tensile failure
						failureMode = TENSILE_FAILURE;
					}
				}
			}
			
			else if(region==4)
			{	// region IV is always failed in tension
				failureMode = TENSILE_FAILURE;
			}
			
			// final region (or above converted to it)
			if(region==1)
			{	if(ds<1.)
				{	// Region I shear failure
					failureMode = SHEAR_FAILURE;
				}
			}
#endif
		}
		
		// If failed, get rotation terms
		if(failureMode!=NO_FAILURE)
		{	// Get eigenvectors (or rotation matrix) to principle stress
			Matrix3 R = strmat.Eigenvectors(s);
			
			// Note the matrix class might reorder eigenvalues
			// Six options - convert to x>y>z
			if(s.y<s.z)
			{	// swap y and z
				R.SwapColumns(1,2);
				double temp = s.y;
				s.y = s.z;
				s.z = temp;
			}
			// Now have x>y>z, y>x>z, or y>z>x
			if(s.x<s.z)
			{	// rotate y->x,z->y,x->z
				R.SwapColumns(1,2);
				R.SwapColumns(0,2);
				Vector t = s;
				s.x = t.y;
				s.y = t.z;
				s.z = t.x;
			}
			else if(s.x<s.y)
			{	// have y>x>z, swap x and y
				R.SwapColumns(0,1);
				double temp = s.y;
				s.y = s.x;
				s.x = temp;
			}

			// flip z if determinant is < 0
			if(R.determinant()<0.)
			{	R(0,2) = -R(0,2);
				R(1,2) = -R(1,2);
				R(2,2) = -R(2,2);
			}
			
			// if not tensile failure, most rotate unit x in principle plane
			// to the calculated normal vector
			// For details see https://en.wikipedia.org/wiki/Rotation_matrix
			if(failureMode!=TENSILE_FAILURE)
			{
#ifndef NO3DTRANSITIONS
				// add rotation of normal to x axis
				if(ny==0.)
				{	// rotate ccw by acos(nx/|n|) about principle y axis
					// This rotates princple x axis to the normal
					double mag = sqrt(nx*nx+nz*nz);
					double c45 = nx/mag;
					double s45 = -nz/mag;
					double ux = R(0,1);
					double uy = R(1,1);
					double uz = R(2,1);
					Matrix3 R45 = Matrix3(c45+ux*ux*(1.-c45),ux*uy*(1.-c45)-uz*s45,ux*uz*(1.-c45)+uy*s45,
										  uy*ux*(1.-c45)+uz*s45,c45+uy*uy*(1.-c45),uy*uz*(1.-c45)-ux*s45,
										  uz*ux*(1.-c45)-uy*s45,uz*uy*(1.-c45)+ux*s45,c45+uz*uz*(1.-c45));
					R = R45*R;
				}
				else if(nz==0.)
				{	// rotate ccw by acos(nx/|n|) about principle z axis
					// This rotates princple x axis to the normal
					double mag = sqrt(nx*nx+ny*ny);
					double c45 = nx/mag;
					double s45 = ny/mag;
					double ux = R(0,2);
					double uy = R(1,2);
					double uz = R(2,2);
					Matrix3 R45 = Matrix3(c45+ux*ux*(1.-c45),ux*uy*(1.-c45)-uz*s45,ux*uz*(1.-c45)+uy*s45,
										  uy*ux*(1.-c45)+uz*s45,c45+uy*uy*(1.-c45),uy*uz*(1.-c45)-ux*s45,
										  uz*ux*(1.-c45)-uy*s45,uz*uy*(1.-c45)+ux*s45,c45+uz*uz*(1.-c45));
					R = R45*R;
				}
				else
				{	// rotate ccw by acos(nx/|n|) about the principle y axis
					double mag = sqrt(nx*nx+ny*ny+nz*nz);
					double ctheta = nx/mag;
					double arg = sqrt(1-ctheta*ctheta);
					double stheta = -sqrt(1.-ctheta*ctheta);
					double ux = R(0,1);
					double uy = R(1,1);
					double uz = R(2,1);
					Matrix3 Rtheta = Matrix3(ctheta+ux*ux*(1.-ctheta),ux*uy*(1.-ctheta)-uz*stheta,ux*uz*(1.-ctheta)+uy*stheta,
										  uy*ux*(1.-ctheta)+uz*stheta,ctheta+uy*uy*(1.-ctheta),uy*uz*(1.-ctheta)-ux*stheta,
										  uz*ux*(1.-ctheta)-uy*stheta,uz*uy*(1.-ctheta)+ux*stheta,ctheta+uz*uz*(1.-ctheta));
					R = Rtheta*R;

					// rotate ccw by dihedral angle atan(-ny/nz) about principle x axis
					arg *= mag;
					ctheta = nz/arg;
					stheta = -ny/arg;
					ux = R(0,0);
					uy = R(1,0);
					uz = R(2,0);
					Rtheta = Matrix3(ctheta+ux*ux*(1.-ctheta),ux*uy*(1.-ctheta)-uz*stheta,ux*uz*(1.-ctheta)+uy*stheta,
											 uy*ux*(1.-ctheta)+uz*stheta,ctheta+uy*uy*(1.-ctheta),uy*uz*(1.-ctheta)-ux*stheta,
											 uz*ux*(1.-ctheta)-uy*stheta,uz*uy*(1.-ctheta)+ux*stheta,ctheta+uz*uz*(1.-ctheta));
					R = Rtheta*R;
				}
#else
				// rotate ccw by 45 degrees about principle y axis
				// This rotates princple x axis to the normal
				double c45 = 0.5*sqrt(2.);				// cos(45deg)
				double s45 = c45;						// sin(45deg)
				double ux = R(0,1);
				double uy = R(1,1);
				double uz = R(2,1);
				Matrix3 R45 = Matrix3(c45+ux*ux*(1.-c45),ux*uy*(1.-c45)-uz*s45,ux*uz*(1.-c45)+uy*s45,
									  uy*ux*(1.-c45)+uz*s45,c45+uy*uy*(1.-c45),uy*uz*(1.-c45)-ux*s45,
									  uz*ux*(1.-c45)-uy*s45,uz*uy*(1.-c45)+ux*s45,c45+uz*uz*(1.-c45));
				R = R45*R;
#endif
			}
			
			// decode to Euler angles
			// see https://en.wikipedia.org/wiki/Euler_angles - but it has the equation wrong
			// There are two solutions (where second angle is one way or the other) for ZYZ Euler angles
			//   and either can be used for softening material to reconstruct the ZYZ R matrix
			// normal to crack will be first column on reconstructed R matrix
			if(R(2,2)<-1.)
				normal->y = acos(-1.);
			else if(R(2,2)>1.)
				normal->y = acos(1.);
			else
				normal->y = acos(R(2,2));				// beta  = arccos(Z3)
			if(R(1,2)==0. && R(0,2)==0.)
				normal->x = 0.;
			else
				normal->x = atan2(R(1,2),R(0,2));		// alpha = atan2(Z2,Z1)
			if(R(2,1)==0. && R(2,0)==0.)
				normal->z = 0.;
			else
				normal->z = atan2(R(2,1),-R(2,0));		// gamma = atan2(Y3,-X3)
#ifdef CHECK_NAN
			if(normal->x!=normal->x || normal->y!=normal->y || normal->z!=normal->z)
			{
#pragma omp critical (output)
				{	cout << "\n# FailureSurface::ShouldInitiateFailure: bad 3D normals" << endl;
					cout << "# Normals = (" << normal->x << "," << normal->y << "," << normal->z << ")" << endl;
					cout << "# From R = " << R << ", Mode = " << failureMode << endl;
				}
			}
#endif
			return failureMode;
		}
	}
	
	else
	{	// 2D initiation
		
		// get principle stress and swap if needed to get s.x>s.y
		Vector s = TensorEigenvalues(str,true,true);
		if(s.x<s.y)
		{	double temp = s.x;
			s.x = s.y;
			s.y = temp;
		}

		// A line through (s1,s2) in direction (1,-1) intersects diagonal at (-P2,-P2)
		// Distance to (s1,s2) from diagonal is r = (s1-s2)/2 = s1+P2
		// A line through (-P2,P2) in units of distance r to (s1,s2) has the equation
		//     (x,y) = (-P2,-P2)+dt*r(1,-1) = (-P2+dt*(s1+P2),-P2-dt*(s1+P2))
		// project failure lines to this line in terms of d
		// See JAN019-98 for details
		double P2 = -0.5*(s.x+s.y);
		double nx=1.,ny=0.;
		double sigi = relStrength*sigmacRed;
		double taui = relStrength*taucRed;
 		
		if(failureSurface==OVOID_SURFACE && sigi>taui)
		{	// Handle this surface using Mohr's stress plot
			// always possible when sigi>taui
			
			// Mohr's circle radius and center (note that r = s.x+P2)
			double r = 0.5*(s.x-s.y);
			
			// Three cases depending on r and ellipse properties
			double tipRcrit = taui*taui/sigi;
			if(r<=tipRcrit)
			{	// Small radius curvation fails in tension when p>sigi-r
				// In principle stress terms fails when -P2>sigi-r OR sigi+P2<r
				// catch tensile-tensile failure too
				double siPlusP2 = sigi+P2;
				if(siPlusP2<0. && s.y>sigi)
				{	// tension-tension failure
					nx = s.x-sigi;
					ny = s.y-sigi;
					failureMode = TENSILE_TENSILE_FAILURE;
				}
				else if(r>siPlusP2)
				{	failureMode = TENSILE_FAILURE;
				}
			}
			
			else if(r>=taui)
			{	// shear failure
				ny = -1.;
				failureMode = SHEAR_FAILURE;
			}

			else
			{	// ovoid failure
				double c = sqrt(sigi*sigi-taui*taui);
				double arg = sqrt(taui*taui-r*r);
				double pc = c*arg/taui;
				if(-P2>pc)
				{	double theta = 0.5*acos(taui*arg/(c*r));
					nx = cos(theta);
					ny = -sin(theta);
					failureMode = TENSILE_SHEAR_FAILURE;
				}
			}
		}
		
        else
		{	// using s.x and taumax alone
            // This code only works for P_INDEPENDENT or P_LINEAR

			// handle all cases
			double siPlusP2 = sigi+P2;
			if(siPlusP2<0.)
			{	if(s.y>sigi)
				{	// tension-tension failure
					nx = s.x-sigi;
					ny = s.y-sigi;
					failureMode = TENSILE_TENSILE_FAILURE;
				}
				else
					failureMode = TENSILE_FAILURE;
			}
			
			else if(siPlusP2<taui)
			{	// tension or shear failure
				double dt = siPlusP2/(s.x+P2);
				if(dt<1.)
				{	failureMode = TENSILE_FAILURE;
					double sigmaStar = sigmacRed-2.*taui ;
					if(s.y<sigmaStar)
					{	// check if closer to shear faaiure
						double vx = s.x-sigi;
						double minusvy = sigmaStar-s.y;
						if(SQRTTWOMINUSONE*vx<minusvy)
						{	ny = -1.;
							failureMode = SHEAR_FAILURE;
						}
					}
				}
			}
			
			else
			{	// shear failure
				double ds = 2.*taui/(s.x-s.y);
				if(ds<1.)
				{	ny = -1.;
					failureMode = SHEAR_FAILURE;
				}
			}
		}
	
		// get normal if needed
		if(failureMode!=NO_FAILURE)
		{	// find angle normal to s1 axis plus s1 to global axis
			// atan2() is ccw angle from max principle direction to crack norm
			// second term is ccw angle from global x axis to max principle
			// sum is ccw angle from global x axis to crack normal
			double theta = atan2(ny,nx) + 0.5*atan2(2.*str->xy,str->xx-str->yy);
			normal->x = cos(theta);
			normal->y = sin(theta);
			normal->z = 0.;
			return failureMode;
		}
	}
	
	// has not failed yet
	return failureMode;
}

// In gap between tension and shear. Find angles to the borders
// Return true if closer to shear than to tension
bool FailureSurface::CloserToShear(Vector *s,double sigi,double P,double taus) const
{	// in gap between tension and shear - pick one
	Vector v = MakeVector(sigi+P,2.*(taus-sigi-P),sigi-2.*taus+P);
	double mag = sqrt(DotVectors(&v,&v));
	ScaleVector(&v,1./mag);
	double v2 = 1/sqrt(2.);
	v.x -= v2;
	v.y += v2;
	Vector smi = MakeVector(s->x-sigi,s->y-2.*(taus-sigi)+3.*P,s->z-sigi+2.*taus);
	if(DotVectors(&smi,&v)>0.) return true;
	return false;
}

#pragma mark FailureSurface::Handle Pressure Dependent Failure

// Methods to implement pressure dependent properties
void FailureSurface::InitPDTerms(void)
{   pdModel = P_INDEPENDENT;
}
char *FailureSurface::InputPDProperty(char *xName,int &input,double &gScaling) { return NULL; }
const char *FailureSurface::VerifyAndLoadPDTerms(double rho) { return NULL; }
void FailureSurface::PrintPDTerms(void) const { cout << endl; }
double FailureSurface::GetPDStrength(double taui,double P,double relStrength) const { return taui; }
double FailureSurface::sigmaIStabilityScale(void) const { return 1.; }
double FailureSurface::sigmaIIStabilityScale(void) const { return 1.; }
bool FailureSurface::IsPDModelConstant(double pressure,double dP,double &constantStrength) const { return true; }

// scaled mode II strength; failure surfaces that depended on other properties should change relStrength
double FailureSurface::sigmaII(double relStrength,double &sigmaIIAlpha,
                                     Tensor &str,double dPe,double C11,double nuterm,double decxx,int np) const
{   sigmaIIAlpha = -1.;
    return relStrength*taucRed;
}

// does this law have pressure dependence
bool FailureSurface::IsPressureDependent(void) const
{   return pdModel==P_INDEPENDENT ? false : true;
}

#pragma mark FailureSurface::Accessors

// reduced mode I stength
double FailureSurface::sigmaI(void) const { return sigmacRed; }

// reduced mode II stength
double FailureSurface::sigmaII(void) const { return taucRed; }

// set failure surface
void  FailureSurface::SetFailureSurface(int surface) { failureSurface = surface; }

// initiation law name
const char *FailureSurface::GetInitiationLawName(void) const { return "Isotropic material failure"; }

