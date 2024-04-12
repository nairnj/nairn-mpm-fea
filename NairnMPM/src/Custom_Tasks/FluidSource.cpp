/********************************************************************************
	FluidSource.cpp
		This task allows for fluid particle to be injected at a given location

	Created by Chad Hammerquist October 2017
	Revised to use Reservoir by John Nairn, 6/14/2021
 
	Parameters:
	1. mat (number) or matname (by name) and mat checked first
	2. FlowRate (number) or FlowFunction (expression) in volume per unit time (mm^3/s in Legacy)
	3. (flow_x,flow_y,flow_z) for vector in flow direction. 2D only needs x and y.
			3D vector must be along x, y, or z axis. Need not be normalized.
	4. (source_x,source_y,source_z) center if inlet
    5. inlet_width,inlet_depth: inlet area. 2D only uses width. 3D are dimensions along
			the other two axes.
	6. particle_size: can adjust particle size relative to default particle size in
			the simulatinos
    7. Reservoir must have that material type and size, otherwise nothinh injected
 
	Ideas:
	1. Tartan grid: It might be close to working, but may need to be careful about
		particle size an about the region with velocity control
	2. Could setting pressure help and may need to track current pressure
	3. 3D pipes and off-axis options
*********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/FluidSource.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/Reservoir.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"

//#define LOCAL_TAIT_C 0.0894

// Constructors
FluidSource::FluidSource() : CustomTask()
{
	//  Initial values
	material = -1;
	matName = NULL;
	FlowRate = -1.0;
	ZeroVector(&FlowNormal);
	width = 0.0;
	depth = 0.0;
	ZeroVector(&sourcePt);
	FlowRate_Function = NULL;

	// private
	inlet_area = 0.0;
	particle_size = 1.;
	particles_per_width = 1;
	particles_per_depth = 1;
	scale_particle_size = 1;
}

// Return name of this task
const char *FluidSource::TaskName(void) { return "Inject Fluid Particles at given location"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *FluidSource::InputParam(char *pName, int &input, double &gScaling){

	if (strcmp(pName, "material") == 0)
	{	input = INT_NUM;
		return (char *)&material;
	}
	else if(strcmp(pName,"matname") == 0)
	{	input = TEXT_PARAMETER;
		return (char *)&matName;
	}
	else if (strcmp(pName, "FlowRate") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&FlowRate;
	}
	else if (strcmp(pName, "FlowFunction") == 0)
	{	input = TEXT_PARAMETER;
		return (char *)&FlowRate_Function;
	}
	else if (strcmp(pName, "flow_x") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&FlowNormal.x;
	}
	else if (strcmp(pName, "flow_y") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&FlowNormal.y;
	}
	else if (strcmp(pName, "flow_z") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&FlowNormal.z;
	}
	else if (strcmp(pName, "source_x") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&sourcePt.x;
	}
	else if (strcmp(pName, "source_y") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&sourcePt.y;
	}
	else if (strcmp(pName, "source_z") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&sourcePt.z;
	}
	else if (strcmp(pName, "inlet_width") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&width;
	}
	else if (strcmp(pName, "inlet_depth") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&depth;
	}
	else if (strcmp(pName, "particle_size") == 0)
	{	input = DOUBLE_NUM;
		return (char *)&scale_particle_size;
	}
	
	return CustomTask::InputParam(pName, input,gScaling);
}

// Actually sets flow rate function
// throws std::bad_alloc, SAXException()
void FluidSource::SetTextParameter(char *fxn, char *ptr)
{
	// FlowRate_Function
	if(ptr == (char *)&FlowRate_Function)
	{	if(FlowRate_Function != NULL)
			ThrowSAXException("Duplicate flow rate function was supplied");
		if(fxn == NULL)
			ThrowSAXException("Flow rate function is missing");
		if(strlen(fxn) == 0)
			ThrowSAXException("Flow rate function is empty");

		FlowRate_Function = Expression::CreateExpression(fxn,"Flow rate function is not valid");
	}
	else if(ptr == (char *)&matName)
	{	// material name needed
		if(matName!=NULL)
			ThrowSAXException("Duplicate material name supplied to task");
		if(fxn==NULL)
			ThrowSAXException("Material name is missing");
		if(strlen(fxn)==0)
			ThrowSAXException("Material name is missing");
		
		// save is
		matName = new char[strlen(fxn)+1];
		strcpy(matName,fxn);
	}
	else
		CustomTask::SetTextParameter(fxn, ptr);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *FluidSource::Initialize(void)
{
	// requires reservoir
	if(mpmReservoir==NULL)
		throw CommonException("FluidSource custom task requires an available reservoir (i.e., a structured grid)",
							  "FluidSource::Initialize");
	
	// Use mat input, but if not provided or invalid, looks for name
	if(material <= 0 || material > nmat)
	{	if(matName == NULL)
			throw CommonException("No material defined", "DeleteDamaged::Initialize");
		
		for(int i=0;i<nmat;i++)
		{	if(strcmp(matName,theMaterials[i]->name)==0)
			{	if(material>0)
					throw CommonException("More than one material has the requested name", "DeleteDamaged::Initialize");
				material = i+1;
			}
		}
		
		if(material <= 0 || material > nmat)
			throw CommonException("No material matching the specified name", "DeleteDamaged::Initialize");
	}
	
	// Make sure things are in multiples of particle sizes
	if (depth == 0) depth = width;
	double inlet_width = width;
	double inlet_depth = depth;
	particle_size = scale_particle_size*mpmgrid.GetMinCellDimension()/((double)fmobj->ptsPerSide);
	if(fmobj->IsThreeD())
	{	// for x flow, yz or width/depth, for y flow, xz are width/depth, for z flow xy are width depth
		particles_per_width = (width <= particle_size) ? 1 : int(round(width / particle_size));
		particles_per_depth = (depth <= particle_size) ? 1 : int(round(depth / particle_size));
		inlet_width = particle_size*double(particles_per_width);
		inlet_depth = particle_size*double(particles_per_depth);
	}
	else
	{	particles_per_width = (width <= particle_size) ? 1 : int(round(width / particle_size));
		inlet_width = particle_size*double(particles_per_width);
		inlet_depth = mpmgrid.GetThickness();
	}
	lpart = MakeVector(particle_size,particle_size,particle_size);
	inlet_area = inlet_depth*inlet_width;

	// normalize flow vector
	if(!fmobj->IsThreeD()) FlowNormal.z = 0.;
	double flowMag =sqrt(DotVectors(&FlowNormal,&FlowNormal));
	if(DbleEqual(flowMag,0.0))
		throw CommonException("The magnitude of the flow vector is zero", "FluidSource::Initialize");
	ScaleVector(&FlowNormal,1./flowMag);
	if(fmobj->IsThreeD())
	{	// currently limited to on axis
		if(!DbleEqual(FlowNormal.x,0.))
		{	if(!DbleEqual(FlowNormal.y,0.) || !DbleEqual(FlowNormal.z,0.))
				throw CommonException("3D flow is limited to x, y, or z direction", "FluidSource::Initialize");
			tangent1 = MakeVector(0.,-FlowNormal.x,0.);
			tangent2 = MakeVector(0.,0.,-FlowNormal.x);
		}
		else if(!DbleEqual(FlowNormal.y,0.))
		{	if(!DbleEqual(FlowNormal.z,0.))
				throw CommonException("3D flow is limited to x, y, or z direction", "FluidSource::Initialize");
			tangent1 = MakeVector(-FlowNormal.y,0.,0.);
			tangent2 = MakeVector(0.,0.,-FlowNormal.y);
		}
		else
		{	tangent1 = MakeVector(-FlowNormal.z,0.,0.);
			tangent2 = MakeVector(0.,-FlowNormal.z,0.);
		}
	}
	else
	{	tangent1 = MakeVector(FlowNormal.y,-FlowNormal.x,0.);
	}

	// get corners of inlet plane
	inlet_plane[0] = sourcePt;
	inlet_plane[1] = sourcePt;
	double halfWidth = 0.5*inlet_width;
	if(fmobj->IsThreeD())
	{	// corners at s +/- (t1*width/2 + t2*width/2)
		double halfDepth = 0.5*inlet_depth;
		AddScaledVector(&inlet_plane[0],&tangent1,halfWidth);
		AddScaledVector(&inlet_plane[0],&tangent2,halfDepth);
		AddScaledVector(&inlet_plane[1],&tangent1,-halfWidth);
		AddScaledVector(&inlet_plane[1],&tangent2,-halfDepth);
	}
	else
	{	// plane from s + t*width/2 to s - t*width/2 where t = (ny,-nx)
		AddScaledVector(&inlet_plane[0],&tangent1,halfWidth);
		AddScaledVector(&inlet_plane[1],&tangent1,-halfWidth);
	}

	// get velocity BC region of interest
	// it is block including 1 cell ahead of inlet plane. The extra 1/4 cell make sure gets all  nodes
	//		in that region (region in front could likely be cell instead of 1.25*cell)
	double region_margin = 0.25*mpmgrid.GetMinCellDimension();
	double region_in_front_of_inlet = 1.25*mpmgrid.GetMinCellDimension();
	double regionHeight = region_margin + region_in_front_of_inlet;
	double regionWidth = inlet_width + 2.*region_margin;
	mvrBox[0] = inlet_plane[0];

	if(fmobj->IsThreeD())
	{	// lower-left
		double regionDepth = inlet_depth + 2.*region_margin;
		AddScaledVector(&mvrBox[0],&tangent1,region_margin);
		AddScaledVector(&mvrBox[0],&tangent2,region_margin);
		AddScaledVector(&mvrBox[0],&FlowNormal,-region_margin);
		mvrBox[1] = mvrBox[0];
		AddScaledVector(&mvrBox[1],&tangent1,-regionWidth);
		AddScaledVector(&mvrBox[1],&tangent2,-regionDepth);
		AddScaledVector(&mvrBox[1],&FlowNormal,regionHeight);
	}
	else
	{	// box around velocity region
		AddScaledVector(&mvrBox[0],&tangent1,region_margin);
		AddScaledVector(&mvrBox[0],&FlowNormal,-region_margin);
		mvrBox[1] = mvrBox[0];
		AddScaledVector(&mvrBox[1],&FlowNormal,regionHeight);
		mvrBox[2] = mvrBox[1];
		AddScaledVector(&mvrBox[2],&tangent1,-regionWidth);
		mvrBox[3] = mvrBox[2];
		AddScaledVector(&mvrBox[3],&FlowNormal,-regionHeight);
	}

	cout << "Fluid Source custom task created" << endl;
	PrintVector("   Source point = ",&sourcePt);
	cout << endl;
	cout << "   Inlet = " << inlet_width << " X " << inlet_depth << " = "
			<< inlet_area << " " << UnitsController::Label(CULENGTH_UNITS) << "^2 " << endl;


	// flow function	
	if(FlowRate<0 && FlowRate_Function == NULL)
		throw CommonException("No Flowrate defined", "FluidSource::Initialize");
	if(FlowRate_Function == NULL)
	{	cout << "   Flowrate = " << FlowRate << " " << UnitsController::Label(CULENGTH_UNITS)
				<< "^3/" << UnitsController::Label(TIME_UNITS) << endl;
	}
	else
	{	cout << "   Flow rate function = " << FlowRate_Function->GetString() << "in " << UnitsController::Label(CULENGTH_UNITS)
				<< "^3/" << UnitsController::Label(TIME_UNITS) << endl;
	}
	PrintVector("   Flow direction = ",&FlowNormal);
	cout << endl;

	cout << "   Injected particle size = " << scale_particle_size/((double)fmobj->ptsPerSide) << " cells" << endl;
	cout << "   Injected material number = " << material << endl;

	// initialize to > particle so that first time step will need particles
	previous_set_displacement = MakeVector(2.*particle_size, 2.*particle_size, 2.*particle_size);
	ZeroVector(&previous_velocity);

	return nextTask;
}

// do custom calculation
CustomTask *FluidSource::StepCalculation(void)
{
	// Get flow rate:
	double flowrate = 0;
	if (FlowRate_Function == NULL)
		flowrate = FlowRate;
	else
		flowrate = FlowRate_Function->TValue(mtime*UnitsController::Scaling(1000.));

	// if there is no flowrate, then we are done here
	if(flowrate == 0.)  return nextTask;

	// update displacement since last injection
	AddScaledVector(&previous_set_displacement, &previous_velocity, timestep);

	// update velocities
	double velocity = flowrate/inlet_area;
	Vector inlet_velocity = MakeVector(0, 0, 0);
	AddScaledVector(&inlet_velocity, &FlowNormal, velocity);
	previous_velocity = inlet_velocity;

	// do we need new particles? (note displacement assumed in one direction only)
	double magDisp = DotVectors(&previous_set_displacement,&previous_set_displacement);
	bool NeedMoreParticles =  magDisp >= particle_size*particle_size;
	
	// when injecting, reset previous displacement to zero
	if(NeedMoreParticles) ZeroVector(&previous_set_displacement);

	// loop through all particles in region of interest update velocities
#pragma omp for
	for(int p = 0; p<nmpmsNR; p++)
	{	MPMBase *point = mpm[p];
		if(point->InReservoir()) continue;
		
		// only look a particles of the injecting material
		if((point->MatID()+1) != material) continue;

		// set velocity boundary conditions for particles in this region
		if(PointInHere(point, mvrBox))
		{	// position update relative to current velocity
			double dx = point->vel.x - inlet_velocity.x;
			double dy = point->vel.y - inlet_velocity.y;
			double dz = point->vel.z - inlet_velocity.z;
			point->pos.x -= dx*timestep;
			point->pos.y -= dy*timestep;
			point->pos.z -= dz*timestep;

			// overwrite velocity
			point->SetVelocity(&inlet_velocity); 
		}
	}

	// if we don't need new particles, exit here
	if(!NeedMoreParticles) return nextTask;

	// Inject new particles
	if(fmobj->IsThreeD())
	{	for(int i=0; i<particles_per_width; i++)
		{	Vector Inlet_location = inlet_plane[0];
			AddScaledVector(&Inlet_location,&tangent1,-(i+0.5)*particle_size);
			AddScaledVector(&Inlet_location,&tangent2,-0.5*particle_size);
			for(int j=0; j<particles_per_depth; j++)
			{	if(!InjectParticle(&Inlet_location, &inlet_velocity)) break;
				AddScaledVector(&Inlet_location,&tangent2,-particle_size);
			}
		}
	}
	else
	{	// 2D injection
		Vector Inlet_location = inlet_plane[0];
		AddScaledVector(&Inlet_location,&tangent1,-0.5*particle_size);
		for(int newp=0;newp<particles_per_width;newp++)
		{	if(!InjectParticle(&Inlet_location,&inlet_velocity)) break;
			AddScaledVector(&Inlet_location,&tangent1,-particle_size);
		}
	}

	// next
	return nextTask;

}

// In region of interest
bool FluidSource::PointInHere(MPMBase *mptr, Vector rect[])
{
	Vector pt = mptr->pos;
	
	if(fmobj->IsThreeD())
	{	if(pt.x<fmin(rect[0].x,rect[1].x) || pt.x>fmax(rect[0].x,rect[1].x)) return false;
		if(pt.y<fmin(rect[0].y,rect[1].y) || pt.y>fmax(rect[0].y,rect[1].y)) return false;
		if(pt.z<fmin(rect[0].z,rect[1].z) || pt.z>fmax(rect[0].z,rect[1].z)) return false;
		return true;
	}
	
	// 2D uses line cross on box with corners in rect[0] to rect[4]
	unsigned i,crossings=0;
	double x1,x2,y1,y2,d;
	
	x1=rect[0].x;
	y1=rect[0].y;
	for(i=1;i<4;i++)
	{	x2=rect[i].x;
		y2=rect[i].y;
		d=(pt.y-y1)*(x2-x1)-(pt.x-x1)*(y2-y1);
		
		// get crossing unless both y's on same side of edge
		if((y1>=pt.y) != (y2>=pt.y))
			crossings += (y2-y1>=0.) ? d>=0. : d<=0. ;
		
		// if d is 0, check if point is on line (and thus in polygon)
		if(!d && fmin(x1,x2)<=pt.x && pt.x<=fmax(x1,x2) &&
							fmin(y1,y2)<=pt.y && pt.y<=fmax(y1,y2))
		{	return true;
		}
		x1=x2;
		y1=y2;
	}
	return (crossings & 0x01);
}

// Function to handle injection particles
bool FluidSource::InjectParticle(Vector *location, Vector *velocity)
{
	MPMBase *mptr = mpmReservoir->InjectParticle(location,&lpart,material);
	if(mptr==NULL) return false;
	
	mptr->SetVelocity(velocity);
	
	if(!fmobj->IsThreeD())
	{	// rotate to inlet plane
		Matrix3 rot = Matrix3::Identity();
		rot.set(0,0,FlowNormal.x);
		rot.set(1,0,FlowNormal.y);
		rot.set(0,1,-FlowNormal.y);
		rot.set(1,1,FlowNormal.x);
		mptr->SetDeformationGradientMatrix(rot);
	}
	
	return true;
}
