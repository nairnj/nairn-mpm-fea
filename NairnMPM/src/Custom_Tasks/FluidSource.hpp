//
//  FluidSource.hpp
//  NairnMPM
//
//  Created by Chad Hammerquist, 2017
//
//

#ifndef _FLUIDSOURCE_

#define _FLUIDSOURCE_

#include "Custom_Tasks/CustomTask.hpp"
class Expression;

class FluidSource : public CustomTask
{
    public:
        // constructors and destructors
        FluidSource();
        // Constants

        // standard methods
        virtual const char *TaskName(void);
        virtual char *InputParam(char *, int &, double &);
        virtual CustomTask *Initialize(void);
        virtual CustomTask *StepCalculation(void);
        virtual void SetTextParameter(char *, char *);


    private:

        // input values
        char *matName;
        int material;
        double FlowRate, width,depth;
        Vector sourcePt;
        double inlet_area, particle_size;
        int particles_per_width,particles_per_depth;
        double scale_particle_size;
        Vector lpart;
        Expression *FlowRate_Function;
        Vector mvrBox[4];
        Vector previous_set_displacement;
        Vector previous_velocity;
        Vector inlet_plane[2];
        Vector FlowNormal;
        Vector tangent1,tangent2;
        bool PointInHere(MPMBase * mptr, Vector rect[]);
        bool InjectParticle(Vector *, Vector *);
};

#endif
