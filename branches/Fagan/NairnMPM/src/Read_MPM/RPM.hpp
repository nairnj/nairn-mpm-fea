/********************************************************************************
    RPM.hpp
    NairnMPM
    
    Created by Tim Fagan on Thu Dec 1 2011.	
********************************************************************************/

#ifndef _RPM_

#define _RPM_

#include "MPM_Classes/MPMBase.hpp" //modiftf #rigidbodyrotation

class Rpm
{
    public:
        double thick;
		static bool rpmApplied;
		double phi, sx, sy, rotx, roty;
		double xcentre, ycentre, zcentre, zDepth;
		double rpm, materialID;
		Vector velSet;
		double simTime;
				
        
        // constructors and destructors
        Rpm();
        //RPM(int,int,double,double);
        
        // methods
		void SetRPM(double, double); // if needed may be a virtual function
		bool CheckRPM(int);
		void SetCentre(double, double, double);
		void SetRotationAngle(double);
		void AddRPM(double &, double &, double &, Vector);
		void AddRPM2(double &, double &, double &, Vector &, Vector *, MPMBase *, double, double);
		void AddRPM3(double &, double &, Vector *, MPMBase *, double, double);
		void SetVel(Vector);
		void UpdateCentre(double);
		void LogTime(double);
		double getxCentre(void);
		double getyCentre(void);
		double getzDepth(void);

};

extern Rpm *rotator;

#endif