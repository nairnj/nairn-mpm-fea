/********************************************************************************
	PolyhedronController.cpp
	nairn-mpm-fea

	Created by John Nairn on 1/7/11.
	Copyright (c) 2011 John A. Nairn, All rights reserved.
********************************************************************************/

#include "PolyhedronController.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_MPM/PolyTriangle.hpp"

#pragma mark PolyhedronController: initializers

// constructor
PolyhedronController::PolyhedronController(int block) : ShapeController(block)
{
}

// destructor
PolyhedronController::~PolyhedronController()
{	unsigned i;
	for(i=0;i<faces.size();i++)
		delete faces[i];
	faces.clear();
}

// set from  arbitrary string data
// not thread safe due to push_back()
// throws std::bad_alloc, SAXException()
void PolyhedronController::SetProperty(char *bData,CommonReadHandler *reader)
{
	int i,numfaces;
	vector<double> pts;
	Vector box[8],xn,yn,zn,origin;
	
	// read the data into a vector
	if(!CommonReadHandler::GetFreeFormatNumbers(bData,pts,distScaling))
        throw SAXException("Invalid data block passed to a <Polyhedron> object - bad number formatting");
	//cout << -180.*atan((pts[5]-pts[2])/(pts[3]-pts[0]))/3.14159 << endl;
	
	switch(style)
	{	case TRICLINIC_POINTS:
		case TRICLINIC_VECTORS:
		{	// requires 12 points
			if(pts.size()!=12)
                throw SAXException("Invalid data block passed to a <Polyhedron> object - wrong number of points");
			
			// Find origin, vectors and make sure first four are points
			origin = MakeVector(pts[0],pts[1],pts[2]);
			if(style==TRICLINIC_POINTS)
			{	xn = MakeVector(pts[3]-pts[0],pts[4]-pts[1],pts[5]-pts[2]);
				yn = MakeVector(pts[6]-pts[0],pts[7]-pts[1],pts[8]-pts[2]);
				zn = MakeVector(pts[9]-pts[0],pts[10]-pts[1],pts[11]-pts[2]);
			}
			else
			{	xn = MakeVector(pts[3],pts[4],pts[5]);
				yn = MakeVector(pts[6],pts[7],pts[8]);
				zn = MakeVector(pts[9],pts[10],pts[11]);
				for(i=0;i<3;i++)
				{	pts[3+i]+=pts[i];
					pts[6+i]+=pts[i];
					pts[9+i]+=pts[i];
				}
			}
			
			// get remaining 4 points to get 24, then fall through to BOX_CORNERS
			Vector v;
			AddVector(AddVector(CopyVector(&v,&origin),&xn),&yn);
			pts.push_back(v.x);
			pts.push_back(v.y);
			pts.push_back(v.z);
			AddVector(&origin,&zn);
			AddVector(CopyVector(&v,&origin),&xn);
			pts.push_back(v.x);
			pts.push_back(v.y);
			pts.push_back(v.z);
			AddVector(&v,&yn);
			pts.push_back(v.x);
			pts.push_back(v.y);
			pts.push_back(v.z);
			SubVector(&v,&xn);
			pts.push_back(v.x);
			pts.push_back(v.y);
			pts.push_back(v.z);
			
			strcpy(order,"12453678");
		}
		case BOX_CORNERS:
			if(pts.size()!=24)
                throw SAXException("Invalid data block passed to a <Polyhedron> object - wrong number of points");
			
			// get box corners in standard order (first 4 one face, second four opposite face)
			for(i=0;i<8;i++)
				box[order[i]-'1'] = MakeVector(pts[3*i],pts[3*i+1],pts[3*i+2]);
			
			// Twelve Triangles
			// 0123 face
			faces.push_back(new PolyTriangle(box[0],box[1],box[2]));
			faces.push_back(new PolyTriangle(box[0],box[2],box[3]));
			// 4567 face
			faces.push_back(new PolyTriangle(box[4],box[5],box[6]));
			faces.push_back(new PolyTriangle(box[4],box[6],box[7]));
			// 0154 face
			faces.push_back(new PolyTriangle(box[0],box[1],box[5]));
			faces.push_back(new PolyTriangle(box[0],box[5],box[4]));
			// 1265 face
			faces.push_back(new PolyTriangle(box[1],box[2],box[6]));
			faces.push_back(new PolyTriangle(box[1],box[6],box[5]));
			// 2376 face
			faces.push_back(new PolyTriangle(box[2],box[3],box[7]));
			faces.push_back(new PolyTriangle(box[2],box[7],box[6]));
			// 3047 face
			faces.push_back(new PolyTriangle(box[3],box[0],box[4]));
			faces.push_back(new PolyTriangle(box[3],box[4],box[7]));
			
			break;
		
		case PYRAMID:
			// 12 is tetrahdron, 15 is square-bottom pyramid
			if(pts.size()!=12 && pts.size()!=15)
                throw SAXException("Invalid data block passed to a <Polyhedron> object - wrong number of points");
			
			// get corner vectors
			// 0 is apex, rest are ccw around the base (not sure if ccw or cw matters)
			numfaces = pts.size()==12 ? 4 : 5 ;
			for(i=0;i<numfaces;i++)
				box[i] = MakeVector(pts[3*i],pts[3*i+1],pts[3*i+2]);
			
			// face from apex
			faces.push_back(new PolyTriangle(box[0],box[1],box[2]));
			faces.push_back(new PolyTriangle(box[0],box[2],box[3]));
			if(numfaces==4)
			{	faces.push_back(new PolyTriangle(box[0],box[3],box[1]));
				faces.push_back(new PolyTriangle(box[1],box[2],box[3]));
			}
			else
			{	faces.push_back(new PolyTriangle(box[0],box[3],box[4]));
				faces.push_back(new PolyTriangle(box[0],box[4],box[1]));
				faces.push_back(new PolyTriangle(box[1],box[2],box[3]));
				faces.push_back(new PolyTriangle(box[3],box[4],box[1]));
			}
			break;
			
		default:
			// invalid style
			throw SAXException("Invalid data block passed to a <Polyhedron> object - unrecognized style");
	}
}

// set a property - only style
void PolyhedronController::SetParameter(const char *aName,const char *value)
{	
	if(strcmp(aName,"style")==0)
	{	if(strcmp(value,"tripts")==0)
            style = TRICLINIC_POINTS;
        else if(strcmp(value,"trivectors")==0)
            style = TRICLINIC_VECTORS;
		else if(strcmp(value,"pyramid")==0)
			style = PYRAMID;
        else if(strlen(value)==8)
        {	int i;
            for(i=0;i<8;i++)
            {	if(value[i]<'1' || value[i]>'8')
                {	style = NO_FACES;
                    break;
                }
            }
            style = BOX_CORNERS;
            strcpy(order,value);
        }
        else
            style = NO_FACES;
	}
}

// Return if has enough parameters to use.
bool PolyhedronController::HasAllParameters(void)
{	if(faces.size()<4) return FALSE;
	
	unsigned i;
	Vector fmin,fmax;
	faces[0]->GetExtents(&pmin,&pmax);
	for(i=1;i<faces.size();i++)
	{	faces[i]->GetExtents(&fmin,&fmax);
		if(fmin.x < pmin.x) pmin.x=fmin.x;
		if(fmin.y < pmin.y) pmin.y=fmin.y;
		if(fmin.z < pmin.z) pmin.z=fmin.z;
		if(fmax.x > pmax.x) pmax.x=fmax.x;
		if(fmax.y > pmax.y) pmax.y=fmax.y;
		if(fmax.z > pmax.z) pmax.z=fmax.z;
	}
	return TRUE;
}

// called read attributes to verify all attributes were there
bool PolyhedronController::FinishParameter(void) { return style!=NO_FACES; }

// called after initialization is done, need to wait for pt parameters
bool PolyhedronController::FinishSetup(void) {	return FALSE; }

#pragma mark PolyhedronController: methods

// return true if point is in this body
// To handle edges, this assumes the faces enclose and area with no open
// space.
bool PolyhedronController::ContainsPoint(Vector& pt)
{	
	unsigned i,crossings=0,edges=0;
	int cross;
	
	// screen extent first
	if(pt.x<pmin.x || pt.x>pmax.x || pt.y < pmin.y || pt.y > pmax.y || pt.z < pmin.z || pt.z > pmax.z) return FALSE;
	
	for(i=0;i<faces.size();i++)
	{	// get parameteric intersection point
		cross=faces[i]->PointCrossesFace(&pt,&edges);
		if(cross==0) return TRUE;
		if(cross>0) crossings+=1;
	}
	// See if totals crossings (- half edges that were double counted) is odd
	return ((crossings+(edges>>1)) & 0x01);
}

#pragma mark PolyhedronController: accessors

// override for 3D objects
bool PolyhedronController::Is2DShape(void) { return FALSE; }

// type of object
const char *PolyhedronController::GetShapeName(void) { return "Polyhedron"; }

