/********************************************************************************
    CommonUtilities.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 1/12/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

// local globals
static int section=1;


#pragma mark Miscellaneous Functions

/*********************************************************************
    Start next section in results window
*********************************************************************/

void PrintSection(const char *text)
{
    char sline[200];
    
    sprintf(sline,"*****%3d. %s\n\n",section,text);
    cout << sline;
	section++;
}

/************************************************************
	See if two double precision numbers are sufficiently
		equal (currently 100 ppm)
	If both are less than 1e-12, they are assumed equal
************************************************************/

int DbleEqual(double db1,double db2)
{
    double ab1=fabs(db1);
    double ab2=fabs(db2);
    double change;
    
    if(db1==db2)
        return(TRUE);
    else if(ab1>ab2)
        change=fabs(db1-db2)/ab1;
    else
        change=fabs(db1-db2)/ab2;
            
    // Equal if different by less than 100 ppm
    if(change<1.e-4)
        return TRUE;
    else
    {	if(ab1<1e-12 && ab2<1e-12)
            return TRUE;
        else
            return FALSE;
    }
}

/************************************************************
    Reverse block of bytes to convert between Intel
	processor and Mac processor
************************************************************/

int Reverse(char *nptr,int nlen)
{
    int i,imax;
    char *rptr,temp;
    
    rptr=nptr+nlen-1;
    imax=nlen>>1;
    for(i=0;i<imax;i++)
    {	temp=*nptr;
        *nptr++=*rptr;
        *rptr--=temp;
    }
    return nlen;
}

/************************************************************
	Factor integer num into prime factors (not including 1)
************************************************************/

void PrimeFactors(int num,vector<int> &factors)
{	int maxFactor = (int)(sqrt((double)num)+.001);
	for(int i=2;i<=maxFactor;i++)
	{	// does i go into num evenly?
		if((num % i) == 0)
		{	factors.push_back(i);
			PrimeFactors(num/i, factors);
			return;
		}
	}
	// if did not find any factors, it is a prime
	factors.push_back(num);
}

#pragma mark Vector Functions

// return vector from components
Vector MakeVector(double x,double y,double z)
{	Vector v;
	v.x=x;
	v.y=y;
	v.z=z;
	return v;
}

// zero a vector and return pointer to zeroed v
Vector *ZeroVector(Vector *v)
{	v->x=v->y=v->z=0.;
	return v;
}

// copy vector and return pointer to new one
Vector *CopyVector(Vector *v,const Vector *oldv)
{	v->x=oldv->x;
	v->y=oldv->y;
	v->z=oldv->z;
	return v;
}

// scale a vector and return pointer to changed vector
Vector *ScaleVector(Vector *v,const double scale)
{	v->x*=scale;
	v->y*=scale;
	v->z*=scale;
	return v;
}

// scale a vector and copy to new one, then return pointer to new v
Vector *CopyScaleVector(Vector *v,const Vector *oldv,const double scale)
{	v->x=scale*oldv->x;
	v->y=scale*oldv->y;
	v->z=scale*oldv->z;
	return v;
}

// add vector toadd to vector v and return pointer to changed v
Vector *AddVector(Vector *v,const Vector *toadd)
{	v->x+=toadd->x;
	v->y+=toadd->y;
	v->z+=toadd->z;
	return v;
}

// sub vector toadd from vector v and return pointer to changed v
Vector *SubVector(Vector *v,const Vector *toadd)
{	v->x-=toadd->x;
	v->y-=toadd->y;
	v->z-=toadd->z;
	return v;
}

// add vector toadd scaled by scale to vector v and return pointer to changed v
Vector *AddScaledVector(Vector *v,const Vector *toadd,const double scale)
{	v->x+=toadd->x*scale;
	v->y+=toadd->y*scale;
	v->z+=toadd->z*scale;
	return v;
}

// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotVectors(const Vector *v1,const Vector *v2)
{	return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z ;
}

// Dot prodiuct of two 2D vectors
double DotVectors2D(const Vector *v1,const Vector *v2)
{	return v1->x*v2->x + v1->y*v2->y ;
}

// Vector cross product of two vectors in a third vector
Vector *CrossProduct(Vector *cp,const Vector *a,const Vector *b)
{	cp->x = a->y*b->z - b->y*a->z ;
	cp->y = b->x*a->z - a->x*b->z ;
	cp->z = a->x*b->y - b->x*a->y ;
	return cp;
}

// Print vector to cout when debugging
void PrintVector(const char *label,const Vector *v)
{	if(strlen(label)>0) cout << label;
	cout << "(" << v->x << ", " << v->y << ", " << v->z << ") ";
}

#pragma mark Tensor Functions

// zero a vector and return pointer to zeroed v
Tensor *ZeroTensor(Tensor *t)
{	t->xx=0.;
	t->yy=0.;
	t->zz=0.;
#ifdef MPM_CODE
	t->yz=0.;
	t->xz=0.;
#endif
	t->xy=0.;
	return t;
}

// add tensor toadd to tensor t and return pointer to changed t
Tensor *AddTensor(Tensor *t,Tensor *toadd)
{	t->xx+=toadd->xx;
	t->yy+=toadd->yy;
	t->zz+=toadd->zz;
#ifdef MPM_CODE
	t->yz+=toadd->yz;
	t->xz+=toadd->xz;
#endif
	t->xy+=toadd->xy;
	return t;
}

// scale tensor t and return pointer to changed t
Tensor *ScaleTensor(Tensor *t,double scale)
{	t->xx*=scale;
	t->yy*=scale;
	t->zz*=scale;
#ifdef MPM_CODE
	t->yz*=scale;
	t->xz*=scale;
#endif
	t->xy*=scale;
	return t;
}

// get element of contracted tensor by ID
double Tensor_i(Tensor *t,int contractedID)
{
	switch(contractedID)
	{	case XX:
			return t->xx;
		case YY:
			return t->yy;
		case ZZ:
			return t->zz;
#ifdef MPM_CODE
		case YZ:
			return t->yz;
		case XZ:
			return t->xz;
#endif
		case XY:
			return t->xy;
		default:
			return 0.;
	}
}

// get element of tensor by row and colume
// For FEA can only get xx or 1,1, xy or 1,2 or 2,1, yy or 2,2, anbd zz or 3,3
double Tensor_ij(Tensor *t,int row,int col)
{
	switch(row)
	{	case 1:
			switch(col)
			{	case 1:
					return t->xx;
				case 2:
					return t->xy;
#ifdef MPM_CODE
				case 3:
					return t->xz;
#endif
				default:
					return 0.;
			}
		case 2:
			switch(col)
			{	case 1:
					return t->xy;
				case 2:
					return t->yy;
#ifdef MPM_CODE
				case 3:
					return t->yz;
#endif
				default:
					return 0.;
			}
		case 3:
			switch(col)
			{
#ifdef MPM_CODE
				case 1:
					return t->xz;
				case 2:
					return t->yz;
#endif
				case 3:
					return t->zz;
				default:
					return 0.;
			}
		default:
			return 0.;
	}
}

// Print vector to cout when debugging
void PrintTensor(const char *label,Tensor *t)
{
	int i,lead=strlen(label);
	for(i=0;i<lead;i++) cout << " ";
#ifdef MPM_CODE
	cout << "(" << t->xx << ", " << t->xy << "," << t->xz << ")" << endl;
	cout << label << "(" << t->xy << ", " << t->yy << "," << t->yz << ")" << endl;
#else
	cout << "(" << t->xx << ", " << t->xy << ", 0.0)" << endl;
	cout << label << "(" << t->xy << ", " << t->yy << ", 0.0)" << endl;
#endif
	for(i=0;i<lead;i++) cout << " ";
#ifdef MPM_CODE
	cout << "(" << t->xz << ", " << t->yz << "," << t->zz << ")" << endl;
#else
	cout << "( 0.0, 0.0," << t->zz << ")" << endl;
#endif
}

#ifdef MPM_CODE

TensorAntisym *ZeroTensorAntisym(TensorAntisym *a)
{	a->xy=a->xz=a->yz=0.;
	return a;
}

#endif




