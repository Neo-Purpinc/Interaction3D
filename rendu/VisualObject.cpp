
/*================================================*\
  Mardi 19 novembre 2019
  Arash Habibi
  VisualObject.cpp

  #define VO_POINT 0
  #define VO_SPHERE 1
  #define VO_LINE 2

  char _type; // VO_POINT, VO_SPHERE, VO_LINE
  int _pt1, _pt2;
  float _radius;

\*================================================*/

#include "VisualObject.h"

static void _trigo(float *tabcos, float *tabsin, int n);
static void _drawSphere(Point pcentre, Point pobs, float rayon);
static void _drawPoint(Point pcentre);
static void _drawLine(Point p1, Point p2);

//---------------------------------------------------------------

VisualObject::VisualObject()
{
	_init(VO_POINT, -1, -1, 0);
}

//---------------------------------------------------------------
// point

VisualObject::VisualObject(char type, int pt)
{
	if(type==VO_POINT)
		_init(VO_POINT, pt, -1, 0);
	else
		cerr << "VisualObject::Visualobject : sure that this is not a point ?" << endl;
}

//---------------------------------------------------------------
// sphere

VisualObject::VisualObject(char type, int pt, float radius)
{
	if(type==VO_SPHERE)
		_init(VO_SPHERE, pt, -1, radius);
	else if(type==VO_LINE)
		_init(VO_LINE, pt, (int)radius, 0);
	else
		cerr << "VisualObject::Visualobject : sure that this is not a line or a sphere ?" << endl;
}

//---------------------------------------------------------------
// line

VisualObject::VisualObject(char type, int pt1, int pt2)
{
	if(type==VO_SPHERE)
		_init(VO_SPHERE, pt1, -1, pt2);
	else if(type==VO_LINE)
		_init(VO_LINE, pt1, pt2, 0);
	else
		cerr << "VisualObject::Visualobject : sure that this is not a line or a sphere ?" << endl;
}

//---------------------------------------------------------------

VisualObject::~VisualObject()
{
}

//---------------------------------------------------------------

void VisualObject::check(char *label)
{
	printf("-------- Check visual object %s --------------\n",label);
	switch(_type)
	{
	case VO_POINT :
		cerr << "type : points" << endl;
		cerr << "indice : " << _pt1 << endl;
		break;
	case VO_LINE :
		cerr << "type : lignes" << endl;
		cerr << "indice1 : " << _pt1 << endl;
		cerr << "indice2 : " << _pt2 << endl;
		break;
	case VO_SPHERE :
		cerr << "type : spheres" << endl;
		cerr << "indice : " << _pt1 << endl;
		cerr << "rayon : " << _radius << endl;
		break;
	}
}

//---------------------------------------------------------------

void VisualObject::display(Point *nuage, int np)
{
	Point origin(0,0,0);
	display(nuage, np, origin);
}

//---------------------------------------------------------------

void VisualObject::display(Point *nuage, int np, Point pobs)
{
	Point p1, p2;

	switch(_type)
	{
	case VO_POINT : p1 = nuage[_pt1]; _drawPoint(p1); break;
	case VO_SPHERE : p1 = nuage[_pt1]; _drawSphere(p1,pobs,_radius); break;
	case VO_LINE : p1 = nuage[_pt1]; p2 = nuage[_pt2]; _drawLine(p1,p2); break;
	default : cerr << "VisualObject::display : unknown visual object type : " << _type << endl;
	}
}

//---------------------------------------------------------------

void VisualObject::_init(char type, int pt1, int pt2, float radius)
{
	_type = type;
	_pt1 = pt1;
	_pt2 = pt2;
	_radius = radius;
}

//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------

static void _trigo(float *tabcos, float *tabsin, int n)
{
	double dtheta = 2*M_PI/n;

	for(int i=0; i<n; i++)
	{
		tabcos[i]=cos(dtheta/2 + i*dtheta);
		tabsin[i]=sin(dtheta/2 + i*dtheta);
	}
}

//------------------------------------------------------------------------

static void _drawSphere(Point pcentre, Point pobs, float rayon)
{
	static float tabcos[VO_SPHERE_RES], tabsin[VO_SPHERE_RES];
	static int nb_calls=0;

	if (nb_calls==0)
		_trigo(tabcos,tabsin,VO_SPHERE_RES);

	Vector u_obs2centre = (pcentre - pobs).unit();
	Vector u_y;
	u_y.base_y();
	Vector u_right = (u_obs2centre ^ u_y).unit();
	Vector u_up = (u_right ^ u_obs2centre).unit();

	double x,y,z;
	glBegin(GL_LINE_LOOP);
	for(int i=0;i<VO_SPHERE_RES;i++)
	{
		Point m = pcentre + u_right*(rayon*tabcos[i])  + u_up*(rayon*tabsin[i]);
		m.get(&x, &y, &z);
		glVertex3f(x,y,z);
	}
	glEnd();

	nb_calls++;
}

//------------------------------------------------------------------------

static void _drawPoint(Point pcentre)
{
	double x,y,z;
	pcentre.get(&x, &y, &z);
	glBegin(GL_POINTS);
	glVertex3f(x,y,z);
	glVertex3f(x,y-0.01,z);
	glVertex3f(x,y+0.01,z);
	glVertex3f(x-0.01,y,z);
	glVertex3f(x+0.01,y,z);
	glVertex3f(x,y,z-0.01);
	glVertex3f(x,y,z+0.01);
	glEnd();
}


//------------------------------------------------------------------------

static void _drawLine(Point p1, Point p2)
{
	double x,y,z;
	glBegin(GL_LINES);
	p1.get(&x, &y, &z);
	glVertex3f(x,y,z);
	p2.get(&x, &y, &z);
	glVertex3f(x,y,z);
	glEnd();
}
