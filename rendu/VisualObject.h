
/*================================================*\
  Mardi 19 novembre 2019
  Arash Habibi
  VisualObject.h
\*================================================*/

#ifndef __VISUAL_OBJECT_H__
#define __VISUAL_OBJECT_H__

#include <iostream>
#include <stdlib.h>
#include <GL/gl.h>

#include "Vector.h"

#define VO_POINT 0
#define VO_SPHERE 1
#define VO_LINE 2

#define VO_SPHERE_RES 32

class VisualObject
{
public:
	VisualObject();
	VisualObject(char type, int pt); // point
	VisualObject(char type, int pt, float radius); // sphere
	VisualObject(char type, int pt1, int pt2); // line
	~VisualObject();

	void check(char *);

	void display(Point *nuage, int np);
	void display(Point *nuage, int np, Point pobs);

protected:

	char _type; // VO_POINT, VO_SPHERE, VO_LINE
	int _pt1, _pt2;
	float _radius;

	void _init(char type, int pt1, int pt2, float radius);
};


#endif // __VISUAL_OBJECT_H__
