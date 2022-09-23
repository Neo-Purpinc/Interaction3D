/*----------------------------------------*\
  Tuesday Novembre 26th 2019
  Arash Habibi
  Vector.h
\*----------------------------------------*/

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
	double x,y,z;
} Vector;


Vector V_new(double x, double y, double z);
void V_get(Vector v, double *x, double *y, double *z);
void V_show(Vector v, char *label);
Vector V_mult(double c, Vector v);
Vector V_add(double c1, Vector v1, double c2, Vector v2);
double V_dot(Vector v1, Vector v2);
double V_module(Vector v);
Vector V_unit(Vector v);

#endif // __VECTOR_H__
