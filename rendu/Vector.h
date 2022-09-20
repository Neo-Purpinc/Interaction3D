/*======================================================*\
  Friday November the 2012
  Arash Habibi
  Vector.h
\*======================================================*/


#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <fstream>
#include <iostream>
#include <stdbool.h>
#include <stdlib.h>		// pour rand
#include <math.h>

using namespace std;

// #define V_NB_COMP (3)
#define V_NB_COMP (4)

#define V_X 0
#define V_Y 1
#define V_Z 2
#define V_W 3

#define V_POINT (1)
#define V_VECTOR (0)

class Vector
{
public :
	Vector();
	Vector(double x, double y, double z, double w);
	Vector(double x, double y, double z);
	~Vector();

	void check(string message);
	void check(ofstream, string message);
	void get(double *x, double *y, double *z);
	void get(double *x, double *y, double *z, double *w);
	void get(double *value, int comp);
	double get(int comp);
	void set(const double x_ref, const double y_ref, const double z_ref);
	void set(const double x_ref, const double y_ref, const double z_ref, const double w_ref);
	void set(int comp, const double value);
	void setRandom(double rayon);

	void base_x();
	void base_y();
	void base_z();
	void base_w();

	double length();
	double length2();
	double maxXyz();
	Vector unit();
	void normalize();
	Vector shuffle();
	bool isZero();
	bool orthotrope();
	int closestAxis();

	void writeVect(ofstream filmFile);

	double operator[](int comp);      // Read component
	Vector operator=(const Vector );  // Assignment
	bool   operator==(const Vector ); // Equality in value
	bool   operator!=(const Vector ); // Inequality in value
	Vector operator+(const Vector );  // Vector addition
	Vector operator-(const Vector );  // Vector substraction
	double operator*(const Vector );  // Dot product
	Vector operator*(const double);   // Vector-Scalar product
	Vector operator^(const Vector );  // Cross product

	void getPointer(double **x, double **y, double **z);

protected :
	double _data[V_NB_COMP];
};

Vector operator * (double k, Vector V);

//-------------------------------------------------------

class Point : public Vector
{
public :
	Point();
	Point(double x, double y, double z, double w);
	Point(double x, double y, double z);
	~Point();

	Point operator=(const Vector );  // Assignment
	Point operator=(const Point  );  // Assignment
	Point operator+(Vector );  // point + vector = point
	Point operator-(Vector );  // point - vector = point
	Vector operator-(Point );  // point - point = vector
	Point  operator*(double);  // point * scalar = poinit
	double operator|(Point );  // distance
	double operator||(Point ); // distance square
	double operator&(Point ); // max distance
	Point operator%(Point );   // random

	double blendingFactor(Point A, Point B);

	double distToPlane(Point P, Vector n);
	double distToPlane(Point A, Point B, Point C);
	double distToLine(Point A, Point B);
	double distToLineSegment(Point A, Point B);
	double distToCircle(Point C, double radius, Vector normal);

	Point orthoProjectionOnPlane(Point P, Vector n);
	Point orthoProjectionOnPlane(Point A, Point B, Point C);
	Vector orthoProjectVectorOnLine(Point A, Point B);
	Point orthoProjectionOnLine(Point A, Point B);

	Point projectionOnSphere(Vector u_ray, Point C, double radius);
	Point projectionOnTrackBall(Vector u_ray, Point C, double radius);
	bool projectionOnPlaneExists(Vector u, Point P, Vector normal);
	bool projectionOnPlaneExists(Vector u, Point A, Point B, Point C);
	Point projectOnPlane(Vector direction, Point P, Vector n);
	Point projectOnPlane(Vector direction, Point A, Point B, Point C);
	Point linePlaneIntersect(Vector u, Point P, Vector normal);
	int lineSphereIntersect(Vector u, Point C, double radius, Point *p1, Point *p2);

	Point closestPointOnLineSegment(Point A, Point B);
	Point closestPointOnCircle(Point C, double radius, Vector normal);

	Point turnAround(Point A, Point B, double angle);
	Point turnAroundCam(Point A, Point B, double angle);
	Point turnAround(Point O, Point A, Point B);

	bool isOnTheSamePlaneAs(Point A, Point B, Point C);

	void placeAtRandomInTriangle(Point A, Point B, Point C);

	bool isInBoundingBox(Point bbmin, Point bbmax);
	double maxAxisDistFromBoundingBox(Point bbmin, Point bbmax);
	double maxAxisDistFromBoundingSphere(Point center, double radius);

	bool isOnTheRight(Point A, Point B, Vector up);
	bool isOnTheFirstEndSide(Point A, Point B);
	bool isOnTheSecondEndSide(Point A, Point B);
};

Point operator * (double k, Point V);

//-------------------------------------------------------

extern double positiveTriangleArea(Point A, Point B, Point C);
extern double triangleArea(Point A, Point B, Point C, Vector normal);
extern double tetrahedronVolume(Point A, Point B, Point C, Point D);

extern void V_oppositePointsOnParallelLines(Point p1, Point p2, Point q1, Point q2,	double *t_p1_onq, double *t_p2_onq, double *t_q1_onp, double *t_q2_onp);
extern int V_oppositePointsOnParallelSegments(Point p1, Point p2, Point q1, Point q2, double *t_p1_onq, double *t_p2_onq, double *t_q1_onp, double *t_q2_onp);
extern void V_lineLineClosestPoints(Point p1, Point p2, Point q1, Point q2, double *t_p, double *t_q);
extern void V_segmentSegmentClosestPoints(Point p1, Point p2, Point q1, Point q2, double *t_p, double *t_q);

extern int V_lineProjectionOnCircle(Point p1, Point p2, Point p_center, double radius, Vector normal, Point *pc1, Point *pc2);
extern int V_segmentCircleClosestPoints(Point p1, Point p2, Point p_center, double radius, Vector normal, Point *pc1, Point *pc2, Point *pl1, Point *pl2, double epsilon);

extern int V_segmentCircleClosestPoints(Point p1, Point p2, Point p_center, double radius, Vector normal, Point *pc1, Point *pc2, Point *pl1, Point *pl2);

extern int V_segmentSphereIntersect(Point A, Point B, Point Center, double radius, Point *P1, Point *P2);
extern bool V_lineSegmentPlaneIntersect(Point p1, Point p2, Point P, Vector normal, double *alpha);

extern void V_circleCircleIntersection(double R1, double R2, double D, double *r1, double *r2, double h);

extern double V_rotationAngle(Point A, Point B, Vector vA, Vector vB, Vector normal, Vector *u_axis);


#endif // _VECTOR_H_
