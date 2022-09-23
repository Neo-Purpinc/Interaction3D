
/*===============================================================*\
  Thursday June the 20th 2013
  Arash HABIBI
  Camera.h
\*===============================================================*/

#ifndef _CAMERA_H_
#define _CAMERA_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Vector.h"

#define CAM_SPACE_INCREMENT 1
#define CAM_ANGLE_INCREMENT 5.1

class Camera
{
public :
	Camera();
	Camera(int width, int height);
	~Camera();

	void initPosit();
	void initSize();

	void check(char *message);

	void setPixelDimensions(int width, int height);

	Point lookfrom();
	Point lookat();
	Vector upVector();
	float fovy();
	float ratio();
	float zNear();
	float zFar();

	void zoom(float dtheta);

	void move(Vector translate);

	void turnAroundX(float dtheta);
	void turnAroundY(float dtheta);
	void turnAroundZ(float dtheta);
	void turnAroundWorldY(float dtheta);

	void turnAroundTrackBall(int previous_x, int previous_y, int x, int y, double sensitivity);

	void frameBoundingSphere(Point p_center, double radius);

	Point pixelPositionIn3D(int x, int y);
	Point pixelPositionInXYPlane(int x, int y);
	void windowHalfSizeIn3D(float *half_width, float *half_height);

	Point* getLookFromPointer();
	Point* getLookAtPointer();

private:

	Point _lookfrom, _lookat;
	Vector _aim, _up_vector, _right_vector;
	int _pixel_width, _pixel_height;
	float _3D_half_width, _3D_half_height; // the width of the window at _lookat
	float _znear, _zfar, _dist_from_at;
	float _half_focal, _focal_length;
	float _trackball_radius, _ratio;

	float _half_focal_to_focal_length(float half_focal);
	float _focal_length_to_half_focal(float focal_length);
	void _updateFrame();
	Point _pixelTo3D(int x, int y);
};

#endif //_CAMERA_H_
