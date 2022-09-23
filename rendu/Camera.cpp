
/*===============================================================*\
  Thursday June the 20th 2013
  Arash HABIBI
  Camera.c
\*===============================================================*/

#include "Camera.h"

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------

Camera::Camera()
{
	_pixel_width = 768;
	_pixel_height = 576;
	initPosit();
	initSize();
}

//---------------------------------------------------------

Camera::Camera(int width, int height)
{
	_pixel_width = width;
	_pixel_height = height;

	initPosit();
	initSize();
}

//---------------------------------------------------------

void Camera::initPosit()
{
	_half_focal = 30;
	_focal_length = _half_focal_to_focal_length(_half_focal);

	_dist_from_at = 10;
	_lookfrom.set     (0,0,_dist_from_at);
	_lookat.set       (0,0,0);
	_up_vector.set    (0,1,0);

	_aim = (_lookat-_lookfrom).unit();
	_right_vector = (_aim ^ _up_vector).unit();

	// turn the camera a bit

	Point px(1,0,0);
	Point py(0,1,0);
	_lookfrom = _lookfrom.turnAroundCam(_lookat,px,-0.5);
	_lookfrom = _lookfrom.turnAroundCam(_lookat,py,0.5);

	// recalculate aim and right

	_aim = (_lookat-_lookfrom).unit();
	_dist_from_at = (_lookat-_lookfrom).length();
	_right_vector = (_aim ^ _up_vector).unit();
}

//---------------------------------------------------------

void Camera::initSize()
{
	_ratio = 1.0*_pixel_width/_pixel_height;

	_3D_half_height = _dist_from_at * tan(M_PI*_half_focal/180);
	_3D_half_width = _3D_half_height * _ratio;

	_trackball_radius = sqrt(_3D_half_height*_3D_half_height + _3D_half_width*_3D_half_width);

	_znear = 0.1;
	_zfar  = 1000;
}

//---------------------------------------------------------

void Camera::check(char *message)
{
	cerr << "--------------- check camera " << message << " ------------------\n" << endl;
	cerr << "Resolution " << _pixel_width << " x " << _pixel_height << endl;
	_lookfrom.check("Position");
	_lookat.check("Aim");
	_up_vector.check("up vector");
	_right_vector.check("R vector");

	cerr << "Clipping planes : [" << _znear << "-" << _zfar << endl;

	cerr << "Focal angle : " << 2*_half_focal << endl;
	cerr << "Focal length : " << _focal_length << endl;
}

//---------------------------------------------------------

void Camera::setPixelDimensions(int width, int height)
{
	_pixel_width = width;
	_pixel_height = height;
	initSize();
}

//---------------------------------------------------------

Point Camera::lookfrom()
{
	return(_lookfrom);
}

//---------------------------------------------------------

Point Camera::lookat()
{
	return(_lookat);
}

//---------------------------------------------------------

Vector Camera::upVector()
{
	return(_up_vector);
}

//---------------------------------------------------------

float Camera::fovy()
{
	return(2*_half_focal);
}

//---------------------------------------------------------

float Camera::ratio()
{
	return _ratio;
}

//---------------------------------------------------------

float Camera::zNear()
{
	return(_znear);
}

//---------------------------------------------------------

float Camera::zFar()
{
	return(_zfar);
}

//---------------------------------------------------------

void Camera::zoom(float dtheta)
{
	float candidate = _half_focal + dtheta;
	if((candidate>0)&&(candidate<90))
	{
		_half_focal = candidate;
		_focal_length = _half_focal_to_focal_length(_half_focal);

		_3D_half_height = _dist_from_at * tan(M_PI*_half_focal/180);
		_3D_half_width = _3D_half_height * _ratio;
	}
}

//---------------------------------------------------------

void Camera::move(Vector translate)
{
	Vector aim, up, right, movelf, movela;
	double tx, ty, tz;
	translate.get(&tx,&ty,&tz);
	right = _right_vector * tx;
	up    = _up_vector * ty;
	aim   = _aim * tz;

	movelf = right + up;
	movela = movelf;
	movelf = movelf + aim;

	_lookfrom = _lookfrom + movelf;
	_lookat   = _lookat + movela;

	_updateFrame();
}

//---------------------------------------------------------

void Camera::turnAroundX(float dtheta)
{
	Point lookat_right, lookat_up;
	Vector u_y;

	u_y.base_y();
	lookat_right = _lookat + _right_vector;

	_lookfrom = _lookfrom.turnAround(lookat_right,_lookat,dtheta);
	_updateFrame();
	_right_vector = (_aim ^ u_y).unit();
	_up_vector = (_right_vector ^ _aim).unit();
}

//---------------------------------------------------------

void Camera::turnAroundY(float dtheta)
{
	Vector right, back, move;

	right = _right_vector * _dist_from_at*sin(M_PI*dtheta/180);
	back  = _aim          * _dist_from_at*(1-cos(M_PI*dtheta/180));
	move  = right + back;
	_lookfrom = _lookfrom + move;

	right= _right_vector * -(1-cos(M_PI*dtheta/180));
	back = _aim          * sin(M_PI*dtheta/180);
	move = right + back;
	_right_vector = _right_vector + move;

	_updateFrame();
}

//---------------------------------------------------------

void Camera::turnAroundZ(float dtheta)
{
}

//---------------------------------------------------------

void Camera::turnAroundWorldY(float dtheta)
{
	Point lookat_up;
	Vector u_y;

	u_y.base_y();
	lookat_up    = _lookat + u_y;

	_lookfrom = _lookfrom.turnAround(_lookat,lookat_up,dtheta);
	_updateFrame();
	_right_vector = (_aim ^ u_y).unit();
	_up_vector = (_right_vector ^ _aim).unit();
}

//---------------------------------------------------------
// previous_x, previous_y, x and y are in the [0,_pixel_width] or
// in the [0,_pixel_height] range.

void Camera::turnAroundTrackBall(int previous_x, int previous_y, int x, int y, double sensitivity)
{
	Point p_begin = _pixelTo3D(previous_x, previous_y);
	Point p_end   = _pixelTo3D(x,y);
	Vector v_begin_end = p_end - p_begin;
	p_end = p_begin + sensitivity*v_begin_end;

	// Project p_begin and p_end on the trackball
	Point proj_p_begin = p_begin.projectionOnTrackBall(_aim,_lookat,_trackball_radius);
	Point proj_p_end   = p_end  .projectionOnTrackBall(_aim,_lookat,_trackball_radius);

	Point p_right = _lookat + _right_vector;
	Point p_up    = _lookat + _up_vector;

	// begin and end are inversed because we are rotating the camera not the environment.
	Point lookfrom_turnaround = _lookfrom.turnAround(_lookat,proj_p_end,proj_p_begin);
	Point p_right_turnaround  = p_right.  turnAround(_lookat,proj_p_end,proj_p_begin);
	Point p_up_turnaround     = p_up.     turnAround(_lookat,proj_p_end,proj_p_begin);

	_lookfrom = lookfrom_turnaround;
	_aim = (_lookat - _lookfrom).unit();
	_right_vector = (p_right_turnaround - _lookat).unit();
	_up_vector    = (p_up_turnaround    - _lookat).unit();
}

//---------------------------------------------------------
// Put the lookat point on p_center, keep the same
// viewing direction and set lookfrom so that a sphere
// of radius radius placed at p_center is entirely visible
// and fills up the screen.

void Camera::frameBoundingSphere(Point p_center, double radius)
{
	Vector u_aim = (_lookat-_lookfrom).unit();

	double _half_focal_in_radians = M_PI*_half_focal/180;
	_lookat = p_center;
	_lookfrom = _lookat - radius*(1+1.0/tan(_half_focal_in_radians))*u_aim;
}

//---------------------------------------------------------
// Typically for the mouse position. (x,y)=(0,0) for the
// top left corner of the image.
// This function returns the position of the (x,y) pixel
// in the 3D space.

Point Camera::pixelPositionIn3D(int x, int y)
{
	Point res;

	_3D_half_height = (_lookfrom | _lookat ) * tan(M_PI*_half_focal/180);
	_3D_half_width = _3D_half_height * _ratio;
	Point top_left_corner = _lookat - (_3D_half_width*_right_vector) + (_3D_half_height*_up_vector);
	res = top_left_corner + (x*2.0*_3D_half_width/_pixel_width)*_right_vector - (y*2.0*_3D_half_height/_pixel_height)*_up_vector;
	return res;
}

//---------------------------------------------------------

void Camera::windowHalfSizeIn3D(float *half_width, float *half_height)
{
	*half_width  = _3D_half_width;
	*half_height = _3D_half_height;
}

//---------------------------------------------------------
// Typically for the mouse position. (x,y)=(0,0) for the
// top left corner of the image.
// This function projects the (x,y) pixel on the XY plane
// and returns the position of that projection.

Point Camera::pixelPositionInXYPlane(int x, int y)
{
	Point res;

	Point p = pixelPositionIn3D(x,y);
	Point origin(0,0,0);
	Vector normal(0,0,1);

	res = _lookfrom.linePlaneIntersect((p-_lookfrom).unit(), origin, normal);

	return res;
}

//---------------------------------------------------------

Point* Camera::getLookFromPointer()
{
	return &_lookfrom;
}

//---------------------------------------------------------

Point* Camera::getLookAtPointer()
{
	return &_lookat;
}


//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------

float Camera::_half_focal_to_focal_length(float half_focal)
{
	return(1.0/tan(M_PI*half_focal/180));
}

//---------------------------------------------------------

float Camera::_focal_length_to_half_focal(float focal_length)
{
	return((180.0/M_PI)*(1.0/focal_length));
}

//---------------------------------------------------------

void Camera::_updateFrame()
{
	_aim = (_lookat-_lookfrom).unit();
	_dist_from_at = (_lookat-_lookfrom).length();
}

//---------------------------------------------------------
// For a given pixel coordinate (x,y), this functon returns
// the corresponding point on the projection plane (the 3D
// plane orthogonal to _aim and containing _lookat)

Point Camera::_pixelTo3D(int x, int y)
{
	// ax and ay are in the [-1,1] range and
	// represent the position of the mouse click (x,y)
	// in the image 2D frame

	double ax = -1.0+2.0*x/_pixel_width;
	double ay =  1.0-2.0*y/_pixel_width;

	return _lookat + ax * _3D_half_width*_right_vector + ay * _3D_half_height*_up_vector;
}
