/*======================================================*\
  Monday November 11th 2012

  Arash Habibi

  Vector functions.
\*======================================================*/


#include "Vector.h"

static double _myCos(double x)
{	return(1-x*x);  }

static double _mySin(double x)
{	return(x);     }


//------------------------------------------------------------------
// Defaut : V_VECTOR

Vector::Vector()
{
	_data[0]=0;
	_data[1]=0;
	_data[2]=0;
	_data[3]=0;
}

//------------------------------------------------------------------

Vector::Vector(double x, double y, double z, double w)
{
	_data[0] = x;
	_data[1] = y;
	_data[2] = z;
	_data[3] = w;
}

//------------------------------------------------------------------
// DEFAULT : V_VECTOR

Vector::Vector(double x, double y, double z)
{
	_data[0] = x;
	_data[1] = y;
	_data[2] = z;
	_data[3] = 0;
}

//------------------------------------------------------------------

Vector::~Vector()
{
}

//------------------------------------------------------------------

void Vector::check(string message)
{
	int i;

	cerr << message << " : ";

	for(i=0;i<V_NB_COMP;i++)
		cerr << "\t" << _data[i];

	cerr << endl;
}

//------------------------------------------------------------------

void Vector::check(ofstream dump, string message)
{
	int i;

	dump << message << " : ";
	for(i=0;i<V_NB_COMP;i++)
		dump << "\t" << _data[i];
	dump << endl;
}

//------------------------------------------------------------------

void Vector::get(double *x, double *y, double *z)
{
	*x = _data[0];
	*y = _data[1];
	*z = _data[2];
}

//------------------------------------------------------------------

void Vector::get(double *x, double *y, double *z, double *w)
{
	get(x,y,z);
	if(V_NB_COMP > 3)
		*w = _data[3];
}

//------------------------------------------------------------------

void Vector::get(double *value, int comp)
{
	if(comp<V_NB_COMP)
		*value = _data[comp];
}

//------------------------------------------------------------------

double Vector::get(int comp)
{
	return(_data[comp]);
}

//------------------------------------------------------------------

void Vector::set(const double x_ref, const double y_ref, const double z_ref)
{
	_data[0] = x_ref;
	_data[1] = y_ref;
	_data[2] = z_ref;
}

//------------------------------------------------------------------

void Vector::set(const double x_ref, const double y_ref, const double z_ref, const double w_ref)
{
	set(x_ref,y_ref,z_ref);

	if(V_NB_COMP > 3)
		_data[3] = w_ref;
}

//------------------------------------------------------------------
// comp is either V_X, V_Y or V_Z

void Vector::set(int comp, const double value)
{
	if(comp<V_NB_COMP)
		_data[comp] = value;
}

//------------------------------------------------------------------
// From a given radius, this function sets the value of the vector
//

void Vector::setRandom(double radius)
{
    double dist, theta, phi;

    dist  = (rand()/32768.) * radius;
    theta = (rand()/32768.) * M_PI;
    phi   = (rand()/32768.) * M_PI * 2;

    _data[0] = dist * cos(theta);
    _data[1] = dist * sin(theta) * sin(phi);
    _data[2] = dist * sin(theta) * cos(phi);
}

//------------------------------------------------------------------

void Vector::base_x()
{
	_data[0] = 1;
	_data[1] = 0;
	_data[2] = 0;
	_data[3] = 0;
}

//------------------------------------------------------------------

void Vector::base_y()
{
	_data[0] = 0;
	_data[1] = 1;
	_data[2] = 0;
	_data[3] = 0;
}

//------------------------------------------------------------------

void Vector::base_z()
{
	_data[0] = 0;
	_data[1] = 0;
	_data[2] = 1;
	_data[3] = 0;
}

//------------------------------------------------------------------

void Vector::base_w()
{
	_data[0] = 0;
	_data[1] = 0;
	_data[2] = 0;
	_data[3] = 1;
}

//------------------------------------------------------------------

double Vector::length()
{
	return(sqrt(length2()));
}

//------------------------------------------------------------------

double Vector::length2()
{
	return(_data[0]*_data[0] + _data[1]*_data[1] + _data[2]*_data[2]);
}

//------------------------------------------------------------------

double Vector::maxXyz()
{
  double x,y,z,max;

  if(_data[0]>0) x = _data[0];   else    x = -_data[0];
  if(_data[1]>0) y = _data[1];   else    y = -_data[1];
  if(_data[2]>0) z = _data[2];   else    z = -_data[2];

  max = x;
  if(y>max) max = y;
  if(z>max) max = z;

  return(max);
}

//------------------------------------------------------------------

Vector Vector::unit()
{
	Vector unit(1,0,0,V_VECTOR);
	double len, OneOverLength;

	if((V_NB_COMP>3)&&(_data[3]!=0))
	{
		cerr << "Vector::unit : error ! the homogeneoous component must be 0 != " << _data[3] << endl;
		return(unit);
	}
	else
	{
		len = length();

		if(len!=0)
		{
			OneOverLength = 1.0/len;
			unit.set(_data[0] * OneOverLength,
					 _data[1] * OneOverLength,
					 _data[2] * OneOverLength);
		}
		/*else
		  cerr << "------------- Vector::unit : zero length vector" << endl;*/
		return(unit);
	}
}

//------------------------------------------------------------------

void Vector::normalize()
{
	Vector v = unit();
	(*this) = v;
}

//------------------------------------------------------------------
// Returns another vector (any other) not colinear to the current
// vector.

Vector Vector::shuffle()
{
    Vector v, candidate, res;
	Vector u_x(1,0,0);
	Vector u_y(0,1,0);
	Vector u_z(0,0,1);
	double l_cross_prod, max=0.0;

	v = (*this);

	int i;
	for(i=0;i<3;i++)
	{
		candidate = v;
		candidate._data[i]++;
		l_cross_prod = (v ^ candidate).length();
		if(l_cross_prod > max)
		{
			max = l_cross_prod;
			res = candidate;
		}
	}
    return(res);
}

//------------------------------------------------------------------

bool Vector::isZero()
{
	int i;
	bool is_zero = true;

	for(i=0;(i<3)&&(is_zero);i++)
		// if(_data[i]!=0)
		if(fabs(_data[i]) > 1e-10)
			is_zero = false;

	return(is_zero);
}

//------------------------------------------------------------------
// Returns true if the vector is parallele to one of the axes and
// false otherwise.

bool Vector::orthotrope()
{
	int nb_none_zero_components=0;
	int i;
	for(i=0;(i<3)&&(nb_none_zero_components<=1);i++)
		if(_data[i]!=0)
			nb_none_zero_components++;

	return(nb_none_zero_components==1);
}

//------------------------------------------------------------------

int Vector::closestAxis()
{
	double xx = fabs(_data[0]);
	double yy = fabs(_data[1]);
	double zz = fabs(_data[2]);
	if(xx>yy)
	{
		if(zz>xx)	return V_Z;
		else		return V_X;
	}
	else
	{
		if(zz>yy)	return V_Z;
		else		return V_Y;
	}
}

//------------------------------------------------------------------
// filmFile is supposed to be opened in binary mode.

void Vector::writeVect(ofstream filmFile)
{
	int i;
	for(i=0;i<3;i++)
		filmFile << _data[i];
}

//------------------------------------------------------------------

double Vector::operator[](int comp)   // Read component
{
	return(_data[comp]);
}

//------------------------------------------------------------------

Vector Vector::operator=(const Vector right)  // Assignment
{
  int i;
  for(i=0;i<V_NB_COMP;i++)
    _data[i] = right._data[i];

  return(*this);
}

//------------------------------------------------------------------

bool Vector::operator==(const Vector right)  // Equality in value
{
	bool are_equal=true;

	int i;
	for(i=0;(i<V_NB_COMP)&&(are_equal);i++)
		if(_data[i] != right._data[i])
			are_equal = false;

	return(are_equal);
}

//------------------------------------------------------------------

bool Vector::operator!=(const Vector right)  // Inequality in value
{
	bool are_equal = ((*this)==right);
	return(!are_equal);
}

//------------------------------------------------------------------

Vector Vector::operator+(const Vector right)  // Vector addition
{
  Vector v_sum;

  int i;
  for(i=0;i<V_NB_COMP;i++)
	  v_sum._data[i] = _data[i] + right._data[i];

  return(v_sum);
}

//------------------------------------------------------------------

Vector Vector::operator-(const Vector right)  // Vector substraction
{
  Vector v_diff;

  int i;
  for(i=0;i<V_NB_COMP;i++)
	  v_diff._data[i] = _data[i] - right._data[i];

  return(v_diff);
}

//------------------------------------------------------------------

double Vector::operator*(const Vector right)    // Dot product
{
	double dot_product=0;

	if((V_NB_COMP>3)&&(_data[3]!=0))
	{
		cerr << "Vector::operator* : error ! the homogeneoous component must be 0 != " << _data[3] << endl;
		return(0);
	}
	else
	{
		int i;
		for(i=0;i<V_NB_COMP;i++)
			dot_product += _data[i] * right._data[i];

		return(dot_product);
	}
}

//------------------------------------------------------------------

Vector Vector::operator*(const double scale)     // Vector-Scalar product
{
	Vector scaled_vector;

	int i;
	for(i=0;i<V_NB_COMP;i++)
		scaled_vector._data[i] = scale * _data[i];

	// scaled_vector._data[3] = 0;
	// This line must be uncommented
	// In the beginning this function used to ignore the homogeneous coordinate.
	// Indeed this operations concerns mainly vectors (not points).
	// In this situation multiplying the homogeneous coordinate or not has no
	// importance. There is at least one situation where this function operates
	// on points and where multiplying the homogeneous coordinate is important
	// This case is the computation of a middle point.

	return(scaled_vector);
}

//------------------------------------------------------------------

Vector Vector::operator^(const Vector right)  // Cross product
{
  Vector v_prod;

  v_prod._data[0] =   _data[1] * right._data[2]  -  _data[2] * right._data[1];
  v_prod._data[1] = - _data[0] * right._data[2]  +  _data[2] * right._data[0];
  v_prod._data[2] =   _data[0] * right._data[1]  -  _data[1] * right._data[0];
  v_prod._data[3] = 0;

  return(v_prod);

}

//------------------------------------------------------------------

void Vector::getPointer(double **x, double **y, double **z)
{
	*x = _data;
	*y = _data+1;
	*z = _data+2;
}

//------------------------------------------------------------------
// So that the scalar-vector product can be commutative.

Vector operator * (const double k, Vector V)
{
	return (V*k);
}


//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Defaut : V_POINT

Point::Point()
{
	_data[0]=0;
	_data[1]=0;
	_data[2]=0;
	_data[3]=1;
}

//------------------------------------------------------------------

Point::Point(double x, double y, double z, double w)
{
	_data[0] = x;
	_data[1] = y;
	_data[2] = z;
	_data[3] = w;
}

//------------------------------------------------------------------
// DEFAULT : V_POINT

Point::Point(double x, double y, double z)
{
	_data[0] = x;
	_data[1] = y;
	_data[2] = z;
	_data[3] = 1.0;
}

//------------------------------------------------------------------

Point::~Point()
{
}

//------------------------------------------------------------------

Point Point::operator=(Vector right)  // assignment
{
  int i;
  for(i=0;i<V_NB_COMP;i++)
	  _data[i] = right.get(i);

  return(*this);
}

//------------------------------------------------------------------

Point Point::operator=(Point right)  // assignment
{
  int i;
  for(i=0;i<V_NB_COMP;i++)
    _data[i] = right._data[i];

  return(*this);
}

//------------------------------------------------------------------

Point Point::operator+(Vector right)  // point + vector = point
{
	Point p_sum;

	int i;
	for(i=0;i<V_NB_COMP;i++)
		p_sum._data[i] = _data[i] + right.get(i);

	return(p_sum);
}

//------------------------------------------------------------------

Point Point::operator-(Vector right)  // point - vector = point
{
	Point p_diff;

	int i;
	for(i=0;i<V_NB_COMP;i++)
		p_diff._data[i] = _data[i] - right.get(i);

	return(p_diff);
}

//------------------------------------------------------------------

Vector Point::operator-(Point right)  // point - point = vector
{
	Vector v_diff;

	int i;
	for(i=0;i<V_NB_COMP;i++)
		v_diff.set(i,_data[i] - right.get(i));

	return(v_diff);
}

//------------------------------------------------------------------
// is relevant in particular for averaging positions.

Point Point::operator*(double right)   // point * scalar = point
{
	Point scaled_point;

	int i;
	for(i=0;i<V_NB_COMP;i++)
		scaled_point._data[i] = right * _data[i];

	return(scaled_point);
}

//------------------------------------------------------------------

double Point::operator|(Point right)  // distance
{
	double distance2=0;

	int i;

	for(i=0;i<3;i++)
		distance2 +=
			(right._data[i]-_data[i]) *
			(right._data[i]-_data[i]);

	return(sqrt(distance2));
}

//------------------------------------------------------------------

double Point::operator|| (Point right)  // distance square
{
	double distance2=0;

	int i;

	for(i=0;i<3;i++)
		distance2 +=
			(right._data[i]-_data[i]) *
			(right._data[i]-_data[i]);

	return(distance2);
}

//------------------------------------------------------------------

double Point::operator & (Point right)  // max coordinate distance
{
	return ((*this) - right).maxXyz();
}

//------------------------------------------------------------------
// So that the scalar-point product can be commutative.

Point operator * (const double k, Point P)
{
	return (P*k);
}

//------------------------------------------------------------------
// p1 % p2 returns a random vector in a block whose opposite vertices
// are p1 and p2.

Point Point::operator%(Point right)  // random
{
	double delta, min, max;
    Point res;

	if(_data[3] != right._data[3])
	{
		cerr << "Point::operator% : ERROR. the operands must have the same homogeneous coordinate." << endl;
		exit(1);
	}
	else
	{
		int i;
		for(i=0;i<3;i++)
		{
			if(_data[i] < right._data[i])
			{	min = _data[i];	max = right._data[i];	}
			else
			{	min = right._data[i];	max = _data[i];	}

			delta = max - min;
			res._data[i] = min + delta*drand48();
		}
	}
	res._data[3] = _data[3];
	return(res);
}





//------------------------------------------------------------------
// Precondition : A and B are disjoint and the current point M is
// supposed to lie on the (AB) line. This function returns alpha
// such that current point M
// M = (1-alpha)*A + alpha*B;
// M = A + alpha*(B-A)

double Point::blendingFactor(Point A, Point B)
{
	Point M = *this;
	Vector v = B - A;
	if(v.isZero())
	{
		cerr << "Point::blendingFactor : A and B are the same point." << endl;
		A.check("A");
		B.check("B");
		return -1;
	}
	else
	{
		int comp = v.closestAxis();
		return (M.get(comp) - A.get(comp)) / v.get(comp);
	}
}

//------------------------------------------------------------------
//  Returns the distance d from the current point (M) to
//  the plane defined by point P and normal vector normal
//  if n points towards A and -d otherwise.
//
//  CAUTION !!! The result is not really a distance since
//  it can be negative.

double Point::distToPlane(Point P, Vector normal)
{
	Vector MP = P - (*this);
	return(MP * (normal.unit()));
}

//------------------------------------------------------------------
// Same thing as before but the plane is defined by three points
// instead of one point and a normal vector.

double Point::distToPlane(Point A, Point B, Point C)
{
  Vector MA, AB, AC;
  Vector normal;

  MA = A - *this;
  AB = B - A;
  AC = C - A;

  normal = AB ^ AC;

  return(distToPlane(A,normal));
}

//------------------------------------------------------------------
// Returns the distance of the current point M to the (AB) line.

double Point::distToLine(Point A, Point B)
{
	Vector MH = orthoProjectVectorOnLine(A,B);
	return(MH.length());
}

//------------------------------------------------------------------

double Point::distToLineSegment(Point A, Point B)
{
	Point H = closestPointOnLineSegment(A,B);
	return (H|(*this));
}

//------------------------------------------------------------------

double Point::distToCircle(Point C, double radius, Vector normal)
{
    Point p_closest = closestPointOnCircle(C,radius,normal);
	return ((*this)|p_closest);
}


//------------------------------------------------------------------

Point Point::orthoProjectionOnPlane(Point P, Vector normal)
{
	Vector MP = P - (*this);
	normal = normal.unit();
	double h = MP * normal;
	return((*this) + h*normal);
}

//------------------------------------------------------------------

Point Point::orthoProjectionOnPlane(Point A, Point B, Point C)
{
  Vector MA, AB, AC;
  Vector normal;

  MA = A - *this;
  AB = B - A;
  AC = C - A;

  normal = AB ^ AC;

  return(orthoProjectionOnPlane(A,normal));
}

//------------------------------------------------------------------
// Let M be the current point. M can be projected on the (AB) line.
// Let H be the projection point. This function returns the MH vector.

Vector Point::orthoProjectVectorOnLine(Point A, Point B)
{
  Vector MA, AB;
  Vector u_projection;

  MA = A - *this;
  AB = B - A;

  u_projection = (AB ^ (MA ^ AB));
  if(u_projection.length()==0)
	  return u_projection;
  else
  {
	  u_projection.normalize();
	  return(u_projection * (MA * u_projection));
  }
}

//------------------------------------------------------------------
// This function projects the current vector on the (AB) line.

Point Point::orthoProjectionOnLine(Point A, Point B)
{
	Vector MH = orthoProjectVectorOnLine(A,B);
	return((*this)+MH);
}

//------------------------------------------------------------------
// This function projects the current point on the sphere of center
// C and radius R along vector u_ray.

Point Point::projectionOnSphere(Vector u_ray, Point C, double radius)
{
	Point M = (*this);
	Vector CH = C.orthoProjectVectorOnLine(M,M+u_ray);
	double d = CH.length();
	Point H = C + CH;
	if(d<=radius)
	{
		double h = sqrt( radius*radius - d*d );
		return(H - u_ray*h);
	}
	else
	{
		fprintf(stderr,"Point::projectionOnSphere : no projection on sphere");
		return(*this);
	}
}

//------------------------------------------------------------------
// Same function as above but when there is no intersection, the
// return position is the position of closest point on the sphere
// to the (M,u_ray) ray.

Point Point::projectionOnTrackBall(Vector u_ray, Point C, double radius)
{
	Point M = (*this);
	Vector CH = C.orthoProjectVectorOnLine(M,M+u_ray);
	double d = CH.length();
	Point H = C + CH;

	if(d>radius)  d = radius;
	double h = sqrt( radius*radius - d*d );
	return(H - u_ray*h);
}

//------------------------------------------------------------------
// In the following set of functions, "projectionOnPlane" is
// viewed in opposition to "orthoProjectionOnPlane" written above.
// Unlike those orthogonal projections, the one described here
// is a projection on plane containing point P and orthogonal to
// vector normal and along the u_ray direction.

bool Point::projectionOnPlaneExists(Vector u_ray, Point P, Vector normal)
{
	return( u_ray * (((P-(*this))*normal)*normal) > 0);
}

//------------------------------------------------------------------

bool Point::projectionOnPlaneExists(Vector u_ray, Point A, Point B, Point C)
{
  Vector MA, AB, AC;
  Vector normal;

  MA = A - *this;
  AB = B - A;
  AC = C - A;

  normal = (AB ^ AC);

  return(projectionOnPlaneExists(u_ray,A,normal));
}

//------------------------------------------------------------------
// Comment : let H be the orthogonal lprojection of the current
// point M on the (P,normal) plane. Let X be the sought projection point.
// This function returns the intersection point between the plane
// containing point P and with normal vector normal and the
// ray stemming from the current point and with direction u_ray.

Point Point::projectOnPlane(Vector u_ray, Point P, Vector normal)
{
	Point X, M=(*this);

	u_ray.normalize();
	normal.normalize();

	Point H = orthoProjectionOnPlane(P,normal);

	if( (H-M) * u_ray <= 0) // no intersection
		X = M;
	else
	{
		double l_MH = (H-M) * normal;

		Vector u_tangent = ( normal ^ ( u_ray ^ normal) ).unit();
		double tan_XH_XM = (u_ray*normal) / (u_ray*u_tangent);

		double l_HX = l_MH / tan_XH_XM;
		Vector HX = l_HX * u_tangent;
		X = H + HX;
	}
	return X;
}

//------------------------------------------------------------------
// Same function as above, but the plane is not defined by P and n
// but by three points A, B and C.

Point Point::projectOnPlane(Vector u, Point A, Point B, Point C)
{
  Vector MA, AB, AC;
  Vector normal;

  MA = A - *this;
  AB = B - A;
  AC = C - A;

  normal = (AB ^ AC);

  return(projectOnPlane(u,A,normal));
}

//------------------------------------------------------------------
// Same thing as the previous function, but instead of calculating
// the intersection between a plane and a RAY, it produces the
// intersection between a plane and a whole line. If there is no
// intersection, then the current point is returned.

Point Point::linePlaneIntersect(Vector u_ray, Point P, Vector normal)
{
	Point X, M=(*this);

	u_ray.normalize();
	normal.normalize();

	Point H = orthoProjectionOnPlane(P,normal);

	if( (H-M) * u_ray == 0) // no intersection
		X = M;
	else
	{
		double l_MH = (H-M) * normal;

		Vector u_tangent = ( normal ^ ( u_ray ^ normal) ).unit();
		double tan_XH_XM = (u_ray*normal) / (u_ray*u_tangent);

		double l_HX = l_MH / tan_XH_XM;
		Vector HX = l_HX * u_tangent;
		X = H + HX;
	}
	return X;
}

//------------------------------------------------------------------

int Point::lineSphereIntersect(Vector u_ray, Point C, double radius, Point *P1, Point *P2)
{
	Point M = (*this);
	int nb_intersections=0;
	Vector CH = C.orthoProjectVectorOnLine(M,M+u_ray);
	Point H = C + CH;

	double d = CH.length();
	if(d==radius)
	{
		*P1 = H;
		*P2 = H;
		nb_intersections = 1;
	}
	if(d<=radius)
	{
		double h = sqrt(radius*radius - d*d);
		*P1 = H - u_ray*h;
		*P2 = H + u_ray*h;
		nb_intersections = 2;
	}
	else
		nb_intersections = 0;

	return nb_intersections;
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------

Point Point::closestPointOnLineSegment(Point A, Point B)
{
	Point M = (*this);
	Point result(0,0,0);

	Vector AM = M - A;
	Vector BM = M - B;
	Vector AB = B - A;

	AM.set(V_W,0);
	BM.set(V_W,0);

	if(M.isOnTheFirstEndSide(A,B))
		result = A;

	else if(M.isOnTheSecondEndSide(A,B))
		result = B;

	else
	{
		Vector MH = M.orthoProjectVectorOnLine(A,B);
		Point H = M + MH;
		result = H;
	}

	return result;
}

//------------------------------------------------------------------

Point Point::closestPointOnCircle(Point C, double radius, Vector normal)
{
	Vector CM, v_right, u_closest;
	Point M, p_closest;

	M = *this;
	CM = M - C;
	if(CM * normal < 0)
		normal = (-1) * normal;

	v_right = CM ^ normal;
	u_closest = (normal ^ v_right).unit();
	p_closest = C + (radius * u_closest);
	return ( p_closest );
}

//------------------------------------------------------------------
// This function returns the result of the rotation of current point
// around the (AB) axis about an angle of angle degrees.

Point Point::turnAround(Point A, Point B, double angle)
{
	angle = M_PI * angle / 180;
	Point H = orthoProjectionOnLine(A,B);
	Vector u_x = ((*this)-H).unit();
	Vector u_z = (B-A).unit();
	Vector u_y = u_z ^ u_x;
	double R = ((*this)-H).length();
	// Vector res = H + (u_x*cos(angle) + u_y*sin(angle))*R;
	Point res = H + (u_x*_myCos(angle) + u_y*_mySin(angle))*R;
	return(res);
}

//------------------------------------------------------------------
// Same thing but the object is assumed to be on the same plane
// as A.

Point Point::turnAroundCam(Point A, Point B, double angle)
{
	Vector u_x = ((*this)-A).unit();
	Vector u_z = (B-A).unit();
	Vector u_y = u_z ^ u_x;
	double R = ((*this)-A).length();
	Point res = A + (u_x*_myCos(angle) + u_y*_mySin(angle))*R;
	return(res);
}

//------------------------------------------------------------------
// This function calculates the rotation that brings OA to OB,
// applies this rotation to the current point and returns the
// resulting point. If OA and OB are colinear, then there is no
// rotation.

Point Point::turnAround(Point O, Point A, Point B)
{
	Vector OA = A - O;
	Vector OB = B - O;
	Vector OAOB = OA ^ OB;
	if(OAOB.length()==0)
		return(*this);
	else
	{
		Vector u_axis = OAOB.unit();
		double sin_angle = OAOB.length() / (OA.length()*OB.length());
		double cos_angle = sqrt(1-sin_angle*sin_angle);

		Point H = this->orthoProjectionOnLine(O,O+u_axis);
		Vector HM = *this-H;

		Vector u_x = HM.unit();
		Vector u_z = u_axis;
		Vector u_y = u_z ^ u_x;

		double R = HM.length();
		Point res = H + (u_x*cos_angle + u_y*sin_angle)*R;
		return(res);
	}
}

//------------------------------------------------------------------

bool Point::isOnTheSamePlaneAs(Point A, Point B, Point C)
{
	Vector AB = B - A;
	Vector AC = C - A;
	Vector normal_ABC_plane = AB^AC;
	Vector AM = *this - A;
	return ( AM * normal_ABC_plane == 0 );
}

//------------------------------------------------------------------

void Point::placeAtRandomInTriangle(Point A, Point B, Point C)
{
	double u = drand48();
	double v = drand48();
	if( u+v > 1)
	{
		double tmp = u;
		u = 1-v;
		v = 1-tmp;
	}
	(*this) = A + u*(B-A) + v*(C-A);
}

//------------------------------------------------------------------
// calculates the area of the ABC triangle. The result is always
// positive and does not depend on the ordre of the vertices.

double positiveTriangleArea(Point A, Point B, Point C)
{
	return( 0.5 * ((B-A)^(C-A)).length() );
}

//------------------------------------------------------------------
// Computes the area of the ABC triangle. This result is positive
// if the ABC triangle is oriented in the same direction as normal.
// otherwise it is negative. Thus the result depends greatly on
// the order or the vertices.

double triangleArea(Point A, Point B, Point C, Vector normal)
{
	normal.normalize();
	return (0.5 * (((B-A)^(C-A)) * normal)) ;
}

//------------------------------------------------------------------
// Computes the volume of the ABCD tetrahedron. The result is
// positive only if D is on the same side of ABC as the AB^AC
// cross product.

double tetrahedronVolume(Point A, Point B, Point C, Point D)
{
	return (((B-D)^(C-D)) * (A-D)) / 6.0;
}

//------------------------------------------------------------------

bool Point::isInBoundingBox(Point bbmin, Point bbmax)
{
	return
		(_data[0]>=bbmin._data[0]) && (_data[0]<=bbmax._data[0]) &&
		(_data[1]>=bbmin._data[1]) && (_data[1]<=bbmax._data[1]) &&
		(_data[2]>=bbmin._data[2]) && (_data[2]<=bbmax._data[2]);
}

//------------------------------------------------------------------
// This function returns the maximum of the following numbers :
// fabs(x-xmin), fabs(y-ymin), fabs(z-zmin),
// fabs(x-xmax), fabs(y-ymax), fabs(z-zmax)

double Point::maxAxisDistFromBoundingBox(Point bbmin, Point bbmax)
{
	double dist1, dist2, dx, dy, dz, dist_max = 0.0;

	dist1 = fabs(_data[0]-bbmin._data[0]);
	dist2 = fabs(_data[1]-bbmin._data[1]);
	if(dist1<dist2) dx = dist1;
	else 			dx = dist2;

	dist1 = fabs(_data[2]-bbmin._data[2]);
	dist2 = fabs(_data[0]-bbmax._data[0]);
	if(dist1<dist2) dy = dist1;
	else 			dy = dist2;

	dist1 = fabs(_data[1]-bbmax._data[1]);
	dist2 = fabs(_data[2]-bbmax._data[2]);
	if(dist1<dist2) dz = dist1;
	else 			dz = dist2;

	if(dx > dy)
		if(dz > dx) dist_max = dz;
		else		dist_max = dx;
	else
		if(dz > dy) dist_max = dz;
		else		dist_max = dy;

	return dist_max;
}

//------------------------------------------------------------------

double Point::maxAxisDistFromBoundingSphere(Point center, double radius)
{
	double dist1, dist2, dx, dy, dz, dist_max = 0.0;

	dist1 = fabs(_data[0]-center._data[0]-radius);
	dist2 = fabs(_data[0]-center._data[0]+radius);
	if(dist1 < dist2)	dx = dist1;
	else 				dx = dist2;

	dist1 = fabs(_data[1]-center._data[1]-radius);
	dist2 = fabs(_data[1]-center._data[1]+radius);
	if(dist1 < dist2)	dy = dist1;
	else 				dy = dist2;

	dist1 = fabs(_data[2]-center._data[2]-radius);
	dist2 = fabs(_data[2]-center._data[2]+radius);
	if(dist1 < dist2)	dz = dist1;
	else 				dz = dist2;

	if(dx > dy)
		if(dz > dx) dist_max = dz;
		else		dist_max = dx;
	else
		if(dz > dy) dist_max = dz;
		else		dist_max = dy;

	return dist_max;
}

//------------------------------------------------------------------
// If an observer is looking at the scene from up, and the current
// point is on the right of line (AB) then the function returns
// true. Otherwise it returns false;

bool Point::isOnTheRight(Point A, Point B, Vector up)
{
	Point M = (*this);
	Vector AM = M - A;
	Vector AB = B - A;
	return ( (AM ^ AB) * up > 0);
}

//------------------------------------------------------------------

bool Point::isOnTheFirstEndSide(Point A, Point B)
{
	Point M = (*this);
	Vector MA = A - M;
	Vector AB = B - A;
	return ( MA * AB > 0 );
}

//------------------------------------------------------------------

bool Point::isOnTheSecondEndSide(Point A, Point B)
{
	Point M = (*this);
	Vector MB = B - M;
	Vector BA = A - B;
	return ( MB * BA > 0 );
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Let there be two parallel lines (p1p2) and (q1q2). If these
// lines are not parallel, then the function exits and returns 0.
// If they are parallel it finds for each of the 4 points,
// another point on the other line which is facing it.
// qA = q1 + t_p1_onq*(q2-q1);
// qB = q1 + t_p2_onq*(q2-q1);
// the pair on p1p2 is characterised by :
// pA = p1 + t_q1_onp*(p2-p1);
// pB = p1 + t_q2_onp*(p2-p1);

extern void V_oppositePointsOnParallelLines(Point p1, Point p2, Point q1, Point q2,
											double *t_p1_onq, double *t_p2_onq, double *t_q1_onp, double *t_q2_onp)
{
	Vector u_p = p2-p1;
	Vector u_q = q2-q1;
	if((u_p^u_q).isZero())
	{
		Point p1_proj = p1.orthoProjectionOnLine(q1,q2);
		Point p2_proj = p2.orthoProjectionOnLine(q1,q2);
		Point q1_proj = q1.orthoProjectionOnLine(p1,p2);
		Point q2_proj = q2.orthoProjectionOnLine(p1,p2);

		if((q1||q2)==0) { *t_p1_onq = 2; *t_p2_onq = 2; }
		else
		{
			*t_p1_onq = p1_proj.blendingFactor(q1,q2);
			*t_p2_onq = p2_proj.blendingFactor(q1,q2);
		}
		if((p1||p2)==0) { *t_q1_onp = 2; *t_q2_onp = 2; }
		else
		{
			*t_q1_onp = q1_proj.blendingFactor(p1,p2);
			*t_q2_onp = q2_proj.blendingFactor(p1,p2);
		}
	}
}

//------------------------------------------------------------------
// Same thing as above but for line segments [p1p2] and [q1q2].
// If the segments are parallel it attempts to find two pairs of points
// on each segment, as far appart as possible, and facing the pair
// on the other segment. If there can be no points facing each
// other, the function exits and returns 1.
// Otherwise, the pair on q1q2 is characterised by :
// qA = q1 + t_p1_onq*(q2-q1);
// qB = q1 + t_p2_onq*(q2-q1);
// the pair on p1p2 is characterised by :
// pA = p1 + t_q1_onp*(p2-p1);
// pB = p1 + t_q2_onp*(p2-p1);

extern int V_oppositePointsOnParallelSegments(Point p1, Point p2, Point q1, Point q2,
											   double *t_p1_onq, double *t_p2_onq, double *t_q1_onp, double *t_q2_onp)
{	int res;
	Vector u_p = p2-p1;
	Vector u_q = q2-q1;
	if((u_p^u_q).isZero())
	{
		Point p1_proj = p1.orthoProjectionOnLine(q1,q2);
		Point p2_proj = p2.orthoProjectionOnLine(q1,q2);
		Point q1_proj = q1.orthoProjectionOnLine(p1,p2);
		Point q2_proj = q2.orthoProjectionOnLine(p1,p2);

		double xp1_onq1q2, xp2_onq1q2, xq1_onp1p2, xq2_onp1p2;
		if((q1||q2)==0) { xp1_onq1q2 = 2; xp2_onq1q2 = 2; }
		else
		{
			xp1_onq1q2 = p1_proj.blendingFactor(q1,q2);
			xp2_onq1q2 = p2_proj.blendingFactor(q1,q2);
		}
		if((p1||p2)==0) { xq1_onp1p2 = 2; xq2_onp1p2 = 2; }
		else
		{
			xq1_onp1p2 = q1_proj.blendingFactor(p1,p2);
			xq2_onp1p2 = q2_proj.blendingFactor(p1,p2);
		}

		// if q1 and q2 are both on the p1 side
		if( (xq1_onp1p2<=0) && (xq2_onp1p2<=0) )
		{
			*t_q1_onp = 0;	*t_q2_onp = 0;

			if(xq1_onp1p2<=xq2_onp1p2) // q1 is closer to p1p2
			{	*t_p1_onq = 1;	*t_p2_onq = 1;	}
			else                      // q2 is closer to p1p2
			{	*t_p1_onq = 0;	*t_p2_onq = 0;	}
			res=1;
		}

		// if q1 and q2 are both on the p2 side
		else if( (xq1_onp1p2>=1) && (xq2_onp1p2>=1) )
		{
			*t_q1_onp = 1;	*t_q2_onp = 1;

			if(xq1_onp1p2<=xq2_onp1p2) // q1 is closer to p1p2
			{	*t_p1_onq = 0;	*t_p2_onq = 0;	}
			else                      // q2 is closer to p1p2
			{	*t_p1_onq = 1;	*t_p2_onq = 1;	}
			res=1;
		}

		// if p1p2 is within the q1q2 range
		else if( ( (xq1_onp1p2<=0) && (xq2_onp1p2>=1) ) ||
				 ( (xq1_onp1p2>=1) && (xq2_onp1p2<=0) ) )
		{
			*t_p1_onq = xp1_onq1q2;	*t_p2_onq = xp2_onq1q2;

			if ( (xq1_onp1p2<=0) && (xq2_onp1p2>=1) )
			{	*t_q1_onp = 0;	*t_q2_onp = 1;	}
			else // if( (xq1_onp1p2>=1) && (xq2_onp1p2<=0) )
			{	*t_q1_onp = 1;	*t_q2_onp = 0;	}
			res=2;
		}

		// if q1q2 is within the p1p2 range
		else if( ( (xp1_onq1q2<=0) && (xp2_onq1q2>=1) ) ||
				 ( (xp1_onq1q2>=1) && (xp2_onq1q2<=0) ) )
		{
			*t_q1_onp = xq1_onp1p2;	*t_q2_onp = xq2_onp1p2;

			if ( (xp1_onq1q2<=0) && (xp2_onq1q2>=1) )
			{	*t_p1_onq = 0;	*t_p2_onq = 1;	}
			else // if( (xp1_onq1q2>=1) && (xp2_onq1q2<=0) )
			{	*t_p1_onq = 1;	*t_p2_onq = 0;	}
			res=2;
		}

		// Only one of the vertices is in the other segment's range
		else if(xp1_onq1q2 <= 0)
		{
			*t_p1_onq = 0;
			*t_p2_onq = xp2_onq1q2;
			*t_q1_onp = xq1_onp1p2;
			*t_q2_onp = 1;
			res=2;
		}
		else if(xp1_onq1q2 >= 1)
		{
			*t_p1_onq = 1;
			*t_p2_onq = xp2_onq1q2;
			*t_q1_onp = 1;
			*t_q2_onp = xq2_onp1p2;
			res=2;
		}
		else if(xp2_onq1q2 <= 0)
		{
			*t_p1_onq = xp1_onq1q2;
			*t_p2_onq = 0;
			*t_q1_onp = xq1_onp1p2;
			*t_q2_onp = 0;
			res=2;
		}
		else if(xp2_onq1q2 >= 1)
		{
			*t_p1_onq = xp1_onq1q2;
			*t_p2_onq = 1;
			*t_q1_onp = 0;
			*t_q2_onp = xq2_onp1p2;
			res=2;
		}
		else // all equal
		{
			cerr << " ALL EQUAL" << endl;
			res=2;
		}
	}
	else
		res = 0;
	return res;
}

//------------------------------------------------------------------
// http://geomalgorithms.com/a07-_distance.html
// Let there be two lines (p1p2) and (q1q2) in space. This function
// finds on each line the closest point to the other line. This
// closest point is represented by the parameter t_p and t_q such that:
// p_inters = p1 + t_p*u_p    where u_p=p2-p1
// q_inters = q1 + t_q*u_q    where u_q=q2-q1

void V_lineLineClosestPoints(Point p1, Point p2, Point q1, Point q2, double *t_p, double *t_q)
{
	Vector u_p = p2-p1;
	Vector u_q = q2-q1;
	Vector q1p1 = p1 - q1;

	double a = u_p*u_p;
	double b = u_p*u_q;
	double c = u_q*u_q;
	double d = u_p*q1p1;
	double e = u_q*q1p1;

	if(a*c-b*b==0) // parallel segments
	{
		*t_p = 0;
		Point H = p1.orthoProjectionOnLine(q1,q2);
		*t_q = H.blendingFactor(q1,q2);
	}
	else
	{
		*t_p = (b*e - c*d)/(a*c-b*b);
		*t_q = (a*e - b*d)/(a*c-b*b);
	}
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------
// Distance between
// p1 + t_p*(p2-p1) and
// q1 + t_q*(q2-q1)


static double distance2BetweenTwoSegmenPoints(Point p1, Point p2, Point q1, Point q2, double t_p, double t_q)
{
	Point p = p1 + t_p*(p2-p1);
	Point q = q1 + t_q*(q2-q1);
	return ( p || q );
}

//------------------------------------------------------------------
// http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
// The [p1p2] and [q1q2] linesegments are defined by :
// p(t_p) = p1 + t_p*u_p where u_p = p2-p1;
// q(t_q) = q1 + t_q*u_q where u_q = q2-q1;
// We wish to find t_p and t_q which minimize the norm of
// w(t_p,t_q) = p(t_p) - q(t_q) =
//            = p1-q1 + t_p*u_p - t_q*u_q
//            = w0 + t_p*u_p - t_q*u_q
// In other words we wish to minimize
// WW(t_p,t_q) = (w0 + t_p*u_p - t_q*u_q) * (w0 + t_p*u_p - t_q*u)
// which is a quadratic function of t_p and t_q

void V_segmentSegmentClosestPoints(Point p1, Point p2, Point q1, Point q2, double *t_p, double *t_q)
{
	Vector u_p = p2 - p1;
	Vector u_q = q2 - q1;
	Vector w0 = p1-q1;
	Vector v_ortho_two_segments = u_p ^ u_q;

	//------------------

	if((u_p.isZero())&&(u_q.isZero())) // p1==p2 && q1==q2
	{	*t_p = 0; 	*t_q = 0;  }

	//------------------

	else if(u_p.isZero()) // p1==p2
	{
		*t_p = 0;
		Point q = p1.closestPointOnLineSegment(q1,q2);
		*t_q = q.blendingFactor(q1,q2);
	}

	//------------------

	else if(u_q.isZero()) // q1==q2
	{
		*t_q = 0;
		Point p = q1.closestPointOnLineSegment(p1,p2);
		*t_p = p.blendingFactor(p1,p2);
	}

	//------------------

	// p1 != p2 and q1 != q2
	else if(v_ortho_two_segments.isZero()) // parallel segments
	{
		double t_p1_onq, t_p2_onq, t_q1_onp, t_q2_onp;
		V_oppositePointsOnParallelSegments(p1, p2, q1, q2, &t_p1_onq, &t_p2_onq, &t_q1_onp, &t_q2_onp);
		*t_p = 0.5*(t_p1_onq+t_p2_onq);
		*t_q = 0.5*(t_q1_onp+t_q2_onp);
	}

	//------------------

	else // segments not parallel
	{
		double t_p_line, t_q_line;
		V_lineLineClosestPoints(p1, p2, q1, q2, &t_p_line, &t_q_line);

		// If the closest point between lines lie inside both
		// line segments, then we have found the answer.
		if( (t_p_line >= 0)&&(t_p_line <= 1)&&
			(t_q_line >= 0)&&(t_q_line <= 1) )
		{
			*t_p = t_p_line;
			*t_q = t_q_line;
		}

		// Otherwise, we study w(t_p,t_q) in the (t_p,t_q) plane
		// We look for the minimum of this function on the ([0,1],[0,1])
		// region of the (t_p,t_q) plane. If this minimum is not inside
		// the region (previous case) then it is on the border of this
		// region, i.e. on one of the edges. It has been shown that this
		// minimum is on a side "SEEN" by C(t_p_line,t_q_line).
		// This leads to a limited number of cases (one edge)
		// t_p=0, t_q=0, t_p=1, t_q=1 (or two edges cf further)
		// For each case, setting the first derivative of w to 0
		// gives the proper answer.

		else
		{
			double up2 = u_p * u_p;
			double uq2 = u_q * u_q;
			double upuq = u_p * u_q;
			double w0up = w0 * u_p;
			double w0uq = w0 * u_q;
			double tmp_tp, tmp_tq;
			double dist2, dist2_min = INFINITY;

			if(t_q_line >= 1)
			{
				tmp_tp = (upuq - w0up) / up2;
				tmp_tq = 1;
				if(tmp_tp<0) tmp_tp=0;
				if(tmp_tp>1) tmp_tp=1;
				dist2 = distance2BetweenTwoSegmenPoints(p1,p2,q1,q2,tmp_tp,tmp_tq);
				if(dist2 < dist2_min)
				{
					*t_p = tmp_tp;
					*t_q = tmp_tq;
					dist2_min = dist2;
				}
			}
			if(t_q_line <= 0 )
			{
				tmp_tp = -w0up / up2;
				tmp_tq = 0;
				if(tmp_tp<0) tmp_tp=0;
				if(tmp_tp>1) tmp_tp=1;
				dist2 = distance2BetweenTwoSegmenPoints(p1,p2,q1,q2,tmp_tp,tmp_tq);
				if(dist2 < dist2_min)
				{
					*t_p = tmp_tp;
					*t_q = tmp_tq;
					dist2_min = dist2;
				}
			}
			if(t_p_line >= 1)
			{
				tmp_tp = 1;
				tmp_tq = (upuq + w0uq) / uq2;
				if(tmp_tq<0) tmp_tq=0;
				if(tmp_tq>1) tmp_tq=1;
				dist2 = distance2BetweenTwoSegmenPoints(p1,p2,q1,q2,tmp_tp,tmp_tq);
				if(dist2 < dist2_min)
				{
					*t_p = tmp_tp;
					*t_q = tmp_tq;
					dist2_min = dist2;
				}
			}
			if(t_p_line <= 0)
			{
				tmp_tp = 0;
				tmp_tq = w0uq / uq2;
				if(tmp_tq<0) tmp_tq=0;
				if(tmp_tq>1) tmp_tq=1;
				dist2 = distance2BetweenTwoSegmenPoints(p1,p2,q1,q2,tmp_tp,tmp_tq);
				if(dist2 < dist2_min)
				{
					*t_p = tmp_tp;
					*t_q = tmp_tq;
					dist2_min = dist2;
				}
			}
		}
	}
}

//------------------------------------------------------------------
// Let there be a circle described by p_center, radius and normal.
// Let there be a list of nb_points points. This function fills
// points with nb_points points distributed on the circle.

extern void _circlePoints(Point p_center, double radius, Vector normal, int nb_points, Point *points)
{
	Vector u_x(1,0,0);
	Vector u_y(0,1,0);
	Vector u_z(0,0,1);

	Vector u_x_circle, u_z_circle;

	if((u_x^normal).isZero()==false)
		u_x_circle = (u_x^normal).unit();

	else if((u_y^normal).isZero()==false)
		u_x_circle = (u_y^normal).unit();

	else u_x_circle = (u_z^normal).unit();

	u_z_circle = (u_x_circle^normal).unit();

	double dtheta = 2*M_PI / nb_points;
	int n;
	for(n=0;n<nb_points;n++)
		points[n] = p_center
			+ (radius*cos(n*dtheta)*u_x_circle)
			+ (radius*sin(n*dtheta)*u_z_circle);
}

//------------------------------------------------------------------
// Let there be the (p1p2) line and the circle C of center p_center,
// radius radius and normal normal. Let us project p1 and p2 on the
// plane containing the circle. Let p1_proj and p2_proj be the
// the result of the projections. This function calculates the
// intersection points between the (p1_proj,p2_proj) line and
// the C circle.

int V_lineProjectionOnCircle(Point p1, Point p2, Point p_center, double radius, Vector normal, Point *pc1, Point *pc2)
{
	Point p1_proj = p1.orthoProjectionOnPlane(p_center, normal);
	Point p2_proj = p2.orthoProjectionOnPlane(p_center, normal);
	Vector u_p1_p2_proj = (p2_proj-p1_proj).unit();

	Vector v = p_center.orthoProjectVectorOnLine(p1_proj, p2_proj);
	double dist2 = v.length2();
	double R2 = radius*radius;

	if(dist2 <= R2)
	{
		Point H = p_center + v;
		double dt = sqrt(radius*radius - dist2);
		*pc1 = H + dt*u_p1_p2_proj;
		*pc2 = H - dt*u_p1_p2_proj;
		return 2;
	}
	else
		return 0;
}

//------------------------------------------------------------------

Point V_closestCircleStepToLine(Point p1, Point p2, Point p_center, double radius, Vector normal, Point p_start, double dtheta)
{
	double cos_dtheta = cos(dtheta);
	double sin_dtheta = sin(dtheta);

	Vector u_x = (p_start-p_center).unit();
	Vector u_y = (normal^u_x).unit();
	Point p_higher_theta, p_lower_theta;
	// p_start = p_center + radius*u_x;
	p_higher_theta = p_start + radius*(  sin_dtheta*u_y + (cos_dtheta-1)*u_x);
	p_lower_theta  = p_start + radius*(- sin_dtheta*u_y + (cos_dtheta-1)*u_x);
	double d_current = p_start.distToLine(p1,p2);
	double d_higher_theta = p_higher_theta.distToLine(p1,p2);
	double d_lower_theta   = p_lower_theta.distToLine(p1,p2);
	int sign;
	double d_next;
	Point p_next;
	if(d_higher_theta < d_lower_theta)	{ sign= 1; d_next = d_higher_theta; p_next = p_higher_theta; }
	else								{ sign=-1; d_next = d_lower_theta;  p_next = p_lower_theta;  }

	while(d_next < d_current)
	{
		p_start = p_next;
		d_current = d_next;

		u_x = (p_start-p_center).unit();
		u_y = (normal^u_x).unit();

		p_next = p_start + radius*(sign*sin_dtheta*u_y + (cos_dtheta-1)*u_x);
		d_next = p_next.distToLine(p1,p2);
	}

	return p_start;
}

//------------------------------------------------------------------
// If the (p1p2) line does not intersect the infinite cylinder
// corresponding to the circle, then the problem is simple.
// There is only one closest point between the line and the
// circle. This intersection always lies on a plane containing
// p_center and its normal. It also lies on a perpendicular
// to the line. So all we have to do is to find the plane
// containing such a perpendicular and to return 1 (the number
// of closest points).
// If the (p1p2) line DOES intersect the infinite cylinder, it
// becomes a much more complicated question which implies
// solving a 4th order polynomial equation. We don't want to do
// that. Instead, we consider a circle as a set of N line segments
// And we search for the closest point between the (p1p2) line
// and each of the small line segments. We should always find
// two cases where the closest point lies within the small line
// segment. We return those 4 points (two circle points pc1 and pc2)
// and two line points (pl1 and pl2) and return 2.
// The closest of the two are pc1 and pl1.

int V_lineCircleClosestPoints(Point p1, Point p2, Point p_center, double radius, Vector normal, Point *pc1, Point *pc2, Point *pl1, Point *pl2, double epsilon)
{
	int nb_closest;

	// Calculate the intersection of line with the infinite cylinder associated with the circle
	double t_on_line, t_on_cylinder_axis;
	V_lineLineClosestPoints(p1, p2, p_center, p_center+normal, &t_on_line, &t_on_cylinder_axis);
	Point p_on_line = p1 + t_on_line*(p2-p1);
	Point p_on_cylinder_axis = p_center + t_on_cylinder_axis*normal;
	double dist2_line_cylinder_axis = p_on_line || p_on_cylinder_axis;

	Vector v_right = normal^(p2-p1);

	// if line perpendicular to the circle plane.
	if(v_right.isZero())
	{
		*pl1 = p_center.orthoProjectionOnLine(p1,p2);
		*pl2 = *pl1;
		Vector u_r = (p_on_line-p_on_cylinder_axis).unit();
		*pc1 = p_center + radius*u_r;

		if(dist2_line_cylinder_axis > radius*radius) // only one closest point
		{	*pc2 = *pc1;	nb_closest = 1;	}
		else
		{	*pc2 = p_center - radius*u_r; nb_closest = 2; }
	}

	else // line non-perpendicular to the circle plane
	{
		// if no intersection between line and infinite cylinder (simple part)
		if((p_on_line || p_on_cylinder_axis) >= radius*radius )
		{
			Vector n_intersection_plane = (normal^v_right).unit();
			*pl1 = p1.linePlaneIntersect(p2-p1, p_center, n_intersection_plane);
			*pl2 = *pl1;
			*pc1 = pl1->closestPointOnCircle(p_center, radius, normal);
			*pc2 = *pc1;
			nb_closest = 1;
		}
		else // more complicated : implies several closest points
		{
			Point p_start1, p_start2;
			V_lineProjectionOnCircle(p1, p2, p_center, radius, normal, &p_start1, &p_start2);

			double dtheta = (p_start1 | p_start2) / (4*radius);

			do
			{
				*pc1 = V_closestCircleStepToLine(p1, p2, p_center, radius, normal, p_start1, dtheta);
				*pc2 = V_closestCircleStepToLine(p1, p2, p_center, radius, normal, p_start2, dtheta);
				dtheta /= 10;
			} while(dtheta > epsilon);

			if ( (*pc1 || *pc2) < epsilon )	nb_closest = 1;
			else							nb_closest = 2;

			Vector v;
			v = pc1->orthoProjectVectorOnLine(p1,p2);
			*pl1 = *pc1 + v;
			v = pc2->orthoProjectVectorOnLine(p1,p2);
			*pl2 = *pc2 + v;
		}
	}
	return nb_closest;
}

//------------------------------------------------------------------

int V_segmentCircleClosestPoints(Point p1, Point p2, Point p_center, double radius, Vector normal, Point *pc1, Point *pc2, Point *pl1, Point *pl2)
{
	int nb_closest = V_lineCircleClosestPoints(p1, p2, p_center, radius, normal, pc1, pc2, pl1, pl2, 1e-5);

	double xpl1 = pl1->blendingFactor(p1,p2);
	double xpl2 = pl2->blendingFactor(p1,p2);

	if((xpl1<0)||(xpl1>1))
	{
		if(xpl1<0) { xpl1=0; *pl1=p1; }
		if(xpl1>1) { xpl1=1; *pl1=p2; }
		*pc1 = pl1->closestPointOnCircle(p_center,radius,normal);
	}

	if((xpl2<0)||(xpl2>1))
	{
		if(xpl2<0) { xpl2=0; *pl2=p1; }
		if(xpl2>1) { xpl2=1; *pl2=p2; }
		*pc2 = pl2->closestPointOnCircle(p_center,radius,normal);
	}

	if(nb_closest==2)
		if(xpl1-xpl2 < 1e-5)
			nb_closest = 1;

	return nb_closest;
}

//------------------------------------------------------------------
// If there is only one intersection then P1 points to that point.

int V_segmentSphereIntersect(Point A, Point B, Point Center, double radius, Point *P1, Point *P2)
{
	int nb_intersections=0;
	Vector CenterH = Center.orthoProjectVectorOnLine(A,B);
	Point H = Center + CenterH;

	double d2 = CenterH.length2();
	double R2 = radius*radius;

	if( d2 > R2 )
		nb_intersections=0;

	else // d2 <= R2
	{
		Vector u_ray = (B-A).unit();
		double h = sqrt( R2 - d2);
		Point p1 = H - u_ray*h;
		Point p2 = H + u_ray*h;

		double alpha1 = p1.blendingFactor(A,B);
		double alpha2 = p2.blendingFactor(A,B);

		if((alpha1>=0)&&(alpha1<=1))
		{
			*P1 = p1;
			if((alpha2>=0)&&(alpha2<=1))
			{
				*P2 = p2;
				nb_intersections=2;
			}
			nb_intersections=1;
		}
		else if((alpha2>=0)&&(alpha2<=1))
		{
			nb_intersections=1;
			*P1 = p2;
		}
	}
	return nb_intersections;
}

//------------------------------------------------------------------

extern bool V_lineSegmentPlaneIntersect(Point p1, Point p2, Point P, Vector normal, double *alpha)
{
	Vector u_p1p2 = (p2-p1).unit();
	Point X = p1.linePlaneIntersect(u_p1p2, P, normal);
	*alpha = X.blendingFactor(p1,p2);

	if((*alpha>=0)&&(*alpha<=1))	return true;
	else							return false;
}

//------------------------------------------------------------------
// Let there be two circles of centers C1 and C2 and of radii R1 and R2.
// The distance between the centers is D=C1C2. This function returns
// the number of intersections between both circles (0, 1, 2 or inf).
// Both intersections X1 and X2 lie on a line perpendicular to C1C2.
// C1C2 and X1X2 intersect at H. At the return of the function r1
// points at C1H and r2 points at C2H. If C1=C2 and R1=R2 there are
// an infinite number of intersections. The returned value is -1.
// The intersections are found by resolving three equations :
// R1^2 = r1^2 + h^2
// R2^2 = r2^2 + h^2
// D = r1 + r2

extern int V_circleCircleIntersection(double R1, double R2, double D, double *r1, double *r2, double *h)
{
	int nb_intersects=-1;
	double Rmin, Rmax;
	double *rmin, *rmax;

	if(R2>R1)
	{
		Rmax=R2;		Rmin=R1;
		rmax=r2;		rmin=r1;
	}
	else
	{
		Rmax=R1;		Rmin=R2;
		rmax=r1;		rmin=r2;
	}

	double D2 = D*D;
	double Rmax2 = Rmax*Rmax;
	double Rmin2 = Rmin*Rmin;

	if((D==0)&&(Rmax==Rmin))
		nb_intersects = -1;

	else if(D>Rmax+Rmin) // disjoint circles
		nb_intersects=0;

	else if(D<Rmax-Rmin) // one circle is inside the other.
		nb_intersects=0;

	else if(D2 > Rmax2 - Rmin2)
	{
		*rmax = (D2+Rmax2-Rmin2)/(2*D);
		*rmin = (D2-Rmax2+Rmin2)/(2*D);
		*h = Rmax2 - (*rmax)*(*rmax);
		nb_intersects = 2;
		if(D2==Rmax+Rmin) 	nb_intersects = 1;
	}
	else // if(D2 <= Rmax2 - Rmin2)
	{
		*rmax = (Rmax2-Rmin2+D2)/(2*D);
		*rmin = (Rmax2-Rmin2-D2)/(2*D);
		*h = Rmax2 - (*rmax)*(*rmax);
		nb_intersects = 2;
		if(D2==Rmax-Rmin) 	nb_intersects = 1;
	}
	return nb_intersects;
}

//------------------------------------------------------------------
// Let there be two points A and B that form the AB segment.
// Let vA and vB be the linear velocities of points A and B respectively.
// If vA and vB are different, then the segment is in rotation around
// a specific axis defined by point p_axis and unit vector u_axis.
// The return value is the estimated rotation angle if the segment
// kept the same length.
// Expressing an angle implies a "normal" vector representing the side of
// the plane at which the observer measures the angle.

extern double V_rotationAngle(Point A, Point B, Vector vA, Vector vB, Vector normal, Vector *u_axis)
{
	Vector v_average = 0.5 * (vA + vB);
	Vector v_diffA = vA - v_average;
	Vector v_diffB = vB - v_average; // v_diffB is -v_diffA

	*u_axis = ((B-A) ^ vB).unit();

	double rotation_angle = 2 * vB.length() / (B-A).length();
	bool rotation_is_clockwise = ( (*u_axis * normal) < 0);
	if(rotation_is_clockwise) rotation_angle = -rotation_angle;

	return rotation_angle;
}
