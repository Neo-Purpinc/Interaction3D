/*----------------------------------------*\
  Tuesday Novembre 26th 2019
  Arash Habibi
  Vector.c
\*----------------------------------------*/

#include "Vector.h"

//------------------------------------------
// param : x,y,z : doubles
// valeur de retour : un vecteur
// A partir de x,y,z, cette fonction retourne
// le vecteur (x,y,z)

Vector V_new(double x, double y, double z)
{
	Vector v;
	v.x = x;
	v.y = y;
	v.z = z;
	return v;
}

//------------------------------------------
// param s : v : Vecteur
// param x,y,z, pointeurs vers doubles
// valeur de retour : rien
// Cette fonction retourne par adresse les
// composantes x,y,z du vecteur v.

void V_get(Vector v, double *x, double *y, double *z)
{
	*x = v.x;
	*y = v.y;
	*z = v.z;
}

//------------------------------------------
// param : v : Vecteur
// param label : chaine de caractères
// valeur de retour : rien
// Cette fonction affiche les coordonnées du vecteur
// (à des fins de debug)

void V_show(Vector v, char *label)
{
	printf("%s : (%f,%f,%f)\n",label,v.x,v.y,v.z);
}

//------------------------------------------
// param : c : double
// param : v : Vecteur
// valeur de retour Vecteur
// Cette fonction retourne le vecteur c * v
// où toutes les composantes ont été
// multipliées par c

Vector V_mult(double c, Vector v)
{
	v.x *= c;
	v.y *= c;
	v.z *= c;
	return v;
}

//------------------------------------------
// params : c1 et c2 : double
// params : v1 et v2 : Vecteurs
// valeur de retour Vecteur
// Cette fonction retourne le vecteur c1*v1 + c2*v2

Vector V_add(double c1, Vector v1, double c2, Vector v2)
{
	Vector res;
	res.x = c1 * v1.x + c2 * v2.x;
	res.y = c1 * v1.y + c2 * v2.y;
	res.z = c1 * v1.z + c2 * v2.z;
	return res;
}

//------------------------------------------
// params : v1 et v2 : vecteurs
// valeur de retour :  double
// cette fonction retourne le produit
// scalaire de v1 et de v2

double V_dot(Vector v1, Vector v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//------------------------------------------
// param : v : vecteur
// valeur de retour : double
// Cette fonction retourne le module (ou norme)
// du vecteur v.

double V_module(Vector v)
{
	return sqrt(V_dot(v,v));
}

//------------------------------------------
// param : v : vecteur
// valeur de retour : vecteur
// Cette fonction retourne un vecteur colinéaire
// à v mais de norme 1.

Vector V_unit(Vector v)
{
	double m = V_module(v);
	if(m==0)
		return V_new(1,0,0);
	else
		return V_mult((1.0/m),v);
}

//------------------------------------------
