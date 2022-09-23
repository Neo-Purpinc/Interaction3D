
/*===============================================================*\
  Thursday June the 20th 2013
  Arash HABIBI
  Interface.h
\*===============================================================*/

#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include <iostream>
//#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <GL/glut.h>
#include <gl.h>

#include "readAnim.h"
#include "Camera.h"
#include "VisualObject.h"
#include "Synch.h"

#define I_PERSPECTIVE    0
#define I_ORTHO_PLUS_Z   1
#define I_ORTHO_PLUS_X   2
#define I_ORTHO_MINUS_Z  3
#define I_ORTHO_MINUS_X  4
#define I_ORTHO_PLUS_Y   5
#define I_ORTHO_MINUS_Y  6

#define I_TRANSLATE 0
#define I_ROTATE 1
#define I_SCALE 2

using namespace std;

class Interface
{
public :
	Interface(char *animfile, int width, int height, bool synch, int argc, char **argv);
	~Interface();

	void exec();
	void check(char *message);
	void _init(int argc, char **argv);
	void _initShade();

	int _width, _height;
	Camera *_cam;
	Point _obs;

	ReadAnim *_ra;
	Synch *_sync;
	int _nb_vo;
	VisualObject **_visual_objects;
	int _nb_points;
	Point *_nuage;

	int _mode; // I_TRANSLATE, I_ROTATE, I_SCALE
	int _smooth;
	int _render_type;  // I_WIREFRAME, I_FLAT, I_SHADE
	int _projection_type; // I_ORTHO_... or I_PERSPECTIVE
	bool _show_normals;
	bool _show_bounding_boxes;
	bool _show_grid;
	bool _show_frame;

	int _mouse_x, _mouse_y;
	int _previous_mouse_x, _previous_mouse_y;
	Point _p_mouse, _p_mouse_xy;
	Vector _v_mouse, _v_mouse_xy;

	int _left_mouse_down;
	int _right_mouse_down;
	int _middle_mouse_down;

	int _transform_mode;

	Point _p_light;
	// GLfloat p_light[4];

	bool _simulation_is_running;
};

void drawLocator(Vector p, Vector color, float size);
void drawFrame(void);
void drawGrid(void);
void display_CB();
void reshape_CB(int largeur, int hauteur);
void mouse_CB(int button, int state, int x, int y);
void mouse_move_CB(int x, int y);
void keyboard_CB(unsigned char key, int x, int y);
void special_CB(int key, int x, int y);
void idle_CB();

void dessineBoule(float x, float y, float z, float rayon);

#endif //_INTERFACE_H_



// setitimer
