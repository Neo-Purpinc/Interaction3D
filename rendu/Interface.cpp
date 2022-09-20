
/*===============================================================*\
  Thursday June the 20th 2013
  Arash HABIBI
  Interface.c
\*===============================================================*/

#include "Interface.h"

Interface *I;

//------------------------------------------------------------------

Interface::Interface(char *animfilename, int width, int height, bool synch, int argc, char **argv)
{
	_width  = width;
	_height = height;
	_cam = new Camera(_width,_height);
	_obs = _cam->lookfrom();

	_smooth = 0;
	_show_grid = true;
	_show_frame = true;

	_projection_type = I_PERSPECTIVE;

	_left_mouse_down=0;
	_right_mouse_down=0;
	_middle_mouse_down=0;

	_mouse_x = _mouse_y = 0;
	_previous_mouse_x = _previous_mouse_y = 0;

	_init(argc,argv);
	_simulation_is_running = false;

	_mode = I_ROTATE;

	_ra = new ReadAnim(animfilename);
	// _ra->check("interface");
	_nb_vo = _ra->nbVisualObjects();
	_visual_objects = new VisualObject*[_nb_vo];
	_ra->visualObjects(_visual_objects,_nb_vo);

	_nb_points = _ra->nbPoints();
	_nuage = new Point[_nb_points];

	_ra->getNextPoints(_nuage,_nb_points);

	_sync=NULL;
	if(synch)
		_sync = new Synch(_ra->timeStep());
	/*
	for(int i=0;i<_nb_vo;i++)
		_visual_objects[i]->check("interface");
	for(int i=0;i<_nb_points;i++)
	_nuage[i].check("interface");
	*/

	I = this;
}

//------------------------------------------------------------------

Interface::~Interface()
{
}

//------------------------------------------------------------------

void Interface::exec()
{
    glutMainLoop();
}

//------------------------------------------------------------------

void Interface::_init(int argc, char **argv)
{
    int windowPosX = 100, windowPosY = 100;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
    glutInitWindowSize(_width,_height);
    glutInitWindowPosition(windowPosX,windowPosY);
    glutCreateWindow(argv[0]);
    glViewport(0, 0, _width, _height);

	glClearColor(0,0,0,0);

    glutDisplayFunc(display_CB);
	glutReshapeFunc(reshape_CB);
	glutKeyboardFunc(keyboard_CB);
	glutSpecialFunc(special_CB);
	glutMouseFunc(mouse_CB);
	glutMotionFunc(mouse_move_CB);
	glutIdleFunc(idle_CB);
	// glutPassiveMotionFunc(passive_mouse_move_CB);
}

//------------------------------------------------------------------

void drawLocator(Point p, float size, float r, float g, float b)
{
	glColor3f(r,g,b);

	double x,y,z;
	p.get(&x,&y,&z);

	glBegin(GL_LINES);
	glVertex3f( x-size, y, z);
	glVertex3f( x+size, y, z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f( x, y-size, z);
	glVertex3f( x, y+size, z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f( x, y, z-size);
	glVertex3f( x, y, z+size);
	glEnd();
}

//---------------------------------------------------------------

void drawFrame(void)
{
	float size = 1;

	glColor3f(1.0, 0.0, 0.0);

	glBegin(GL_LINES);
	glVertex3f( 0.0, 0.0, 0.0);
	glVertex3f( size, 0.0, 0.0);
	glEnd();

	glColor3f(0.0, 1.0, 0.0);

	glBegin(GL_LINES);
	glVertex3f( 0.0, 0.0, 0.0);
	glVertex3f( 0.0, size, 0.0);
	glEnd();

	glColor3f(0.0, 0.0, 1.0);

	glBegin(GL_LINES);
	glVertex3f( 0.0, 0.0, 0.0);
	glVertex3f( 0.0, 0.0, size);
	glEnd();
}

//---------------------------------------------------------------

void drawGrid(void)
{
	int grid_size = 10;
	int i;

	glColor3f(0.25, 0.25, 0.25);
	for(i=-grid_size;i<=grid_size;i++)
	{
		glBegin(GL_LINES);
		glVertex3f(-grid_size, 0.0, i);
		glVertex3f( grid_size, 0.0, i);
		glEnd();
	}

	for(i=-grid_size;i<=grid_size;i++)
	{
		glBegin(GL_LINES);
		glVertex3f(i, 0.0,-grid_size);
		glVertex3f(i, 0.0, grid_size);
		glEnd();
	}
}



//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------


//------------------------------------------------------------------
//	C'est le display callback. A chaque fois qu'il faut
//	redessiner l'image, c'est cette fonction qui est
//	appelee. Tous les dessins doivent etre faits a partir
//	de cette fonction.
//------------------------------------------------------------------

void display_CB()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	float window_half_width, window_half_height;
	I->_cam->windowHalfSizeIn3D(&window_half_width, &window_half_height);
	window_half_width *= 2;
	window_half_height *= 2;

	Point lf = I->_cam->lookfrom();
	Point la = I->_cam->lookat();

	if(I->_projection_type == I_PERSPECTIVE)
	{
		gluPerspective(I->_cam->fovy(),I->_cam->ratio(),I->_cam->zNear(),I->_cam->zFar());

		double lfx, lfy, lfz, lax, lay, laz, uvx, uvy, uvz;
		Vector uv = I->_cam->upVector();
		lf.get(&lfx,&lfy,&lfz);
		la.get(&lax,&lay,&laz);
		uv.get(&uvx,&uvy,&uvz);
		gluLookAt( lfx,lfy,lfz, lax,lay,laz, uvx,uvy,uvz);
	}
	else
	{
		glOrtho(-window_half_width, +window_half_width, -window_half_height, +window_half_height, I->_cam->zNear(),I->_cam->zFar());

		glTranslated(-lf.get(V_X), -lf.get(V_Y), -lf.get(V_Z));
		switch(I->_projection_type)
		{
		case I_ORTHO_PLUS_X : glRotatef(-90,0,1,0); I->_obs.set(1000,0,0); break;
		case I_ORTHO_MINUS_X : glRotatef(90,0,1,0); I->_obs.set(-1000,0,0); break;
		case I_ORTHO_PLUS_Z : glRotatef(0,0,1,0); I->_obs.set(0,0,1000); break;
		case I_ORTHO_MINUS_Z : glRotatef(180,0,1,0); I->_obs.set(0,0,-1000); break;
		case I_ORTHO_PLUS_Y : glRotatef(90,1,0,0); I->_obs.set(0,1000,0); break;
		case I_ORTHO_MINUS_Y : glRotatef(90,-1,0,0); I->_obs.set(0,-1000,0); break;
		default : cerr << "diplay_CB : unknown projection type : " << I->_projection_type << endl;
		}
	}

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//-----------------------

	if(I->_show_grid)	drawGrid(); // dessiner la grille
	if(I->_show_frame)	drawFrame();// dessiner repere

	glColor3f(1,1,1);

	for(int i=0; i<I->_nb_vo; i++)
		I->_visual_objects[i]->display(I->_nuage,I->_nb_points,I->_obs);

	//---------------------------

	glFlush();
    glutSwapBuffers();
}

//------------------------------------------------------------------

void reshape_CB(int width, int height)
{
	I->_width = width;
	I->_height = height;
	I->_cam->setPixelDimensions(width,height);
    glViewport(0, 0,I->_width,I->_height);
}

//------------------------------------------------------------------
// Cette fonction permet de réagir en fonction de la position de
// la souris (x,y), en fonction du bouton de la souris qui a été
// pressé ou relaché.
//------------------------------------------------------------------

void mouse_CB(int button, int state, int x, int y)
{
	I->_mouse_x = x;
	I->_mouse_y = y;

	if(state==GLUT_DOWN)
	{
		switch(button)
		{
		case GLUT_LEFT_BUTTON : I->_left_mouse_down = 1; break;
		case GLUT_RIGHT_BUTTON : I->_right_mouse_down = 1; break;
		case GLUT_MIDDLE_BUTTON : I->_middle_mouse_down = 1; break;
		}
		I->_previous_mouse_x = x;
		I->_previous_mouse_y = y;
	}
	else
	{
		switch(button)
		{
		case GLUT_LEFT_BUTTON : I->_left_mouse_down = 0; break;
		case GLUT_RIGHT_BUTTON : I->_right_mouse_down = 0; break;
		case GLUT_MIDDLE_BUTTON : I->_middle_mouse_down = 0; break;
		}
		I->_v_mouse.set(0,0,0);
 		I->_v_mouse_xy.set(0,0,0);
	}
	glutPostRedisplay();
}

//------------------------------------------------------------------
// Cette fonction permet de réagir aux mouvements de la souris
// quand l'un des boutons de la souris est pressé.
//------------------------------------------------------------------

void mouse_move_CB(int x, int y)
{
	I->_mouse_x = x;
	I->_mouse_y = y;

	int dx = I->_mouse_x - I->_previous_mouse_x;
	int dy = I->_mouse_y - I->_previous_mouse_y;

	float translate_sensitivity = 0.05;
	float rotate_sensitivity = 0.5;

	if(I->_left_mouse_down)
	{
		I->_cam->turnAroundX(rotate_sensitivity * dy);
		I->_cam->turnAroundWorldY(rotate_sensitivity * -dx);
	}

	if(I->_middle_mouse_down)
	{
		Vector move(translate_sensitivity*(-dx),translate_sensitivity*dy,0);
		I->_cam->move(move);
	}

	if(I->_right_mouse_down)
	{
		Vector move(0,0,translate_sensitivity*dx);
		I->_cam->move(move);
	}

	I->_p_mouse = I->_cam->pixelPositionIn3D(x,y);
	I->_p_mouse_xy = I->_cam->pixelPositionInXYPlane(x,y);

	Point previous = I->_cam->pixelPositionIn3D(I->_previous_mouse_x,I->_previous_mouse_y);
	Point previous_xy = I->_cam->pixelPositionInXYPlane(I->_previous_mouse_x,I->_previous_mouse_y);

	if(I->_left_mouse_down)
	{
		I->_v_mouse = I->_p_mouse - previous;
		I->_v_mouse_xy = I->_p_mouse_xy - previous_xy;
	}
	else
	{
		I->_v_mouse.set(0,0,0);
		I->_v_mouse_xy.set(0,0,0);
	}

	I->_previous_mouse_x = x;
	I->_previous_mouse_y = y;

	I->_obs = I->_cam->lookfrom();

	glutPostRedisplay();
}

//------------------------------------------------------------------
// Cette fonction permet de réagir au fait que l'utilisateur
// presse une touche (non-spéciale) du clavier.
//------------------------------------------------------------------

void keyboard_CB(unsigned char key, int x, int y)
{
	char projection_type[50];
	Vector u_devant   (0,0,1);
	Vector u_derriere (0,0,-1);

	switch(key)
	{
	case 27 :  exit(0); break;
	case 's' : I->_ra->getNextPoints(I->_nuage,I->_nb_points); break; // Une image
	case 'S' : /*simulateOneStep();*/ break; // un échantillon
	case 'g' : I->_show_grid = ! I->_show_grid; break;
	case 'i' :
		I->_cam->initPosit(); I->_obs = I->_cam->lookfrom(); break;
	case 'o':
	{	I->_ra->goToFrame(0); I->_ra->getNextPoints(I->_nuage,I->_nb_points); I->_simulation_is_running=false; }
		break;

	case 't' : I->_mode = I_TRANSLATE; printf("Translate mode\n"); break;
	case 'r' : I->_mode = I_ROTATE; printf("Rotate mode\n"); break;
	case 'z' : I->_mode = I_SCALE; printf("Zoom mode \n"); break;

	case 'f' : I->_show_frame = ! I->_show_frame; break;

	case '+' : I->_cam->move(u_devant); break;
	case '-' : I->_cam->move(u_derriere); break;

	case ' ' :
			I->_projection_type = (I->_projection_type+1)%7;
			switch(I->_projection_type)
			{
			case I_PERSPECTIVE : strcpy(projection_type,"Projection perspective"); break;
			case I_ORTHO_PLUS_X : strcpy(projection_type,"Vue de droite (+x)"); break;
			case I_ORTHO_MINUS_X :strcpy(projection_type,"Vue de gauche (-x)"); break;
			case I_ORTHO_PLUS_Z : strcpy(projection_type,"Vue de face (+z)"); break;
			case I_ORTHO_MINUS_Z : strcpy(projection_type,"Vue de dos (-z)"); break;
			case I_ORTHO_PLUS_Y : strcpy(projection_type,"Vue de dessus (+y) "); break;
			case I_ORTHO_MINUS_Y : strcpy(projection_type,"Vue de dessous (-y)"); break;
			default : strcpy(projection_type,"Point de vue inconnu"); break;
			}

			cerr << "projection_type = " << projection_type << endl; break;

	case 13 : I->_simulation_is_running = ! I->_simulation_is_running ; break; // lancer la simulation
		/*
		  case 45 : printf ("6\n"); break;
		  case 232 : printf ("7\n"); break;
		  case 95 : printf ("8\n"); break;
		  case 231 : printf ("9\n"); break;
		  case 224 : printf ("0\n"); break;
		*/
	}
	// printf("key=%c %d\n",key,key);
	glutPostRedisplay();
}

//------------------------------------------------------------------
// Cette fonction permet de réagir au fait que l'utilisateur
// presse une touche spéciale (F1, F2 ... F12, home, end, insert,
// haut, bas, droite, gauche etc).
//------------------------------------------------------------------

void special_CB(int key, int x, int y)
{
	Vector u_droite   (1,0,0);
	Vector u_gauche   (-1,0,0);
	Vector u_haut     (0,1,0);
	Vector u_bas      (0,-1,0);
	Vector u_devant   (0,0,1);
	Vector u_derriere (0,0,-1);

	switch(key)
	{
 	case GLUT_KEY_LEFT :
		if (I->_mode==I_TRANSLATE)
			I->_cam->move(u_droite);
		else if (I->_mode==I_ROTATE)
			I->_cam->turnAroundWorldY(5.1);
		else
			I->_cam->zoom(5.0);
		break;

	case GLUT_KEY_RIGHT:
		if (I->_mode==I_TRANSLATE)
			I->_cam->move(u_gauche);
		else if (I->_mode==I_ROTATE)
			I->_cam->turnAroundWorldY(-5.1);
		else
			I->_cam->zoom(-5.0);
		break;

 	case GLUT_KEY_UP   :
		if (I->_mode==I_TRANSLATE)
			I->_cam->move(u_bas);
		else if (I->_mode==I_ROTATE)
			I->_cam->turnAroundX(-5.1);
		else
			I->_cam->zoom(-5.0);
		break;

 	case GLUT_KEY_DOWN :
		if (I->_mode==I_TRANSLATE)
			I->_cam->move(u_haut);
		else if (I->_mode==I_ROTATE)
			I->_cam->turnAroundX(5.1);
		else
			I->_cam->zoom(5.0);
		break;

 	case GLUT_KEY_PAGE_UP   :
		if (I->_mode==I_TRANSLATE)
			I->_cam->move(u_derriere);
		else if (I->_mode==I_ROTATE)
			I->_cam->turnAroundZ(5.1);
		else
			I->_cam->zoom(5.0);
		break;

 	case GLUT_KEY_PAGE_DOWN :
		if (I->_mode==I_TRANSLATE)
			I->_cam->move(u_devant);
		else if (I->_mode==I_ROTATE)
			I->_cam->turnAroundZ(-5.1);
		else
			I->_cam->zoom(-5.0);
		break;

 	case GLUT_KEY_F1 :
		break;

	default : fprintf(stderr,"touche inconnue.");
	}

	I->_obs = I->_cam->lookfrom();

	glutPostRedisplay();
}

//------------------------------------------------------------------------

void idle_CB()
{

	if(I->_simulation_is_running== true)
	{
		int nb_skips=0;

		if(I->_sync!=NULL)
			I->_sync->attendSignal();

		int fr = I->_ra->getNextPoints(I->_nuage,I->_nb_points);
		if (fr >= I->_ra->nbFrames())
			I->_simulation_is_running=false;

		if(I->_sync!=NULL)
			nb_skips = I->_sync->nbSkips();

		for(int i=0; i<nb_skips; i++)
		{
			fr = I->_ra->getNextPoints(I->_nuage,I->_nb_points);
			if (fr >= I->_ra->nbFrames())
				I->_simulation_is_running=false;
		}
	}

	I->_v_mouse.set(0,0,0);
	I->_v_mouse_xy.set(0,0,0);

	glutPostRedisplay();
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
