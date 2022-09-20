/*======================================================*\
  Jeudi 21 novembre 2019
  Arash Habibi
  readAnim.cpp

  Lire les fichiers d'animation dont voici le format :
  #PV==
  0.04 // pas de temps : 25 images par seconde
  4  // nombre de points
  1  // nombre d'images
  points 0:1
  spheres 0:2 0.25
  cline 0:3
  oline 3-4
  ===
  -1 0 1  (point d'indice 0)
  -1 0 -1 (point d'indice 1)
  1 0 -1  (point d'indice 2)
  1 0 1   (idem : 3)
  -1 0 1  (point d'indice 0 encore et ainsi de suite)

\*======================================================*/


#include "readAnim.h"

//-----------------------------------------------

ReadAnim::ReadAnim()
{
	strcpy(_anim_file_name,"");
	_anim_file = NULL;

	_nb_visual_objects=0;
	_nb_points = 0;
	_nb_frames = 0;
	_time_step = 0;
	_current_frame = 0;
	_header_read = false;
}

//-----------------------------------------------

ReadAnim::ReadAnim(char *animfilename)
{
	strcpy(_anim_file_name,animfilename);
	_anim_file = fopen(animfilename,"r");
	if(_anim_file==NULL)
	{
		perror(animfilename);
		exit(1);
	}
	else
	{
		_nb_visual_objects=0;
		_nb_points = 0;
		_nb_frames = 0;
		_time_step = 0;
		_current_frame = 0;
		_header_read = false;

		_readHeader();
	}
}

//-----------------------------------------------

ReadAnim::~ReadAnim()
{
	if(_anim_file!=NULL)
		fclose(_anim_file);
}

//-----------------------------------------------

void ReadAnim::check(char *label)
{
	printf("-------- Check anim file %s --------------\n",label);
	printf("Fichier : %s\n",_anim_file_name);
	printf("_nb_visual_objects=%d\n",_nb_visual_objects);
	printf("_nb_points=%d\n",_nb_points);
	printf("_nb_frames=%d\n",_nb_frames);
	printf("_time_step=%f\n",_time_step);
	printf("_current_frame=%d\n",_current_frame);
}

//-----------------------------------------------

int ReadAnim::nbVisualObjects(void)
{
	if(!_header_read)
		_readHeader();

	return _nb_visual_objects;
}

//-----------------------------------------------

int ReadAnim::visualObjects(VisualObject **vo, int nvo)
{
	char buffer[RA_LINE_MAX], keyword[50];
	int pt1, pt2, count=0;
	float radius;

	rewind(_anim_file);
	for(int i=0; i<4; i++)
		fgets(buffer,RA_LINE_MAX,_anim_file);

	count=0;
	while(count<nvo)
	{
		bool indiv;
		fgets(buffer,RA_LINE_MAX,_anim_file);
		if(sscanf(buffer,"%s %d:%d",keyword,&pt1,&pt2)==3)
			indiv=false;
		else if(sscanf(buffer,"%s %d-%d",keyword,&pt1,&pt2)==3)
		{
			indiv=true;
			if(strncasecmp(keyword,"oline",5)!=0)
			{
				printf("ReadAnim::visualObjects : the individual linking can only be used with open lines.");
				exit(1);
			}
		}
		else
		{	printf("---"); exit(1); }


		if ( strncasecmp(keyword,"point",5)==0 )
			for(int ind=pt1;ind<=pt2;ind++)
				vo[count++] = new VisualObject(VO_POINT,ind);

		else if ( strncasecmp(keyword,"oline",5)==0 )
		{
			if(indiv)
				vo[count++] = new VisualObject(VO_LINE,pt1,pt2);
			else
				for(int ind=pt1;ind<pt2;ind++)
					vo[count++] = new VisualObject(VO_LINE,ind,ind+1);
		}

		else if ( strncasecmp(keyword,"cline",5)==0 )
		{
			for(int ind=pt1;ind<pt2;ind++)
				vo[count++] = new VisualObject(VO_LINE,ind,ind+1);
			vo[count++] = new VisualObject(VO_LINE,pt2,pt1);
		}

		else if ( strncasecmp(keyword,"sphere",6)==0 )
		{
			sscanf(buffer,"%s %d:%d %f",keyword,&pt1,&pt2,&radius);
			for(int ind=pt1;ind<=pt2;ind++)
				vo[count++] = new VisualObject(VO_SPHERE,ind,radius);
		}
		else
		{
			printf("ReadAnim::visualObjects : Unknown visible object specifier %s\n",keyword);
			return -1;
		}
	}

	fgets(buffer,RA_LINE_MAX,_anim_file);
	sscanf(buffer,"%s %d:%d",keyword,&pt1,&pt2);
	// printf("buffer=%s keyword=%s\n",buffer,keyword);
	if(keyword[0]!='=')
	{
		printf("ReadAnim::visualObjects : Header should end with '==='\n");
		return -1;
	}

	return 0;
}

//-----------------------------------------------

float ReadAnim::timeStep(void)
{
	if(!_header_read)
		_readHeader();

	return _time_step;
}

//-----------------------------------------------

int ReadAnim::nbPoints(void)
{
	if(!_header_read)
		_readHeader();

	return _nb_points;
}

//-----------------------------------------------

int ReadAnim::nbFrames(void)
{
	if(!_header_read)
		_readHeader();

	return _nb_frames;
}

//-----------------------------------------------

int ReadAnim::getNextPoints(Point *nuage, int np)
{
	char buffer[RA_LINE_MAX];
	float x,y,z;

	if (_current_frame < _nb_frames)
	{
		for(int i=0; i<np; i++)
		{
			if(fgets(buffer,RA_LINE_MAX,_anim_file)!=NULL)
			{
				sscanf(buffer,"%f %f %f",&x, &y, &z);
				nuage[i].set(x,y,z);
			}
			else
			{
				printf("ReadAnim::getNextPoints : unable to read further. Stopped at frame %d point %d\n",_current_frame,i);
				exit(1);
			}
		}
		_current_frame++;
	}
	return _current_frame;
}

//-----------------------------------------------
// On commence avec le frame 1

int ReadAnim::goToFrame(int fr)
{
	char buffer[RA_LINE_MAX];

	if(fr<0)
	{
		printf("ReadAnim::goToFrame : invalid frame number. First frame is frame number 0.\n");
		exit(1);
	}
	else if(fr >= _nb_frames)
	{
		printf("ReadAnim::goToFrame : invalid frame number. Last frame is frame number %d.\n",_nb_frames-1);
		exit(1);
	}
	else
	{
		_readHeader();

		for(int f=0;f<fr-1;f++)
			for(int np=0;np<_nb_points;np++)
				if(fgets(buffer,RA_LINE_MAX,_anim_file)==NULL)
				{
					printf("ReadAnim::goToFrame : unable to reach frame %d. Stopped at frame %d.\n",fr,f);
					exit(1);
				}
	}
	return 0;
}

//-----------------------------------------------

int ReadAnim::nextFrame(void)
{
	return 0;
}

//-----------------------------------------------
//-----------------------------------------------
//-----------------------------------------------

int ReadAnim::_errorMessage(int error_id, int line, char *buffer)
{
	printf("ReadAnim::_readHeadder : file %s ligne %d error %d : \n",_anim_file_name,line,error_id);
	printf("%s",buffer);
	switch(error_id)
	{
	case 0 : printf("File does not begin with #PV== and therefore is not a valid particle view file.\n"); break;
	case 1 : printf("File ended before header was complete.\n"); break;
	case 2 : printf("Unable to read time step.\n"); break;
	case 3 : printf("Unable to read number of points.\n"); break;
	case 4 : printf("Unable to read visual object specifications.\n"); break;
	case 5 : printf("Point start index out of range.\n"); break;
	case 6 : printf("Point end index out of range.\n"); break;
	case 7 : printf("Point end index smaller than point start index.\n"); break;
	case 8 : printf("Unrecognized visual object specification (possible keywords : oline, cline, points, spheres).\n"); break;
	case 9 : printf("Unable to read number of frames.\n"); break;
	case 10 : printf("The individual linking syntaxe (n-m) can only be used with open lines.\n"); break;
	default : printf("Unrecognized error identifier.\n");break;
	}
	exit(-1);
}

//-----------------------------------------------
// reads the magic title, the number of  points
// and the time step and returns the number of
// lines read.

int ReadAnim::_readHeader(void)
{
	char buffer[RA_LINE_MAX], keyword[50];
	int pt1, pt2, line=1, count=0;
	bool header=true;

	rewind(_anim_file);

	if (fgets(buffer,RA_LINE_MAX,_anim_file)==NULL) _errorMessage(1,line,buffer);
	if (strncmp(buffer,"#PV==",5)!=0) _errorMessage(0, line,buffer);

	line++;
	if (fgets(buffer,RA_LINE_MAX,_anim_file)==NULL) _errorMessage(1,line,buffer);
	if(sscanf(buffer,"%f",&_time_step)!=1)   		_errorMessage(2,line,buffer);

	line++;
	if (fgets(buffer,RA_LINE_MAX,_anim_file)==NULL) _errorMessage(1,line,buffer);
	if(sscanf(buffer,"%d",&_nb_points)!=1)   		_errorMessage(3,line,buffer);

	line++;
	if (fgets(buffer,RA_LINE_MAX,_anim_file)==NULL) _errorMessage(1,line,buffer);
	if(sscanf(buffer,"%d",&_nb_frames)!=1)   		_errorMessage(3,line,buffer);

	while(header)
	{
		line++;
		if (fgets(buffer,RA_LINE_MAX,_anim_file)==NULL)  _errorMessage(1,line,buffer);

		if (sscanf(buffer,"%s",keyword)>=1)
		{
			if(keyword[0]=='=')
				header = false;
			else
			{
				bool indiv;
				if (sscanf(buffer,"%s %d:%d",keyword,&pt1,&pt2)==3)
					indiv=false;
				else if(sscanf(buffer,"%s %d-%d",keyword,&pt1,&pt2)==3)
				{
					indiv=true;
					if(strncasecmp(keyword,"oline",5)!=0)
						_errorMessage(10,line,buffer);
				}
				else
					_errorMessage(4,line,buffer);

				if (pt1<0 || pt1>=_nb_points) _errorMessage(5,line,buffer);
				if (pt2<0 || pt2>=_nb_points) _errorMessage(6,line,buffer);

				if (pt1 > pt2) _errorMessage(7,line,buffer);

				if ( strncasecmp(keyword,"point",5)==0 )
					count += pt2-pt1+1;
				else if ( strncasecmp(keyword,"cline",5)==0 )
					count += pt2-pt1+1;
				else if ( strncasecmp(keyword,"oline",5)==0 )
					if (indiv)	count++;
					else	count += pt2-pt1;
				else if ( strncasecmp(keyword,"sphere",6)==0 )
					count += pt2-pt1+1;
				else
					_errorMessage(8,line,buffer);
			}
		}
	}
	_nb_visual_objects=count;
	_current_frame=0;
	_header_read = true;
	return line;
}
