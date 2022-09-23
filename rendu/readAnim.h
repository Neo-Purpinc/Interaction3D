/*======================================================*\
  Jeudi 21 novembre 2019
  Arash Habibi
  readAnim.h
\*======================================================*/

#ifndef __READ_ANIM_H__
#define __READ_ANIM_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>

#include "Vector.h"
#include "VisualObject.h"

#define RA_NAME_MAX 50
#define RA_LINE_MAX 300

class ReadAnim
{
public:
	ReadAnim();
	ReadAnim(char *animfilename);
	~ReadAnim();

	void check(char *);

	int nbVisualObjects(void);
	int visualObjects(VisualObject **vo, int nvo);

	float timeStep(void);
	int nbPoints(void);
	int nbFrames(void);
	int getNextPoints(Point *nuage, int np);
	int goToFrame(int fr);
	int nextFrame(void);


protected:
	char _anim_file_name[RA_NAME_MAX];
	FILE *_anim_file;
	int _nb_visual_objects;
	int _nb_points;
	int _nb_frames;
	int _current_frame;
	float _time_step;
	bool _header_read;

	int _errorMessage(int error_id, int line, char *buffer);
	int _readHeader(void);
};

#endif // __READ_ANIM_H__
