/*----------------------------------------*\
  Tuesday Novembre 26th 2019
  Arash Habibi
  01_battement.c
\*----------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Vector.h"

int main()
{
	FILE *f = fopen("01_battement.pv","w");
	if(f==NULL)	{	perror("01_battement.pv"); exit(1);	}

	int nb_points=4;
	int nb_frames=100;
	float time_step=0.1;

	// placement des 4 points
	Vector tous_les_points[nb_points];
	tous_les_points[0] = V_new(-1,0,-1);
	tous_les_points[1] = V_new(-1,0, 1);
	tous_les_points[2] = V_new( 1,0, 1);
	tous_les_points[3] = V_new( 1,0,-1);

	// ecriture de l'entete
	fprintf(f,"#PV==\n%f\n%d\n%d\n",time_step,nb_points,nb_frames);
	fprintf(f,"cLine 0:3\n");
	fprintf(f,"=====\n");

	// -------------- animation
	float pulsation=0.5;
	float amplitude=0.3;

	// parcours des frames (images)
	for(int fr=0;fr<nb_frames;fr++)
	{
		// parcours des points
		for(int n=0;n<nb_points;n++)
		{
			double x,y,z;
			V_get(tous_les_points[n],&x,&y,&z);
			y=amplitude * sin(pulsation*fr-0.5*(x-z));
			fprintf(f,"%f %f %f\n",x,y,z);
		}
	}
	fclose(f);

	return 0;
}
