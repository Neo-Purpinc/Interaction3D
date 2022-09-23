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
	FILE *f = fopen("01_vagues.pv","w");
	if(f==NULL)	{	perror("01_vagues.pv"); exit(1);	}

	float time_step = 0.1;
	int nb_frames = 150;
	int nb_points = 100;

	// Placement de tous les points
	int count=0;
	Vector tous_les_points[nb_points];
	for(int lgn=0; lgn<10; lgn++)
		for(int col=0; col<10; col++)
			tous_les_points[count++] = V_new(col-5,0,lgn-5);

	// ecriture de l'entete
	fprintf(f,"#PV==\n");
	fprintf(f,"%f\n",time_step);
	fprintf(f,"%d\n",nb_points);
	fprintf(f,"%d\n",nb_frames);
	for(int lgn=0; lgn<10; lgn++)
		fprintf(f,"oLine %d:%d\n",10*lgn,10*(lgn+1)-1);
	fprintf(f,"====\n");

	// -------------- animation
	float pulsation=0.15;
	float amplitude=0.3;

	double x,y,z;
	// parcours des frames (images)
	for(int fr=0; fr<nb_frames; fr++)
		// parcours des points
		for(int n=0; n<nb_points; n++)
		{
			V_get(tous_les_points[n],&x,&y,&z);
			y=amplitude*sin(pulsation*fr - 0.5*(x-z));
			fprintf(f,"%f %f %f\n",x,y,z);
		}

	fclose(f);
	return 0;
}
