#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "Vector.h"

int main()
{
	FILE *f = fopen("00_noanim.pv","w");
	if(f==NULL)	{ perror("00_nomanim.pv");	exit(1); }

	int nb_points = 4;
	float time_step = 1;
	int nb_frames=2;

	// placement des 4 points
	Vector tous_les_points[nb_points];
	tous_les_points[0] = V_new(-1,0,-1);
	tous_les_points[1] = V_new(-1,0, 1);
	tous_les_points[2] = V_new( 1,0, 1);
	tous_les_points[3] = V_new( 1,0,-1);

	// ecriture de l'entete
	fprintf(f,"#PV==\n%f\n%d\n%d\n",time_step,nb_points,nb_frames);
	fprintf(f,"oLine 0:2\n");
	fprintf(f,"Sphere 2:2\n");
	fprintf(f,"Point 2:3\n");
	fprintf(f,"====\n");

	// parcours des frames (images)
	for(int fr=0; fr<nb_frames; fr++)
	{
		// parcours des points
		double x,y,z;
		for(int n=0;n<nb_points;n++)
		{
			V_get(tous_les_points[n],&x,&y,&z);
			fprintf(f,"%f %f %f\n",x,y,z);
		}
	}
	fclose(f);

	return 0;
}
