
/*========================================================================*\
  Wednesday November the 7th 2012
  Main.cpp
  Arash HABIBI
\*========================================================================*/

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "Interface.h"

//---------------------------------------------------------------

void readArgs(int argc, char **argv, char *animfile, bool *synch, int *largeur, int *hauteur)
{
	if(argc!=2 && argc!=3 && argc!=5 && argc!=6)
	{
		printf("\nUsage : %s <animfile> [-s] [-x larg haut]\n\n",argv[0]);
		exit(1);
	}
	else
	{
		strcpy(animfile,argv[1]);

		*synch=false;
		for(int i=1; i<argc; i++)
			if(strncmp(argv[i],"-s",2)==0)
				*synch=true;

		*largeur = *hauteur = 800;
		for(int i=1; i<argc-2; i++)
			if(strncmp(argv[i],"-x",2)==0)
			{
				*largeur=atoi(argv[i+1]);
				*hauteur=atoi(argv[i+2]);
			}
	}
}

//---------------------------------------------------------------

int main(int argc, char **argv)
{
	char animfilename[100];
	bool synch=false;
	int largeur, hauteur;

	readArgs(argc, argv, animfilename, &synch, &largeur, &hauteur);

	Interface I(animfilename, largeur, hauteur, synch, argc, argv);
	I.exec();

	return 0;
}
