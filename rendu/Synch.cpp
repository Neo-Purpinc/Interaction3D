
/*================================================*\
  Mardi 19 novembre 2019
  Arash Habibi
  Synch.cpp
  Un objet permettant d'attendre jusqu'au prochain
  top d'horloge. (La version précédente, permettait
  d'appeler une fonction successivement).
\*================================================*/

#include "Synch.h"

Synch *GB_synch;
sigset_t masque;
sigset_t demasque;

//------------------------------------------------------

static void _my_handler(int bidon)
{
	GB_synch->countOneMoreTick();
}

//------------------------------------------------------

Synch::Synch(float periode)
{
	_init(periode);
	GB_synch = 	this;
}

//------------------------------------------------------

Synch::Synch()
{
	_init(0.04);
	GB_synch = 	this;
}

//------------------------------------------------------

Synch::~Synch()
{
}

//------------------------------------------------------

void Synch::attendSignal(void)
{
	_nb_tops=0;
    sigprocmask(SIG_SETMASK, &demasque, NULL);
	pause();
}

//------------------------------------------------------

void Synch::countOneMoreTick(void)
{
	_nb_tops++;
}

//------------------------------------------------------

int Synch::nbSkips(void)
{
    sigprocmask(SIG_SETMASK, &masque, NULL);

	if(_nb_tops>1)
		printf("Out of synch\n");

	return _nb_tops-1;
}

//------------------------------------------------------

void Synch::_init(float periode)
{
	_period = periode;

    sigemptyset(&demasque);
	sigemptyset(&masque);
    sigaddset(&masque, SIGALRM);
	// _masque et _demasque sont des ensembles de signaux
	// _demasque est nul et _masque contient juste SIGALRM
	// Quand on veut recevoir les signaux : on utilise demasque
	// Quand on veut être tranquille, on utilise masque.

    sigprocmask(SIG_SETMASK, &masque, NULL);
	// Pour l'instant on veut être tranquille
    // Apres setitimer, les signaux d'horloge vont
    // pleuvoir. Heureusement qu'il y a le masque.

	int periode_sec = (int)(periode);
	int periode_usec = (int)((periode - periode_sec)*1000000);

    _timer.it_interval.tv_sec = periode_sec;
    _timer.it_interval.tv_usec = periode_usec;

    _timer.it_value.tv_sec = 0;
    _timer.it_value.tv_usec = 1;

    setitimer(ITIMER_REAL, &_timer, &_old_timer);
	// A partir de maintenant, on va recevoir plein de signaux SIGALRM

	_p_action.sa_handler = &_my_handler;
    sigaction(SIGALRM, &_p_action, &_p_action_anc);
	// Associer p_action à SIGALRM
}
