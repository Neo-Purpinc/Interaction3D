
/*================================================*\
  Mardi 19 novembre 2019
  Arash Habibi
  Synch.h
  Un objet permettant d'attendre jusqu'au prochain
  top d'horloge. (La version précédente, permettait
  d'appeler une fonction successivement).
\*================================================*/

#ifndef __SYNCH_H__
#define __SYNCH_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/time.h>
#include <signal.h>
#include <sys/signal.h>

class Synch
{
public:
	Synch();
	Synch(float periode);
	~Synch();

	void attendSignal(void);
	void countOneMoreTick(void);
	int nbSkips(void);

protected:
	struct itimerval _timer;
	struct itimerval _old_timer;
	struct sigaction _p_action;
	struct sigaction _p_action_anc;

	float _period;
	int _nb_tops;

	void _init(float periode);

};

static void _my_handler(int);

#endif // __SYNCH_H__
