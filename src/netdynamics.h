#ifndef STD_H
#define STD_H
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <time.h>
        #include <unistd.h> // needed for pipes
        #include <fcntl.h> // needed only for the nonblocking read
        #include <sys/types.h>
        #include <sys/wait.h>
	#include <signal.h>
#endif

#ifndef GSL_H
#define GSL_H
        #include <gsl/gsl_rng.h>
	#include <gsl/gsl_rstat.h>
#endif

#ifndef SYSTEMSL_H
#define SYSTEMSL_H
        #include "system.h"
#endif

#ifndef BOOL_H
#define BOOL_H
	typedef enum { false, true } bool;
#endif

#ifndef OPENMP_H
#define OPENMP_H
        #include <omp.h>
#endif

#define NORMA(x,y) (sqrt((x)*(x)+(y)*(y)))

int parent_thread(int argc, char **argv, int *pfd);
void print_field(System *sys, FILE *fout);
void snapshot(double *x, int N, char *name);
void print_pattern(System *sys, FILE *fout);
