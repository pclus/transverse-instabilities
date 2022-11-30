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

#ifndef OPENMP_H
#define OPENMP_H
	#include <omp.h>
#endif

#ifndef GSL_H
#define GSL_H
        #include <gsl/gsl_rng.h>
        #include <gsl/gsl_rstat.h>
#endif

#define RAND (gsl_rng_uniform(rrr))//((double)rand())/RAND_MAX 
#define RANDI(n) (gsl_rng_uniform_int(rrr,n))
#define RANDIFF (-1+2.*RAND)

///////////// GSL RNG ////////////////////
const gsl_rng_type * T;
gsl_rng * rrr;
////////////////////////////////////

#define SIGM_E0 (2.5)
#define SIGM_R (0.56)
#define SIGM_V0 (6.0)

#define PAR_A (3.25)
#define PAR_a (100.0)
#define PAR_B (22.0)
#define PAR_b (50.0)
#define PAR_C1 (135.0)
#define PAR_C2 (0.8*PAR_C1)
#define PAR_C3 (0.25*PAR_C1)
#define PAR_C4 (0.25*PAR_C1)

#define MATRIX(A,i,j) A[(i)*6+j]
#define TMAT(A,node,iLE,var) A[ (iLE)*N*6 + (node)*6 + var ]

typedef struct graph{
        double *input;
        double *weights;
        double *sigmoid;
        double *dsigmoid;
        double *weightsum;
        int *neighbour;
        int *cumulative;
        int *degree;
}Graph;

typedef struct sys{
	int N;
	int nLE;
	double *x;
	double *x0;
	double *xf;
	double *k;
	double t;
	double dt;
	double p;
	double eps;
}System;

int system_initialize(int argc, char **argv, System *sym, Graph *g );
void system_rk4(System *sys , Graph *g);
void initial_conditions(double *x, int N, double p);
void velocity_fields(System *sym, Graph *g);
double Sigm(double y);
int initialize_network(char *namenet, Graph *g);
void clear_workspace(System *sys, Graph *g);
int read_ic(double *x, int N, char *name);
