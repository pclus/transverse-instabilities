#ifndef STD_H
#define STD_H
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <time.h>
#endif

#ifndef BOOL_H
#define BOOL_H
        typedef enum { false, true } bool;
#endif

#ifndef GSL_H
#define GSL_H
        #include <gsl/gsl_rng.h>
        #include <gsl/gsl_rstat.h>
        #include <gsl/gsl_matrix.h>
        #include <gsl/gsl_vector.h>
        #include <gsl/gsl_linalg.h>
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

#define MATRIX(A,i,j) A[(i)*SYS_M+j]
#define TMAT(A,node,iLE,var) A[ (N*SYS_M) + (node*SYS_M + var)*nLE + (iLE) ] 

#define SYS_M (6)

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
	int m;
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

typedef struct lyapunov_exponents
{
        size_t rows;
        size_t cols;
        gsl_matrix_view m;
        gsl_matrix *M;
        gsl_matrix *T;
        gsl_matrix *Q;
        gsl_matrix *R;
        double *expo;
        int count;
}Lyapunov;


int system_initialize(int argc, char **argv, System *sym, Graph *g );
void system_rk4(System *sys , Graph *g);
void initial_conditions(System *sys);
void velocity_fields(System *sym, Graph *g);
double Sigm(double y);
double dSigm(double y);
int initialize_network(char *namenet, Graph *g, int nLE);
void clear_workspace(System *sys, Graph *g);

void lyap_inititalize(Lyapunov *lyap, size_t rows, size_t cols, double *x);
void lyap_finish(Lyapunov *lyap);
void lyap_exp(Lyapunov *lyap, bool flag);

