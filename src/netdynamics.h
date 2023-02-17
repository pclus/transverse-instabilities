#ifndef STD_H
#define STD_H
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <time.h>
#endif

#ifndef GSL_H
#define GSL_H
        #include <gsl/gsl_rng.h>
        #include <gsl/gsl_rstat.h>
        #include <gsl/gsl_matrix.h>
        #include <gsl/gsl_vector.h>
        #include <gsl/gsl_linalg.h>
#endif

#ifndef SYSTEMSL_H
#define SYSTEMSL_H
	#include "system.h"
#endif

#ifndef BOOL_H
#define BOOL_H
	typedef enum { false, true } bool;
#endif

void print_field(System *sys, FILE *fout);
void print_pattern(System *sys, FILE *fout);
