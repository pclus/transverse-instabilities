/* ---- COMPUTATION OF LARGEST LYAPUNOV EXPONENT IN THE JANSEN-RIT MODEL ----
 * Code written by P. Clusella (pau.clusella@upf.edu)
 * 11/03/2021
 * Free to use, modify, and distribute. Non-comercial.
 * The Jansen-Rit model is a Neural Mass Model composed of 6 ODE.
 * In order to study the stability of limit-cycle solutions in a coupled
 * network of JR-NMM one needs to study the Master Stability Function
 * which relates the Floquet exponents (or Lyapunov exponents) of the
 * syncrhonized with that of the connectivity matrix.
 * 
 * Given parameter values of the system and and eigenvalue of the connectivity
 * matrix this code returns the largest Lyapunov exponent. Applying
 * this routine to all eigenvalues of the connectivity matrix provides
 * the dispersion relation of the system.
 *
 * The period and an initial condition on the limit cycle must be provided.
 *
 * The code uses Gnu Scientific Library (GSL) and BLAS for readability.
 * Code not otpimized by readability, not efficiency.
 *
 * Compilation example:
 * gcc floquet.c -o ../floquet -lm -lgsl -lgslcblas -lblas -O3
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


#define VECTOR(v,i) (v->data[i*v->stride])
#define MATRIX_SET(a,i,j,x) gsl_matrix_set(a,i,j,x)

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

double Sigm(double y);
double dSigm(double y);
void velocity_fields(gsl_vector *x, gsl_vector *k, double p, double eps, double lam);
void system_rk4(gsl_vector *x, gsl_vector *xf, gsl_vector *x0, gsl_vector *k, double dt, double p, double eps, double lam);
void jacobian(gsl_matrix *j, gsl_vector *x, double eps, double lambda );
double normalize(gsl_vector *v );

int main(int argc, char **argv)
{
	double p=atof(argv[1]); 
	double eps=atof(argv[2]); 
	double lambda=atof(argv[3]);
	double period; 

	//FILE *fic=fopen(argv[1],"r");
        //if(fic==NULL)
        //{
        //        fprintf(stderr,"File %s not found!\n",argv[1]);
        //        return -1;
        //}

	int n=6;
	int m=2*n;
	gsl_vector *x=gsl_vector_alloc(m);
	gsl_vector *x0=gsl_vector_alloc(m);
	gsl_vector *xf=gsl_vector_alloc(m);
	gsl_vector *k=gsl_vector_alloc(m);

        //int err=fscanf(fic,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&p,&eps,&period,&VECTOR(x,0), &VECTOR(x,1), &VECTOR(x,2), &VECTOR(x,3), &VECTOR(x,4), &VECTOR(x,5) );
        //fclose(fic);

	double p0,eps0;
	FILE *fin=fopen("7_results/StabilityAnalysis/FlIcs/xtenflic_all.dat","r");
        while(fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                                &p0,&eps0,&period,
                                &VECTOR(x,0),
                                &VECTOR(x,1),
                                &VECTOR(x,2),
                                &VECTOR(x,3),
                                &VECTOR(x,4),
                                &VECTOR(x,5)
                                )!=EOF ){ if(p0==p && eps0==eps) break; }
        fclose(fin);
	if(eps0!=eps || p0!=p)
	{
		fprintf(stderr,"WARNING! Initial values not found: %lf %lf vs  %lf %lf\n",p,eps,p0,eps0);
	}


	gsl_vector_view xt = gsl_vector_subvector(x, n, n);
	for(int i=0;i<6;i++)
		VECTOR(x,n+i)=1.0;
	normalize((gsl_vector *) &xt);

	fprintf(stderr,"\e[92mLargest Lyapunov Exponent of the coupled Jansen-Rit model\n");
	fprintf(stderr,"\e[96mConstant input p = %lf\n",p);
	fprintf(stderr,"\e[96mCoupling strength eps = %lf\n",eps);
	fprintf(stderr,"\e[96mConnectivity matrix eigenvalue = %lf\n",lambda);
	fprintf(stderr,"\e[96mPeriod of limit cycle per = %lf\n",period);


	int nn=1000;				// number of cycles
	double dt=period/(nn*1.0);		// time step
	double l=0;				// sum of the LE
	double t=0;				// time
	int count=0, transient=10*nn, T=10000*nn;	// counter, transient, total time
	for(int i=0;i<T;i++)
	{
		system_rk4(x,xf,x0,k, dt,p,eps,lambda);
		t+=dt;
		count++;
		//fprintf(stdout,"%lf %lf\n",t,VECTOR(x,1)-VECTOR(x,2));
		if(count%nn==0)
		{
			if( i > transient )
			{
				double l0=log(normalize((gsl_vector *) &xt));
				l+=l0;
				//fprintf(stdout,"%.16g %.16g %.16g\n",l0,l,VECTOR(x,1)-VECTOR(x,2));
			}else
			{
				normalize((gsl_vector *) &xt);
			}
		}
	}
	double le = l/(dt*(T-transient));
	fprintf(stderr,"\e[92mLE is %.16g\e[0m\n",le);
	fprintf(stdout,"%lf %lf %lf %.16g\n",p,eps,lambda,le);
	
        gsl_vector_free(x);
        gsl_vector_free(x0);
        gsl_vector_free(xf);
        gsl_vector_free(k);

	return 0;
}

double normalize(gsl_vector *v )
{
	double nv=gsl_blas_dnrm2( v );
	int n=v->size;
	for(int i=0;i<n;i++)
		VECTOR(v,i)=VECTOR(v,i)/nv;
	return nv;
}

double Sigm(double y)
{
        return 2*SIGM_E0/(1+exp(SIGM_R*(SIGM_V0-y)));
}

double dSigm(double y)
{
        double aux_exp=exp(SIGM_R*(SIGM_V0-y));
        return 2*SIGM_E0*SIGM_R*aux_exp/((1+aux_exp)*(1+aux_exp));
}

void velocity_fields(gsl_vector *x, gsl_vector *k, double p, double eps, double lam)
{
	double y0,y1,y2,y3,y4,y5;
	int n=6;

       	y0=VECTOR(x,0);
       	y1=VECTOR(x,1);
       	y2=VECTOR(x,2);
       	y3=VECTOR(x,3);
       	y4=VECTOR(x,4);
       	y5=VECTOR(x,5);

       	double s=Sigm(y1-y2);
       	VECTOR(k,0)=y3;
       	VECTOR(k,1)=y4;
       	VECTOR(k,2)=y5;
       	VECTOR(k,3)=PAR_a*( PAR_A*s-2*y3-PAR_a*y0);
       	VECTOR(k,4)=PAR_A*PAR_a*( p + PAR_C2*Sigm(PAR_C1*y0) + eps*s ) - PAR_a*(2*y4+PAR_a*y1);
       	VECTOR(k,5)=PAR_b*( PAR_B*PAR_C4*Sigm(PAR_C3*y0)-2*y5-PAR_b*y2);

       	gsl_vector_view xt = gsl_vector_subvector(x, n, n);
       	gsl_vector_view kt = gsl_vector_subvector(k, n, n);
	gsl_matrix *j = gsl_matrix_alloc(n, n);
        gsl_matrix_set_zero( j );
	jacobian( j, x, eps, lam );
	gsl_blas_dgemv (CblasNoTrans,  1.0, j , (gsl_vector *) &xt, 0.0,(gsl_vector *) &kt );
        gsl_matrix_free(j);

        return ;
}

void system_rk4(gsl_vector *x, gsl_vector *xf, gsl_vector *x0, gsl_vector *k, double dt, double p, double eps, double lam)
{
        double dtsixth=dt/6.;

	gsl_vector_memcpy(x0,x);
	gsl_vector_memcpy(xf,x);

	velocity_fields(x,k,p,eps, lam);
	gsl_vector_axpby(dtsixth, k, 1.0,xf);
	gsl_vector_axpby( 0.5*dt, k, 1.0,x);
	
	velocity_fields(x,k,p,eps, lam);
	gsl_vector_axpby(2.0*dtsixth, k, 1.0,xf);
	gsl_vector_memcpy(x,x0);
	gsl_vector_axpby( 0.5*dt, k, 1.0,x);
	
	velocity_fields(x,k,p,eps, lam);
	gsl_vector_axpby(2.0*dtsixth, k, 1.0,xf);
	gsl_vector_memcpy(x,x0);
	gsl_vector_axpby( dt, k, 1.0,x);
	
	velocity_fields(x,k,p,eps, lam);
	gsl_vector_axpby( dtsixth, k, 1.0,xf);
	gsl_vector_memcpy(x,xf);

        return ;
}

void jacobian(gsl_matrix *j, gsl_vector *x, double eps, double lambda )
{

        double y0=VECTOR(x,0);
        double y=VECTOR(x,1)-VECTOR(x,2);

        MATRIX_SET(j,0,3, 1.0);
        MATRIX_SET(j,1,4, 1.0);
        MATRIX_SET(j,2,5, 1.0);

        MATRIX_SET(j,3,3, -2*PAR_a);
        MATRIX_SET(j,4,4, -2*PAR_a);
        MATRIX_SET(j,5,5, -2*PAR_b);

        MATRIX_SET(j,3,0, -PAR_a*PAR_a); 
        MATRIX_SET(j,4,1, -PAR_a*PAR_a + PAR_A*PAR_a*lambda*eps*dSigm(y) );
        MATRIX_SET(j,4,2, -PAR_A*PAR_a*lambda*eps*dSigm(y) );
        MATRIX_SET(j,5,2, -PAR_b*PAR_b);

        MATRIX_SET(j,3,1, PAR_a*PAR_A*dSigm(y) );
        MATRIX_SET(j,3,2, -PAR_a*PAR_A*dSigm(y) );
        MATRIX_SET(j,4,0, PAR_a*PAR_A*PAR_C1*PAR_C2*dSigm(PAR_C1*y0) );
        MATRIX_SET(j,5,0, PAR_b*PAR_B*PAR_C3*PAR_C4*dSigm(PAR_C3*y0) );

        return ;
}

