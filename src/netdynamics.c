#include "netdynamics.h"
#define FLAG_PATTERN true

int parent_thread(int argc, char **argv, int *pfd)
{
	System sys;
	Graph g;
	int i;
	int count=0; 
	if(system_initialize(argc, argv, &sys, &g)==-1) return -1;

	//------------------ SET OUTPUT FILES ------------
	FILE *fout;
	FILE *fout2;
	char *namm=malloc(sizeof(char)*200);
	sprintf(namm,"4_outputs/mean_%s.dat",argv[1]);
	fout=fopen(namm,"w");
	if(FLAG_PATTERN)
	{
		sprintf(namm,"4_outputs/pattern_%s.bin",argv[1]);
		fout2=fopen(namm,"wb");
	}
	//-----------------------------------------------
	
	//-----------------------------------------------
	double sec=1/sys.dt;
	double T=2e3;
	double transient=T-1e3;
	//double transient=10e4-1e3;
	while( sys.t < T ) 
	{
		if(count%10000==0) fprintf(stderr,"Integrating t=%.0f\r",sys.t);

		system_rk4( &sys, &g );
		count++;

                if( count%((int)sec/1000)==0 && sys.t >= transient ) 
		{
                	print_field(&sys,fout);
			if(FLAG_PATTERN) print_pattern(&sys,fout2);
		}
		if(count%10000==0)
		{
			fflush(fout);
			if(FLAG_PATTERN) fflush(fout2);
		}

	}
	fprintf(stderr,"\033[1;00m\n");
	//snapshot(sys.x,sys.N,"final_state.dat");

	free(namm);
	fclose(fout);
	if(FLAG_PATTERN) fclose(fout2);
	clear_workspace(&sys, &g);

	fprintf(stderr,"Ending parent process\n");
	return 1;
}


void print_field(System *sys, FILE *fout)
{
	gsl_rstat_workspace *rstat_p = gsl_rstat_alloc();
	for(int i=0; i<sys->N;i++)
		gsl_rstat_add( MATRIX( sys->x,i,1 )-MATRIX( sys->x,i,2 ), rstat_p);

        double mean = gsl_rstat_mean(rstat_p);
        double var = gsl_rstat_variance(rstat_p);
	double ex1=MATRIX(sys->x,0,1)-MATRIX(sys->x,0,2);
	double ex2=MATRIX(sys->x,1,1)-MATRIX(sys->x,1,2);

        fprintf(fout,"%lf %lf %lf %lf %lf\n",sys->t,mean, var, ex1 , ex2 );

	gsl_rstat_free(rstat_p);

        return ;
}


void snapshot(double *x, int N, char *name)
{
	FILE *fout=fopen(name,"w");
	fprintf(stderr,"State snapshot at %s\n",name);
	for(int i=0;i<N*6;i++)
		fprintf(fout,"%.16g\n",x[i]);
	fclose(fout);
	return ;
}

void print_pattern(System *sys, FILE *fout)
{
	for(int i=0;i<sys->N;i++)
	{
		double ex1=MATRIX(sys->x,i,1)-MATRIX(sys->x,i,2);
                fwrite(&ex1,sizeof(double),1,fout);
	}
	return ;
}


int main(int argc, char **argv )
{
        int pfd[2];
        parent_thread(argc, argv, pfd);

        return -1;
}

