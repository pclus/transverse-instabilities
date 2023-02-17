#include "netdynamics.h"

int main(int argc, char **argv )
{
	// Initialize system
	System sys;
	Graph g;
	if(system_initialize(argc, argv, &sys, &g)==-1) return -1;
	Lyapunov lyap;
        lyap_inititalize(&lyap, sys.N*sys.m, sys.nLE, sys.x+sys.m*sys.N);
        lyap_exp(&lyap, false);

	// Set output files
	FILE *fout;
	FILE *fout2;
	char *namm=malloc(sizeof(char)*200);
	sprintf(namm,"outputs/mean_%s.dat",argv[1]);
	fout=fopen(namm,"w");
	sprintf(namm,"outputs/pattern_%s.bin",argv[1]);
	fout2=fopen(namm,"wb");
	
	// Main integration
	int count=0; 
	int nn=100;
	double sec=1/sys.dt;
	double T=2e3;
	double transient=T-1e3;
	double transientLE=1e3;
	fprintf(stdout,"%lf\n",MATRIX(sys.x,0,3));
	while( sys.t < T ) 
	{
		if(count%10000==0) fprintf(stderr,"Integrating t=%.0f\r",sys.t);

		system_rk4( &sys, &g );
		count++;

                if( count%nn==0  )
                {
                        lyap_exp(&lyap, sys.t>= transientLE );
                }
                if( count%((int)sec/1000)==0 && sys.t >= transient ) 
		{
        	   	print_field(&sys,fout);
			print_pattern(&sys,fout2);
		}
		if(count%10000==0)
		{
			fflush(fout);
			fflush(fout2);
		}
	}
	fprintf(stderr,"\033[1;00m\n");

        fprintf(stdout,"%lf %lf ",sys.p,sys.eps);
        for(int i=0;i<sys.nLE;i++)
                fprintf(stdout,"%lf ",lyap.expo[i]/(sys.t-transientLE));
        fprintf(stdout,"\n");


	// Freeing memory
	free(namm);
	fclose(fout);
	fclose(fout2);
	lyap_finish(&lyap);
	clear_workspace(&sys, &g);

	fprintf(stderr,"Ending parent process\n");
	return 1;
}

// Calculate and save the mean-field
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

// Print the spatiotemporal dynamics in a binary file
void print_pattern(System *sys, FILE *fout)
{
	for(int i=0;i<sys->N;i++)
	{
		double ex1=MATRIX(sys->x,i,1)-MATRIX(sys->x,i,2);
                fwrite(&ex1,sizeof(double),1,fout);
	}
	return ;
}


