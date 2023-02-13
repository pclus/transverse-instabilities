#ifndef SYSTEMSL_H
#define SYSTEMSL_H
        #include "system.h"
#endif

double newpar_alpha;

int system_initialize(int argc, char **argv, System *sys , Graph *g)
{

	newpar_alpha=1.0;
        if(argc<2)
        {
                //fprintf(stderr,"USAGE: ./....");
                return -1;
        }
        sys->t=0;
	sys->dt=0.001;
        int N=initialize_network(argv[2], g);
        if(N==-1)
        {
                fprintf(stderr,"Errors in system.c, program breaking.\n");
                return -1;
        }
	sys->p=atof(argv[3]);
        sys->eps=atof(argv[4]);
	sys->N=N;
	int n=8;
	sys->n=n;

	fprintf(stderr,"\033[1;96mCoupled Jansen-Rit NMM\n");
        fprintf(stderr,"p = %lf\n",sys->p);
        fprintf(stderr,"eps = %lf\n",sys->eps);
        fprintf(stderr,"\033[1;31m\n");
        ////////////////////////////////////////

        ////////////////////////////////////////
        gsl_rng_env_setup();

        T = gsl_rng_default;
        rrr = gsl_rng_alloc(T);

        // Different perturbation       
        srand(time(0));
        gsl_rng_set(rrr,time(0));

        // Same perturbation always
        //srand(0);
        //gsl_rng_set(rrr,0);
        ////////////////////////////////////////
        //omp_set_dynamic(0);
        //omp_set_num_threads(4);

        sys->x=(double *)malloc(sizeof(double)*N*n);
        sys->x0=(double *)malloc(sizeof(double)*N*n);
        sys->xf=(double *)malloc(sizeof(double)*N*n);
        sys->k=(double *)malloc(sizeof(double)*N*n);

	if(read_ic(sys->x,sys->N, "initial_conditions.dat")==-1)
        	initial_conditions(sys->x,sys->N,sys->p);

        return 1;
}

void system_rk4(System *sys , Graph *g)
{
        double dtsixth=sys->dt/6.;
        double w[4]={1,2,2,1};
        double dw[4]={0.5,0.5,1,1};
	double dww,ww;

	int N=sys->N;
	double dt=sys->dt;

	{
                velocity_fields(sys, g);

                ww=w[0];
                dww=dw[0];

	        //#pragma omp parallel for default(shared)
                for(int i=0;i<N;i++)
                for(int j=0;j<sys->n;j++)
                {
			MATRIX(sys->x0,i,j)=MATRIX(sys->x,i,j);
                        MATRIX(sys->xf,i,j)=MATRIX(sys->x,i,j)+ww*dtsixth*MATRIX(sys->k,i,j);
                        MATRIX(sys->x ,i,j)=MATRIX(sys->x0,i,j)+dww*dt*MATRIX(sys->k,i,j);
                }
	}

        for(int s=1;s<4;s++)
        {
                velocity_fields(sys, g);

                ww=w[s];
                dww=dw[s];
	        //#pragma omp parallel for default(shared)
		for(int i=0;i<N;i++)
		for(int j=0;j<sys->n;j++)
		{
			MATRIX(sys->xf,i,j)=MATRIX(sys->xf,i,j)+ww*dtsixth*MATRIX(sys->k,i,j);
			MATRIX(sys->x ,i,j)=MATRIX(sys->x0,i,j)+dww*dt*MATRIX(sys->k,i,j);
		}
        }
	sys->t+=dt;
        double *temp=sys->x; sys->x=sys->xf; sys->xf=temp;

        return ;
}

void initial_conditions(double *x, int N, double p)
{
	double a=1.0,b=0.001;
	double y[sys->n];
	for(int i=0;i<sys->n;i++) y[i]=0;
	
	y[0]=0.01;
	y[1]=6.82;
	y[2]=2.156;
        y[3]=1.3
       	y[4]=0.589;
        y[5]=0.5
        y[6]=0.0
        y[7]=0.01;

	double r;
	for(int i=0; i<N;i++)
	{
		r=RANDIFF;
		pert=r*r;
		MATRIX(x,i,0)=a*y[0]+b*r;
		r=RANDIFF;
		MATRIX(x,i,1)=a*y[1]+b*r;
		r=RANDIFF;
		MATRIX(x,i,2)=a*y[2]+b*r;
		r=RANDIFF;
		MATRIX(x,i,3)=a*y[3]+b*r;
		r=RANDIFF;
		MATRIX(x,i,4)=a*y[4]+b*r;
		r=RANDIFF;
		MATRIX(x,i,5)=a*y[5]+b*r;
		r=RANDIFF;
		MATRIX(x,i,6)=a*y[5]+b*r;
		r=RANDIFF;
		MATRIX(x,i,7)=a*y[5]+b*r;
	}

        return ;
}

int read_ic(double *x, int N, char *name)
{
	FILE *fin=fopen(name,"r");
	if(fin==NULL)
	{
		fprintf(stderr,"\e[42m\e[39mFile %s does not exist, new IC will be generated\e[m\n",name);
		return -1;
	}
	fprintf(stderr,"\e[92mInitial conditions from '%s' file\n",name);
	fprintf(stderr,"\033[1;31m\n");

	int i=0;
	while(fscanf(fin,"%lf\n",x+i)!=EOF)
	{
		x[i]+=0.001*RAND;
		i++;
		if(i>6*N)
		{
			fprintf(stderr,"\e[1;31mWARNING! Reading initial conditions:");
			fprintf(stderr," file larger than expected\n");
			break;
		}

	}
	if(i!=6*N)
	{
		fprintf(stderr,"\e[1;31mSizes missmatch in read_ic: %d %d\n",i,N);
	}
	fprintf(stderr,"\033[1;00m\n");
	fclose(fin);

	return 1;
}

void velocity_fields(System *sys, Graph *g)
{
	int N=sys->N;
	//double y0,y1,y2,y3,y4,y5;

	//#pragma omp parallel for default(shared)
	for(int i=0;i<N;i++)
	{
		double y1=MATRIX(sys->x,i,1);
		double y2=MATRIX(sys->x,i,2);
		g->sigmoid[i]=Sigm(y1-y2);
	}

	//#pragma omp parallel for default(shared)
	for(int i=0;i<N;i++)
        {
		g->input[i]=0;
                for(int j=0;j<g->degree[i];j++)
                {
                        int jjj=g->cumulative[i]+j;
                        int jj=g->neighbour[jjj];
                        g->input[i]+=g->weights[jjj]*g->sigmoid[jj];
                }
        }

	//#pragma omp parallel for default(shared)
	for(int i=0;i<N;i++)
	{
		double y0=MATRIX(sys->x,i,0);
		double y1=MATRIX(sys->x,i,1);
		double y2=MATRIX(sys->x,i,2);
		double y3=MATRIX(sys->x,i,3);
		double y4=MATRIX(sys->x,i,4);
		double y5=MATRIX(sys->x,i,5);
		double x1=MATRIX(sys->x,i,6);
		double x2=MATRIX(sys->x,i,7);

		MATRIX(sys->k,i,0)=y3;
		MATRIX(sys->k,i,1)=y4;
		MATRIX(sys->k,i,2)=y5;
		MATRIX(sys->k,i,3)=PAR_a*( PAR_A*Sigm(y1-y2 +x1)-2*y3-PAR_a*y0);
		MATRIX(sys->k,i,4)=PAR_a* (PAR_A*( p + PAR_C2*Sigm(PAR_C1*y0) ) - 2*y4 - PAR_a*y1);
		MATRIX(sys->k,i,5)=PAR_b*( PAR_B*PAR_C4*Sigm(PAR_C3*y0+x1)-2*y5-PAR_b*y2);
		MATRIX(sys->k,i,6)=x2;
		MATRIX(sys->k,i,7)=PAR_a*(PAR_A*sys->eps*g->input[i] - 2*x2-PAR_a*x1);


	}

	return ;
}

double Sigm(double y)
{
        return 2*SIGM_E0/(1+exp(SIGM_R*(SIGM_V0-y)));
}

int initialize_network(char *namenet, Graph *g)
{
        int i,cumul=0,j=0;
        int a,b;
        double c;
	int N;

        FILE *fnet=fopen(namenet,"r");

        if(fnet==NULL)
        {
                fprintf(stderr,"ERROR (initialize_network()): File %s not found!\n",namenet);
                return -1;
        }

        fscanf(fnet,"%d\n",&N);
        fprintf(stderr,"Reading %d nodes\n",N);
        g->degree=(int *) malloc(sizeof(int)*N);
        g->weightsum=(double *) malloc(sizeof(double)*N);
        for(i=0;i<N;i++)
	{
		g->degree[i]=0;
		g->weightsum[i]=0;
	}
        while(fscanf(fnet,"%d %d %lf\n",&a,&b,&c)!=EOF)
        {
                if(a<0 || a>N || b<0 || b>N)
                {
                        fprintf(stderr,"ERROR (initialize_network()): Index out of range\n a=%d\nb=%d\n",a,b);
                        return -1;
                }
                g->degree[a]++;
                g->weightsum[a]+=c;
                cumul++;
        }
	double mean_strength=0;
	for(i=0;i<N;i++)
		mean_strength+=g->weightsum[i];
	mean_strength=mean_strength/N;
	fprintf(stderr,"Mean node strength %lf\n",mean_strength);

        fclose(fnet);
        g->cumulative=(int *) malloc(sizeof(int)*N);
        g->cumulative[0]=0;
        for(i=1;i<N;i++)
                g->cumulative[i]=g->cumulative[i-1]+g->degree[i-1];

        g->weights=(double *) malloc(sizeof(double)*cumul);
        g->neighbour=(int *) malloc(sizeof(int)*cumul);
        fnet=fopen(namenet,"r");
        fscanf(fnet,"%*d\n");
        while(fscanf(fnet,"%d %d %lf\n",&a,&b,&c)!=EOF)
        {
                g->neighbour[j]=b;
                //g->weights[j]=c;				// raw setup
                //g->weights[j]=c/mean_strength;             	// realistic setup
                //g->weights[j]=c/(1.0*g->degree[a]);        	// binary homogneous
                g->weights[j]=c/(1.0*g->weightsum[a]);     	// weighted homogeneous    
                j++;
        }
        fclose(fnet);
        for(i=0;i<N;i++)
        {
                if(g->degree[i]<=0)
                {
                        fprintf(stderr,"WARNING: node %d is isolated\n",i);
                }
        }
        if(j!=cumul)
        {
                fprintf(stderr,"ERROR (initialize_network()): Missmatch between j=%d and cumulative=%d\n",j,cumul);
                return -1;
        }
        g->input=(double *) malloc(sizeof(double)*N);
        g->sigmoid=(double *) malloc(sizeof(double)*N);
        return N;
}

void clear_workspace(System *sys, Graph *g)
{
	fprintf(stderr,"Freeing workspace\n");
	free(g->input);
	free(g->weights);
	free(g->neighbour);
	free(g->cumulative);
	free(g->degree);
	free(g->weightsum);
	free(g->sigmoid);

	free(sys->xf);
	free(sys->x);
	free(sys->x0);
	free(sys->k);
	fprintf(stderr,"Done\n");
	return ;
}
