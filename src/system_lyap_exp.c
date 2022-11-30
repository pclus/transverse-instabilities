#ifndef SYSTEMSL_H
#define SYSTEMSL_H
        #include "system.h"
#endif

int system_initialize(int argc, char **argv, System *sys , Graph *g)
{
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

        sys->x=(double *)malloc(sizeof(double)*N*6);
        sys->x0=(double *)malloc(sizeof(double)*N*6);
        sys->xf=(double *)malloc(sizeof(double)*N*6);
        sys->k=(double *)malloc(sizeof(double)*N*6);

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
                for(int j=0;j<6;j++)
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
		for(int j=0;j<6;j++)
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
	double a=0,b=0.1;
	double y[6];
	for(int i=0;i<6;i++) y[i]=0;
	//double x0= -0.738561; // p=30 eps=50
	double x0= 7.0; 
	y[0]=Sigm(x0)*PAR_A/PAR_a;
        y[1]=(p+PAR_C2*Sigm(PAR_C1*y[0]))*PAR_A/PAR_a;
        y[2]=PAR_C4*Sigm(PAR_C3*y[0])*PAR_B/PAR_b;
        y[3]=0;
        y[4]=0;
        y[5]=0;


	for(int i=0; i<N;i++)
	{

		MATRIX(x,i,0)=a*y[0]+b*RANDIFF;
		MATRIX(x,i,1)=a*y[1]+b*RANDIFF;
		MATRIX(x,i,2)=a*y[2]+b*RANDIFF;
		MATRIX(x,i,3)=a*y[3]+b*RANDIFF;
		MATRIX(x,i,4)=a*y[4]+b*RANDIFF;
		MATRIX(x,i,5)=a*y[5]+b*RANDIFF;
	}
//	MATRIX(x,0,0)=a+0.1*RAND;
//	MATRIX(x,0,1)=b+0.1*RAND;
//	MATRIX(x,0,2)=c+0.1*RAND;
//	MATRIX(x,1,0)=d+0.1*RAND;
//	MATRIX(x,1,1)=e+0.1*RAND;
//	MATRIX(x,1,2)=f+0.1*RAND;

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
		i++;
		if(i>N)
		{
			fprintf(stderr,"\e[1;31mWARNING! Reading initial conditions:");
			fprintf(stderr,"file larger than expected\n");
			break;
		}

	}
	if(i!=N)
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

	for(int iLE=0;iLE<nLE;iLE++)
	for(int a=0;a<N;a++)
	{
		jacobian(sys->jac,sys->x+6*a);
		for(int i=0;i<6;i++)
        	{
			TMAT(sys->k,a,iLE,i)=0;
               		for(int j=0;j<6;j++)
			{
				TMAT(sys->k,a, iLE, i)+=TMAT(sys->x,a, iLE, j)*MATRIX(sys->jac,i,j);
			}
        	}
	}	

	for(int i=0;i<N;i++)
	{
		double y1=MATRIX(sys->x,i,1);
		double y2=MATRIX(sys->x,i,2);
		g->sigmoid[i]=Sigm(y1-y2);
		g->dsigmoid[i]=dSigm(y1-y2);
	}

	for(int i=0;i<N;i++)
        {
		g->input[i]=0;
                for(int j=0;j<g->degree[i];j++)
                {
                        int jjj=g->cumulative[i]+j;
                        int jj=g->neighbour[jjj];
                        g->input[i]+=g->weights[jjj]*g->sigmoid[jj];
			for(int iLE=0;iLE<nLE;i++)
			{
				TMAT(sys->k,i,iLE,4)+=g->weights[jjj]*g->dsigmoid[jj]*\
						      (TMAT(sys->x,jj,iLE,1)-TMAT(sys->x,jj,iLE,2));
			}
                }
        }

	for(int i=0;i<N;i++)
	{
		double y0=MATRIX(sys->x,i,0);
		double y1=MATRIX(sys->x,i,1);
		double y2=MATRIX(sys->x,i,2);
		double y3=MATRIX(sys->x,i,3);
		double y4=MATRIX(sys->x,i,4);
		double y5=MATRIX(sys->x,i,5);

		MATRIX(sys->k,i,0)=y3;
		MATRIX(sys->k,i,1)=y4;
		MATRIX(sys->k,i,2)=y5;
		MATRIX(sys->k,i,3)=PAR_a*( PAR_A*g->sigmoid[i]-2*y3-PAR_a*y0);
		MATRIX(sys->k,i,4)=PAR_A*PAR_a*( sys->p + PAR_C2*Sigm(PAR_C1*y0) + sys->eps*g->input[i] ) - PAR_a*(2*y4+PAR_a*y1);
		MATRIX(sys->k,i,5)=PAR_b*( PAR_B*PAR_C4*Sigm(PAR_C3*y0)-2*y5-PAR_b*y2);

	}


	return ;
}

void jacobian(double *j, double *x)
{
        int n=6;
        double y0=x[0];
        double y=x[1]-x[2];
	
		
        MATRIX(j,0,3) = (1.0);
        MATRIX(j,1,4) = (1.0);
        MATRIX(j,2,5) = (1.0);

        MATRIX(j,3,3) = (-2*PAR_a);
        MATRIX(j,4,4) = (-2*PAR_a);
        MATRIX(j,5,5) = (-2*PAR_b);

        MATRIX(j,3,0) = (-PAR_a*PAR_a);
        MATRIX(j,4,1) = (-PAR_a*PAR_a);
        MATRIX(j,5,2) = (-PAR_b*PAR_b);

        MATRIX(j,3,1) = (PAR_a*PAR_A*dSigm(y) );
        MATRIX(j,3,2) = (-PAR_a*PAR_A*dSigm(y) );
        MATRIX(j,4,0) = (PAR_a*PAR_A*PAR_C1*PAR_C2*dSigm(PAR_C1*y0) );
        MATRIX(j,5,0) = (PAR_b*PAR_B*PAR_C3*PAR_C4*dSigm(PAR_C3*y0) );

        return ;
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
                g->weights[j]=c/mean_strength;         
                //g->weights[j]=c/(1.0*g->degree[a]);         
                //g->weights[j]=c/(1.0*g->weightsum[a]);         
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
