#ifndef SYSTEMSL_H
#define SYSTEMSL_H
        #include "system.h"
#endif

const double PAR_2a = 2*PAR_a ;
const double PAR_2b = 2*PAR_b ;
const double PAR_a2 = PAR_a*PAR_a ;
const double PAR_b2 = PAR_b*PAR_b ;
const double PAR_aA = PAR_a*PAR_A ;
const double PAR_bB = PAR_b*PAR_B ;
const double PAR_bBC3C4 = PAR_b*PAR_B*PAR_C3*PAR_C4;
const double PAR_bBC4 = PAR_b*PAR_B*PAR_C4;
const double PAR_aAC1C2 = PAR_a*PAR_A*PAR_C1*PAR_C2;

int system_initialize(int argc, char **argv, System *sys , Graph *g)
{
        if(argc!=6)
        {
                fprintf(stderr,
		"Run a simulation with: ./netdyn <output files id> <network file> <p> <eps>\n");
                return -1;
        }
	int nLE=atoi(argv[5]);
        int N=initialize_network(argv[2], g,nLE);
	int m=SYS_M;
        if(N==-1)
        {
                fprintf(stderr,"Errors in system.c, program breaking.\n");
                return -1;
        }
        sys->t=0;
	sys->dt=0.001;
	sys->p=atof(argv[3]);
        sys->eps=atof(argv[4]);
	sys->N=N;
	sys->m=m;
	sys->nLE=nLE;

	fprintf(stderr,"\033[1;96mCoupled Jansen-Rit NMM\n");
        fprintf(stderr,"p = %lf\n",sys->p);
        fprintf(stderr,"eps = %lf\n",sys->eps);
        fprintf(stderr,"\033[1;31m\n");

        gsl_rng_env_setup();
        T = gsl_rng_default;
        rrr = gsl_rng_alloc(T);
        srand(time(0));
        gsl_rng_set(rrr,time(0));

        sys->x=(double *)malloc(sizeof(double)*(nLE+1)*N*m);
        sys->x0=(double *)malloc(sizeof(double)*(nLE+1)*N*m);
        sys->xf=(double *)malloc(sizeof(double)*(nLE+1)*N*m);
        sys->k=(double *)malloc(sizeof(double)*(nLE+1)*N*m);

        initial_conditions(sys);

        return 1;
}

void initial_conditions(System *sys)
{
	for(int i=0;i<sys->m*sys->N*(sys->nLE+1);i++)
		sys->x[i]=RANDIFF;

        return ;
}

void system_rk4(System *sys , Graph *g)
{
        double dtsixth=sys->dt/6.;
        double w[4]={1,2,2,1};
        double dw[4]={0.5,0.5,1,1};
	double dww,ww;

	double dt=sys->dt;
	int M=sys->N*(1+sys->nLE);
	int m=sys->m;

	{
                velocity_fields(sys, g);

                ww=w[0];
                dww=dw[0];

                for(int i=0;i<M;i++)
                for(int j=0;j<m;j++)
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
		for(int i=0;i<M;i++)
		for(int j=0;j<m;j++)
		{
			MATRIX(sys->xf,i,j)=MATRIX(sys->xf,i,j)+ww*dtsixth*MATRIX(sys->k,i,j);
			MATRIX(sys->x ,i,j)=MATRIX(sys->x0,i,j)+dww*dt*MATRIX(sys->k,i,j);
		}
        }
	sys->t+=dt;
        double *temp=sys->x; sys->x=sys->xf; sys->xf=temp;

        return ;
}

void velocity_fields(System *sys, Graph *g)
{
        int N=sys->N;
        int nLE=sys->nLE;

        for(int i=0;i<N;i++)
        {
                double y1=MATRIX(sys->x,i,1);
                double y2=MATRIX(sys->x,i,2);
                double y=y1-y2;
                g->sigmoid[i]=Sigm(y);
                double dS=dSigm(y);
                for(int j=0;j<nLE;j++)
                {
                        g->dsigmoid[i*nLE+j]=dS*(TMAT(sys->x,i,j,1)-TMAT(sys->x,i,j,2));
                        TMAT(sys->k,i,j,4)=0;
                }
        }

        for(int i=0;i<N;i++)
        {
                g->input[i]=0;
                for(int j=0;j<g->degree[i];j++)
                {
                        int jjj=g->cumulative[i]+j;
                        int jj=g->neighbour[jjj];
                        g->input[i]+=g->weights[jjj]*g->sigmoid[jj];
                        for(int iLE=0;iLE<nLE;iLE++)
                        {
                                TMAT(sys->k,i,iLE,4)+=g->weights[jjj]*g->dsigmoid[jj*nLE+iLE];
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
                MATRIX(sys->k,i,3)= PAR_aA*g->sigmoid[i] - PAR_2a*y3 - PAR_a2*y0;
                MATRIX(sys->k,i,4)=PAR_aA*( sys->p +
                                PAR_C2*Sigm(PAR_C1*y0) + sys->eps*g->input[i] )
                                - PAR_2a*y4 - PAR_a2*y1; 
                MATRIX(sys->k,i,5)=PAR_bBC4*Sigm(PAR_C3*y0) - PAR_2b*y5 - PAR_b2*y2;

                double dS=dSigm(y1-y2);
                for(int iLE=0;iLE<nLE;iLE++)
                {
                        double dy0=TMAT(sys->x,i,iLE,0);
                        double dy1=TMAT(sys->x,i,iLE,1);
                        double dy2=TMAT(sys->x,i,iLE,2);
                        double dy3=TMAT(sys->x,i,iLE,3);
                        double dy4=TMAT(sys->x,i,iLE,4);
                        double dy5=TMAT(sys->x,i,iLE,5);

                        TMAT(sys->k,i,iLE,0)=dy3;
                        TMAT(sys->k,i,iLE,1)=dy4;
                        TMAT(sys->k,i,iLE,2)=dy5;

                        TMAT(sys->k,i,iLE,3)=-PAR_2a*dy3 - PAR_a2*dy0 
                                             + PAR_aA*dS*(dy1-dy2); 

                        TMAT(sys->k,i,iLE,4)*=sys->eps*PAR_aA;
                        TMAT(sys->k,i,iLE,4)+=-PAR_2a*dy4 - PAR_a2*dy1 
                                             + PAR_aAC1C2*dSigm(PAR_C1*y0)*dy0;

                        TMAT(sys->k,i,iLE,5)= -PAR_2b*dy5 - PAR_b2*dy2 
                                            + PAR_bBC3C4*dSigm(PAR_C3*y0)*dy0;
                }
        }

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

int initialize_network(char *namenet, Graph *g, int nLE)
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
                        fprintf(stderr,
			"ERROR (initialize_network()): Index out of range\n a=%d\nb=%d\n",a,b);
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
                fprintf(stderr,"ERROR (initialize_network()):\
			       	Missmatch between j=%d and cumulative=%d\n",j,cumul);
                return -1;
        }
        g->input=(double *) malloc(sizeof(double)*N);
        g->sigmoid=(double *) malloc(sizeof(double)*N);
        g->dsigmoid=(double *) malloc(sizeof(double)*N*nLE);
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
	free(g->dsigmoid);

	free(sys->xf);
	free(sys->x);
	free(sys->x0);
	free(sys->k);
	fprintf(stderr,"Done\n");
	return ;
}

void lyap_inititalize(Lyapunov *lyap, size_t rows, size_t cols, double *x)  // cols=nLE ; rows=6*N;
{
        lyap->rows=rows;
        lyap->cols=cols;
        lyap->T=gsl_matrix_alloc(cols,cols);
        lyap->Q=gsl_matrix_alloc(rows,rows);
        lyap->R=gsl_matrix_alloc(cols,cols);    // only upper diagonal
        lyap->expo=malloc(sizeof(double)*cols);
        for(int i=0;i<cols;i++) lyap->expo[i]=0;
        lyap->m = gsl_matrix_view_array(x, rows, cols);
        lyap->M = (gsl_matrix *) &lyap->m;
        lyap->count=0;
        return ;
}

void lyap_finish(Lyapunov *lyap)
{
        free(lyap->expo);
        gsl_matrix_free(lyap->T);
        gsl_matrix_free(lyap->Q);
        gsl_matrix_free(lyap->R);
        return ;
}

void lyap_exp(Lyapunov *lyap, bool flag)
{
        int ncols=lyap->cols;
        int nrows=lyap->rows;
        gsl_linalg_QR_decomp_r( lyap->M, lyap->T);
        gsl_linalg_QR_unpack_r( lyap->M, lyap->T, lyap->Q, lyap->R);
        for(int i=0;i<ncols;i++)
        {
                gsl_vector_view c=gsl_matrix_column( lyap->Q , i );
                gsl_matrix_set_col( lyap->M, i, (gsl_vector *) &c);
        }

        if(flag)
        {
                double l;
                for(int i=0;i<ncols;i++)
                {
                        l=log(fabs(gsl_matrix_get(lyap->R,i,i)));
                        lyap->expo[i]+=l;
                }
                lyap->count++;
        }
        return ;
}

