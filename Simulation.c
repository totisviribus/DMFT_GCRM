/*
   2022. 6. 16.

   Simulation Code for 'Generalized Consumer-Resource Model' with self-regulation.
   Reference -- arXiv:2305.12341


   Packages below are needed
   -GSL
   -OpenMP
 */
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<omp.h>

static double delta_t;
static double th_conv;
static double zero_threshold;
static double aetol;

static double aveM=1.;
static double aveD=1.;
static double aver=0.1;

int omp_thrd=5;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//							Butcher Tableau of Runge-Kutta-Fehlberg method
static double a[5][5]	 ={ {1./4.,			0., 			0., 			0., 			0.	 	},
	{3./32., 		9./32., 		0., 			0., 			0. 		},
	{1932./2197., 	-7200./2197., 	7296./2197., 	0., 			0. 		},
	{439./216., 	-8., 			2680./513., 	-845./4104., 	0. 		},
	{-8./27., 		2., 			-3544./2565., 	1859./4104., 	-11./40.} };
//
static double b5[6]	 =	{16./135., 		0., 			6656./12825., 	28561./56430., 	-9./50., 	2./55.};
static double b4[6]	 =	{25./216., 		0., 			1408./2565., 	2197./4104.,	-1./5., 	0.	  };
static double delta_b[6]=	{1./360.,		0.,				-128./4275.,	-2197./75240.,	1./50.,		2./55.};
//
static double c[7]=		{0.,			1./4., 			3./8., 			12./13.,  		1., 		1./2. };
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void N_Dot(int S, double *N, int M, double *R, double **C, double *m, double *r, double **kN, int Nth);
void R_Dot(int S, double *N, int M, double *R, double **C, double *K, double *D, double **kR, int Nth);
double RKF(int S, double *N, int M, double *R, double **C, double *m, double *r, double *K, double *D, double **kN, double **kR, double dt, double *tempN, double *tempR);
void ConvTest(int S,int M,double *dN,double *dR,double *K,int *ConvN,int *ConvR,double dt);
//
int main(int argc, char *argv[])
{
	delta_t		   =		pow(10., -3.);
	th_conv		   =		pow(10., -8.);
	zero_threshold =		pow(10., -9.);
	aetol		   =	    pow(10., -8.); // absolute error tolerance of RKF

	if(argc!=10){
		printf("Arg. Error!\n1. S\n2. M\n3. ave C\n4. sig C\n5. sig M\n6. aveK\n7. sig K\n8. sig D\n9. Ens\n");
		return 0;
	}
	int S = atoi(argv[1]);
	int M = atoi(argv[2]);

	double aveC = atof(argv[3]);
	double sigmaC = atof(argv[4]);

	double sigmaM = atof(argv[5]);
	double aveK = atof(argv[6]);
	double sigmaK = atof(argv[7]);
	double sigmaD = atof(argv[8]);

	int Ensemble= atoi(argv[9]);

	long int Seed;

	///	gsl package to make random number
	const gsl_rng_type *MT;
	gsl_rng *random;
	gsl_rng_env_setup();
	MT=gsl_rng_mt19937;
	random=gsl_rng_alloc(MT);
	gsl_rng_set(random,Ensemble);
	////
	int i,j;

	double* N = (double*)calloc(S,sizeof(double));
	double* R = (double*)calloc(M,sizeof(double));

	double* dN = (double*)calloc(S,sizeof(double));
	double* dR = (double*)calloc(M,sizeof(double));

	double** kN = (double**)calloc(6,sizeof(double*));
	double** kR = (double**)calloc(6,sizeof(double*));
	for(i=0;i<6;i++){
		kN[i] = (double*)calloc(S,sizeof(double));
		kR[i] = (double*)calloc(M,sizeof(double));
	}

	double* m = (double*)calloc(S,sizeof(double));
	double* r = (double*)calloc(S,sizeof(double));
	double** C = (double**)calloc(S,sizeof(double*));

	for(i=0;i<S;i++)
		C[i] = (double*)calloc(M,sizeof(double));

	double* K = (double*)calloc(M,sizeof(double));
	double* D = (double*)calloc(M,sizeof(double));


	double* tempN = (double*)calloc(S,sizeof(double));
	double* tempR = (double*)calloc(M,sizeof(double));

	double aveN2,aveR2;

	//	Files
	char FILE_Average_Dynamics[100];
	sprintf(FILE_Average_Dynamics,"./Results/%d_ConvTest_sigC%.2lf.txt",Ensemble,sigmaC);
	FILE *fp_ave_dynamics;

	char FILE_Results_N[100],FILE_Results_R[100];
	sprintf(FILE_Results_N,"./Results/%d_Final_N_sigC%.2lf.txt",Ensemble,sigmaC);
	sprintf(FILE_Results_R,"./Results/%d_Final_R_sigC%.2lf.txt",Ensemble,sigmaC);
	FILE *fp_Nresults,*fp_Rresults;

	char FILE_Results[100];
	sprintf(FILE_Results,"./Results/%d_Results.txt",Ensemble);
	FILE *fp_results;
	//
	clock_t begin,inter,end;
	double time_cost;

	int cntN,cntR;
	double aveN,aveR;

	int Step=pow(10,6);
	long int cnt=1;
	double t=0.;
	double dt=delta_t;
	int ConvN=0,ConvR=0;
	//
	//
	FILE *fp_params;
	char FILE_Params[100];

	int trial=0;
	// Time evolving
Re: // When the system is stucked in limit cycle...
	trial++;

	Seed=time(NULL);
	
	//
	printf("Ensemble: %d\t\tSeed: %ld\t\tsigC: %.3lf\t\t",Ensemble,Seed,sigmaC);
	gsl_rng_set(random,Seed);
	//
	begin=clock();
	//
	time_cost=0.;
	cnt=1;
	t=0.;
	dt=delta_t;
	ConvN=ConvR=0;
	//
	for(i=0;i<6;i++){
		memset(kN[i],0.,S*sizeof(double));
		memset(kR[i],0.,M*sizeof(double));
	}
	memset(tempN,0.,S*sizeof(double));
	memset(tempR,0.,M*sizeof(double));
	//
	//
	for(i=0;i<S;i++){
		N[i]=gsl_rng_uniform_pos(random) + 0.5;

		r[i]=aver;
		m[i]=gsl_ran_gaussian( random, sigmaM ) + aveM;
		for(j=0;j<M;j++)
			C[i][j]=gsl_ran_gaussian( random, sigmaC/sqrt(M) ) + aveC/M;
	}
	for(i=0;i<M;i++){
		R[i]=gsl_rng_uniform_pos(random) + 0.5;
		K[i]=gsl_ran_gaussian( random, sigmaK ) + aveK; // Gaussian dist.
		D[i]=gsl_ran_gaussian( random, sigmaD ) + aveD; // Gaussian dist.
	}
	//
	//	//
	while(1){
		dt=RKF(S,N,M,R,C,m,r,K,D,kN,kR,dt,tempN,tempR);
		//		Print average values
		if(cnt % Step == 0){
			//		Convergence test
			ConvTest(S,M,dN,dR,K,&ConvN,&ConvR,dt);
			if(ConvN == S && ConvR == M){ break; }
			else{
				inter=clock();
				time_cost=(double)(inter-begin)/CLOCKS_PER_SEC;
				//
				if(cnt>=2e8)
					goto Re;
				else{
					// N
					for(aveN=0.,i=0;i<S;i++)
						aveN+=N[i];
					aveN/=S;
					// R
					for(aveR=0.,i=0;i<M;i++)
						aveR+=R[i];
					aveR/=M;
					//
					fp_ave_dynamics = fopen(FILE_Average_Dynamics,"w");
					fprintf(fp_ave_dynamics,"%.5le\t%.5le\t\t\t%d\t\t%d\n\nCount: %ld\t\tS: %.5le\tM: %.5le\tH: %.5le\n\nTrial: %d\nSeed : %ld",t,dt,ConvN,ConvR,cnt,time_cost,time_cost/60.,time_cost/3600.,trial,Seed);
					fclose(fp_ave_dynamics);
					//
				}
			}
		}

		memset(dN,0.,S*sizeof(double));
		for(i=0;i<S;i++){
			if(N[i]<zero_threshold){N[i]=0.;}
			else{
				for(j=0;j<5;j++){dN[i]+=b4[j]*kN[j][i]*dt;}
				N[i]+=dN[i];
			}
		}
		memset(dR,0.,M*sizeof(double));
		for(i=0;i<M;i++){
			if(R[i]<zero_threshold){R[i]=0.;}
			else{
				for(j=0;j<5;j++){dR[i]+=b4[j]*kR[j][i]*dt;}
				R[i]+=dR[i];
			}
		}
		////////////////////////////////////
		t+=dt;
		cnt++;
	}

	end=clock();
	time_cost=(double)(end-begin)/CLOCKS_PER_SEC;


	fp_Nresults = fopen(FILE_Results_N,"w");
	for(i=0;i<S;i++)
		fprintf(fp_Nresults,"%.12le\t\t%.12le\n",N[i],kN[0][i]);
	fclose(fp_Nresults);

	fp_Rresults = fopen(FILE_Results_R,"w");
	for(i=0;i<M;i++)
		fprintf(fp_Rresults,"%.12le\t\t%.12le\n",R[i],kR[0][i]);
	fclose(fp_Rresults);

	//	N
	for(cntN=0,aveN=aveN2=0.,i=0;i<S;i++){
		if(N[i]!=0.0){
			aveN+=N[i];
			aveN2+=N[i]*N[i];
			cntN++;
		}
	}
	//	R
	for(cntR=0,aveR=aveR2=0.,i=0;i<M;i++){
		if(R[i]!=0.0){
			aveR+=R[i];
			aveR2+=R[i]*R[i];
			cntR++;
		}
	}
	//
	fp_results = fopen(FILE_Results,"a");
	fprintf(fp_results,"%.3le\t%d\t%.16le\t%.16le\t\t%d\t%.16le\t%.16le\n",sigmaC,cntN,aveN,aveN2,cntR,aveR,aveR2);
	fclose(fp_results);


	//
	gsl_rng_free(random);

	free(N);
	free(R);

	free(m);
	free(K);
	free(D);

	for(i=0;i<6;i++){
		free(kN[i]);
		free(kR[i]);
	}

	free(kN);
	free(kR);

	for(i=0;i<S;i++)
		free(C[i]);
	free(C);

	free(tempN);
	free(tempR);

	free(dN);
	free(dR);

	
	sprintf(FILE_Params,"./Results/%d_Params.txt",Ensemble);
	fp_params = fopen(FILE_Params,"a");
	fprintf(fp_params,"Ensemble: %d\t\tSigC: %.3lf\t\tSeed: %ld\t\t\tTrial: %d\n",Ensemble,sigmaC,Seed,trial);
	fclose(fp_params);

	return 0;
}
///////////
void N_Dot(int S, double *N, int M, double *R, double **C, double *m, double *r, double **kN, int Nth){
	int i,j;
	double x;

#pragma omp parallel for schedule(dynamic) collapse(1) private(i,j) shared(N,R,C) num_threads(omp_thrd)
	for(i=0;i<S;i++){
		for(x=0.,j=0;j<M;j++)
			x+=C[i][j]*R[j];
		kN[Nth][i] = N[i]*( x - m[i] - aver*N[i] );
		//		kN[Nth][i] = N[i]*( x - m[i] - r[i]*N[i] );
	}
}
void R_Dot(int S, double *N, int M, double *R, double **C, double *K, double *D, double **kR, int Nth){
	int i,j;
	double x;

#pragma omp parallel for schedule(dynamic) collapse(1) private(i,j) shared(N,R,C) num_threads(omp_thrd)
	for(j=0;j<M;j++){
		x=0.;
		for(i=0;i<S;i++){
			x += N[i]*C[i][j];
		}
		kR[Nth][j] = K[j] - (D[j] + x)*R[j];
	}
}
double RKF(int S, double *N, int M, double *R, double **C, double *m, double *r, double *K, double *D, double **kN, double **kR, double dt, double *tempN, double *tempR){
	int i,j;
	int Nth;
	double sum;
	double delta;
	double err_N,err_R;
	double s=0.;
	double Error=1.;


	N_Dot(S,N,M,R,C,m,r,kN,0);
	R_Dot(S,N,M,R,C,K,D,kR,0);

	while(Error>aetol){
		//Reject_and_Recalc:
		for(Nth=1;Nth<6;Nth++){
			//	N
#pragma omp parallel for schedule(dynamic) collapse(1) private(i,j) shared(a,kN) num_threads(omp_thrd)
			for(i=0;i<S;i++){
				sum=0.;
				for(j=0;j<Nth;j++)
					sum+=a[Nth-1][j]*kN[j][i];
				tempN[i]=N[i] + dt*sum;
			}
			//	R
#pragma omp parallel for schedule(dynamic) collapse(1) private(i,j) shared(a,kR) num_threads(omp_thrd)
			for(i=0;i<M;i++){
				sum=0.;
				for(j=0;j<Nth;j++)
					sum+=a[Nth-1][j]*kR[j][i];
				tempR[i]=R[i] + dt*sum;
			}

			N_Dot(S,tempN,M,tempR,C,m,r,kN,Nth);
			R_Dot(S,tempN,M,tempR,C,K,D,kR,Nth);
		}

		err_N=0.;
		for(i=0;i<S;i++){
			delta=0.;
			for(j=0;j<6;j++)
				delta+=delta_b[j]*kN[j][i];
			err_N += delta*delta;
		}

		err_R=0.;
		for(i=0;i<M;i++){
			delta=0.;
			for(j=0;j<6;j++)
				delta+=delta_b[j]*kR[j][i];
			err_R += delta*delta;
		}

		Error=sqrt( (err_N + err_R) / (S+M) )*dt;
		if(Error!=0.){s=0.87*pow(aetol/Error,0.2);}
		else{s=1.;}
		dt=s*dt;
	}

	return dt;
}
void ConvTest(int S,int M,double *dN,double *dR,double *K,int *ConvN,int *ConvR,double dt){
	int i;
	double speed_conv=th_conv*dt;
	for(*ConvN=0,i=0;i<S;i++)
		*ConvN+=(fabs(dN[i]) < speed_conv);
	for(*ConvR=0,i=0;i<M;i++)
		*ConvR+=(fabs(dR[i]) < speed_conv);
}


