// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
double PE[]    = {1, 0, 0.0}; // {PE1, PE2, PE3}
int N[]        = {4, 8};
double B[]     = {8, 16, 32};
double Alpha[] = {1.0};
double Beta[]  = {1.0};
double C[]     = {0.5, 1};
double Delta[] = {0.5, 1};
double Gamma[] = {1};
double X0[]    = {1};
int    K[]     = {2, 5, 10, 20};
double Eta[]   = {1}; 
double Omega[] = {0, 0.25, 0.5, 0.75, 1};
double Discount = 0.2;
/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 9;       
int      G         = 1000;   
int      PROTOCOL  = 0;
double   T         = 10000;   
double   F0        = 1; 
double   Mu        = 0.25; 
double   Sigma     = 0.25; 
double   Sigma_B   = 0.0; 
double   Sigma_dxi = 0.25; 
double   Sigma_dsi = 0.25; 
double   S0        = 1;
double   Phi       = 5; 
double   E         = 1; 
double   e         = 0.1;

/* end of other basic parameters */



typedef struct config{ 
  int n;   
  double b;
  double c;
  double alpha;
  double beta;
  double gma;  
  double delta;  
  double x0;
  double eta;
  double omega;
  int k;
} Config;

void set_config(Config *cg, int n, double b, double c, double alpha, double beta, double delta, double gma, double x0, double eta, double omega, int k);
void create_config(Config *, int);
void create_script(int i, int n, double b, double c, double alpha, double beta, double delta, double gma, double x0, double eta, double omega, int k);

int main()
{
  int ni, bi, ci, bti, ali, dli, gi, xi, ei, ki, oi, sc;
  int ns, bs, bts, cs, als, dls, gs, xs, es, ks, os;
  char scall[50];
  ns  = sizeof(N)/sizeof(int);
  bs  = sizeof(B)/sizeof(double);
  cs  = sizeof(C)/sizeof(double);
  bts = sizeof(Beta)/sizeof(double);
  als = sizeof(Alpha)/sizeof(double);
  dls = sizeof(Delta)/sizeof(double);
  gs  = sizeof(Gamma)/sizeof(double);
  xs  = sizeof(X0)/sizeof(double);
  es  = sizeof(Eta)/sizeof(double);
  ks  = sizeof(K)/sizeof(int);
  os  = sizeof(Omega)/sizeof(double);
  int i = 1;
  Config cg;
  for(ni = 0; ni < ns; ni++)
    for(bi = 0; bi < bs; bi++)
      for(ci = 0; ci < cs; ci++)
	for(bti = 0; bti < bts; bti++)
	  for(ali = 0; ali < als; ali++)
	    for(dli = 0; dli < dls; dli++)
	      for(gi = 0; gi < gs; gi++)
		for(xi = 0; xi < xs; xi++)
		  for(ei = 0; ei < es; ei++)
		    for(ki = 0; ki < ks; ki++)
		      for(oi = 0; oi <os; oi++, i++){
			set_config(&cg, N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[dli], Gamma[gi], X0[xi], Eta[ei], Omega[oi], K[ki]);	      
			create_config(&cg, i);
			create_script(i, N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[dli], Gamma[gi], X0[xi], Eta[ei], Omega[oi], K[ki]);	   
			sprintf(scall, "qsub script.%d.sh", i);
			sc = system(scall);
			if(sc) printf("Error submitting jobs!!\n");
		      }
  printf("\n%d\n", i);
  return 0;
}

void set_config(Config *cg, int n, double b, double c, double alpha, double beta, double delta, double gma, double x0, double eta, double omega, int k)
{  
  cg->gma  = gma;  
  cg->n    = n;
  cg->b    = b;
  cg->c    = c;
  cg->alpha= alpha;
  cg->beta = beta;
  cg->delta= delta;
  cg->x0   = x0;
  cg->eta  = eta;
  cg->omega = omega;
  cg->k    = k;
}

void create_config(Config *con, int i)
{
  FILE *fp;
  char cfile[20];
  sprintf(cfile, "sa%d.config", i);
  if(!(fp = fopen(cfile, "w"))){
    printf("Error!! %s couldn't be created.\n", cfile);
    return;
  }
  fprintf(fp, "unsigned Seed      = %u; \n", Seed);
  fprintf(fp, "int      Runs      = %d; \n", Runs);
  fprintf(fp, "int      N         = %d; \n", con->n);
  fprintf(fp, "int      G         = %d; \n", G);
  fprintf(fp, "double   T         = %lf; \n", T);
  fprintf(fp, "\n");
  fprintf(fp, "double   PE[0]     = %lf; \n", PE[0]);
  fprintf(fp, "double   PE[1]     = %lf; \n", PE[1]);
  fprintf(fp, "double   PE[2]     = %lf; \n", PE[2]);
  fprintf(fp, "\n");
  fprintf(fp, "double   F0        = %lf; \n", F0);
  fprintf(fp, "double   B         = %lf; \n", con->b);
  fprintf(fp, "double   C         = %lf; \n", con->c);
  fprintf(fp, "double   X0        = %lf; \n", con->x0);
  fprintf(fp, "double   Discount  = %lf; \n", Discount);
  fprintf(fp, "double   Omega     = %lf; \n", con->omega);
  fprintf(fp, "\n");
  fprintf(fp, "double   Alpha     = %lf; \n", con->alpha);
  fprintf(fp, "double   Beta      = %lf; \n", con->beta);
  fprintf(fp, "double   Gamma     = %lf; \n", con->gma);
  fprintf(fp, "\n");
  fprintf(fp, "double   Delta     = %lf; \n", con->delta);
  fprintf(fp, "\n");
  fprintf(fp, "int      PROTOCOL  = %d; \n", PROTOCOL);
  fprintf(fp, "int      K         = %d; \n", con->k);
  fprintf(fp, "double   Eta       = %lf; \n", con->eta);  
  fprintf(fp, "\n");
  fprintf(fp, "double   E         = %lf; \n", E);
  fprintf(fp, "double   S0        = %lf; \n", S0);
  fprintf(fp, "double   Phi       = %lf; \n", Phi);
  fprintf(fp, "double   e         = %lf; \n", e);
  fprintf(fp, "\n");
  fprintf(fp, "double   Mu        = %lf; \n", Mu);
  fprintf(fp, "double   Sigma     = %lf; \n", Sigma);
  fprintf(fp, "double   Sigma_B   = %lf; \n", Sigma_B);
  fprintf(fp, "double   Sigma_dxi = %lf; \n", Sigma_dxi);
  fprintf(fp, "double   Sigma_dsi = %lf; \n", Sigma_dsi);  

  fclose(fp);
}

void create_script(int i, int n, double b, double c, double alpha, double beta, double delta, double gma, double x0, double eta, double omega, int k )
{
  char sfile[20];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "script.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N sa.%d.%.2lf.%.2lf.%.2lf.%.2lf%.2lf.%.2lf.%.2lf%.3lf%.2lf%d\n", n, b, c, alpha, beta, gma, delta, x0, eta, omega, k);
  fprintf(fp, "#PBS -l walltime=100:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./pun sa%d.config\n", i);
  fclose(fp);  
}


/** Usage:
    compile : gcc -o sajob sajob.c
    run     : ./sajob
**/

