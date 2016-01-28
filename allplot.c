#include<stdio.h>
#include<stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

struct stat st = {0};

/* sets of parameters to be simulated on */
int Runs       = 9;
int N[]        = {4, 8};
double B[]     = {8, 16, 32};
double Alpha[] = {1};
double Beta[]  = {1.0};
double C[]     = {0.5, 1};
double Delta[] = {0.5, 1};
double X0[]    = {1};
double Gamma[] = {1};
double Eta[]   = {1};
double Omega[] = {0, 0.25, 0.5, 0.75, 1};
int    K[]     = {2, 5, 10, 20};
/* end of sets of parameters to be simulated on */


int main()
{
  int ni, bi, ci, bti, ali, dli, xi, gi, ki, ei, oi;
  int ns, bs, bts, cs, als, dls, xs, gs, ks, es, os;
  char scall[100];
  //char str[100];
  ns  = sizeof(N)/sizeof(int);
  bs  = sizeof(B)/sizeof(double);
  cs  = sizeof(C)/sizeof(double);
  bts = sizeof(Beta)/sizeof(double);
  als = sizeof(Alpha)/sizeof(double);
  dls = sizeof(Delta)/sizeof(double);
  xs  = sizeof(X0)/sizeof(double);
  ks  = sizeof(K)/sizeof(int);
  es  = sizeof(Eta)/sizeof(double);
  gs  = sizeof(Gamma)/sizeof(double);
  os  = sizeof(Omega)/sizeof(double);
	
  int i = 1;
  for(ni = 0; ni < ns; ni++)
    for(bi = 0; bi < bs; bi++)
      for(ci = 0; ci < cs; ci++)
	for(bti = 0; bti < bts; bti++)
	  for(ali = 0; ali < als; ali++)
	    for(dli = 0; dli < dls; dli++)
	      for(xi = 0; xi < xs; xi++)
		for(gi = 0; gi < gs; gi++)
		  for(ei = 0; ei < es; ei++)
		    for(oi = 0; oi < os; oi++)
		      for(ki = 0; ki < ks; ki++, i++){
		      //sprintf(str,"n%d.b%.2lf.c%.2lf.be%.2lfd%.2lfx0%.2lf", N[ni], B[bi], C[ci], Beta[bti], Delta[dli], X0[xi]);
		      /*if (stat(str, &st) == -1) {
			mkdir(str, 0700);
		      } */ 
			sprintf(scall, "./fplot %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", N[ni], C[ci], B[bi], Alpha[ali], Beta[bti], Delta[dli], Gamma[gi], X0[xi], Eta[ei], Omega[oi], K[ki], Runs);   
			if(!system(scall)) printf("Error while calling fplot!!\n");
		      //sprintf(scall, "mv *.png %s", str);  
		      //system(scall);          
		      }
  printf("\n%d\n", i);
  return 0;

}


/** Usage:
    compile : gcc -o allplot allplot.c
    run     : ./allplot
**/



