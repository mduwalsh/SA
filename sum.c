#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h> 
#include<stdbool.h>

#define DEBUG 1
#if DEBUG
#define __USE_GNU
#include <fenv.h>       /* enable floating exceptions */
unsigned long Where;    // debugging counter
#endif

/* Parameters for plotting */

double N[]     = {4, 8};
double B[]     = {8, 16, 32};
double C[]     = {0.5};
double Alpha[] = {1};
double Beta[]  = {1.0};
double Delta[] = {0.5};
double Gamma[] = {1};
double Eta[]   = {1};
double K[]     = {2, 5, 10, 20};
double X0[]    = {1};
double Omega[] = {0.25, 0.5, 0.75};

#define XAXIS "B"                   // horizontal axis in a plot in a graph
#define YAXIS "cost" // options are: effort, payoff, threshold, aggressiveness, cost, punishment
#define VBAR  "N"
#define PLOTS "K"                   // variation in plots in a graph
#define MULTIPLOT_LAYOUT "2,3"      // layout of multiple plots in a graph

double *P, *Q, *R, *S, *T, *U, *V, *W, *X, *Y, *Z;
char Pa[3], Qa[3], Ra[3], Sa[3], Ta[3], Ua[3], Va[3], Wa[3], Xa[3], Ya[3], Za[3];
int Ps, Qss, Rs, Ss, Ts, Us, Vs, Ws, Xs, Ys, Zs;

void set_param();
void prep_infile(char *fname, char *appnd, double ri, double si, double ti, double ui, double vi, double wi, double xi, double yi, double zi, double qi, double pi);  


/* Parameters for plotting */
#define EPS 0                     // 1: eps image; 0: png image

#define HIST_TYPE "cluster"
#define BOXWIDTH 0.75


#define PREC "7.4" /* printing precision */
  
// write data in d array to file fp in n columns
void write_delta_ef(FILE *fp, double *d, char *ba, double bi, int n, double gma)
{
  int i;
  double X = 0.0;
  fprintf(fp, "%s=%.0f ", ba, bi);
  for(i = 0; i < n; i++){
    X += pow(d[i], 1.0/gma);
    fprintf(fp, "%"PREC"lf ", d[i]);
  }
  //fprintf(fp, "%"PREC"lf ", pow(X, gma));     // total effort
  fprintf(fp, "\n");
}

// calc averages of data in columns 
int calc_data(FILE *fp, double *d, int n, bool toNormal)
{
  int pos, nc = 0; 
  double val, sum = 0.0;
  char line[1024], *str;
  while( fgets(line, 1024, fp) !=NULL ) {
      // Just search for the latest line, do nothing in the loop
  } 
  // just read first line and do nothing
  str = line;
  while(sscanf(str, "%lf%n", &val, &pos)==1){
    d[nc] = val;
    if(toNormal){
      sum += d[nc];
    }
    str += pos;
    nc++;   
  }
  if(nc == n){
    if(toNormal){          // normalization of values
      while(nc--){
	d[nc] /= sum;
      }
    }
    return 1;
  }
  else 
    return 0;
}

void hist_sum(FILE *gp, char *edata, int n, double gma, char *title, char *outputfile)
{
  int i;
  char data[300]; char str[20];
  //int bs  = sizeof(B)/sizeof(double); // no. of Bs in array corresponds to no. of clusters
  
  if(EPS)
    fprintf(gp, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
  else
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,12\" \n");

  fprintf(gp, "set output '%s' \n", outputfile);
  fprintf(gp, "set key autotitle columnheader\n");
  //fprintf(gp, "unset xtics \n");
  fprintf(gp, "set key outside \n");
  fprintf(gp, "set style data histogram \n");
  fprintf(gp, "set style histogram %s title offset 0, -1 \n", HIST_TYPE);
  fprintf(gp, "set style fill solid border -1\n");
  fprintf(gp, "set boxwidth %lf relative\n", BOXWIDTH);
  //fprintf(gp, "set tmargin %lf \n", T_MARGIN);
  fprintf(gp, "set bmargin 2.0 \n");
  //fprintf(gp, "set lmargin %lf \n", L_MARGIN);
  //fprintf(gp, "set rmargin %lf \n", R_MARGIN);
  fprintf(gp, "set title ''\n");  
  
  if(strcmp(PLOTS, "") == 0){
    fprintf(gp, "set multiplot layout %s title '%s' \n", MULTIPLOT_LAYOUT, title);    // set subplots layout
    //fprintf(gp, "set xlabel '' \n");
    fprintf(gp, "set ylabel 'efforts' \n"); 
    fprintf(gp, "set auto x \n");  
    
    if((int)(gma*1000) != 1000){
      fprintf(gp, "set ytics %lf nomirror \n", 0.5);    
    }
    else{ 
      //fprintf(gp, "set ytics %lf mirror \n", E_YTICS);   
    }
    if((int)(gma*1000) != 1000){
      fprintf(gp, "set logscale y2 \n");
      fprintf(gp, "set y2tics nomirror \n");    
      fprintf(gp, "set y2range %s\n", "[0:]");
    }
    
    fprintf(gp, "set style arrow 1 nohead lw 2 \n");
    sprintf(data, "%sehist.dat", edata);
    fprintf(gp, "plot \\\n");  
    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle  \\\n", n+2, data);   
    fprintf(gp, "\n");
    
    fprintf(gp, "set ylabel 'threshold' \n");
    fprintf(gp, "unset y2tics \n");  
    //fprintf(gp, "unset autoscale y \n");
    fprintf(gp, "set yrange [0:] \n");
    fprintf(gp, "set ytics auto mirror\n");  
    sprintf(data, "%sdxihist.dat", edata);
    fprintf(gp, "plot \\\n");   
    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+2, data);    
    fprintf(gp, "\n");  
    
    fprintf(gp, "set ylabel 'aggressiveness' \n");
    fprintf(gp, "unset y2tics \n");  
    //fprintf(gp, "unset autoscale y \n");
    fprintf(gp, "set yrange [0:] \n");
    fprintf(gp, "set ytics auto mirror\n"); 
    sprintf(data, "%sdsihist.dat", edata);
    fprintf(gp, "plot \\\n");   
    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) title columnheader \\\n", n+2, data);    
    fprintf(gp, "\n"); 
    
    fprintf(gp, "set ylabel 'payoffs' \n");
    fprintf(gp, "unset y2tics \n");  
    //fprintf(gp, "unset autoscale y \n");
    fprintf(gp, "set yrange [0:] \n");
    fprintf(gp, "set ytics 0.1 mirror\n");  
    sprintf(data, "%sfhist.dat", edata);
    fprintf(gp, "plot \\\n");   
    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+2, data);    
    fprintf(gp, "\n"); 
    
    fprintf(gp, "set ylabel 'cost' \n");
    fprintf(gp, "unset y2tics \n");  
    //fprintf(gp, "set autoscale y \n");
    fprintf(gp, "set yrange [0:] \n");
    fprintf(gp, "set ytics auto mirror\n");  
    sprintf(data, "%spunihist.dat", edata);
    fprintf(gp, "plot \\\n");   
    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+2, data);    
    fprintf(gp, "\n"); 
    
    fprintf(gp, "set ylabel 'punishment' \n");
    fprintf(gp, "unset y2tics \n");  
    //fprintf(gp, "set autoscale y \n");
    fprintf(gp, "set yrange [0:] \n");
    fprintf(gp, "set ytics auto mirror\n");  
    sprintf(data, "%spunjhist.dat", edata);
    fprintf(gp, "plot \\\n");   
    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+2, data);    
    fprintf(gp, "\n"); 
  }
  else{
    
    fprintf(gp, "set multiplot layout 1, %d title '%s' \n", Qss, title);    // set subplots layout
    // if n is set as xaxis, set n to be the highest number of individuals in a group used in simulation
    if(strcmp(Pa, "n") == 0){
      for(n = 0, i = 0; i < Ps; i++){
	if(n < P[i]) n = P[i];
      }
    }
    
    fprintf(gp, "set auto x \n");
    
    if((int)(gma*1000) != 1000){
      fprintf(gp, "set ytics %lf nomirror \n", 0.5);    
    }
    else{ 
      //fprintf(gp, "set ytics %lf mirror \n", E_YTICS);   
    }
    if((int)(gma*1000) != 1000){
      fprintf(gp, "set logscale y2 \n");
      fprintf(gp, "set y2tics nomirror \n");    
      fprintf(gp, "set y2range %s\n", "[0:]");
    }
    
    fprintf(gp, "set style arrow 1 nohead lw 2 \n");
    
    if(strcmp(YAXIS, "effort") == 0) sprintf(str, "ehist.dat");
    else if(strcmp(YAXIS, "payoff") == 0) sprintf(str, "fhist.dat");
    else if(strcmp(YAXIS, "threshold") == 0) sprintf(str, "dxihist.dat");
    else if(strcmp(YAXIS, "aggressiveness") == 0) sprintf(str, "dsihist.dat");
    else if(strcmp(YAXIS, "cost") == 0) sprintf(str, "punihist.dat");
    else if(strcmp(YAXIS, "punishment") == 0) sprintf(str, "punjhist.dat");
    
    for(i = 0; i < Qss; i++){
      if( strcmp(Qa, "n") == 0){
	n = Q[i];
      }
      if(strcmp(Qa, "ga") == 0) gma = Q[i];
      fprintf(gp, "set ylabel '%s, %s=%g' \n", YAXIS, Qa, Q[i]); 
    
      sprintf(data, "%s%s%.2f%s", edata, Qa, Q[i], str);
      fprintf(gp, "plot \\\n");   
      fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+2, data);    
      fprintf(gp, "\n"); 
    }
  }
  
  fprintf(gp, "unset multiplot \n");
}

int main(int argc, char **argv)
{  
  char infile[200], edata[200], fdata[200], pi_gdata[200], dxidata[200], dsidata[200], punjdata[200], punidata[200];
  double *d, val;  
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8;
  char str[100], title[100];    

  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode  
  
  int j;
  int pi, qi, ri, si, ti, ui, vi, wi, xi, yi, zi, ni; 
  double ga; 
  
  set_param();
  for(ri = 0; ri < Rs; ri++){
    if( strcmp(Ra, "n") == 0){
      d = malloc(R[ri]*sizeof(double));  // to hold data after each read of columns in data files      
    }    
      for(si = 0; si < Ss; si++){
	for(ti = 0; ti < Ts; ti++){
	  for(ui = 0; ui < Us; ui++){    
	    for(vi = 0; vi < Vs; vi++){
	      if(strcmp(Va, "ga") == 0) ga = V[vi];
	      for(wi = 0; wi < Ws; wi++){
		if(strcmp(Wa, "ga") == 0) ga = W[wi];
		for(xi = 0; xi < Xs; xi++){
		  if(strcmp(Wa, "ga") == 0) ga = W[wi];
		  for(yi = 0; yi < Ys; yi++){
		    for(zi = 0; zi < Zs; zi++){
		      for(qi = 0; qi < Qss; qi++){
			if(strcmp(Qa, "ga") == 0) ga = Q[qi];
			sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2fehist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);   // data for histogram of efforts
			sprintf(fdata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2ffhist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);   // data for histogram of payoff
			sprintf(pi_gdata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2fpi_ghist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);   // data for histogram of payoff
			sprintf(dxidata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2fdxihist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);   // data for threshold
			sprintf(dsidata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2fdsihist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);   // data for aggressiveness
			sprintf(punidata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2fpunihist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);
			sprintf(punjdata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2fpunjhist.dat", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);
			fp2 = fopen(edata, "w");
			fp3 = fopen(fdata, "w");
			fp4 = fopen(pi_gdata, "w");
			fp5 = fopen(dxidata, "w");
			fp6 = fopen(dsidata, "w");
			fp7 = fopen(punidata, "w");
			fp8 = fopen(punjdata, "w");
			fprintf(fp2, "%s's ", Pa);
			fprintf(fp3, "%s's ", Pa);
			fprintf(fp4, "%s's ", Pa);
			fprintf(fp5, "%s's ", Pa);
			fprintf(fp6, "%s's ", Pa);
			fprintf(fp7, "%s's ", Pa);
			fprintf(fp8, "%s's ", Pa);
			if( strcmp(Qa, "n") == 0){
			  d = malloc(Q[qi]*sizeof(double));  // to hold data after each read of columns in data files
			  ni = (int)Q[qi];
			  for(j = 0; j < Q[qi]; j++){
			    fprintf(fp2, "%d  ", j+1); fprintf(fp3, "%d  ", j+1); fprintf(fp5, "%d  ", j+1); fprintf(fp6, "%d  ", j+1); fprintf(fp7, "%d  ", j+1); fprintf(fp8, "%d  ", j+1);  // adding headers
			  }
			} 
			if( strcmp(Ra, "n") == 0){
			  ni = (int)R[ri];
			  for(j = 0; j < R[ri]; j++){
			    fprintf(fp2, "%d  ", j+1); fprintf(fp3, "%d  ", j+1); fprintf(fp5, "%d  ", j+1); fprintf(fp6, "%d  ", j+1); fprintf(fp7, "%d  ", j+1); fprintf(fp8, "%d  ", j+1);  // adding headers
			  }
			}
			
			for(pi = 0; pi < Ps; pi++){
			  if(strcmp(Pa, "ga") == 0) ga = P[pi];
			  if( strcmp(Pa, "n") == 0 ){
			    d = malloc(P[pi]*sizeof(double));  // to hold data after each read of columns in data files
			    ni = (int)P[pi];
			    if(pi == 0){
			      for(j = 0; j < P[pi]; j++){
				fprintf(fp2, "%d  ", j+1); fprintf(fp3, "%d  ", j+1); fprintf(fp5, "%d  ", j+1); fprintf(fp6, "%d  ", j+1); fprintf(fp7, "%d  ", j+1); fprintf(fp8, "%d  ", j+1);  // adding headers
			      }
			      fprintf(fp2, "X \n"); fprintf(fp3, "X \n"); fprintf(fp4, "pi_g \n"); fprintf(fp5, "X \n"); fprintf(fp6, "X \n"); fprintf(fp7, "X\n"); fprintf(fp8, "X\n");
			    }   
			  }
			  // files to be read for efforts for set of paramters
			  prep_infile(infile, "esum.dat", R[ri], S[si], T[ti], U[ui], V[vi], W[wi], X[xi], Y[yi], Z[zi], Q[qi], P[pi]);
			  //sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fesum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);  
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize			    
			    fclose(fp1);
			    write_delta_ef(fp2, d, Pa, P[pi], ni, ga);   // write data to histogram data file for efforts // j is passed as delta for linear spacing    
			  }
			  else{
			    fclose(fp1);
			    printf("error in data read '%s': column count does not match!!\n", infile);
			    exit(1);
			  }
			  
			// for fertilities
			// files to be read for fertilities for set of paramters
			  prep_infile(infile, "fsum.dat", R[ri], S[si], T[ti], U[ui], V[vi], W[wi], X[xi], Y[yi], Z[zi], Q[qi], P[pi]);
			  //sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2ffsum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(calc_data(fp1, d, ni, 1)){          // calculate data and normalize fertitilties
			    fclose(fp1);
			    write_delta_ef(fp3, d, Pa, P[pi], ni, ga);   // write data to histogram data file for fertilities // j is passed as delta for linear spacing
			  }
			  else{
			    fclose(fp1);
			    printf("error in data read '%s': column count does not match!!\n", infile);
			    exit(1);
			  }
			  
			  // for group payoff
			  /*sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fpi_gsum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi]);
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(fscanf(fp1, "%lf", &val) == 1){
			    fprintf(fp4, "b=%.2f %.4lf\n", B[bi], val);
			  }*/
			  
			  // files to be read for threshold for set of paramters
			  prep_infile(infile, "dxisum.dat", R[ri], S[si], T[ti], U[ui], V[vi], W[wi], X[xi], Y[yi], Z[zi], Q[qi], P[pi]);
			  //sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fdxisum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);  
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize
			    fclose(fp1);
			    write_delta_ef(fp5, d, Pa, P[pi], ni, ga);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			  }
			  else{
			    fclose(fp1);
			    printf("error in data read '%s': column count does not match!!\n", infile);
			    exit(1);
			  }
			
			  // files to be read for aggressiveness for set of paramters
			  prep_infile(infile, "dsisum.dat", R[ri], S[si], T[ti], U[ui], V[vi], W[wi], X[xi], Y[yi], Z[zi], Q[qi], P[pi]);
			  //sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fdsisum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);  
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize
			    fclose(fp1);
			    write_delta_ef(fp6, d, Pa, P[pi], ni, ga);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			  }
			  else{
			    fclose(fp1);
			    printf("error in data read '%s': column count does not match!!\n", infile);
			    exit(1);
			  }
			  
			  // for group effort
			  /*sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fx_gsum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi]);
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(fscanf(fp1, "%lf", &val) == 1){
			    fprintf(fp7, "b=%.2f %.4lf\n", B[bi], val);
			  }*/
			  
			  // for punishment
			  prep_infile(infile, "punisum.dat", R[ri], S[si], T[ti], U[ui], V[vi], W[wi], X[xi], Y[yi], Z[zi], Q[qi], P[pi]);
			  //sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunisum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize
			    fclose(fp1);
			    write_delta_ef(fp7, d, Pa, P[pi], ni, ga);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			  }
			  else{
			    fclose(fp1);
			    printf("error in data read '%s': column count does not match!!\n", infile);
			    exit(1);
			  }
			  // for punishment
			  prep_infile(infile, "punjsum.dat", R[ri], S[si], T[ti], U[ui], V[vi], W[wi], X[xi], Y[yi], Z[zi], Q[qi], P[pi]);
			  //sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunjsum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
			  if(!(fp1 = fopen(infile, "r"))){
			    printf("Could not open file '%s'",infile);
			    exit(1);
			  }
			  if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize
			    fclose(fp1);
			    write_delta_ef(fp8, d, Pa, P[pi], ni, ga);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			  }
			  else{
			    fclose(fp1);
			    printf("error in data read '%s': column count does not match!!\n", infile);
			    exit(1);
			  }
			  
			  if( strcmp(Pa, "n") == 0 ){
			    free(d);
			  }
			}    
			fclose(fp2); fclose(fp3); fclose(fp4); fclose(fp5); fclose(fp6); fclose(fp7); fclose(fp8);
			if( strcmp(Qa, "n") == 0){
			  free(d);
			}  
		      }
		    }
		}
	      }
	    }
	  }	 
	}
      }
    }
    if( strcmp(Ra, "n") == 0){
      free(d);
    }     
  }  
   // reading data and prepare histogram data file for multiple values of b and delta but same other parameters completed
   
  // plotting summary histograms all at a time // the one with cluster base is removed from for loop (like for B in original case)
  ni = 0; ga = 1;
    for(ri = 0; ri < Rs; ri++){      
      if( strcmp(Ra, "n") == 0){
	ni = R[ri];
      }
      for(si = 0; si < Ss; si++){
	for(ti = 0; ti < Ts; ti++){
	  for(ui = 0; ui < Us; ui++){    
	    for(vi = 0; vi < Vs; vi++){
	      if(strcmp(Va, "ga") == 0) ga = V[vi];
	      for(wi = 0; wi < Ws; wi++){
		if(strcmp(Wa, "ga") == 0) ga = W[wi];
		for(xi = 0; xi < Xs; xi++){
		  if(strcmp(Xa, "ga") == 0) ga = X[xi];
		  for(yi = 0; yi < Ys; yi++){
		    for(zi = 0; zi < Zs; zi++){
		      
		      if(strcmp(PLOTS, "") == 0){
			for(qi = 0; qi < Qss; qi++){
			  if( strcmp(Qa, "n") == 0){
			    ni = Q[qi];
			  }
			  if(strcmp(Qa, "ga") == 0) ga = Q[qi];
			  sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);   // data for histogram of efforts
			  // plot histograms      
			  sprintf(title, "%s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);
			  if(EPS)
			    sprintf(str, "sum_%s%g%s%g%s%.g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.eps",Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);
			  else
			    sprintf(str, "sum_%s%g%s%g%s%.g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.png",Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi], Qa, Q[qi]);
			  hist_sum(gp, edata, ni, ga, title, str);
			  // histograms for efforts and payoffs for one set of n, c and X0 compeleted
			}
		      }
		      else{
			sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi]);   // data for histogram of efforts
			// plot histograms      
			sprintf(title, "%s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g", Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi]);
			if(EPS)
			  sprintf(str, "sum_%s%g%s%g%s%.g%s%g%s%g%s%g%s%g%s%g%s%g.eps",Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi]);
			else
			  sprintf(str, "sum_%s%g%s%g%s%.g%s%g%s%g%s%g%s%g%s%g%s%g.png",Ra, R[ri], Sa, S[si], Ta, T[ti], Ua, U[ui], Va, V[vi], Wa, W[wi], Xa, X[xi], Ya, Y[yi], Za, Z[zi]);
			hist_sum(gp, edata, ni, ga, title, str);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
  fflush(gp); 
  pclose(gp);
  
  return 0;
}

void set_param(){
  int i, h, p;
  double *Param[11];
  char ParamStr[11][10];
  char ParamA[11][3];
  int ParamS[11];
  int ns, bs, bts, cs, als, dls, xs, gs, ks, es, os;
  
  ns  = sizeof(N)/sizeof(double);
  bs  = sizeof(B)/sizeof(double);
  cs  = sizeof(C)/sizeof(double);
  bts = sizeof(Beta)/sizeof(double);
  als = sizeof(Alpha)/sizeof(double);
  dls = sizeof(Delta)/sizeof(double);
  xs  = sizeof(X0)/sizeof(double);
  ks  = sizeof(K)/sizeof(double);
  es  = sizeof(Eta)/sizeof(double);
  gs  = sizeof(Gamma)/sizeof(double);	
  os  = sizeof(Omega)/sizeof(double);
  
  Param[0] = N; 
  Param[1] = B; 
  Param[2] = C;
  Param[3] = Alpha;
  Param[4] = Beta;
  Param[5] = Delta;
  Param[6] = Gamma;
  Param[7] = Eta;
  Param[8] = K;
  Param[9] = X0;
  Param[10] = Omega;  
  ParamS[0] = ns; 
  ParamS[1] = bs; 
  ParamS[2] = cs;
  ParamS[3] = als;
  ParamS[4] = bts;
  ParamS[5] = dls;
  ParamS[6] = gs;
  ParamS[7] = es;
  ParamS[8] = ks;
  ParamS[9] = xs;
  ParamS[10] = os;  
  sprintf(ParamStr[0], "N");
  sprintf(ParamStr[1], "B");
  sprintf(ParamStr[2], "C");
  sprintf(ParamStr[3], "Alpha");
  sprintf(ParamStr[4], "Beta");
  sprintf(ParamStr[5], "Delta");
  sprintf(ParamStr[6], "Gamma");
  sprintf(ParamStr[7], "Eta");
  sprintf(ParamStr[8], "K");
  sprintf(ParamStr[9], "X0");
  sprintf(ParamStr[10], "Omega");
  sprintf(ParamA[0], "n");
  sprintf(ParamA[1], "b");
  sprintf(ParamA[2], "c");
  sprintf(ParamA[3], "a");
  sprintf(ParamA[4], "be");
  sprintf(ParamA[5], "dl");
  sprintf(ParamA[6], "ga");
  sprintf(ParamA[7], "e");
  sprintf(ParamA[8], "k");
  sprintf(ParamA[9], "xx");
  sprintf(ParamA[10], "w");
  
    
  for(i = 0; i < 11; i ++){
    if( strcmp(ParamStr[i], XAXIS) == 0){
      h = i;
      P = Param[i];
      sprintf(Pa,"%s", ParamA[i]);
      Ps = ParamS[i];
    }
    else if( strcmp(ParamStr[i], PLOTS) == 0 ){
      p = i;
      Q = Param[i];
      sprintf(Qa, "%s", ParamA[i]);
      Qss = ParamS[i];      
    }
  }
  if( strcmp(PLOTS, "") == 0 ){
    p = -1;
  }
  
  i = 0;
  if( i != h && i != p){
    R = Param[i];
    sprintf(Ra,"%s", ParamA[i]);
    Rs = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      R = Param[i];
      sprintf(Ra,"%s", ParamA[i]);
      Rs = ParamS[i];
    }
    else{
      i++;
      R = Param[i];
      sprintf(Ra,"%s", ParamA[i]);
      Rs = ParamS[i];      
    }
  }
  i++;
  
  if( i != h && i != p){
    S = Param[i];
    sprintf(Sa,"%s", ParamA[i]);
    Ss = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      S = Param[i];
      sprintf(Sa,"%s", ParamA[i]);
      Ss = ParamS[i];
    }
    else{
      i++;
      S = Param[i];
      sprintf(Sa,"%s", ParamA[i]);
      Ss = ParamS[i];      
    }
  }
  i++;
  
  if( i != h && i != p){
    T = Param[i];
    sprintf(Ta,"%s", ParamA[i]);
    Ts = ParamS[i]; 
  }
  else{
    i++;
    if (i != h && i != p){
      T = Param[i];
      sprintf(Ta,"%s", ParamA[i]);
      Ts = ParamS[i]; 
    }
    else{
      i++;
      T = Param[i];
      sprintf(Ta,"%s", ParamA[i]);    
      Ts = ParamS[i]; 
    }
  } 
  i++;
  
  if( i != h && i != p){
    U = Param[i];
    sprintf(Ua,"%s", ParamA[i]);
    Us = ParamS[i]; 
  }
  else{
    i++;
    if (i != h && i != p){
      U = Param[i];
      sprintf(Ua,"%s", ParamA[i]);
      Us = ParamS[i];
    }
    else{
      i++;
      U = Param[i];
      sprintf(Ua,"%s", ParamA[i]);  
      Us = ParamS[i];
    }
  }
  i++;
  
  if( i != h && i != p){
    V = Param[i];
    sprintf(Va,"%s", ParamA[i]);
    Vs = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      V = Param[i];
      sprintf(Va,"%s", ParamA[i]);
      Vs = ParamS[i];
    }
    else{
      i++;
      V = Param[i];
      sprintf(Va,"%s", ParamA[i]);
      Vs = ParamS[i];     
    }
  }
  i++;
  
  if( i != h && i != p){
    W = Param[i];
    sprintf(Wa,"%s", ParamA[i]);
    Ws = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      W = Param[i];
      sprintf(Wa,"%s", ParamA[i]);
      Ws = ParamS[i];
    }
    else{
      i++;
      W = Param[i];
      sprintf(Wa,"%s", ParamA[i]);      
      Ws = ParamS[i];
    }
  }
  i++;
  
  if( i != h && i != p){
    X = Param[i];
    sprintf(Xa,"%s", ParamA[i]);
    Xs = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      X = Param[i];
      sprintf(Xa,"%s", ParamA[i]);
      Xs = ParamS[i];
    }
    else{
      i++;
      X = Param[i];
      sprintf(Xa,"%s", ParamA[i]);   
      Xs = ParamS[i];
    }
  }
  i++;
  
  if( i != h && i != p){
    Y = Param[i];
    sprintf(Ya,"%s", ParamA[i]);
    Ys = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      Y = Param[i];
      sprintf(Ya,"%s", ParamA[i]);
      Ys = ParamS[i];
    }
    else{
      i++;
      Y = Param[i];
      sprintf(Ya,"%s", ParamA[i]);      
      Ys = ParamS[i];
    }
  }
  i++;
  
  if( i != h && i != p){
    Z = Param[i];
    sprintf(Za,"%s", ParamA[i]);
    Zs = ParamS[i];
  }
  else{
    i++;
    if (i != h && i != p){
      Z = Param[i];
      sprintf(Za,"%s", ParamA[i]);
      Zs = ParamS[i];
    }
    else{
      i++;
      Z = Param[i];
      sprintf(Za,"%s", ParamA[i]);   
      Zs = ParamS[i];
    }
  }
  i++;
  
  if(strcmp(PLOTS, "") == 0){
    printf("entered \n");
    Q = Param[i];
    sprintf(Qa,"%s", ParamA[i]);
    Qss = ParamS[i];
  }
  i++;  
  
}


/** Usage:
    compile : gcc -o hist hist.c
    run     : ./hist
**/
 

void prep_infile(char *fname, char *appnd, double ri, double si, double ti, double ui, double vi, double wi, double xi, double yi, double zi, double qi, double pi)
{
  char str[10];
  strcpy(str, "");
  strcpy(fname, "");
  // for n
  if(strcmp(Ra,"n") == 0){
    sprintf(str,"n%.0f", ri);
  }
  else if(strcmp(Pa,"n") == 0){
    sprintf(str,"n%.0f", pi);
  }
  else if(strcmp(Qa,"n") == 0){
    sprintf(str,"n%.0f", qi);
  }
  strcat(fname, str);
  
  // for b
  if(strcmp(Ra,"b") == 0){
    sprintf(str,"b%.2f", ri);
  }
  else if(strcmp(Sa, "b") == 0){
   sprintf(str,"b%.2f", si);
  }
  else if(strcmp(Pa,"b") == 0){
    sprintf(str,"b%.2f", pi);
  }
  else if(strcmp(Qa,"b") == 0){
    sprintf(str,"b%.2f", qi);
  }
  strcat(fname, str);
  
  // for c
  if(strcmp(Ra,"c") == 0){
    sprintf(str,"c%.2f", ri);
  }
  else if(strcmp(Sa, "c") == 0){
   sprintf(str,"c%.2f", si);
  }
  else if(strcmp(Ta, "c") == 0){
   sprintf(str,"c%.2f", ti);
  }
  else if(strcmp(Pa,"c") == 0){
    sprintf(str,"c%.2f", pi);
  }
  else if(strcmp(Qa,"c") == 0){
    sprintf(str,"c%.2f", qi);
  }
  strcat(fname, str);
  
  // for a
  if(strcmp(Sa, "a") == 0){
   sprintf(str,"a%.2f", si);
  }
  else if(strcmp(Ta, "a") == 0){
   sprintf(str,"a%.2f", ti);
  }
  if(strcmp(Ua,"a") == 0){
    sprintf(str,"a%.2f", ui);
  }
  else if(strcmp(Pa,"a") == 0){
    sprintf(str,"a%.2f", pi);
  }
  else if(strcmp(Qa,"a") == 0){
    sprintf(str,"a%.2f", qi);
  }
  strcat(fname, str);
  
  // for be
  if(strcmp(Ta, "be") == 0){
   sprintf(str,"b%.2f", ti);
  }
  else if(strcmp(Ua,"be") == 0){
    sprintf(str,"b%.2f", ui);
  }
  else if(strcmp(Va, "be") == 0){
   sprintf(str,"b%.2f", vi);
  }  
  else if(strcmp(Pa,"be") == 0){
    sprintf(str,"b%.2f", pi);
  }
  else if(strcmp(Qa,"be") == 0){
    sprintf(str,"b%.2f", qi);
  }
  strcat(fname, str);
  
  // for d
  if(strcmp(Ua,"dl") == 0){
    sprintf(str,"d%.2f", ui);
  }
  else if(strcmp(Va, "dl") == 0){
   sprintf(str,"d%.2f", vi);
  }  
  else if(strcmp(Wa, "dl") == 0){
   sprintf(str,"d%.2f", wi);
  }  
  else if(strcmp(Pa,"dl") == 0){
    sprintf(str,"d%.2f", pi);
  }
  else if(strcmp(Qa,"dl") == 0){
    sprintf(str,"d%.2f", qi);
  }
  strcat(fname, str);
  
  // for g
  if(strcmp(Va, "ga") == 0){
   sprintf(str,"g%.2f", vi);
  }  
  else if(strcmp(Wa, "ga") == 0){
   sprintf(str,"g%.2f", wi);
  }  
  else if(strcmp(Xa,"ga") == 0){
    sprintf(str,"g%.2f", xi);
  }  
  else if(strcmp(Pa,"ga") == 0){
    sprintf(str,"g%.2f", pi);
  }
  else if(strcmp(Qa,"ga") == 0){
    sprintf(str,"g%.2f", qi);
  }
  strcat(fname, str);
  
  // for e
  if(strcmp(Wa, "e") == 0){
   sprintf(str,"e%.3f", wi);
  }  
  else if(strcmp(Xa,"e") == 0){
    sprintf(str,"e%.3f", xi);
  }  
  else if(strcmp(Ya, "e") == 0){
   sprintf(str,"e%.3f", yi);
  }
  else if(strcmp(Pa,"e") == 0){
    sprintf(str,"e%.3f", pi);
  }
  else if(strcmp(Qa,"e") == 0){
    sprintf(str,"e%.3f", qi);
  }
  strcat(fname, str);
  
  // for k
  if(strcmp(Xa,"k") == 0){
    sprintf(str,"k%.0f", xi);
  }  
  else if(strcmp(Ya, "k") == 0){
   sprintf(str,"k%.0f", yi);
  }
  else if(strcmp(Za, "k") == 0){
   sprintf(str,"k%.0f", zi);
  }    
  else if(strcmp(Pa,"k") == 0){
    sprintf(str,"k%.0f", pi);
  }
  else if(strcmp(Qa,"k") == 0){
    sprintf(str,"k%.0f", qi);
  }
  strcat(fname, str);
  
  strcpy(str, "");
  // for x0
  //printf("x0 = %s\n", Ya);
  if(strcmp(Ya, "xx") == 0){
    //printf("entered as correct\n");
   sprintf(str,"x0%.2f", yi);
  }
  else if(strcmp(Za, "xx") == 0){
   sprintf(str,"x0%.2f", zi);
  }    
  else if(strcmp(Pa,"xx") == 0){
    sprintf(str,"x0%.2f", pi);
  }
  else if(strcmp(Qa,"xx") == 0){
    sprintf(str,"x0%.2f", qi);
  }
  strcat(fname, str);
  
  // for w
  if(strcmp(Za, "w") == 0){
   sprintf(str,"w%.2f", zi);
  }    
  else if(strcmp(Pa,"w") == 0){
    sprintf(str,"w%.2f", pi);
  }
  else if(strcmp(Qa,"w") == 0){
    sprintf(str,"w%.2f", qi);
  }
  strcat(fname, str);
  
  strcat(fname, appnd);
} 
  
 
