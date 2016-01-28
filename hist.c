#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h> 
#include<stdbool.h>

/* Parameters for plotting */

int    N[]     = {4, 8};
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

/* Parameters for plotting */
#define EPS 0                     // 1: eps image; 0: png image


#define HIST_TYPE "cluster"
#define BOXWIDTH 0.75
#define MULTIPLOT_LAYOUT "2,3"



#define PREC "7.4" /* printing precision */
  
// write data in d array to file fp in n columns
void write_delta_ef(FILE *fp, double *d, int bi, int n, double gma)
{
  int i;
  double X = 0.0;
  fprintf(fp, "b=%.0f ", B[bi]);
  for(i = 0; i < n; i++){
    X += pow(d[i], 1.0/gma);
    fprintf(fp, "%"PREC"lf ", d[i]);
  }
  fprintf(fp, "%"PREC"lf ", pow(X, gma));     // total effort
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

void hist_sum(FILE *gp, char *edata, char *fdata, char *dxidata, char *dsidata, char *punidata, char *punjdata, int n, double gma, char *title, char *outputfile)
{
  /*int i;
  char data[300];
  int bs  = sizeof(B)/sizeof(double); // no. of Bs in array corresponds to no. of clusters
  */
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
  fprintf(gp, "set multiplot layout %s title '%s' \n", MULTIPLOT_LAYOUT, title);    // set subplots layout
  //fprintf(gp, "set xlabel '' \n");
  fprintf(gp, "set ylabel 'efforts' \n"); 
  fprintf(gp, "set auto x \n");
  
  //fprintf(gp, "set autoscale y \n");
  //fprintf(gp, "set yrange %s \n", E_YRANGE);
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
  fprintf(gp, "plot \\\n");  
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle  \\\n", n+1, edata);   
  fprintf(gp, "\n");
  
  fprintf(gp, "set ylabel 'threshold' \n");
  fprintf(gp, "unset y2tics \n");  
  //fprintf(gp, "unset autoscale y \n");
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ytics auto mirror\n");  
  fprintf(gp, "plot \\\n");   
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n",  n+1, dxidata);    
  fprintf(gp, "\n");  
  
  fprintf(gp, "set ylabel 'aggressiveness' \n");
  fprintf(gp, "unset y2tics \n");  
  //fprintf(gp, "unset autoscale y \n");
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ytics auto mirror\n");  
  fprintf(gp, "plot \\\n");   
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) title columnheader \\\n",  n+1, dsidata);    
  fprintf(gp, "\n"); 
  
  fprintf(gp, "set ylabel 'payoffs' \n");
  fprintf(gp, "unset y2tics \n");  
  //fprintf(gp, "unset autoscale y \n");
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ytics 0.1 mirror\n");  
  fprintf(gp, "plot \\\n");   
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n",  n+1, fdata);    
  fprintf(gp, "\n"); 
  
  fprintf(gp, "set ylabel 'cost' \n");
  fprintf(gp, "unset y2tics \n");  
  //fprintf(gp, "set autoscale y \n");
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ytics auto mirror\n");  
  fprintf(gp, "plot \\\n");   
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n",  n+1, punidata);    
  fprintf(gp, "\n"); 
  
  fprintf(gp, "set ylabel 'punishment' \n");
  fprintf(gp, "unset y2tics \n");  
  //fprintf(gp, "set autoscale y \n");
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ytics auto mirror\n");  
  fprintf(gp, "plot \\\n");   
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n",  n+1, punjdata);    
  fprintf(gp, "\n"); 
  
  fprintf(gp, "unset multiplot \n");
}
/*
void hist_pi_g(FILE *gp, char *pi_gdata, char *x_gdata, char *pun_gdata, char *title, char *outputfile)
{
  int i;   
  
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
  fprintf(gp, "set bmargin %lf \n", B_MARGIN);
  //fprintf(gp, "set lmargin %lf \n", L_MARGIN);
  //fprintf(gp, "set rmargin %lf \n", R_MARGIN);  
  
  fprintf(gp, "set multiplot layout 2,2 title '%s' \n", title);    // set subplots layout
  //fprintf(gp, "set xlabel '' \n");
  fprintf(gp, "set ylabel 'group payoff' \n"); 
  fprintf(gp, "set auto x \n");  
  fprintf(gp, "set yrange %s \n", E_YRANGE);  
  fprintf(gp, "set ytics auto mirror \n");    
  fprintf(gp, "plot \\\n");  
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) title columnheader \\\n", 2, pi_gdata);  
  fprintf(gp, "\n");  
  
  fprintf(gp, "set ylabel 'group effort' \n"); 
  fprintf(gp, "set auto x \n");  
  fprintf(gp, "set yrange %s \n", E_YRANGE);  
  fprintf(gp, "set ytics auto mirror \n");    
  fprintf(gp, "plot \\\n");  
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) title columnheader \\\n", 2, x_gdata);    
  fprintf(gp, "\n"); 
  
  fprintf(gp, "set ylabel 'group punishment' \n"); 
  fprintf(gp, "set auto x \n");  
  fprintf(gp, "set yrange %s \n", E_YRANGE);  
  fprintf(gp, "set ytics auto mirror \n");    
  fprintf(gp, "plot \\\n");  
  fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) title columnheader \\\n", 2, pun_gdata);     
  fprintf(gp, "\n"); 
  
  fprintf(gp, "unset multiplot \n");
}*/


int main(int argc, char **argv)
{  
  char infile[200], edata[200], fdata[200], pi_gdata[200], dxidata[200], dsidata[200], punjdata[200], punidata[200];
  double *d/*, val*/;  
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8;
  char str[100], title[100];    

  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode  
  
  int j;
  int ni, bi, ci, bti, ali, /*dli,*/ xi, gi, ki, ei, oi;
  int ns, bs, bts, cs, als, /*dls,*/ xs, gs, ks, es, os;
  
  ns  = sizeof(N)/sizeof(int);
  bs  = sizeof(B)/sizeof(double);
  cs  = sizeof(C)/sizeof(double);
  bts = sizeof(Beta)/sizeof(double);
  als = sizeof(Alpha)/sizeof(double);
  //dls = sizeof(Delta)/sizeof(double);
  xs  = sizeof(X0)/sizeof(double);
  ks  = sizeof(K)/sizeof(int);
  es  = sizeof(Eta)/sizeof(double);
  gs  = sizeof(Gamma)/sizeof(double);	
  os  = sizeof(Omega)/sizeof(double);
  
  for(ni = 0; ni < ns; ni++){
    d = malloc(N[ni]*sizeof(double));  // to hold data after each read of columns in data files
    
      for(ci = 0; ci < cs; ci++){
	for(bti = 0; bti < bts; bti++){
	  for(ali = 0; ali < als; ali++){    
	    for(xi = 0; xi < xs; xi++){
	      for(gi = 0; gi < gs; gi++){
		for(ei = 0; ei < es; ei++){
		  for(ki = 0; ki < ks; ki++){
		    for(oi = 0; oi < os; oi++){
		    sprintf(edata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fehist.dat",N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of efforts
		    sprintf(fdata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2ffhist.dat",  N[ni],  C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of payoff
		    sprintf(pi_gdata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fpi_ghist.dat",  N[ni],  C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of payoff
		    sprintf(dxidata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fdxihist.dat",N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for threshold
		    sprintf(dsidata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fdsihist.dat",  N[ni],  C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]); // data for aggressiveness
		    sprintf(punidata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunihist.dat",  N[ni],  C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
		    sprintf(punjdata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunjhist.dat",  N[ni],  C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
		    fp2 = fopen(edata, "w");
		    fp3 = fopen(fdata, "w");
		    fp4 = fopen(pi_gdata, "w");
		    fp5 = fopen(dxidata, "w");
		    fp6 = fopen(dsidata, "w");
		    fp7 = fopen(punidata, "w");
		    fp8 = fopen(punjdata, "w");
		    fprintf(fp2, "b's  ");
		    fprintf(fp3, "b's  ");
		    fprintf(fp4, "b's  ");
		    fprintf(fp5, "b's  ");
		    fprintf(fp6, "b's  ");
		    fprintf(fp7, "b's  ");
		    fprintf(fp8, "b's  ");
		    for(j = 0; j < N[ni]; j++){
		      fprintf(fp2, "%d  ", j+1); fprintf(fp3, "%d  ", j+1); fprintf(fp5, "%d  ", j+1); fprintf(fp6, "%d  ", j+1); fprintf(fp7, "%d  ", j+1); fprintf(fp8, "%d  ", j+1);  // adding headers
		    }
		    
		    fprintf(fp2, "X \n"); fprintf(fp3, "X \n"); fprintf(fp4, "pi_g \n"); fprintf(fp5, "X \n"); fprintf(fp6, "X \n"); fprintf(fp7, "X\n"); fprintf(fp8, "X\n");
		    for(bi = 0; bi < bs; bi++){
		      
			// files to be read for efforts for set of paramters
			sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fesum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);  
			if(!(fp1 = fopen(infile, "r"))){
			  printf("Could not open file '%s'",infile);
			  exit(1);
			}
			if(calc_data(fp1, d, N[ni], 0)){   // calculate data and not normalize
			  fclose(fp1);
			  write_delta_ef(fp2, d, bi, N[ni], Gamma[gi]);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			}
			else{
			  fclose(fp1);
			  printf("error in data read '%s': column count does not match!!\n", infile);
			  exit(1);
			}
		     
		      

		      // for fertilities
		       // files to be read for fertilities for set of paramters
			sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2ffsum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
			if(!(fp1 = fopen(infile, "r"))){
			  printf("Could not open file '%s'",infile);
			  exit(1);
			}
			if(calc_data(fp1, d, N[ni], 1)){          // calculate data and normalize fertitilties
			  fclose(fp1);
			  write_delta_ef(fp3, d, bi, N[ni], Gamma[gi]);   // write data to histogram data file for fertilities // j is passed as delta for linear spacing
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
			sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fdxisum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);  
			if(!(fp1 = fopen(infile, "r"))){
			  printf("Could not open file '%s'",infile);
			  exit(1);
			}
			if(calc_data(fp1, d, N[ni], 0)){   // calculate data and not normalize
			  fclose(fp1);
			  write_delta_ef(fp5, d, bi, N[ni], Gamma[gi]);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			}
			else{
			  fclose(fp1);
			  printf("error in data read '%s': column count does not match!!\n", infile);
			  exit(1);
			}
		      
			// files to be read for aggressiveness for set of paramters
			sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fdsisum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);  
			if(!(fp1 = fopen(infile, "r"))){
			  printf("Could not open file '%s'",infile);
			  exit(1);
			}
			if(calc_data(fp1, d, N[ni], 0)){   // calculate data and not normalize
			  fclose(fp1);
			  write_delta_ef(fp6, d, bi, N[ni], Gamma[gi]);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
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
			sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunisum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
			if(!(fp1 = fopen(infile, "r"))){
			  printf("Could not open file '%s'",infile);
			  exit(1);
			}
			if(calc_data(fp1, d, N[ni], 0)){   // calculate data and not normalize
			  fclose(fp1);
			  write_delta_ef(fp7, d, bi, N[ni], Gamma[gi]);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			}
			else{
			  fclose(fp1);
			  printf("error in data read '%s': column count does not match!!\n", infile);
			  exit(1);
			}
			// for punishment
			sprintf(infile, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunjsum.dat",  N[ni], B[bi], C[ci], Alpha[ali], Beta[bti], Delta[0],Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
			if(!(fp1 = fopen(infile, "r"))){
			  printf("Could not open file '%s'",infile);
			  exit(1);
			}
			if(calc_data(fp1, d, N[ni], 0)){   // calculate data and not normalize
			  fclose(fp1);
			  write_delta_ef(fp8, d, bi, N[ni], Gamma[gi]);   // write data to histogram data file for efforts // j is passed as delta for linear spacing
			}
			else{
			  fclose(fp1);
			  printf("error in data read '%s': column count does not match!!\n", infile);
			  exit(1);
			}
		    }		    
		    fclose(fp2); fclose(fp3); fclose(fp4); fclose(fp5); fclose(fp6); fclose(fp7); fclose(fp8); 
		    }
		}
	      }
	    }
	  }	 
	}
      }
    }
    free(d);
  }  
   // reading data and prepare histogram data file for multiple values of b and delta but same other parameters completed
   
// plotting summary histograms all at a time // the one with cluster base is removed from for loop (like for B in original case)
  for(ni = 0; ni < ns; ni++){    
    for(ci = 0; ci < cs; ci++){
      for(bti = 0; bti < bts; bti++){
	for(ali = 0; ali < als; ali++){    
	  for(xi = 0; xi < xs; xi++){
	    for(gi = 0; gi < gs; gi++){
	      for(ei = 0; ei < es; ei++){
		for(ki = 0; ki < ks; ki++){
		  for(oi = 0; oi < os; oi++){
		    sprintf(edata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fehist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of efforts
		    sprintf(fdata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2ffhist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of payoff   
		    sprintf(pi_gdata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fpi_ghist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of payoff
		    sprintf(dxidata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fdxihist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of efforts
		    sprintf(dsidata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fdsihist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);   // data for histogram of payoff
		    sprintf(punidata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunihist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]); 
		    sprintf(punjdata, "n%dc%.2fa%.2fb%.2fg%.2fe%.3fk%dx0%.2fw%.2fpunjhist.dat", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]); 
		    
		    // plot histograms      
		    sprintf(title, "n= %d, c= %.2f, al= %.2f, be= %.2f, ga= %.2f, eta= %.3f, k= %d, X0= %.2f, w=%.2f", N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
		    if(EPS)
		      sprintf(str, "sum_n%02dc%.2fa%.2fb%.2fg%.2fe%.3fk%02dx0%.2lfw%.2f.eps",N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
		    else
		      sprintf(str, "sum_n%02dc%.2fa%.2fb%.2fg%.2fe%.3fk%02dx0%.2lfw%.2f.png",N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
		    hist_sum(gp, edata, fdata, dxidata, dsidata, punidata, punjdata, N[ni], Gamma[gi], title, str);
		    sprintf(str, "sumg_n%02dc%.2fa%.2fb%.2fg%.2fe%.3fk%02dx0%.2lfw%.2f.png",N[ni], C[ci], Alpha[ali], Beta[bti], Gamma[gi], Eta[ei], K[ki], X0[xi], Omega[oi]);
		    //hist_pi_g(gp, pi_gdata, x_gdata, pun_gdata, title, str);
		    // histograms for efforts and payoffs for one set of n, c and X0 compeleted
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



/** Usage:
    compile : gcc -o hist hist.c
    run     : ./hist
**/

 
