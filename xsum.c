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

// Parameters to specify type of plots
#define WINPOPUP 0           // 0 or 1; if 0, png file is created; if 1, window pop up is generated
#define STACKED 0            //  1 or 0; if 1 VBAR should be defined
#define xAXIS "N"            // horizontal axis in a plot in a graph
#define yAXIS "effort"   // options are: effort, payoff, threshold, aggressiveness, cost, punishment
#define XAXIS "B"            // variation in plots in columns in a graph
#define YAXIS "Omega"            // variation in plots in rows in a graph
#define VBAR  "K"        // variation in vertical stacked bars in a cluster of a histogram


double Pn[11], *P[11];        // Pn stores value of parameter currently under multiple loops of parameters for each indexed parameter; P array points to reference of parameter used in plotting
char Pa[11][3];               // stores alphbet strings array used to denote specific parameter in creating data file 
int Ps[11];                   // Ps stores no. of values under one parameter indexed by its index.

// methods declaration
void set_param();
void prep_infile(char *fname, char *appnd);  
void prep_outfile(char *file, char *appnd);
void plot_graphs();
void prepare_data_files();


/* Parameters for plotting */
#define EPS 0                     // 1: eps image; 0: png image

#define BOXWIDTH 0.75


#define PREC "7.4" /* printing precision */
  
// write data in d array to file fp in n columns
/*
 * fp : file pointer
 * d  : pointer to data to be written
 * ba : alphabet representation of parameter for which data is to be written
 * bi : value of parameter for wihic data is to be written
 * n  : no. of columns of data to be written
 * pad: 1 or 0; 1 if no. of columms of data for each row is different. it pads columns with zero to max possible columns; in this case, data is collected for different individuals in a group. group size may vary
 */
void write_data(FILE *fp, double *d, char *ba, double bi, int n, int pad)
{
  int i, k, nx;
  fprintf(fp, "%s=%.0f ", ba, bi);
  for(i = 0; i < n; i++){
    fprintf(fp, "%"PREC"lf ", d[i]);
  }
  if(pad){
    for(nx = 0, i = 0; i < sizeof(N)/sizeof(double); i++){
      if(nx < N[i]) nx = N[i];
    }
    if(n < nx){
      for(k = n; k < nx; k++)
	fprintf(fp, "%.4f ", 0.0);
    }
  }
  fprintf(fp, "\n");
}

// calc averages of data in columns 
/*
 * fp: file pointer
 * d : pointer array in which data to be stored
 * n : no. of columns of data to be read/written
 * toNormal : normalize data or not; 1 or 0
 */
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

void hist_sum(FILE *gp, char *edata, int n, char *title, char *outputfile, int wid)
{
  int i, xN = 6;
  char xlabel[6][200];
  
#if STACKED 
  int pi;
#endif
  char data[300]; char str[20];
  //int bs  = sizeof(B)/sizeof(double); // no. of Bs in array corresponds to no. of clusters
#if WINPOPUP 
  fprintf(gp, "set term x11 %d \n", wid);
#else
  if(EPS)
    fprintf(gp, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
  else
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,10\" \n");
  fprintf(gp, "set output '%s' \n", outputfile);  
#endif
  
  fprintf(gp, "set key autotitle columnheader\n");
  fprintf(gp, "set key outside \n");
  fprintf(gp, "set style data histogram \n");  
#if STACKED
  fprintf(gp, "set style histogram rowstacked gap 1 title offset 0, -1 \n");
#else
  fprintf(gp, "set style histogram cluster title offset 0, -1 \n");
#endif
  fprintf(gp, "set style fill solid border -1\n");
  fprintf(gp, "set boxwidth %lf relative\n", BOXWIDTH);
  
  //fprintf(gp, "set bmargin 2.0 \n");
  fprintf(gp, "set title ''\n");    
  fprintf(gp, "set yrange [0:] \n");
  
  int k;
#if STACKED
  char tstr[30], ext[50];
  sprintf(ext, " %s=", Pa[7]);
  
  for(i = 0; i < Ps[7]; i++){
    if(i)
      sprintf(tstr,",%g", P[7][i]);
    else
      sprintf(tstr, "%g", P[7][i]);
    strcat(ext, tstr);
  }
  if(strcmp(XAXIS, "") == 0){
    fprintf(gp, "set multiplot layout %d, %d title '%s' \n", 2, 3, title);
  }
  else{
    fprintf(gp, "set multiplot layout %d, %d title '%s, %s' \n", Ps[8], Ps[9], title, ext);                        // set subplots layout   
  }
#else
  if(strcmp(XAXIS, "") == 0){
    fprintf(gp, "set multiplot layout %d, %d title '%s' \n", 2, 3, title);
  }
  else{
    fprintf(gp, "set multiplot layout %d, %d title '%s' \n", Ps[8], Ps[9], title);                        // set subplots layout  
  }
#endif
  
  // if n is set as xaxis, set n to be the highest number of individuals in a group used in simulation
  if(strcmp(Pa[10], "n") == 0){
    for(n = 0, i = 0; i < Ps[10]; i++){
      if(n < P[10][i]) n = P[10][i];
    }
  }
  else if(strcmp(Pa[7], "n") == 0){
    for(n = 0, i = 0; i < Ps[7]; i++){
      if(n < P[7][i]) n = P[7][i];
    }
  }     
  
#if STACKED
    fprintf(gp, "unset xtics \n");
#else
    fprintf(gp, "set auto x \n");
#endif  
  
  fprintf(gp, "set style arrow 1 nohead lw 2 \n");
  
  if(strcmp(XAXIS, "") != 0){                                          // when xaxis is specified
    if(strcmp(yAXIS, "effort") == 0) sprintf(str, "ehist.dat");
    else if(strcmp(yAXIS, "payoff") == 0) sprintf(str, "fhist.dat");
    else if(strcmp(yAXIS, "threshold") == 0) sprintf(str, "dxihist.dat");
    else if(strcmp(yAXIS, "aggressiveness") == 0) sprintf(str, "dsihist.dat");
    else if(strcmp(yAXIS, "cost") == 0) sprintf(str, "punihist.dat");
    else if(strcmp(yAXIS, "punishment") == 0) sprintf(str, "punjhist.dat");
    for(k = 0; k < Ps[8]; k++){
      if( strcmp(Pa[8], "n") == 0){
	n = P[8][k];
      }   
      for(i = 0; i < Ps[9]; i++){
	if( strcmp(Pa[9], "n") == 0){
	  n = P[9][i];
	}
	fprintf(gp, "set ylabel '%s, %s=%g, %s=%g' \n", yAXIS, Pa[8], P[8][k], Pa[9], P[9][i]); 
      
	sprintf(data, "%s%s%.2f%s%.2f%s", edata, Pa[8], P[8][k], Pa[9], P[9][i], str);
#if !STACKED
	fprintf(gp, "plot \\\n");   
	fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+1, data);    
	fprintf(gp, "\n"); 
#endif
#if STACKED      
	fprintf(gp, "plot \\\n");   
	for(pi = 0; pi < Ps[10]; pi++){    
	  if(pi != 0) fprintf(gp, ", ");
	  fprintf(gp, "newhistogram \"%s=%g\" lt 1 , for[col=2:%d] '%s' index %d u col:xtic(1) notitle \\\n", Pa[10], P[10][pi], n+1, data, pi);
	}
	fprintf(gp, "\n");
#endif 
      }
    }
  }
  else{                                                                 // when xaxis is not specified; 
    sprintf(xlabel[0], "effort"); sprintf(xlabel[1], "threshold"); sprintf(xlabel[2], "aggressiveness"); sprintf(xlabel[3], "payoff"); sprintf(xlabel[4], "cost"); sprintf(xlabel[5], "punishment"); 
    for(i = 0; i < xN; i++){
      if(i == 0) sprintf(str, "ehist.dat");
      else if(i == 1) sprintf(str, "dxihist.dat");
      else if(i == 2) sprintf(str, "dsihist.dat");
      else if(i == 3) sprintf(str, "fhist.dat");
      else if(i == 4) sprintf(str, "punihist.dat");
      else if(i == 5) sprintf(str, "punjhist.dat");
      fprintf(gp, "set ylabel '%s, %s=%g, %s=%g' \n", xlabel[i], Pa[8], P[8][k], Pa[9], P[9][i]);      
      sprintf(data, "%s%s", edata, str);
      #if !STACKED
	fprintf(gp, "plot \\\n");   
	fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) notitle \\\n", n+1, data);    
	fprintf(gp, "\n"); 
#endif
#if STACKED      
	fprintf(gp, "plot \\\n");   
	for(pi = 0; pi < Ps[10]; pi++){    
	  if(pi != 0) fprintf(gp, ", ");
	  fprintf(gp, "newhistogram \"%s=%g\" lt 1 , for[col=2:%d] '%s' index %d u col:xtic(1) notitle \\\n", Pa[10], P[10][pi], n+1, data, pi);
	}
	fprintf(gp, "\n");
#endif 
    }
    
  }
  
  
  fprintf(gp, "unset multiplot \n");
}

int main(int argc, char **argv)
{      
  if(STACKED){
    if(strcmp(VBAR, "") == 0 || !VBAR){
      printf("Please define VBAR parameter! \n");
      exit(1);
    }
  }  
  set_param();                                            // sets P array with pointers to array of parameters like B, N, Beta,...
  prepare_data_files();                                   // readies data files for plotting graphs  
  plot_graphs();  
  return 0;
}

void set_param(){
  int i, h = -1, p = -1, z = -1, y = -1;
  double *Param[11];
  char ParamStr[11][10];
  char ParamA[11][3];
  int ParamS[11];
  int ns, bs, bts, cs, als, dls, xs, gs, ks, es, os;
  
  // calculating size of arrays
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
  
  // copying array address
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
  // setting array name as strings
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
  // setting alphabets representing parameters
  sprintf(ParamA[0], "n");
  sprintf(ParamA[1], "b");
  sprintf(ParamA[2], "c");
  sprintf(ParamA[3], "a");
  sprintf(ParamA[4], "be");
  sprintf(ParamA[5], "dl");
  sprintf(ParamA[6], "ga");
  sprintf(ParamA[7], "e");
  sprintf(ParamA[8], "k");
  sprintf(ParamA[9], "x0");
  sprintf(ParamA[10], "w");
  
  // assigning array reference according to alphabets for easy location of main axes's parameter
  for(i = 0; i < 11; i++){
    if( strcmp(ParamStr[i], xAXIS) == 0){
      h = i;
      P[10] = Param[i];
      sprintf(Pa[10],"%s", ParamA[i]);
      Ps[10] = ParamS[i];
    }
    else if( strcmp(ParamStr[i], XAXIS) == 0 ){
      p = i;
      P[9] = Param[i];
      sprintf(Pa[9], "%s", ParamA[i]);
      Ps[9] = ParamS[i];         
    }
    else if( strcmp(ParamStr[i], YAXIS) == 0 && strcmp("", XAXIS) != 0){
      z = i;
      P[8] = Param[i];
      sprintf(Pa[8], "%s", ParamA[i]);
      Ps[8] = ParamS[i];
    }
#if STACKED
    else if( strcmp(ParamStr[i], VBAR) == 0){
      y = i;
      P[7] = Param[i];
      sprintf(Pa[7], "%s", ParamA[i]);
      Ps[7] = ParamS[i];
    }    
#endif
  }
  if( strcmp(XAXIS, "") == 0 ){
    p = -1;
  }
  int k;
  for(k = 0, i = 0; i < 11; i++){    
    printf("h:%d p:%d y:%d z:%d i:%d k:%d\n", h, p, y, z, i, k);
    if( i != h && i != p && i != z && i != y){            
      P[k] =  Param[i];
      sprintf(Pa[k],"%s", ParamA[i]);
      Ps[k] = ParamS[i];
      k++;      
    }
  }  
} 

// prepares name of file from which data to be read
void prep_infile(char *fname, char *appnd)
{
  char str[10];
  strcpy(str, "");
  strcpy(fname, "");
  int i;
  // for n
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "n") == 0){
      sprintf(str,"n%.0f", Pn[i]);
    }
  }  
  strcat(fname, str);
  
  // for b
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "b") == 0){
      sprintf(str,"b%.2f", Pn[i]);
    }
  }    
  strcat(fname, str);
  
  // for c
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "c") == 0){
      sprintf(str,"c%.2f", Pn[i]);
    }    
  }     
  strcat(fname, str);
  
  // for a
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "a") == 0){
      sprintf(str,"a%.2f", Pn[i]);
    }    
  }   
  strcat(fname, str);
  
  // for be
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "be") == 0){
      sprintf(str,"b%.2f", Pn[i]);
    }    
  }     
  strcat(fname, str);
  
  // for d
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "dl") == 0){
      sprintf(str,"d%.2f", Pn[i]);
    }    
  }   
  strcat(fname, str);
  
  // for g
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "ga") == 0){
      sprintf(str,"g%.2f", Pn[i]);
    }    
  }  
  strcat(fname, str);
  
  // for e
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "e") == 0){
      sprintf(str,"e%.3f", Pn[i]);
    }    
  }    
  strcat(fname, str);
  
  // for k
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "k") == 0){
      sprintf(str,"k%.0f", Pn[i]);
    }    
  }  
  strcat(fname, str);
  
  strcpy(str, "");
  // for x0
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "x0") == 0){
      sprintf(str,"x0%.2f", Pn[i]);
    }    
  }  
  strcat(fname, str);
  
  // for w
  for(i = 0; i < 11; i++){
    if(strcmp(Pa[i], "w") == 0){
      sprintf(str,"w%.2f", Pn[i]);
    }    
  }  
  strcat(fname, str);
  
  strcat(fname, appnd);
} 

// prepares file to write data for graphs
void prep_outfile(char *file, char *appnd)
{
#if !STACKED
  sprintf(file, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s", Pa[0], Pn[0], Pa[1], Pn[1], Pa[2], Pn[2], Pa[3], Pn[3], Pa[4], Pn[4], Pa[5], Pn[5], Pa[6], Pn[6], Pa[7], Pn[7], Pa[8], Pn[8], Pa[9], Pn[9], appnd);
#else
  sprintf(file, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s", Pa[0], Pn[0], Pa[1], Pn[1], Pa[2], Pn[2], Pa[3], Pn[3], Pa[4], Pn[4], Pa[5], Pn[5], Pa[6], Pn[6], Pa[8], Pn[8], Pa[9], Pn[9], appnd);
#endif
}

void prepare_data_files()
{
  int xN = 6;                                                  // no. of plots in a graph if XAXIS is not defined; By default effort, threshold, aggressiveness, payoff, cost and punishment
  char infile[200], edata[6][200], str[100];
  double *d = NULL, v;
  char a[3];
  FILE *fp1, *fp2[6];
  int i, j, pi, qi, ri, si, ti, ui, vi, wi, xi, yi, zi, ni = 0, pad = 0;   
  
  // looping through all parameter arrays
  for(ri = 0; ri < Ps[0]; ri++){
    if( strcmp(Pa[0], "n") == 0){
      d = malloc(P[0][ri]*sizeof(double));  // to hold data after each read of columns in data files   
      ni = (int)P[0][ri];
    }    
    Pn[0] = P[0][ri];
      for(si = 0; si < Ps[1]; si++){
	Pn[1] = P[1][si];
	for(ti = 0; ti < Ps[2]; ti++){
	  Pn[2] = P[2][ti];
	  for(ui = 0; ui < Ps[3]; ui++){    	    
	    Pn[3] = P[3][ui];
	    for(vi = 0; vi < Ps[4]; vi++){	      
	      Pn[4] = P[4][vi];
	      for(wi = 0; wi < Ps[5]; wi++){
		Pn[5] = P[5][wi];		
		for(xi = 0; xi < Ps[6]; xi++){
		  Pn[6] = P[6][xi];		  
#if !STACKED
		  for(yi = 0; yi < Ps[7]; yi++){
		    Pn[7] = P[7][yi];		        
#endif
		    for(zi = 0; zi < Ps[8]; zi++){
		      Pn[8] = P[8][zi];
		      if( strcmp(Pa[8], "n") == 0){
			d = malloc(P[8][zi]*sizeof(double));  // to hold data after each read of columns in data files
			ni = (int)P[8][zi];
		      }
		      for(qi = 0; qi < Ps[9]; qi++){
			Pn[9] = P[9][qi];
			if( strcmp(Pa[9], "n") == 0){
			  d = malloc(P[9][qi]*sizeof(double));  // to hold data after each read of columns in data files
			  ni = (int)P[9][qi];
			}
			for(i = 0; i < xN; i++){
			  if(i == 0) sprintf(str, "ehist.dat");
			  else if(i == 1) sprintf(str, "fhist.dat");
			  else if(i == 2) sprintf(str, "dxihist.dat");
			  else if(i == 3) sprintf(str, "dsihist.dat");
			  else if(i == 4) sprintf(str, "punihist.dat");
			  else if(i == 5) sprintf(str, "punjhist.dat");
			  
			  prep_outfile(edata[i], str);                                                               // data for histogram of property of interest
			  fp2[i] = fopen(edata[i], "w");
			}
			
#if !STACKED
			// prepare headers
			for( i = 0; i < xN; i++){
			  fprintf(fp2[i], "%s's ", Pa[10]);
			  
			  if( strcmp(Pa[8], "n") == 0){                        // no. of columns will be according to group size n
			    for(j = 0; j < P[8][zi]; j++){
			      fprintf(fp2[i], "%d  ", j+1); 
			    }
			    fprintf(fp2[i], "X \n"); 
			  } 
			  else if( strcmp(Pa[9], "n") == 0){
			    for(j = 0; j < P[9][qi]; j++){
			      fprintf(fp2[i], "%d  ", j+1);  
			    }
			    fprintf(fp2[i], "X \n"); 
			  } 
			  else if( strcmp(Pa[0], "n") == 0){
			    for(j = 0; j < P[0][ri]; j++){
			      fprintf(fp2[i], "%d  ", j+1); 
			    }
			    fprintf(fp2[i], "X \n"); 
			  }
			  else if(strcmp(Pa[10], "n") == 0){    
			    for(j = 0; j < P[10][0]; j++){
			      fprintf(fp2[i], "%d  ", j+1); 
			    }
			    fprintf(fp2[i], "X \n"); 	     
			  }
			}
		      
#endif
			for(pi = 0; pi < Ps[10]; pi++){
			  Pn[10] = P[10][pi];
			  if( strcmp(Pa[10], "n") == 0 ){
			    d = malloc(P[10][pi]*sizeof(double));  // to hold data after each read of columns in data files
			    ni = (int)P[10][pi];    
			  }
#if !STACKED  
			  sprintf(a,"%s", Pa[10]); v = P[10][pi];
#endif
#if STACKED
			  // prepare headers
			  for(i = 0; i < xN; i++){
			    fprintf(fp2[i], "#%s = %g\n", Pa[10], P[10][pi]); fprintf(fp2[i], "%s ", Pa[7]);	  
			    
			    if( strcmp(Pa[0], "n") == 0){
			      for(j = 0; j < P[0][ri]; j++){
				fprintf(fp2{i], "%d  ", j+1); 
			      }
			      fprintf(fp2[i}, "X \n"); 
			    }
			    else if( strcmp(Pa[7], "n") == 0){
			      for(j = 0; j < P[7][0]; j++){
				fprintf(fp2[i], "%d  ", j+1); 
			      }
			      fprintf(fp2[i], "X \n"); 
			    }
			    else if( strcmp(Pa[8], "n") == 0){
			      for(j = 0; j < P[8][zi]; j++){
				fprintf(fp2[i], "%d  ", j+1);
			      }
			      fprintf(fp2[i], "X \n"); 
			    } 
			    else if( strcmp(Pa[9], "n") == 0){
			      for(j = 0; j < P[9][qi]; j++){
				fprintf(fp2[i], "%d  ", j+1); 
			      }
			      fprintf(fp2[i], "X \n"); 
			    } 
			    else if(strcmp(Pa[10], "n") == 0){			    
			      for(j = 0; j < P[10][pi]; j++){
				fprintf(fp2[i], "%d  ", j+1); 
			      }
			      fprintf(fp2[i], "X \n");      
			    }
			  }
			  // headers end
			  
			  for(pad = 0, yi = 0; yi < Ps[7]; yi++){
			    Pn[7] = P[7][yi];
			    if( strcmp(Pa[7], "n") == 0 ){
			      d = malloc(P[7][yi]*sizeof(double));  // to hold data after each read of columns in data files
			      ni = (int)P[7][yi];    
			      pad = 1;
			    }
			    sprintf(a, "%s", Pa[7]); v = P[7][yi];
#endif
			    // files to be read for any property of interest for set of paramters
			    for(i = 0; i < xN;i++){
			      if(strcmp(XAXIS, "") ==0){
				if(i == 0) prep_infile(infile, "esum.dat");
				else if(i == 1) prep_infile(infile, "fsum.dat");
				else if(i == 2) prep_infile(infile, "dxisum.dat");
				else if(i == 3) prep_infile(infile, "dsisum.dat");
				else if(i == 4) prep_infile(infile, "punisum.dat");
				else if(i == 5) prep_infile(infile, "punjsum.dat");
			      }
			      else{
				if(strcmp(yAXIS, "effort") == 0) prep_infile(infile, "esum.dat");
				else if(strcmp(yAXIS, "payoff") == 0) prep_infile(infile, "fsum.dat");
				else if(strcmp(yAXIS, "threshold") == 0) prep_infile(infile, "dxisum.dat");
				else if(strcmp(yAXIS, "aggressiveness") == 0) prep_infile(infile, "dsisum.dat");
				else if(strcmp(yAXIS, "cost") == 0) prep_infile(infile, "punisum.dat");
				else if(strcmp(yAXIS, "punishment") == 0) prep_infile(infile, "punjsum.dat");
			      }
			      if(!(fp1 = fopen(infile, "r"))){
				printf("Could not open file '%s'",infile);
				exit(1);
			      }
			      if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize			    
				fclose(fp1);
				write_data(fp2[i], d, a, v, ni, pad);   // write data to histogram data file for efforts // j is passed as delta for linear spacing    
			      }
			      else{
				fclose(fp1);
				printf("error in data read '%s': column count does not match!!\n", infile);
				exit(1);
			      }    
			    }
#if STACKED
			    if( strcmp(Pa[7], "n") == 0 ){
			      free(d);   
			    }
			  }
			  for( i = 0; i < xN; i++)
			    fprintf(fp2[i],"\n\n"); 
#endif
	  
			  if( strcmp(Pa[10], "n") == 0 ){
			    free(d);
			  }
			  
			}// pi loop ends
			for( i = 0; i < xN; i++)
			  fclose(fp2[i]); 
			if( strcmp(Pa[9], "n") == 0){
			  free(d);
			}  
		      }// qi loop ends		    
		      if( strcmp(Pa[8], "n") == 0){
			free(d);
		      } 
		    } // zi loop ends
#if !STACKED
		}  // yi loop ends
#endif
	      }  // xi loop ends
	    }  // wi loop ends
	  }  // vi loop ends
	}  // ui loop ends
      }  // ti loop ends
    }  // si loop ends
    if( strcmp(Pa[0], "n") == 0){
      free(d);
    }     
  }  // ri loop ends
   // reading data and prepare histogram data file for multiple values of b and delta but same other parameters completed
}

void plot_graphs()
{
  char str[100], title[100], edata[200];   
  int ri, si, ti, ui, vi, wi, xi, zi, qi, ni = 0;    
#if !STACKED
  int yi;
#endif
  int wid = 0;
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode  
  // plotting summary histograms all at a time // the one with cluster base is removed from for loop (like for B in original case)
  
  
  for(ri = 0; ri < Ps[0]; ri++){      
    if( strcmp(Pa[0], "n") == 0){
      ni = P[0][ri];
    }     
    for(si = 0; si < Ps[1]; si++){
     
      for(ti = 0; ti < Ps[2]; ti++){
	for(ui = 0; ui < Ps[3]; ui++){    
	  for(vi = 0; vi < Ps[4]; vi++){
	    for(wi = 0; wi < Ps[5]; wi++){
	      
	      for(xi = 0; xi < Ps[6]; xi++){
#if !STACKED		
		if(strcmp(XAXIS, "") != 0){                      // when xaxis is specified
		  
		  for(yi = 0; yi < Ps[7]; yi++){
		    sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi]);   // data for histogram of efforts
		    // plot histograms      
		    sprintf(title, "%s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi]);
		    if(EPS)
		      sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.eps",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi]);
		    else
		      sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.png",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi]);
		    hist_sum(gp, edata, ni, title, str, ++wid); 
		    
		  }
		}
		else{                                            // when xaxis is not specified, it will plot all properties of a simulation
		  for(yi = 0; yi < Ps[7]; yi++){
		    for(zi = 0; zi < Ps[8]; zi++){
		      for(qi = 0; qi < Ps[9]; qi++){
			sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi], Pa[8], P[8][zi], Pa[9], P[9][qi]);   // data for histogram of efforts
			// plot histograms      
			sprintf(title, "%s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi], Pa[8], P[8][zi], Pa[9], P[9][qi]);
			if(EPS)
			  sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.eps",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi], Pa[8], P[8][zi], Pa[9], P[9][qi]);
			else
			  sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.png",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][yi], Pa[8], P[8][zi], Pa[9], P[9][qi]);
			
			hist_sum(gp, edata, ni, title, str, ++wid); 
			
		      }
		    }
		  }
		}
	      
#endif
#if STACKED
	      if(strcmp(XAXIS, "") != 0){                      // when xaxis is specified
		sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi]);   // data for histogram of efforts
		// plot histograms      
		sprintf(title, "%s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi]);
		if(EPS)
		  sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g.eps",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi]);
		else
		  sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g.png",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi]);
		hist_sum(gp, edata, ni, title, str, ++wid);
	      }
	      else{
		for(zi = 0; zi < Ps[8]; zi++){
		  for(qi = 0; qi < Ps[9]; qi++){
		    sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[8], P[8][yi], Pa[9], P[9][zi]);   // data for histogram of efforts
		    // plot histograms      
		    sprintf(title, "%s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[8], P[8][yi], Pa[9], P[9][zi]);
		    if(EPS)
		      sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.eps",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[8], P[8][yi], Pa[9], P[9][zi]);
		    else
		      sprintf(str, "sum_%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g.png",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[8], P[8][yi], Pa[9], P[9][zi]);
		    hist_sum(gp, edata, ni, title, str, ++wid);
		  }
		  
		}
	      }
#endif
		
	      }
	    }
	  }
	}
      }
    }
  }
    
  fflush(gp); 
  pclose(gp);
}
  
/** Usage:
    compile : gcc -o xsum xsum.c -lm
    run     : ./xsum
**/
 
