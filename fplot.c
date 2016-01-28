#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h> 
#include<stdbool.h>

/* Parameters for plotting */
int      Runs;
int      N;       
double   C;       
double   B;        
double   Alpha;
double   Beta;    
double   Delta;
double   Gamma;
double   X0;
double   Eta;
double   Omega;
int      K;

/* Parameters for plotting */




#define PREC "7.4" /* printing precision */

// prepare filename to write data
void prep_file(char *filevar, char *apnd)
{  
  sprintf(filevar, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2f%s", N,B,C,Alpha,Beta,Delta,Gamma,Eta,K,X0,Omega, apnd);  
}


// plot data using lines for all columns as y axis data
void plotef(FILE *gp, int wid, char *edata, char *fdata, char *dxidata, char *dsidata, char *pun_idata, char *pun_jdata, int datacolumn, char *title, char *outputfile, int ind )
{
  fprintf(gp, "set key outside \n");        
  fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
  fprintf(gp, "set output '%s' \n", outputfile);
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  
  fprintf(gp, "set multiplot layout 2,3 title '%s' \n", title);    // set subplots layout
  
  fprintf(gp, "set ylabel 'efforts' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader \n", datacolumn, edata);
  
  fprintf(gp, "set ylabel 'threshold' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader \n", datacolumn, dxidata);
  
  fprintf(gp, "set ylabel 'aggressiveness' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader \n", datacolumn, dsidata);
  
  fprintf(gp, "set ylabel 'payoff' \n");
  if(!ind){
    fprintf(gp, "set style line 5 lc rgb '#808080' lt 2 lw 2 pi -1 ps 1.0\n");
    fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader, '%s' using 1:%d with line ls 5 title columnheader \n", datacolumn, fdata, fdata, datacolumn+1);  
  }
  else{
    fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader \n", datacolumn, fdata);
  }
   
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel 'cost of punishing' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader \n", datacolumn, pun_idata);
  
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel 'punishment' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 3 title columnheader \n", datacolumn, pun_jdata);
  
  fprintf(gp, "unset multiplot \n");
  
}

void plotpi_g(FILE *gp, int wid, char *data, int datacolumn, char *title, char *ylabel, char *outputfile )
{
  fprintf(gp, "set key outside \n");        
  fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,12\" \n");
  fprintf(gp, "set output '%s' \n", outputfile);
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title '%s' \n", title);    
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel '%s' \n", ylabel);
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines notitle \n", datacolumn, data);  
}

void plotpun(FILE *gp, int wid, char *data, int datacolumn, char *title, char *ylabel, char *outputfile )
{
  fprintf(gp, "set key outside \n");        
  fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,12\" \n");
  fprintf(gp, "set output '%s' \n", outputfile);
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title '%s' \n", title);    
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel '%s' \n", ylabel);
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, data);  
}


/*
// plot data using lines for all columns as y axis data
void plotef(FILE *gp, int wid, char *dxidata, char *dsidata, int datacolumn, char *title, char *outputfile )
{
  fprintf(gp, "set key outside \n");        
  fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,12\" \n");
  fprintf(gp, "set output '%s' \n", outputfile);
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 2,1 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set ylabel 'efforts' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, edata);
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, fdata);  

  fprintf(gp, "unset multiplot \n");
  
}
*/


int main(int argc, char **argv)
{
  if(argc ^ 13){ printf("Usage: ./fplot n c b alpha beta delta gamma x0 eta omega k runs\n"); return 0;}
  
  N = atoi(argv[1]);
  C = atof(argv[2]);
  B = atof(argv[3]);
  Alpha = atof(argv[4]);
  Beta  = atof(argv[5]);
  Delta = atof(argv[6]);
  Gamma = atof(argv[7]);
  X0    = atof(argv[8]);
  Eta   = atof(argv[9]);
  Omega = atof(argv[10]);
  K     = atoi(argv[11]);
  Runs  = atoi(argv[12]);
	 
  //open file for writing average efforts per rank per generation
  int i;
  char xdata[200], fdata[200], pi_gdata[200], tstr[200], dxidata[200], dsidata[200], pundata[200], pun_idata[200], pun_jdata[200];
  prep_file(xdata, "e.dat"); prep_file(fdata,"f.dat"); prep_file(pi_gdata, "pi_g.dat"); 
  prep_file(dxidata, "dxi.dat"); prep_file(dsidata, "dsi.dat");
  prep_file(pundata, "pun.dat"); prep_file(pun_idata, "pun_i.dat"); prep_file(pun_jdata, "pun_j.dat");
  
  char title[200], xpng[200];
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "n= %d, b= %.2f, c= %.2f, al= %.2f, be= %.2f, dl= %.2f, eta= %.3f, k= %d, X0= %.2f, w= %.2f", N, B, C, Alpha, Beta, Delta, Eta, K, X0, Omega);
  sprintf(xpng, "ef_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2fw%.2f.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0,Omega);
  plotef(gp, 0, xdata, fdata, dxidata, dsidata, pun_idata, pun_jdata, N+1, title, xpng, 0 );  
  //sprintf(xpng, "pi_g_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2f.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0);
  //plotpi_g(gp, 0, pi_gdata, 2, title, "group payoff", xpng );
  //sprintf(xpng, "pun_g_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2f.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0);
  //plotpun(gp, 0, pundata, N+1, title, "punishment", xpng );
  
  for(i = 0; i < Runs; i++){
    sprintf(tstr, "e%d.dat", i);
    prep_file(xdata, tstr);
    sprintf(tstr, "f%d.dat", i);
    prep_file(fdata, tstr);
    sprintf(tstr, "pi_g%d.dat", i);
    prep_file(pi_gdata, tstr);
    sprintf(tstr, "dxi%d.dat", i);
    prep_file(dxidata, tstr);
    sprintf(tstr, "dsi%d.dat", i);
    prep_file(dsidata, tstr);
    sprintf(tstr, "pun_i%d.dat", i);
    prep_file(pun_idata, tstr);
    sprintf(tstr, "pun_j%d.dat", i);
    prep_file(pun_jdata, tstr);
    sprintf(title, "n= %d, b= %.2f, c= %.2f, al= %.2f, be= %.2f, dl= %.2f, eta= %.3f, k= %d, X0= %.2f, w= %.2f, %d", N, B, C, Alpha, Beta, Delta, Eta, K, X0, Omega, i+1);
    sprintf(xpng, "ef_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2fw%.2f_%d.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0, Omega, i+1);
    plotef(gp, 0, xdata, fdata, dxidata, dsidata, pun_idata, pun_jdata, N+1, title, xpng, 1 );      
    //sprintf(xpng, "pi_g_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2f_%d.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0, i+1);
    //plotpi_g(gp, 0, pi_gdata, 2, title, "group payoff", xpng );
  }
  fflush(gp); 
  pclose(gp);

  return 0;
}



/** Usage:
    compile : gcc -o fplot fplot.c
    run     : ./fplot n c b alpha beta delta gamma x0 runs
**/

