
/************************* writing of data to file ****************************/

void writefile_effort_fertility(){
  int h, j;
  double pm;  
  char xdata[200], fdata[200], pi_gdata[200];
  FILE *fp1 = NULL, *fp2 = NULL, *fp3 = NULL;   
  prep_file(xdata, "e.dat"); prep_file(fdata,"f.dat"); prep_file(pi_gdata, "pi_g.dat");    
  fp1 = fopen(xdata, "w"); fp2 = fopen(fdata, "w"); fp3 = fopen(pi_gdata, "w");
  
  double *xtu_sum, *pitu_sum;     // holds average of summary of results for last few considered time unit for summary
  
  xtu_sum  = calloc(N,sizeof(double));     // hold average value of efforts for each rank for last few generations considered for summary
  pitu_sum  = calloc(N,sizeof(double));     // hold average value of fertilities for each rank for last few generations considered for summary
  
  
#if ALLDATAFILE
  double STU = (T-SUMT)/SKIP;
  double pi_g_sum;                // holds average of group payoffs for last few consider time units for summary
  double xg_sum;
  pi_g_sum  = 0;                            // hold average value of group payoff for summary time units
  xg_sum = 0;  
  char esumdata[200], fsumdata[200], pi_gsumdata[200], xgsumdata[200];
  FILE *fp4 = NULL, *fp5 = NULL, *fp6 = NULL, *fp7 = NULL;
  prep_file(esumdata, "esum.dat"); prep_file(fsumdata,"fsum.dat"); prep_file(pi_gsumdata, "pi_gsum.dat");
  prep_file(xgsumdata, "x_gsum.dat"); 
  fp4 = fopen(esumdata, "w"); fp5 = fopen(fsumdata, "w"); fp6 = fopen(pi_gsumdata, "w"); fp7 = fopen(xgsumdata, "w");
#endif
  
  // write headers
  fprintf(fp1, "%d\t", 0);
  fprintf(fp2, "%d\t", 0);
  for(h = 0; h < N; h++){
    fprintf(fp1, "%d\t", h + 1);
    fprintf(fp2, "%d\t", h + 1);
  }    
  fprintf(fp1, "\n");    
  fprintf(fp2, "avg \n");  

  // write data
  for(j = 0; j < (int)(T/SKIP); j++){
    fprintf(fp1, "%d\t", (int)(j*SKIP));
    fprintf(fp2, "%d\t", (int)(j*SKIP)); 
    fprintf(fp3, "%d\t", (int)(j*SKIP));    
    pm = 0;
    for(h = 0; h < N; h++){
      xmean[j][h] = xmean[j][h] / (double)Runs; // average effort per rank over all groups over number of runs
      pimean[j][h] = pimean[j][h] / (double)Runs;
      pm += pimean[j][h];
      fprintf(fp1, "%.4lf  ", xmean[j][h]);
      fprintf(fp2, "%.4lf  ", pimean[j][h]);                 
    }       
    fprintf(fp2, "%.4lf  ", pm/(double)N);
    fprintf(fp3, "%.4lf  ", pi_gmean[j]/(double)Runs);
    fprintf(fp1, "\n");        
    fprintf(fp2, "\n");  
    fprintf(fp3, "\n");
#if ALLDATAFILE
    if(j < STU ) continue;  // calculate summary below
    for(h = 0; h < N; h++){
      xtu_sum[h]  += xmean[j][h];
      pitu_sum[h] += pimean[j][h];    
      pi_g_sum    += pimean[j][h];        
    }
#endif
  }   
  fclose(fp1); fclose(fp2); fclose(fp3);  
  
#if ALLDATAFILE
  int i;
  double st = SUMT/(double)SKIP;
  // write summary of efforts and payof to file  
  for(i = 0; i < N; i++){
    xtu_sum[i] = xtu_sum[i]/st;
    pitu_sum[i] = pitu_sum[i]/st;
    pi_g_sum = pi_g_sum/st;   
    xg_sum      += xtu_sum[i];
    fprintf(fp4, "%.6f\t", xtu_sum[i]);
    fprintf(fp5, "%.6f\t", pitu_sum[i]);
    fprintf(fp6, "%.6f\t", pi_g_sum);        
  }
  fprintf(fp7, "%.6f\n", xg_sum);
  fprintf(fp4, "\n");
  fprintf(fp5, "\n");  
  fprintf(fp6, "\n");   
  
  fclose(fp4); fclose(fp5); fclose(fp6); fclose(fp7);
#endif
  free(xtu_sum); free(pitu_sum); 
}

void writefile_threshold_aggressiveness(){
  int h, j;  
  char dxidata[200], dsidata[200], pun_idata[200], pun_jdata[200];
  FILE *fp1 = NULL, *fp2 = NULL, *fp3 = NULL, *fp4 = NULL;   
  prep_file(dxidata, "dxi.dat"); prep_file(dsidata,"dsi.dat"); prep_file(pun_idata, "pun_i.dat"); prep_file(pun_jdata, "pun_j.dat");
  fp1 = fopen(dxidata, "w"); fp2 = fopen(dsidata, "w"); fp3 = fopen(pun_idata, "w"); fp4 = fopen(pun_jdata, "w"); 
  
    
#if ALLDATAFILE
  int i;
  double STU = (T-SUMT)/SKIP;
  FILE *fp5 = NULL, *fp6 = NULL, *fp7 = NULL, *fp8 = NULL;
  char dxisumdata[200], dsisumdata[200], punisumdata[200], punjsumdata[200];
  double *dxitu_sum, *dsitu_sum, *punitu_sum, *punjtu_sum; // holds average of summary of results for last few considered time unit for summary    
  dxitu_sum  = calloc(N,sizeof(double));                   // hold average value of thresholds for each rank for last few generations considered for summary
  dsitu_sum  = calloc(N, sizeof(double));                  // hold average value of aggressiveness for each rank for last few generations considered for summary
  punitu_sum  = calloc(N, sizeof(double));
  punjtu_sum  = calloc(N, sizeof(double));
  
  prep_file(dxisumdata, "dxisum.dat"); prep_file(dsisumdata,"dsisum.dat");
  prep_file(punisumdata, "punisum.dat"); prep_file(punjsumdata, "punjsum.dat"); 
  fp5 = fopen(dxisumdata, "w"); fp6 = fopen(dsisumdata, "w");  fp7 = fopen(punisumdata, "w"); fp8 = fopen(punjsumdata, "w");
#endif
  // write headers
  fprintf(fp1, "%d\t", 0);
  fprintf(fp2, "%d\t", 0);
  fprintf(fp3, "%d\t", 0);
  fprintf(fp4, "%d\t", 0);
  for(h = 0; h < N; h++){
    fprintf(fp1, "%d\t", h + 1);
    fprintf(fp2, "%d\t", h + 1);
    fprintf(fp3, "%d\t", h + 1);
    fprintf(fp4, "%d\t", h + 1);
  }    
  fprintf(fp1, "\n");    
  fprintf(fp2, "\n");  
  fprintf(fp3, "\n");  
  fprintf(fp4, "\n"); 

  // write data
  for(j = 0; j < (int)(T/SKIP); j++){
    fprintf(fp1, "%d\t", (int)(j*SKIP));
    fprintf(fp2, "%d\t", (int)(j*SKIP)); 
    fprintf(fp3, "%d\t", (int)(j*SKIP));
    fprintf(fp4, "%d\t", (int)(j*SKIP));   
    for(h = 0; h < N; h++){
      dximean[j][h] = dximean[j][h] / (double)Runs; // average threshold per rank over all groups over number of runs
      dsimean[j][h] = dsimean[j][h] / (double)Runs;
      pun_imean[j][h] = pun_imean[j][h] / (double)Runs;  
      pun_jmean[j][h] = pun_jmean[j][h] / (double)Runs;      
      fprintf(fp1, "%.4lf  ", dximean[j][h]);
      fprintf(fp2, "%.4lf  ", dsimean[j][h]);  
      fprintf(fp3, "%.4lf  ", pun_imean[j][h]);        // payoff reduction due to punishing others
      fprintf(fp4, "%.4lf  ", pun_jmean[j][h]);        // payoff reduction due to punishment from others               
    }          
    fprintf(fp1, "\n");        
    fprintf(fp2, "\n");   
    fprintf(fp3, "\n");   
    fprintf(fp4, "\n"); 
#if ALLDATAFILE
    if(j < STU) continue;
    for(h = 0; h < N; h++){
      dxitu_sum[h]  += dximean[j][h];
      dsitu_sum[h]  += dsimean[j][h];
      punitu_sum[h] += pun_imean[j][h]; 
      punjtu_sum[h] += pun_jmean[j][h];
    }
#endif
  }   
  fclose(fp1); fclose(fp2); fclose(fp3); fclose(fp4);
#if ALLDATAFILE  
  double st = SUMT/(double)SKIP;
  // write summary of efforts and payof to file  
  for(i = 0; i < N; i++){
    dxitu_sum[i] = dxitu_sum[i]/st;
    dsitu_sum[i] = dsitu_sum[i]/st;
    punitu_sum[i] = punitu_sum[i]/st;
    punjtu_sum[i] = punjtu_sum[i]/st;
   
    fprintf(fp5, "%.6f\t", dxitu_sum[i]);
    fprintf(fp6, "%.6f\t", dsitu_sum[i]);    
    fprintf(fp7, "%.6f\t", punitu_sum[i]);
    fprintf(fp8, "%.6lf  ", punjtu_sum[i]); 
  }  
  fprintf(fp5, "\n");
  fprintf(fp6, "\n");   
  fprintf(fp7, "\n");
  fprintf(fp8, "\n");
  fclose(fp5); fclose(fp6); fclose(fp7); fclose(fp8); 
  free(dxitu_sum); free(dsitu_sum); free(punitu_sum); free(punjtu_sum);
#endif
  
}

void plotall(int m)
{
    // write data to file
  char xdata[200], fdata[200], pi_gdata[200], dxidata[200], dsidata[200], pun_idata[200], pun_jdata[200];
  char title[200], xpng[200];
  int datacolumn = N+1;
  prep_file(xdata, "e.dat"); prep_file(fdata,"f.dat"); prep_file(pi_gdata, "pi_g.dat");  
  prep_file(dxidata, "dxi.dat"); prep_file(dsidata,"dsi.dat"); prep_file(pun_idata, "pun_i.dat"); prep_file(pun_jdata, "pun_j.dat");  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "n= %d, b= %.2f, c= %.2f, al= %.2f, be= %.2f, dl= %.2f, eta= %.3f, k= %d, X0= %.2f, w= %.2f", N, B, C, Alpha, Beta, Delta, Eta, K, X0, Omega);
  
  fprintf(gp, "set key outside \n");        
  if(m){
    sprintf(xpng, "ef_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2fw%.2f.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0, Omega); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 dashed %d \n", 0);
  }
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 2,3 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set ylabel 'efforts' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, xdata);
  fprintf(gp, "set ylabel 'threshold' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, dxidata);
  fprintf(gp, "set ylabel 'aggressiveness' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, dsidata);
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "set style line 5 lc rgb '#808080' lt 2 lw 2 pi -1 ps 1.0\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader, '%s' using 1:%d with line ls 5 title columnheader \n", datacolumn, fdata, fdata, datacolumn+1);  
  
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel 'cost' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, pun_idata);
  
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel 'punishment' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, pun_jdata);
  
  fprintf(gp, "unset multiplot \n");
  
  fflush(gp); 
  pclose(gp);
}

void plotallIndividualRun(int r, int m)
{
    // write data to file
  char xdata[200], fdata[200], pi_gdata[200], dxidata[200], dsidata[200], pun_idata[200], pun_jdata[200];
  char title[200], xpng[200], str[200];
  int datacolumn = N+1;
  sprintf(str, "e%d.dat", r); prep_file(xdata, str);
  sprintf(str, "f%d.dat", r); prep_file(fdata, str); 
  sprintf(str, "pi_g%d.dat", r); prep_file(pi_gdata, str);  
  sprintf(str, "dxi%d.dat", r); prep_file(dxidata, str);
  sprintf(str, "dsi%d.dat", r); prep_file(dsidata,str);
  sprintf(str, "pun_i%d.dat", r); prep_file(pun_idata, str); 
  sprintf(str, "pun_j%d.dat", r); prep_file(pun_jdata, str);  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "Run:%d, n= %d, b= %.2f, c= %.2f, al= %.2f, be= %.2f, dl= %.2f, eta= %.3f, k= %d, X0= %.2f, w= %.2f", r, N, B, C, Alpha, Beta, Delta, Eta, K, X0, Omega);
   
  
  fprintf(gp, "set key outside \n");      

  if(m){
    sprintf(xpng, "ef_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2fw%.2f%d.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0, Omega, r); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 %d \n", 0);
  }
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 2,3 title '%s' \n", title);    // set subplots layout
    
  fprintf(gp, "set ylabel 'efforts' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, xdata);
  fprintf(gp, "set ylabel 'threshold' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, dxidata);
  fprintf(gp, "set ylabel 'aggressiveness' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, dsidata);
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader, '%s' using 1:%d with line lt 3 lc -1 title columnheader \n", datacolumn, fdata, fdata, datacolumn+1);  
  
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel 'cost' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, pun_idata);
  
  fprintf(gp, "set yrange [0:] \n");
  fprintf(gp, "set ylabel 'punishment' \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines title columnheader \n", datacolumn, pun_jdata);
  
  fprintf(gp, "unset multiplot \n"); 

  fflush(gp); 
  pclose(gp);
}

void plotTraits(int r, int m)
{
  char tdata[200], str[100], xpng[200];
  sprintf(str, "traits%d.dat", r); prep_file(tdata, str);  
  int i;
  FILE *gp = popen("gnuplot -persistent", "w");
  fprintf(gp, "set key outside \n");        
  if(m){
    sprintf(xpng, "ef_n%02db%05.2fc%.2fga%.2fbe%.2fal%.2fdl%.2fe%.3fk%02dx0%.2fw%.2ftraits%d.png", N,B,C,Gamma,Beta,Alpha,Delta,Eta,K, X0, Omega, r); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 %d \n", 0);
  }
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set xlabel 'threshold' \n");
  fprintf(gp, "set ylabel 'aggressiveness' \n");
  fprintf(gp, "set multiplot layout 2,%d title 'Run:%d' \n", N/2, r);    // set subplots layout
  for( i = 0; i < N; i++){
    fprintf(gp, "set title 'rank %d' \n", i+1);
    fprintf(gp, "plot '%s' using %d:%d with points lw 1 pt 7 notitle\n", tdata, 2*i+1, 2*i+2);
  }
  fprintf(gp, "unset multiplot \n");

  fflush(gp);
  pclose(gp);
} 
