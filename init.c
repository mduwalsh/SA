
void *Malloc(size_t size)
{
  void *p = malloc(size);

  if (p) return p;
  printf("malloc failed\n");
  _exit(2);
}

#define malloc(x) Malloc(x) 

void *Calloc(size_t nmemb, size_t size)
{
  void *p = calloc(nmemb, size);

  if (p) return p;
  printf("calloc failed\n");
  _exit(3);
}

#define calloc(x,y) Calloc(x,y) 

#define EXPECT(a,b,c) if ((a) != fscanf(f, b"%*[ ^\n]\n",c)){ fclose(f); printf("Error: %s\n",b); return 1; }

int read_config(char *file_name)
{
  FILE *f;

  if (!(f = fopen(file_name,"r"))) return 1;
  EXPECT(1, "unsigned Seed      = %u;", &Seed);
  EXPECT(1, "int      Runs      = %d;", &Runs);
  EXPECT(1, "int      N         = %d;", &N);
  EXPECT(1, "int      G         = %d;", &G);  
  EXPECT(1, "double   T         = %lf;", &T);
  
  EXPECT(1, "double   PE[0]     = %lf;", PE+0);
  EXPECT(1, "double   PE[1]     = %lf;", PE+1);
  EXPECT(1, "double   PE[2]     = %lf;", PE+2);
  
  EXPECT(1, "double   F0        = %lf;", &F0);
  EXPECT(1, "double   B         = %lf;", &B);
  EXPECT(1, "double   C         = %lf;", &C);
  EXPECT(1, "double   X0        = %lf;", &X0);
  EXPECT(1, "double   Discount  = %lf;", &Discount);
  EXPECT(1, "double   Omega     = %lf;", &Omega);
  
  EXPECT(1, "double   Alpha     = %lf;", &Alpha);
  EXPECT(1, "double   Beta      = %lf;", &Beta);
  EXPECT(1, "double   Gamma     = %lf;", &Gamma);
  
  EXPECT(1, "double   Delta     = %lf;", &Delta);
    
  EXPECT(1, "int      PROTOCOL  = %d;", &PROTOCOL);
  EXPECT(1, "int      K         = %d;", &K);
  EXPECT(1, "double   Eta       = %lf;", &Eta);  
  
  EXPECT(1, "double   E         = %lf;", &E);
  EXPECT(1, "double   S0        = %lf;", &S0);
  EXPECT(1, "double   Phi       = %lf;", &Phi);
  EXPECT(1, "double   e         = %lf;", &e);
  
  EXPECT(1, "double   Mu        = %lf;", &Mu);
  EXPECT(1, "double   Sigma     = %lf;", &Sigma);
  EXPECT(1, "double   Sigma_B   = %lf;", &Sigma_B);
  EXPECT(1, "double   Sigma_dxi = %lf;", &Sigma_dxi);
  EXPECT(1, "double   Sigma_dsi = %lf;", &Sigma_dsi);  

  fclose(f); return 0;
}
#undef EXPECT

void init()
{
  int i,g;  
  // set polity size array (no. of groups in each polity) and group size array (no. of individuals in each group)
  //PS = malloc(NP*sizeof(int));
  GS = malloc(G*sizeof(int));
  for(i = G; i--;)                 // through all groups
    GS[i] = N;                // set to initial group size
    
  
  // set strategies and roles for all individuals
  x = malloc(G*sizeof(double *)); // mem allocated along groups for strategies   
  V = malloc(G*sizeof(double *));  // mem allocated along groups for valuations
  X = calloc(G, sizeof(double));
  P = calloc(G, sizeof(double));
  S = malloc(G*sizeof(double *));
  dxi = malloc(G*sizeof(double *));
  dsi = malloc(G*sizeof(double *));
  
  for(g = 0; g < G; g++){     // through each group 
    x[g] = malloc( N*sizeof(double));  // mem allocated along individuals for each group
    V[g] = malloc( N*sizeof(double));    // mem allocated along individuals of each group for valuations  
    S[g] = calloc( N, sizeof(double));
    dxi[g] = calloc( N, sizeof(double));
    dsi[g] = calloc( N, sizeof(double));    
  }   
  
  // allocate mem for payoffs individual, group and polity
  //pp = malloc(NP*sizeof(double));                      // for accumulated payoffs of polity
  pi_g    = calloc(G,sizeof(double));                    // mem along polities for accumulated payoffs of groups
  pi      = malloc(G*sizeof(double*));                   // mem along polities for payoffs of individuals
  pun_i   = malloc(G*sizeof(double*));                   // punishment due to punishing others
  pun_j   = malloc(G*sizeof(double*));                   // punishment due to punishment from others
  pun_g   = calloc(G,sizeof(double));                    // group punishment
  Api     = malloc(G*sizeof(double*));                   // accumulated payoffs
  for(g = 0; g < G; g++){
    pi[g]    = malloc(N*sizeof(double));                 // mem along individuals for payoffs of individuals in each group
    pun_i[g] = calloc(N,sizeof(double));
    pun_j[g] = calloc(N,sizeof(double));
    Api[g]   = malloc(N*sizeof(double));
  }
  
  Event_dist = allocdist(EVENTS);
  
  // allocate mem for candidate strategies
  CS = malloc((K)*sizeof(double *));
  #if FORESIGHT
    FCS = calloc( K, sizeof(double[TRAITS]));
  #endif
    for(i = K; i--;)
    CS[i] = calloc(TRAITS, sizeof(double));     
  
  // for statistical variables
  // statistical variables
  xmean   = malloc( (int)(T/SKIP) * sizeof(double*) ); // hold mean value of effort 
  pimean   = malloc( (int)(T/SKIP) * sizeof(double*) ); // hold mean value of relative fertility 
  dximean   = malloc( (int)(T/SKIP) * sizeof(double*) ); // hold mean value of threshold 
  dsimean   = malloc( (int)(T/SKIP) * sizeof(double*) ); // hold mean value of aggressiveness 
  pi_gmean   = calloc( (int)(T/SKIP), sizeof(double) );
  pun_imean   = malloc( (int)(T/SKIP) * sizeof(double*) ); // hold mean value of punishment
  pun_jmean   = malloc( (int)(T/SKIP) * sizeof(double*) ); // hold mean value of punishment
  pun_gmean   = calloc( (int)(T/SKIP), sizeof(double) );
  
  for(i = 0; i < (int)(T/SKIP); i++){
    xmean[i] = calloc(N,sizeof(double)); // to hold mean value of effort per rank after each event in SKIP gap
    pimean[i] = calloc(N,sizeof(double));  
    dximean[i] = calloc(N,sizeof(double)); // to hold mean value of threshold per rank after each event in SKIP gap
    dsimean[i] = calloc(N,sizeof(double)); // to hold mean value of threshold per rank after each event in SKIP gap    
    pun_imean[i] = calloc(N,sizeof(double));
    pun_jmean[i] = calloc(N,sizeof(double));
  }  
  
  Ge_dist = allocdist(G);       // allocate memory for distribution for group extinction
  Gr_dist = allocdist(G);       // allocate memory for distribution for group replication
  
  Sij = malloc( G * sizeof(double **));          // Sij matrix for all groups
  Delta_ij = malloc( G * sizeof(double **));    // Delta matrix for all groups
  Sij_avg = malloc( N * sizeof(double *));          // Sij matrix avg
  Delta_ij_avg = malloc( N * sizeof(double *));    // Delta matrix avg
  for( g = 0; g < G; g ++){
    Sij[g] = malloc( N * sizeof(double *));          // Sij matrix
    Delta_ij[g] = malloc( N * sizeof(double *));    // Delta matrix
    for(i = 0; i < N; i++){
      Sij[g][i] = calloc(N, sizeof(double));
      Delta_ij[g][i] = calloc(N, sizeof(double));
    }
  }
  for(i = 0; i < N; i++){
    Sij_avg[i] = calloc(N, sizeof(double));
    Delta_ij_avg[i] = calloc(N, sizeof(double));
  }
}

void cleanup()
{
  int i,g;    
  // free strategies and roles for individuals  
  for(g = 0; g < G; g++){    
    free(x[g]);           // free mem allocated along individuals for each group
    free(V[g]);           // free mem allocated along individuals of each group for valuations
    free(S[g]);
    free(dxi[g]);
    free(dsi[g]);   
    free(pi[g]);          // free mem along individuals for payoffs of individuals in each group
    free(pun_i[g]);
    free(pun_j[g]);
    free(Api[g]);
  }  
  free(x); free(V); free(X); free(P); free(S); free(dxi); free(dsi); 
  free(pi_g); free(pi); free(pun_i); free(pun_j); free(pun_g);
  free(Api);
  // free group size array (no. of individuals in each group)  
  free(GS);  
  
  freedist(Event_dist);          // free distribution for events
  
  for(i = K; i--;)               // free candidate strategies array
    free(CS[i]);
  free(CS);    
  #if FORESIGHT
    free(FCS);
  #endif
  
  // statistical variables
  for(i = 0; i < (int)(T/SKIP); i++){
    free(xmean[i]);
    free(pimean[i]);    
    free(dximean[i]);
    free(dsimean[i]);       
    free(pun_imean[i]);  
    free(pun_jmean[i]);  
  }
  free(xmean); free(pimean); free(dximean); free(dsimean); free(pi_gmean); 
  free(pun_imean); free(pun_jmean); free(pun_gmean); 
  
  freedist(Ge_dist);
  freedist(Gr_dist);
  
  for(g = 0; g < G; g++)
  {
    for( i = 0; i < N; i++){
      free(Sij[g][i]);
      free(Delta_ij[g][i]);
    }
    free(Sij[g]); free(Delta_ij[g]);    
  }
  for( i = 0; i < N; i++){
    free(Delta_ij_avg[i]);
    free(Sij_avg[i]);
  }
  free(Sij); free(Delta_ij); free(Sij_avg); free(Delta_ij_avg);
    
}



