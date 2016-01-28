#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define DEBUG 1
#if DEBUG
#define __USE_GNU
#include <fenv.h>       /* enable floating exceptions */
unsigned long Where;    // debugging counter
#endif

#include "rand.c"

#define Accumulatedpayoff 0  // use accumulated payoff or not
#define PUNISH 1
#define AGGR   0
#define FORESIGHT 1

#define CLUSTER 0            // is this simulation on cluster yes or no (1 or 0)
#define EGALITARIAN 0        // simulation for egalitarian yes or no (1 or 0)
#define DISP_MATRIX 0        // show statistics matrix 
#define ALLDATAFILE 0        // if 1, generate all data files for individual runs and summary too 
#define GRAPHS      0

#define SKIP 10.0           // skip time units after which statistics are calculated

#define EVENTS 3             // no. of events
#define TRAITS 3             // no. of strategy traits
#define NGUT 2               // no. of groups playing us vs them game
#define SUMT        2000     // last considered generations for summary
#define MINPAYOFF 0.00001


#define MAX(x,y)         (((x)>(y))? (x): (y))
#define MIN(x,y)         (((x)<(y))? (x): (y))

// Global variables
/** configuration parameters **/
unsigned Seed;
int Runs;
int G;          // no. of groups
int N;          // no. of individuals
int K;          // no. of candidate strategies
int PROTOCOL;   // 0: K candidate strategies, 1: 1 candidate strategies and best of two selection
double T;       // max time of simulation
double C;       // cost parameter
double B;       // expected benefit for a game
double F0;
double Eta;      
double Delta;   // within group inequality
double Alpha;   // cost exponent: -cx^Alpha
double Beta;    // group strength exponent: S^Beta
double Gamma;   // synergicity exponent: S=[sum x_i^(1/Gamma)]^Gamma
double Sigma;            // sigma for distribution of mutation
double Sigma_B;          // sigma for distribution of B
double Sigma_dsi;        // sigma for distribution of aggressiveness dsi
double Sigma_dxi;        // sigma for distribution of threshold
double Mu;               // probability of strategy update
double X0;        // base effort parameter in us vs nature game
double S0;        // cost coefficient of punishing
double Phi;       // exponential cost coefficient of punishing
double E;         // cost coefficient fo punishment
double e;         // error parameter in strength comparison
double Discount;  // discount used in calculating accumulated payoff
double Omega;     // weight used in calculating expected payoff from level 1 and level 2 expected payoff

/** simulation / object attributes **/
double **x;             // contribution strategy for every individual; strategy x[g][i] for individual i of g group 
double **pi;            // payoffs; pi[g][i] for individual i of group g  // it is immediate payoff
double **V;             // valuations of individuals within a group; V[g][i] for individual i for group g in polity p
double **S;             // strength of each individual
double **dxi;           // punishment threshold delta_x_i
double **dsi;           // aggressiveness delta_s_i 
double *X;              // group effort of each group 
double *P;              // probability of success of each group
double *pi_g;           // payoff for group g pi_g[p][g] in polity p
double *pi_p;           // payoff for polity p pi_p[p]
double **Api;           // accumulated payoff

int *PS;                // polity size array
int *GS;                // group size array
int *PL;                // polity label array
double PE[EVENTS];      // Probability of selecting an event; PE[j] for event j

/** supporting variables **/
dist *Event_dist;        // distribution for event
double **CS;             // candidate strategies
dist *Ge_dist, *Gr_dist;
double (*Tx)[TRAITS];    // to hold original state of strategies for a group to be used in update of strategies using quantal response 
double (*FCS)[TRAITS];   // strategy vector with modified strategy for i and candidate foresight strategy vector for j's
double (*Strat)[TRAITS]; // holds strategies prior to update
double *Xmax;            // holds maximum effort allowed for individuals in a group
int (*Lookup)[3];        // holds punishment lookup table to be used in update_stragegy
int *Lidx;               // holds index for lookup table in update_stragey() 
double (*FS)[TRAITS];    // holds foresights for all individuals
int (*Lt)[3];            // holds punishment lookup table to be used in foresight
int *Idx;                // holds indices for lookup table to be used in foresight()       

double Bo;               // original (expected mean) value of B
double Lambda;           // 1/(summation of P[j])
double I_Gamma;          // 1.0/Gamma
double I_Alpha;          // 1.0/Alpha
double X0_Be;            // pow(X0, Beta)
double GaBe;             // Gamma * Beta

// statistical variables
double **xmean;          // hold mean value of effort groups
double **pimean;          // hold mean value of payoff all groups
double **dximean;        // hold mean value of thresholds
double **dsimean;        // hold mean value of aggressiveness
double *pi_gmean;       // hold group payoff of all groups for every skip data point
double **pun_i, **pun_j;            // hold punishment incurred by each individual
double **pun_imean, **pun_jmean;            // average punishment incurred by each rank individual
double *pun_g;           // hold punishment incurred by all individuals in a group
double *pun_gmean;       // average punishment incurred by group

double ***Sij, ***Delta_ij;
double **Sij_avg, **Delta_ij_avg;

// prepare filename to write data
void prep_file(char *filevar, char *apnd)
{  
  sprintf(filevar, "n%db%.2fc%.2fa%.2fb%.2fd%.2fg%.2fe%.3fk%dx0%.2fw%.2f%s", N,Bo,C,Alpha,Beta,Delta,Gamma,Eta,K,X0,Omega, apnd);  
}

#include "init.c"        // variables mallocs and frees; read configuration 
#include "dataplot.c"

// generates random number following exponential distribution
double randexp(double lambda)
/*
 * lambda: rate parameter of exponential distribution
 */
{
  return -log(1.0-U01()) / lambda;  
}

// merge sort
inline void merge (double *a, int n, int m) {
    int i, j, k;
    double *x = malloc(n * sizeof(double));
    for (i = 0, j = m, k = 0; k < n; k++) {
        x[k] = j == n      ? a[i++]
             : i == m      ? a[j++]
             : a[j] < a[i] ? a[j++]
             :               a[i++];
    }
    for (i = n; i--;) {
        a[i] = x[i];
    }
    free(x);
}
 
inline void merge_sort (double *a, int n) {
    if (n < 2)
        return;
    int m = n>>1;   // divide by 2
    merge_sort(a, m);
    merge_sort(a + m, n - m);
    merge(a, n, m);
}

// sets initial strategies, roles, valuations and payoffs
void set_init_xrvpi()
{
  int g, i;
  double v_sum; 

  for(g = 0; g < G; g++){     // through each group                               
    for(pi_g[g] = 0, v_sum = 0, i = 0; i < N; i++){   // through each individual in a group      
      x[g][i] = U01()*0.05;              // initialize strategy vector for all individuals      
      pi[g][i] = 1.0;                           // set initial payoffs for each individual to 1   
      Api[g][i] = 0;
      // start calculation of valuations of each individual
#if !PUNISH  
      V[g][i]  = pow(N - i, Delta);                // rank based valuations
#if EGALITARIAN
      V[g][i] = 1;
#endif
      v_sum += V[g][i];                            // sum of valuations of individuals in a group for normalization 
#endif      
    }              
    pi_g[g] = 1;
#if !PUNISH
    // normalize valuations
    for(i = 0; i < N; i++){  // through each individual 
      V[g][i] /= v_sum;         // normalize valuations      
    }     
#endif
  }
  
#if PUNISH

  int j;
  double *ss = malloc(N*sizeof(double));            // scratch space for strength sorting
    
  for(g = 0; g < G; g++){     // through each group          
    
    for( i = 0; i < N; i++){ ss[i] = U01();  /*ss[i] = 1. - i/(double)N;*/}
    merge_sort(ss, N);   // sort strength                     
    for(v_sum = 0, i = 0; i < N; i++){   // through each individual in a group  
      S[g][i] = ss[N-i-1];                      // assign strength in descending order 
      dxi[g][i] = U01()*0.05;
#if AGGR
      dsi[g][i] = 1 - U01()*0.05;              // setting value near to 1   
#else
      dsi[g][i] = 0;
#endif
      
#if EGALITARIAN
      S[g][i] = 1;
      V[g][i] = 1;                              // same valuations
#else
      V[g][i] = pow(S[g][i], Beta);             // strength based valuations
#endif
      v_sum += V[g][i];
    }
    // normalize valuations
    for(i = 0; i < N; i++){  // through each individual 
      V[g][i] /= v_sum;         // normalize valuations      
    }  
  }

  // computing Sij matrix for all groups
  for(g = 0; g < G; g++){
    for(i = 0; i < N; i++){
      for(j = 0; j < N; j++){
	if( i == j) continue;
	Sij[g][i][j] = S0 *( exp( Phi*(S[g][j] - S[g][i]) ) );
      }   
    }
  }  
  free(ss);
#endif
}

// updates payoff of individuals in a group g in polity p and accumulated payoff of group and polity
void updateIndividualPayoffAfterProduction(int g, int ng, double pbs)
/* g: group index in a polity of a group
 * ng: no. of groups playing game
 * pbs: probability of success of a group
 */
{
  int i;
  double bgp;
  bgp = B*ng*pbs;
  for(i = GS[g]; i--;){
    pi[g][i] = F0*(1. + (bgp*V[g][i]) - (C*pow(x[g][i], Alpha)) );      
#if Accumulatedpayoff
    Api[g][i] = (1-Discount)*Api[g][i] + pi[g][i];   
    Api[g][i] = MAX(Api[g][i], MINPAYOFF);
#else
    pi[g][i] = MAX(pi[g][i], MINPAYOFF);
    Api[g][i] = pi[g][i];
#endif
  }  
}

// update payoff of group
inline void updateGroupPayoff(int g)
{
  int i; double sum = 0;  
  for(i = GS[g]; i--;){
#if Accumulatedpayoff
    sum += Api[g][i];
#else
    sum += pi[g][i];
#endif
  }
  pi_g[g] = sum/GS[g];                       // average group payoff
  //pi_g[g] = sum;
}

// returns total group strategy value for a group
double get_X(int g)
/*
 * g: group index of a group in a polity
 */
{
  double xx = 0;
  int i;
  for(i = 0; i < GS[g]; i++){ // move through all individuals in the group g 
   xx += pow(x[g][i], I_Gamma);           
  }
  return pow(xx, Gamma); 
}

// returns probability of success of group playing CA1 (us vs nature)
inline void update_XP_CA1(int g)
/*
 * g: group index of a group in a polity
 */
{
  double xx;
  X[g] = get_X(g);                                                // update group effort
  xx = pow(X[g], Beta);
  P[g] = xx/( xx + X0_Be );       // update probability of success
}


void update_XP_CA2(int *g, int ng)
/* g: array of group indices within polity of groups chosen
 * ng: no. of groups chosen or selected for playing us vs them game
*/
{
  int i;
  double sum;
  for(sum = 0, i = ng; i--;){
    X[g[i]] = get_X(g[i]);                    // update group effort for group g[i]
    sum += pow(X[g[i]], Beta);
  }
  sum = MAX(sum, 0.001);
  for(i = ng; i--;){
    P[g[i]] = pow(X[g[i]], Beta) / sum;       // update probability of success of     
  }  
}

static inline void perturb(int g, int i, double *c, double xx, double t, double a, int m)
/* 
 * g: group index
 * i: individual index
 * c: candidate strategy vector pointer 
 * x: effort
 * t: threshold
 * a: aggressiveness
 * m: to update threshold or not
 * */
{
  c[0] = normal(xx, Sigma);                                             // perturbation in effort
  if( c[0] > Xmax[i]) c[0] = Xmax[i];
  else if( c[0] < 0.0) c[0] = 0.0;
  
#if PUNISH
  if(AGGR || m){                                                        // if m !=1, it means i does not punish and use old threshold
    c[1] = normal(t, Sigma_dxi);                                        // perturb threshold 
    c[1] = MAX(c[1],0);
  }
  else{
    c[1] = t;
  }
#if AGGR
  c[2] = normal(a, Sigma_dsi);                                          // perturbation in aggressiveness     
  c[2] = MAX(c[2], 0.0);  
#endif
#endif  
}

static inline void exp_initdist(dist *d, double m) // m = max d->p[i]
{
  int i = d->n;  double s = 0;
  while(i--) s += (d->p[i] = exp(d->p[i] - m));
  initdist(d,s);
}

#if PUNISH

#if AGGR
#define strengthcheck(j,i) ((S[g][j] < dsi[g][i])? 0: 1)
#else
#define strengthcheck(j,i) ((S[g][j] < (S[g][i] + e*normal(0,1)))? 0: 1)
#endif

// punishment
void punish(int g)
{
    int i, j, p, q, lsize, (*punLookup)[2], *l_idx;
    double dl, pdl;
    lsize = GS[g]*(GS[g]-1);                                             // punishment lookup table size
    punLookup = malloc(lsize * sizeof(int[2]));
    l_idx = malloc(lsize*sizeof(int));
    
    for(p = 0, i = 0; i < GS[g]; i++){                                   // prepare lookup table of punisher and punishee pair   
      for(j = 0; j < GS[g]; j++){
	if(i != j){
	  punLookup[p][0] = i;
	  punLookup[p][1] = j;
	  p++;
	}
      }
    }      
    for( i = lsize; i--;)                                                 // initialize index of lookup table for punishment pairs
      l_idx[i] = i;
    
    for( i = 0; i < N; i++){                                              // reinitializing values to zero for next event
      pun_i[g][i] = 0;
      pun_j[g][i] = 0;      
    }
    i = lsize;
    while(i > 0){                                                         // random permutation of indices of punishment lookup table 
      p = rnd(i); q = l_idx[--i]; l_idx[i] = l_idx[p]; l_idx[p] = q;
    }
    for( p = 0; p < lsize; p++){
      i = punLookup[l_idx[p]][0];                                          // punisher
      j = punLookup[l_idx[p]][1];                                          // punishee      
      if(!strengthcheck(j,i)){                                             // if j is not strong enough
	dl = dxi[g][i] - x[g][j];
	if( dl > 0){                                                       // if j is not making enough effort
	  dl = MIN(Api[g][j], E*dl);                                       // punishment j can pay
	  pdl = dl*Sij[g][i][j];                                           // cost of punishing for i
	  if(Api[g][i] < pdl){                                             // if i doesn't have enough payoff ot punish, punish partially
	    pdl = Api[g][i];
	    dl = pdl/Sij[g][i][j];                                         // partial punishment to j by i
	  }
	  Api[g][j] -= dl;                                                 // punishment deduction from payoff of j
	  Api[g][i] -= pdl;                                                // cost of punishing deduction from payoff of i
	  Api[g][j] = MAX(Api[g][j], MINPAYOFF);
	  Api[g][i] = MAX(Api[g][i], MINPAYOFF);
	  pun_j[g][j] += dl;                                               // log punishment from others
	  pun_i[g][i] += pdl;                                              // log cost of punishing
	  Delta_ij[g][i][j] = dl;                                          // udpate punishment matrix
	}
    }           
  }
  free(punLookup);  free(l_idx);
}

#undef strengthcheck

#define effort(j)         ((j==i)?     cs[0]: ust[j][0])
#define threshold(j)      ((j==i)?     cs[1]: ust[j][1])
#define aggressiveness(j) ((j==i)?     cs[2]: ust[j][2])
#if AGGR
#define strengthcheck(j,i) ((S[g][j] < aggressiveness(i))? 0: 1)
#else
#define strengthcheck(j,i) (!lt[idx[p]][2])
#endif
inline double expectedPayoffAfterPunishment(int g, int i, double *cs, double f, int (*lt)[3], int *idx, double (*ust)[TRAITS])
/*
 * g: group index
 * i: individual index
 * cs: candidate strategy vector
 * f: new payoff
 * lt:  punishment lookup table
 * idx: indices of lookup table
 * ust: unchanged current strategy array
 */
{
  int j, s, p, q;
  double dl; //, pdl;
  s = 2*(GS[g]-1);                                                               // punishment lookup table size  
  
#if AGGR
  double th;
  th = cs[1];                                                            // temporarily store candidate threshold strategy
  cs[1] = ust[i][1];                                                     // use old threshold unless it is aggressive enough to punish other individuals
#endif
  for(p = 0; p < s; p++){                                                // start punishing
    q = lt[idx[p]][0];                                                   // punisher
    j = lt[idx[p]][1];                                                   // punishee    

    if(!strengthcheck(j,q))                                              // if strength of j (punishee) less than agressiveness of q (punisher)
    {
#if AGGR                                                                 // only needed when AGGR = 1, for AGGR = 0, it is already decided whether to update threshold 
      if(q == i){                                                        // if punisher is i, update threshold
	cs[1] = th;
      }
#endif
      dl = threshold(q) - effort(j);      
      if( dl > 0){                                                       // if j is not making enough effort
	if(q == i){                                                      // if i is punisher	  
	  f -= E*dl*Sij[g][q][j];                                        // cost of punishing paid by q  
	  f = MAX(0.00, f);
	}
	else{                                                            // i is punishee
	  f -= MIN(f, E*dl);                                             // punishment from j
	}
      }
    }
  }
  return f;                                                               // return expected payoff
}

#undef effort
#undef threshold
#undef aggressiveness
#undef strengthcheck
#endif 

inline int prepPunLookupTableUpdateStrategy(int g, int i, int (*lt)[3], int *idx)
{
  static int p, s, j, c, q;
  s = 2*(GS[g]-1);                                                               // punishment lookup table size
  for( c = 0, p = 0, j = 0; j < N; j++){                                         // prepare lookup table for punisher and punishee pair
    if( i == j) continue;
    lt[p][0] = i;                                                                // punisher
    lt[p][1] = j;                                                                // punishee
#if AGGR == 0
    if(S[g][j] < (S[g][i] + e*normal(0,1))){                                     // check if i is punishing any j
      lt[p++][2] = 1;                                                              // i punishes j 
      c = 1;                                                                     // if punisher is i, then mark c = 1 for update of threshold
    }
    else{
      lt[p++][2] = 0;                                                              // i does not punish j
    }
#else
    p++;
#endif   
    lt[p][0] = j;                                                                // punisher
    lt[p][1] = i;                                                                // punishee
#if AGGR == 0
    if(S[g][i] < (S[g][j] + e*normal(0,1))){                                     // check if any j is punishing i
      lt[p++][2] = 1;                                                              // j punishes i 
    }
    else{
      lt[p++][2] = 0;                                                              // j does not punish i
    }
#else
    p++;
#endif    
  }
  j = s;
  while(j--){
    idx[j] = j;                                                                  // initialize indices
  }    
  j = s;
  while(j > 0){                                                                  // random permutation of indices of punishment lookup table 
    p = rnd(j); q = idx[--j]; idx[j] = idx[p]; idx[p] = q;
  }
  return c;
}

double get_payoff(int g, int i, double x1, double x2, double bvg, double api, double *cs, int (*lt)[3], int *idx, double (*ust)[TRAITS])
{
  double f;
  f = F0*(1. + bvg*x1/(x1+x2) - ( C*pow(cs[0], Alpha) ) );                 // payoff after production of new strategy for cs  
#if Accumulatedpayoff
  f   = api + f;    
#endif
  f   = MAX(f, MINPAYOFF);   
#if PUNISH
  f = expectedPayoffAfterPunishment(g, i, cs, f, lt, idx, ust);             // payoff after punishment  
#endif
  f   = MAX(f, MINPAYOFF);   
  return f;
}

#if FORESIGHT
double foresight(int g, int i, int ng, double gx, double x2)
{
  int h, j, c, m;  
  double x1, tp, bvg, f, pgx, wx, api;
               
  // generate foresighted candidate strategies for every individual other than i and select one
  for( wx = 0, j = 0; j < GS[g]; j++){
    if(j != i){
      api = (1-Discount)*Api[g][j];
      bvg = B*ng*V[g][j];      
      pgx = gx - pow(Strat[j][0], I_Gamma);                                         // subtract effort of individual who is generating foresight strategy
      m = prepPunLookupTableUpdateStrategy(g, j, Lt, Idx);                          // prepare lookup for punishment interactions   
      
      FCS[0][0] = Strat[j][0];  FCS[0][1] = Strat[j][1]; FCS[0][2] = Strat[j][2];   // 1st candidate strategy to be present strategy
      for(tp = -1, c = 0, h = 0; h < K; h++){                                       // next K candidate foresight strategies 
	if(h){ // perturbation	
	  perturb(g, j, FCS[h], Strat[j][0], Strat[j][1], Strat[j][2], m); 
	}
	x1 = pgx + pow(FCS[h][0], I_Gamma);                                         // add effort of individual which is generating foresight strategy
	x1 = pow(x1, GaBe);                                                         // new group effort with beta exponent	  
	f = get_payoff(g, j, x1, x2, bvg, api, FCS[h], Lt, Idx, Strat);             // payoff
	if(tp < f){ tp = f; c = h;}                                                 // keep track of strategy with higher payoff   
      }
      // select one strategy using protocol
      FS[j][0] = FCS[c][0]; FS[j][1] = FCS[c][1]; FS[j][2] = FCS[c][2]; 
      wx += pow(FCS[c][0], I_Gamma);
    }
    // complete selection of foresight strategy for each member
  }  
  // 2nd level expected payoff  
  FS[i][0] = Strat[i][0]; FS[i][1] = Strat[i][1]; FS[i][2] = Strat[i][2];          // candidate strategy i
  
  wx += pow(FS[i][0], I_Gamma);
  wx = pow(wx, GaBe); 
  f = get_payoff(g, i, wx, x2, B*ng*V[g][i], (1-Discount)*Api[g][i], FS[i], Lookup, Lidx, FS);
  return f;
}
#endif

// update strategy of an individual using quantal response approach
void updateStrategy(int g, int i, int ng, double x2)
/*
 * g: group index
 * i: individual index
 * ng: no. of groups involved in a game
 * x2: sum of groups contribution except 'g' group with Beta exponent playing the game
 */
{
  int j, c, m;
  double gx, x1, f, tp, bvg;    
  double api = 0;
  
  gx = pow(X[g], I_Gamma);                                               // power of inverse of Gamma of total effort; gx = sum of pow(xi, I_gamma);
  gx = gx - pow(x[g][i], I_Gamma);                                       // subtract contribution of current effort    
  bvg = B*ng*V[g][i];  
  
  #if Accumulatedpayoff
    api = (1-Discount)*Api[g][i];                                    // accumulated payoff of i with discount   
  #endif     
  #if PUNISH
    m = prepPunLookupTableUpdateStrategy(g, i, Lookup, Lidx);                   // prepare lookup for punishment  
    #if FORESIGHT
      double f2;
      for( j = GS[g]; j--;){   
	if( j != i){                                                                // copy current strategies except for i
	  Strat[j][0] = Tx[j][0];                                                   // copy effort
	  Strat[j][1] = Tx[j][1];                                                   // copy threshold
	  Strat[j][2] = Tx[j][2];                                                   // copy aggressiveness  	
	   #if Accumulatedpayoff      
	    Xmax[j] = pow( ( (1-Discount)*Api[g][j] + (1.0 + B*V[g][j]) )/C, I_Alpha );
	  #else
	    Xmax[j] = pow( (1.0 + B*V[g][j])/C, I_Alpha);                                 // upper bound of x at each role 
	  #endif
	}    
      }      
    #endif
  #endif    
  for(tp = -1, c = 0, j = 0; j < K; j++){                                 // calculate expected payoff for next K candidate strategies and use exponent of it as probability of getting selected
    if(!j){
      CS[j][0] = x[g][i];  CS[j][1] = dxi[g][i]; CS[j][2] = dsi[g][i];   // 1st candidate strategy to be present strategy      
    }
    else{
      perturb(g, i, CS[j], x[g][i], dxi[g][i], dsi[g][i], m);              // perturb strategies       
    }
    x1 = gx + pow(CS[j][0], I_Gamma);                                   // group effort with new candidate strategy effort
    x1 = pow(x1, GaBe);                                                 // new group effort with beta power   
    f = get_payoff(g, i, x1, x2, bvg, api, CS[j], Lookup, Lidx, Tx);     // level 1 payoff
    #if PUNISH           
      #if FORESIGHT
	// copy strategies for i
	Strat[i][0] = CS[j][0];                                           // copy effort
	Strat[i][1] = CS[j][1];                                           // copy threshold
	Strat[i][2] = CS[j][2];                                           // copy aggressiveness  
	f2 = foresight(g, i, ng, x1, x2);                          // level 2 payoff 
	f = (1.0-Omega)*f + Omega*f2;                                     // expected payoff
      #endif        
      f = MAX(f, MINPAYOFF); 
    #endif    
    if(tp < f){ 
      tp = f;                                                             // largest payoff
      c = j;                                                              // highest payoff candidate selection
    }
  } 
  x[g][i] = CS[c][0];  dxi[g][i] = CS[c][1];  dsi[g][i] = CS[c][2];       // update strategy with strategy in index c of CS  
}


// play CA1 event (us vs nature game for one group)
void event_CA1()
{
  int i, g;  
  // choose a group randomly
  g = rnd(G); 
  // play us vs nature game
  B = normal(Bo, Sigma_B);                                     // randomly choose Benefit B
  B = MAX(B,0);                                                // check if B is less than zero
  update_XP_CA1(g);                                            // updates group effort and probability of success of the group
  updateIndividualPayoffAfterProduction(g, 1, P[g]);           // update payoff of individuals in the group and of the group 
#if PUNISH
  punish(g);
#endif
  updateGroupPayoff(g);
  
  Tx = malloc( GS[g]*sizeof(double [TRAITS]));  
  Xmax  = malloc(GS[g]*sizeof(double));
    
  // get state of strategies before update
  for(i = 0; i < GS[g]; i++){
    Tx[i][0]   = x[g][i];
    Tx[i][1] = dxi[g][i];
    Tx[i][2] = dsi[g][i];        
#if Accumulatedpayoff      
  Xmax[i] = pow( ( (1-Discount)*Api[g][i] + (1.0 + B*V[g][i]) )/C, I_Alpha );
#else
  Xmax[i] = pow( (1.0 + B*V[g][i])/C, I_Alpha);                                 // upper bound of x at each role 
#endif      
  }  
#if PUNISH
  Lookup = malloc(2*(GS[g]-1)*sizeof(int[3]));                                     // allocate memory for lookup table  
  Lidx   = malloc(2*(GS[g]-1)*sizeof(int));
  #if FORESIGHT
    Strat = malloc(GS[g]*sizeof(double[TRAITS]));
    FS = malloc(GS[g]*sizeof(double[TRAITS]));
    Lt = malloc(2*(GS[g]-1)*sizeof(int[3])); 
    Idx = malloc(2*(GS[g]-1)*sizeof(int));
  #endif
#endif  
  // update strategies of individual with probability Mu
  for(i = 0; i < GS[g]; i++){   // through every individual in the group
    if( U01() < Mu)          // if random number is less than Mu, then update strategy
      updateStrategy(g, i, 1, X0_Be);    // use quantal response approach to update strategy (each individual strategy is updated believing other individuals do not change)    
  }   
  free(Tx); free(Xmax);
#if PUNISH
  free(Lookup); free(Lidx);  
  #if FORESIGHT
    free(Strat); free(FS); free(Lt); free(Idx);    
  #endif
#endif  
}

// play CA2 event (us vs them game for ng no. of groups)
void event_CA2(int ng)
{
  int i, j, k, *g, *sg, preexists;
  double wgx, xb;    // sum of groups contribution with Beta exponent
  
  g = malloc(ng*sizeof(int));       // array for group index of selected group in its polity  
  sg = malloc(ng*sizeof(int));      // array for selected unique group indices
  
  // choose ng no. of groups randomly
  for(k = 0; k < ng; k++){
    do{
      g[k] = rnd(G);            // selection from uniform distribution; // generate random no. less than total no. of groups
      preexists = 0;
      for(j = 0; j < k; j++){      // check if it already exists in the selected group array 'sg' for playing game
	  if(g[k] == sg[j]){
	    preexists = 1; break;
	  }
	}
    }while(preexists == 1);

    sg[k] = g[k];
  }
   
  // play us vs them game
  B = normal(Bo, Sigma_B);                // randomly choose Benefit B
  if(B < 0) B = 0;                     // check if B is less than zero
  update_XP_CA2(g, ng);        // caculate probability of success of ng groups in array pbs    
  
  for(k = 0; k < ng; k++){
    updateIndividualPayoffAfterProduction(g[k], ng, P[g[k]]);          // update payoff of individuals in each group and of each group 
#if PUNISH
    punish(g[k]);
#endif
    updateGroupPayoff(g[k]);
  } 
    
  for(wgx = 0, k = 0; k < ng; k++)
    wgx += pow(X[g[k]], Beta);        // sum of group contribution with Beta exponent
  
  for(k = 0; k < ng; k++){    
    Tx    = malloc( GS[g[k]]*sizeof(double [TRAITS]));
    Xmax  = malloc(GS[g[k]]*sizeof(double));
    // get state of strategies before update
    for(i = 0; i < GS[g[k]]; i++){
      Tx[i][0] = x[g[k]][i];
      Tx[i][1] = dxi[g[k]][i];
      Tx[i][2] = dsi[g[k]][i];      
#if Accumulatedpayoff      
  Xmax[i] = pow( ( (1-Discount)*Api[g[k]][i] + (1.0 + B*V[g[k]][i]) )/C, I_Alpha );         // upper bound of x at each role 
#else
  Xmax[i] = pow( (1.0 + B*V[g[k]][i])/C, I_Alpha);                                 
#endif 
    }    
#if PUNISH
    Lookup = malloc(2*(GS[g[k]]-1)*sizeof(int[3]));                                     // allocate memory for lookup table  
    Lidx   = malloc(2*(GS[g[k]]-1)*sizeof(int));
    #if FORESIGHT
      Strat = malloc(GS[g[k]]*sizeof(double[TRAITS]));      
      FS = malloc(GS[g[k]]*sizeof(double[TRAITS]));
      Lt = malloc(2*(GS[g[k]]-1)*sizeof(int[3])); 
      Idx = malloc(2*(GS[g[k]]-1)*sizeof(int));
    #endif
#endif  
    // update strategies of individual with probability Mu
    for(i = 0; i < GS[g[k]]; i++){   // through every individual in group
      if( U01() < Mu){          // if random number is less than Mu, then update strategy
	xb = wgx - pow( X[g[k]], Beta );
	xb = MAX(xb, 0.001);
	updateStrategy(g[k], i, ng, xb);    // use quantal response approach to update strategy (each individual strategy is updated believing other individuals do not change)
      }
    }
    free(Tx); 
#if PUNISH
    free(Lookup); free(Lidx);  
    #if FORESIGHT
      free(Strat); free(FS); free(Lt); free(Idx);    
    #endif
#endif 
  }
  free(g); free(sg); free(Xmax); 
}

// plays event group extinction and group replication
void event_group_ext_rep()
{    
  int j, g, gr;    
  double m, t;
  m = pi_g[0];
  t = -m;
  // set distribution probability    
  for(j = G; j--; ){    
    Gr_dist->p[j] = pi_g[j];                // probabilith of recolonisation of group proportional to +ve of group payoff
    Ge_dist->p[j] = -pi_g[j];             // probability of extinction of group proporttional to -ve of group payoff
    if(m < Gr_dist->p[j]) m = Gr_dist->p[j];
    if(t < Ge_dist->p[j]) t = Ge_dist->p[j];    
  }  
  exp_initdist(Gr_dist, 0);
  exp_initdist(Ge_dist, 0);
  g = drand(Ge_dist);        // select group to extinction   
  gr = drand(Gr_dist);      // select group to repopulate     
  if( g == gr) return;       // if same, then noop
  // make selected group extinct and repopulate by another
  for(j = N; j--;){    
    // copy strategies from repopulating group individuals to extincting individuals    
    x[g][j] = x[gr][j];    // replace strategy            
    dxi[g][j] = dxi[gr][j];  //replace threshold
    dsi[g][j] = dsi[gr][j];  // replace aggressiveness
    pi[g][j] = pi[gr][j];     // replace payoff     
    Api[g][j] = Api[gr][j];
  }    
  pi_g[g] = pi_g[gr];                          // replace group payoff  
}

// returns index of randomly chosen event with probability associated with event
int select_event()
{
  int j;
  double sum;  
  for(sum = 0.0, j = 0; j < EVENTS; j++){
    Event_dist->p[j] = PE[j];
    sum += PE[j];
  }
  initdist(Event_dist, sum);         // initialise distribution event_dist
  
  return drand(Event_dist);          // selects one event randomly from distribution and returns index  
}

// plays event according to inded of event e
void play_event(int ev)
{
  if(ev == 0)
    event_CA1();
  else if(ev == 1)
    event_CA2(NGUT);
  else if(ev == 2)
    event_group_ext_rep();    
} 

// calculates all the stats
void calc_stat(int k, int r) 
{
  double xm, pm, dxim, dsim, pun_im, pun_jm;
  int h, g;    
  double *xsum, *pisum, *dxisum, *dsisum, *pun_isum, *pun_jsum;  
  
#if ALLDATAFILE
  int j;
  static FILE **fp = NULL;  // file pointers for individual run  
  double spi_g = 0;  // sum of payoffs for individual run
  if(k < 0){
    for (h = 7; h--; fclose(fp[h]));
    free(fp); 
    return;
  }     
  
  if(!k){   
    char ixdata[200], ifdata[200], ipi_gdata[200], idxidata[200], idsidata[200], ipun_idata[200], ipun_jdata[200], tstr[100];
    fp = malloc(7*sizeof(FILE *));    
    sprintf(tstr, "e%d.dat", r); prep_file(ixdata, tstr);
    sprintf(tstr, "f%d.dat", r); prep_file(ifdata, tstr);
    sprintf(tstr, "dxi%d.dat", r); prep_file(idxidata, tstr);
    sprintf(tstr, "dsi%d.dat", r); prep_file(idsidata, tstr);
    sprintf(tstr, "pi_g%d.dat", r); prep_file(ipi_gdata, tstr);    
    sprintf(tstr, "pun_i%d.dat", r); prep_file(ipun_idata, tstr);
    sprintf(tstr, "pun_j%d.dat", r); prep_file(ipun_jdata, tstr);
    fp[0] = fopen(ixdata, "w");
    fp[1] = fopen(ifdata, "w");
    fp[2] = fopen(idxidata, "w");
    fp[3] = fopen(idsidata, "w");
    fp[4] = fopen(ipun_idata, "w");
    fp[5] = fopen(ipun_jdata, "w");
    fp[6] = fopen(ipi_gdata, "w");        
  }    
  // write headers
  if( k == 0){
    for( j = 0; j < 6; j++){
      fprintf(fp[j], "%d\t", 0);   
      for(h = 0; h < N; h++)  fprintf(fp[j], "%d\t", h + 1);  
      if(j == 1) fprintf(fp[j], "avg ");     // fp[1] used for pi and avg needs to be shown
      fprintf(fp[j], "\n"); 
    }     
  }
  for( j = 0; j < 7; j++) fprintf(fp[j], "%d\t", (int)(k*SKIP));        // write time units     
#endif 
  // calculate average effort per rank over all groups in each generation
  xsum = calloc( N,sizeof(double));
  pisum = calloc( N,sizeof(double));
  dxisum = calloc( N,sizeof(double));
  dsisum = calloc( N,sizeof(double));
  pun_isum = calloc( N,sizeof(double));
  pun_jsum = calloc( N,sizeof(double));    
  for(g = G; g--; ){                    // through each group
    pun_g[g] = 0;
    for(h = 0; h < GS[g]; h++){             // through each individual               
      xsum[h] += x[g][h];        // sum of efforts for each rank      
      pisum[h] += Api[g][h];  // sum of payoffs for each rank
      dxisum[h] += dxi[g][h];        // sum of thresholds for each rank      
      dsisum[h] += dsi[g][h];  // sum of aggressivenes for each rank
      pun_isum[h] += pun_i[g][h];
      pun_jsum[h] += pun_j[g][h];
      pun_g[g] += (pun_i[g][h]+pun_j[g][h]);   // sum of punishment accrued by individuals in a group      
    }     
  } 
  
  for(h = 0; h < N; h++){
    pm = pisum[h] /(double)(G);    // average pay off for each role
    xm = xsum[h] /(double)(G);  // average effort for each role
    dxim = dxisum[h] /(double)(G);    // average threshodl for each role    
    dsim = dsisum[h] /(double)(G);  // average aggressiveness for each role  
    pun_im = pun_isum[h] / (double)G;
    pun_jm = pun_jsum[h] / (double)G;
    xmean[k][h] += xm;                // average effort per role over all groups in addition for no. of runs
    pimean[k][h] += pm;
    dximean[k][h] += dxim;                // average threshold per role over all groups in addition for no. of runs
    dsimean[k][h] += dsim;
    pun_imean[k][h] += pun_im;
    pun_jmean[k][h] += pun_jm;
    pi_gmean[k] += pm;              // average group payoffs over runs
#if ALLDATAFILE    
    spi_g += pm;             // average group payoffs over each run    
    // write data for individual runs
    fprintf(fp[0], "%.4lf  ", xm);
    fprintf(fp[1], "%.4lf  ", pm);
    fprintf(fp[2], "%.4lf  ", dxim);
    fprintf(fp[3], "%.4lf  ", dsim);
    fprintf(fp[4], "%.4lf  ", pun_im);
    fprintf(fp[5], "%.4lf  ", pun_jm);    
#endif
  }  
#if ALLDATAFILE  
  fprintf(fp[6], "%.4lf  \n", spi_g/N);
  fprintf(fp[0], "\n");
  fprintf(fp[1], "\n");
  fprintf(fp[2], "\n");
  fprintf(fp[3], "\n");
  fprintf(fp[4], "\n");
  fprintf(fp[5], "\n");
#endif
  
  free(xsum);  free(pisum); free(dxisum);  free(dsisum);  free(pun_isum); free(pun_jsum);  
}

void displayMatrix()
{
  int g, i, j;
  for(g = 0; g < G; g++){
    for(i = 0; i < N; i++){
      for(j = 0; j < N; j++){
	Delta_ij_avg[i][j] += Delta_ij[g][i][j];
	Sij_avg[i][j] += Sij[g][i][j];
      }
    }
  }
  printf("\n S matrix\n");
  for(i = 0; i < N; i++){
      for(j = 0; j < N; j++){	
	printf("%.2lf \t\t", Sij_avg[i][j]/(G)) ;
      }  
      printf("\n");
    }   
  printf("\n");
  printf("\n Delta matrix \n");
  for(i = 0; i < N; i++){
      for(j = 0; j < N; j++){	
	printf("%.4lf \t\t", Delta_ij_avg[i][j]/(G)) ;
      }  
      printf("\n");
    }   
  printf("\n");
  printf("\n strength ratio matrix \n");
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){	
      printf("%.2lf \t\t", S[0][j]/S[0][i]) ;
    }  
    printf("\n");
  }   
  printf("\n");   
}

int main(int argc, char **argv)
{
#if DEBUG
  feenableexcept(FE_DIVBYZERO| FE_INVALID|FE_OVERFLOW); // enable exceptions
#endif
  if(argc ^ 2){
    printf("Usage: ./sa sa.config\n");
    exit(1);
  }
  if(read_config(argv[1])){           // read config
    printf("READDATA: Can't process %s \n", argv[1]);
    return 1;
  }  
  T = T+SKIP;               // increase T by SKIP to go from 0 to T all the way in time units
  
  // precomputing  
  I_Gamma = 1.0/Gamma;
  I_Alpha = 1.0/Alpha;    
  Bo = B;    
  X0_Be = pow(X0, Beta);
  GaBe  = Gamma * Beta;

  // event probabilites scaling by G;
  PE[0] = PE[0]*(double)G;
  PE[1] = PE[1]*(double)G;
  PE[2] = PE[2]*(double)G;
  
  initrand(Seed);
  init();          // note: method is in init.c  // allocate memory at initial  
  
  int i, j, ev, k = 0; // m;
  unsigned long seed;
  double dt, ct, tt;      // tt: time of simulation; dt: delta time after which next event occurs; ct: counter time to check with skip time for stat calculation
  time_t now;   
 
  for( i = 0; i < Runs; i++){
    tt = 0.0, ct = 0.0;
    k = 0; //m = 0;
    if(Seed == 0){      
      now = time(0);
      seed = ((unsigned long)now);     
      initrand(seed);
      printf("\n run: %d rand seed: %lu\n",i+1, seed);
    }
    set_init_xrvpi();   // sets initial strategies vector, roles, valuations and payoffs       
     
    calc_stat(0, i);
    
    for( j = 0; j < N; j++) printf("%.4lf  ", V[0][j]); printf("\n");
    
    // Gillespie's stochastic simulation 
    while(tt < T){     
      // calc lambda; lambda = summation of P[j]s
      for(Lambda = 0, j = EVENTS; j--;)
	Lambda += PE[j];
      // select time for next event
      dt = randexp(Lambda); //printf("dt: %.4lf\n", dt);
      // select next event
      ev = select_event();
      // play the event      
      play_event(ev);
      // update time
      tt += dt;
      // calc stat
      ct += dt;
      if(ct >= SKIP){
        calc_stat(++k, i); // calculates all the stats	
	ct = 0.0;
      }
      // update events associated probabilites if necessary
    }
#if ALLDATAFILE
    calc_stat(-1, -1);    // free file pointers for individual run data files
#if !CLUSTER
    plotallIndividualRun(i, 0);      // note: method is in dataplot.c  
#if GRAPHS
    plotallIndividualRun(i, 1);      // note: method is in dataplot.c
#endif
#endif
    // write traits of final state
    {
      int m, n;
      char tstr[200], str[100];
      sprintf(str, "traits%d.dat", i);
      sprintf(str, "traits%d.dat", i); prep_file(tstr, str);
      FILE *fp = fopen(tstr, "w");
      for(m = 0; m < G; m++){
	for( n = 0; n < GS[m]; n++){
	  fprintf(fp, "%.4lf %.4lf    ", dxi[m][n], dsi[m][n]);
	}
	fprintf(fp, "\n");
      }
      fclose(fp);
    }
#if PUNISH
#if DISP_MATRIX
    plotTraits(i, 0);   // note: method is in dataplot.c
#endif
#if GRAPHS
    plotTraits(i, 1);   // note: method is in dataplot.c
#endif
#endif
#endif
  }  
  // write effort and payoff data
  writefile_effort_fertility();              // note: method is in dataplot.c
  writefile_threshold_aggressiveness();      // note: method is in dataplot.c
  // plot data with graphics using gnuplot
  if(!CLUSTER){    
    plotall(0);                              // note: method is in dataplot.c
#if GRAPHS
    plotall(1);                              // note: method is in dataplot.c
#endif
  }
#if PUNISH
#if DISP_MATRIX
  displayMatrix();                           
#endif  
#endif
  
  cleanup();                                 // note: method is in init.c
  return 0;
}

/* 
gcc -Wall -O2 -march=native -pipe -o pun pun.c -lm
./pun sa.config
valgrind -v --track-origins=yes --leak-check=full --show-leak-kinds=all ./pun sa.config
gcc -g -o pun pun.c -lm
gdb pun
run sa.config

//profiling code
gcc -Wall -O2 -march=native -pipe -pg pun.c -o pun -lm
./pun sa.config
 gprof pun gmon.out > analysis.txt  
*/
