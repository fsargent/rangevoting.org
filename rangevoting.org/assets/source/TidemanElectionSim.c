#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

// Election simulator for Tideman-statistics and random-rank-order-vote statistics
// elections, by Warren D. Smith, Sept. 2011.

// Compile:   gcc -Wall -O6 TidemanElectionSim.c

/***
As Jameson Quinn hinted, perhaps Tideman's election model is not very good.
Tideman based his model on a fit to 87 elections having from 3 to 26 candidates.

However, I do not have a better model than Tideman's
to offer.  I'm also unaware of anybody besides Tideman who attempted
to make a reality-based statistical model of this kind (perhaps
Tideman knows of some examples?).   I made one myself once for IEVS, 
but it has no simple description since it is based on large stored tables of 
election data.
In contrast, Tideman's model is specified by just 4 real parameters.

This software for simulating Tideman and random-vote-elections
can handle any reasonably simple election method based on rank-order
ballots which depends only on the pairwise-margins matrix.
Examples: Borda, Schulze-Condorcet, Nanson-Condorcet, Copeland-Condorcet, 
  basic (least reversal) Condorcet, Simpson-Kramer Condorcet.
Examples of voting methods it cannot handle: Instant runoff,
  plurality, approval, range voting, highest median score.
***/

#define uint unsigned int
#define uint64 unsigned long long
#define real double
#define uchar unsigned char
#define UC (uchar)
#define MAXUINT32 0xFFFFFFFF
#define MAXUINT64 ((((uint64)(MAXUINT32))<<32) | ((uint64)MAXUINT32))
#define TWOPOWER64 18446744073709551616.0
#define TWOPI 6.2831853071795864769252867665590057683943387987502116419498891846

#define MAXCANDDTS 5002

#define TRIES 10000
#define RANDSEED 18140670317

int sign(real x){
  if(x>0.0) return(1);
  if(x<0.0) return(-1);
  return(0);
}

uint64 MLt[256]; 
uchar MLc;
//Simple, fast, ultra-hi-period extremely random uniform uint random generator 
//by George Marsaglia (he called it "LFIB4");
//the MLt[0..255] array must be seeded with random bits initially, 
//but in such a way that so that not all MLt[j] are even numbers.

real RandUniform(){ // returns random uniform real on interval (0,1)
  uint64 x;
  MLc++;
  x = MLt[MLc]+MLt[UC(MLc+58)]+MLt[UC(MLc+119)]+MLt[UC(MLc+178)];
  MLt[MLc] = x;
  return( (x+0.5) / TWOPOWER64 );
}

void InitRandom( uint seed ){
  uint64 j;
  uint i;
  srand48(seed);
  for(j=1; j; j += j){ // fill randgen array with "random" bits
    for(i=0; i<256; i++){ MLt[i] += (drand48() > 0.5)? j:0; }
  }
  // make sure not all odd by ORing 32 of them with 1
  for(i=0; i<32; i++){  j = drand48()*256; MLt[j] |= 1; }
  MLc = 1;
  // "warm up" randgen by calling it several thousand times
  for(i=0; i<3000; i++){ RandUniform(); }
}

real RandExponential(){
  return  (-log(RandUniform()));
}

// Generates a random normal-distributed real variable 
// using Box-Muller polar method:
// [if want mean=mu, variance=sigma^2 then replace output x by mu+x*sigma]:
real RandNormal(){
  real theta, g2rad, z;
  static real y;
  static uint mode = 0;
  if(mode!=0){ mode=0; return(y); }
  theta = RandUniform() * TWOPI;
  g2rad = sqrt( -2.0 * log(1.0 - RandUniform()) );
  z = cos(theta) * g2rad;
  y = sin(theta) * g2rad;
  mode = 99;
  return(z);
}

// ProbabilityDensity(x) = x^(k-1) * exp(-x) / Gamma(k).
// Conditions on the parameter:  k>0.
// Returned values range between 0 and +infinity.
// Hisashi Tanizaki: A simple gamma random number generator for arbitrary shape parameters,
// Economics Bulletin  3,7 (2008) 1-10.
// Another method which looks good (but I have not programmed it) is
// G.Marsaglia & W.W.Tsang: Monty Python method, ACM TOMS 24,3 (1998) 341-350
//   http://www.cparity.com/projects/AcmClassification/samples/292453.pdf
real RandGamma(real k){
  real n, b1, b2, c1, c2, v1, v2, w1, w2, x, y;
  assert(k>0.0);
  // do stuff depending on k only:
  if(k <= 0.4){
    n = 1.0/k;
    b1 = 0.0;  b2 = k+k;  c1 = 0.0;
  }else{
    if(k<=4.0){
      n = ((k - 0.4)/3.6 + 1.0)/k;
    }else{ // k>4.0
      n = 1.0/sqrt(k);
    }
    b1 = k - 1.0/n;  b2 = k + 1.0/n;
    c1 = b1 * 0.5 * (log(b1)-1.0);
  }
  c2 =  b2 * 0.5 * (log(b2)-1.0);
  // now generate the gamma(k) random variate x via rejection ploys:
  for(;;){
    do{
      v1 = RandUniform();
      v2 = RandUniform();
      w1 = c1 + log(v1);
      w2 = c2 + log(v2);
      y = n*(b1*w2 - b2*w1);
    }while(y<0.0);
    x = n*(w2-w1);
    x = exp(x);
    if(y >= x) return(x);
  }
}

// Generates a random beta-distributed variable.
// Conditions on the parameters are alpha > 0 and beta > 0.
// Returned values range between 0 and 1.
// [Knuth Vol 2 Edition 3 page 134 "the beta distribution"]
real RandBeta(real alpha, real beta){
  real x,y;
  assert(alpha>0.0);
  assert(beta>0.0);
  y = RandGamma(alpha);
  if(y==0.0) return( 0.0 );
  x = RandGamma(beta);
  return( y / (y+x) );
}

// NOTE the margin matrix will be correct up to a scaling factor.
// Simple method invented by Warren D. Smith to compute exact distribution
void GenerateRandElectionPairwiseMargins(uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS]){
  uint a,b;
  real x,y;
  for(a=0; a<NumCanddts; a++){
    MarginMat[a][a] = RandNormal();
  }
  for(a=0; a<NumCanddts; a++){
    for(b=a+1; b<NumCanddts; b++){
      x = RandNormal();
      y = MarginMat[a][a] - MarginMat[b][b] - x;
      MarginMat[a][b] =  y;
      MarginMat[b][a] = -y;
    }
  }
  for(a=0; a<NumCanddts; a++){
      MarginMat[a][a] = 0.0; 
  }
}

// Tideman model described here:
//   http://rangevoting.org/TidemanElModel.html
void GenerateTidemanPairwiseMargins( uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS] ){
  const real  a1 = -0.532,   a2 = -0.789,   a3 = -2.486,   a4 = -1.281;
  uint i,j;
  real x,y,r,eta,alpha,beta,u;
  for(j=1; j<NumCanddts; j++){
    MarginMat[j][j] = 0.0; 
    for(i=0; i<j; i++){
      x = (j-i)/(NumCanddts+1.0);
      y = 1.0 - (j+i)/(NumCanddts+1.0);
      r = 0.5*exp((a1 + a2*x + a3*y*y)*x);
      eta = pow( (1-2*r)*r*0.5, a4);
      alpha = r*eta;
      beta = eta-alpha;
      u = RandBeta(alpha, beta) - 0.5;
      MarginMat[i][j] =  u;
      MarginMat[j][i] = -u;
    }
  } 
}

// Returns index of Borda winner
uint CopelandAndBordaScores( uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS], 
			     real Borda[MAXCANDDTS], uint Copeland[MAXCANDDTS],
			     real CopeBord[MAXCANDDTS] ){ 
  uint i,j;
  real mx;
  for(j=0; j<NumCanddts; j++){
    Borda[j] = 0;
    Copeland[j] = NumCanddts;
    for(i=0; i<NumCanddts; i++){
      Borda[j] += MarginMat[j][i];  // can be either sign, max for winner
      Copeland[j] += sign(MarginMat[j][i]);  // 2 points for a pairwise win, 0 for loss, 1 for tie
    }
  }

  j=0;
  mx = Borda[j];
  for(i=1; i<NumCanddts; i++){
    if(Borda[i] > mx){ j = i; mx = Borda[i]; }
  }
  assert( Borda[j] >= Borda[0] );
  for(i=1; i<NumCanddts; i++){
    CopeBord[i] = Copeland[i] * 999999.9 + Borda[i];
  }
  return(j);
}

const int Ciura[] = {0, 1, 4, 10, 23, 57, 132, 301, 701, 1750};

// Uses shell sort to produce CopeOrder[k] which is the index {0..NumCanddts-1} of the candidate
// with the kth greatest Copeland score, k=0..NumCanddts-1:
// outputs number of Copeland co-winners.
uint SortByCopeland( uint NumCanddts, uint Copeland[MAXCANDDTS], 
		     real CopeBord[MAXCANDDTS], uint CopeOrder[MAXCANDDTS] ){
  uint j,i,tempkey,temp,inc,iidx;
  for(j=0; j<NumCanddts; j++){
    CopeOrder[j] = j;
  }

  iidx = 9;
  inc = Ciura[iidx];
  while(inc > 0){
    for(i=inc; i<NumCanddts; i++){
      temp = CopeOrder[i];
      j = i;
      tempkey = Copeland[temp];
      while( j>=inc && Copeland[CopeOrder[j-inc]] < tempkey ){
	CopeOrder[j] = CopeOrder[j-inc];
	j -= inc;
      }
      CopeOrder[j] = temp;
    }
    iidx--;
    inc = Ciura[iidx];
  }
  assert( CopeOrder[0] < NumCanddts );
  assert( CopeOrder[1] < NumCanddts );
  assert( CopeOrder[0] != CopeOrder[1] );
  assert( Copeland[CopeOrder[0]] >= Copeland[CopeOrder[1]] ); 

  i = 1;
  tempkey = Copeland[CopeOrder[0]]; // Copeland score of winner
  for(j=1; j<NumCanddts; j++){
    if( Copeland[CopeOrder[j]] == tempkey ) i++; else break;
  }
  return(i);
}

// Arrange the candidates by the Copeland order breaking ties with Borda order.
// Then the "Tideman hardness" is the max over
// negative components of the upper triangle of (j - i + 1)
// where i and j are row and column respectively.
uint TidemanHardnessHack( uint NumCanddts, uint SmithSize,
			  real MarginMat[MAXCANDDTS][MAXCANDDTS], 
		     real CopeBord[MAXCANDDTS], uint CBOrder[MAXCANDDTS] ){
  uint j,i,temp,inc,iidx,mx;
  real tempkey;
  for(j=0; j<NumCanddts; j++){
    CBOrder[j] = j;
  }

  iidx = 9;
  inc = Ciura[iidx];
  while(inc > 0){
    for(i=inc; i<NumCanddts; i++){
      temp = CBOrder[i];
      j = i;
      tempkey = CopeBord[temp];
      while( j>=inc && CopeBord[CBOrder[j-inc]] < tempkey ){
	CBOrder[j] = CBOrder[j-inc];
	j -= inc;
      }
      CBOrder[j] = temp;
    }
    iidx--;
    inc = Ciura[iidx];
  }
  assert( CBOrder[0] < NumCanddts );
  assert( CBOrder[1] < NumCanddts );
  assert( CBOrder[0] != CBOrder[1] );
  assert( CopeBord[CBOrder[0]] >= CopeBord[CBOrder[1]] ); 

  mx = 0;
  for(j=1; j<SmithSize; j++){
    for(i=0; i<j; i++){
      tempkey = MarginMat[CBOrder[j]][CBOrder[i]];
      if(tempkey>0){
	temp = j-i+1;
	if(mx < temp) mx = temp;
      }
    }
  }
  assert(mx>=0);
  return(mx);
}

uint SmithSetSize(uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS], uint CopeOrder[MAXCANDDTS]){
  int j,i,k;
  k=0;
  for(j=0; j<NumCanddts; j++){
    for(i=NumCanddts-1; i>j; i--){
      if(MarginMat[CopeOrder[j]][CopeOrder[i]] < 0.0) break;
    }
    if(i>k) k=i;
    if(j>=k) return(j+1); // Smith set found as prefix of Copeland ordering
  }
  return(NumCanddts);
}

// Finds index of "Basic Condorcet" winner -- reversing the least number of
// pairwise prefs on ballots will cause him to become a beats-all winner.
int BasicCondorcet(uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS] ){
  int k, j, i;
  real sum, mn;
  mn = HUGE;
  k = -1;
  for(j=0; j<NumCanddts; j++){
    sum = 0;
    for(i=0; i<NumCanddts; i++){
      if(MarginMat[j][i] < 0) sum += MarginMat[j][i];
    }
    sum = -sum;
    assert(sum>=0.0);
    if(mn > sum){ mn=sum; k=j; }
  }
  assert(k>=0);
  return(k);
}

// Finds index of "Simpson Kramer Condorcet" winner (MinMax) -- 
// the candidate with the weakest worst-pairwise-defeat, wins.
int SimpsonKramerCondorcet(uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS] ){
  int k, j, i;
  real mn, mx;
  mn = HUGE;
  k = -1;
  for(j=0; j<NumCanddts; j++){
    mx = 0;
    for(i=0; i<NumCanddts; i++){
      if(MarginMat[j][i] < mx) mx = MarginMat[j][i];
    }
    mx = -mx;
    assert(mx >= 0.0);
    if(mn > mx){ mn=mx; k=j; }
  }
  assert(k>=0);
  return(k);
}

// Finds index of "Nanson Baldwin Condorcet" winner -- 
// repeatedly eliminate the Borda loser.
int NansonBaldwinCondorcet(uint NumCanddts, real MarginMat[MAXCANDDTS][MAXCANDDTS] ){
  int k, j, i;
  real sum, mn;
  real Bord[MAXCANDDTS];
  uint Live[MAXCANDDTS];
  mn = HUGE;
  for(j=0; j<NumCanddts; j++){ // initial Borda scores, all live:
    sum = 0;
    for(i=0; i<NumCanddts; i++){
      sum += MarginMat[j][i];
    }
    Bord[j] = sum;
    Live[j] = 99;
  }
  for(j=1; j<NumCanddts; j++){ // eliminations:
      mn = HUGE; k = -1;
      for(i=0; i<NumCanddts; i++) if(Live[i]){
	  if( mn > Bord[i] ){ mn = Bord[i]; k = i; }
      }
      assert(k>=0);
      Live[k] = 0;  //eliminate k
      for(i=0; i<NumCanddts; i++){ // update Borda scores:
	Bord[i] -= MarginMat[i][k];
      }
  }
  for(i=0; i<NumCanddts; i++) if(Live[i]){ // survivor wins:
      return(i);
  }
  assert(0);
}

void PrintHisto(uint n, uint Hist[] ){ // prints histogram
  uint i;
  for(i=0; i<=n; i++){
    if(Hist[i]>0){
      printf(" %d:%d", i, Hist[i]);
    }
  }
  printf("\n");
  fflush(stdout);
}

// Speed: version of 25 sept 2011 runs 1 million TRIES up to 200 in about 13 hours
// Version of 28 sept runs 10000 tries up to 2000 in 8 hours
main(){
  uint ElSize[] = {2,3,4,5,6, 7,8,9,10,11, 12,13,14,15,20, 30,50,100,150,200,
  	   300,500,1000,2000,5000, 999999999};
  uint j, i, nc, RSm, TSm, TCW, RCW, RBo, TBo, Tcondwin, Rcondwin, RTH, TTH;
  uint Rbcd, Tbcd, TcoNonUq, RcoNonUq, TSmSum, RSmSum;
  uint Rlrc, Rskc, Rnbc, Rblk,     Tlrc, Tskc, Tnbc, Tblk;
  uint Rbld, Rbsd, Rbnd, Rlsd, Rlnd, Rsnd;
  uint Tbld, Tbsd, Tbnd, Tlsd, Tlnd, Tsnd;
  static real RandElnMar[MAXCANDDTS][MAXCANDDTS];
  static real TidemanMar[MAXCANDDTS][MAXCANDDTS];
  real RCopeBord[MAXCANDDTS], TCopeBord[MAXCANDDTS];
  real RBorda[MAXCANDDTS], TBorda[MAXCANDDTS];
  uint RCBOrder[MAXCANDDTS], TCBOrder[MAXCANDDTS];
  uint RCopeland[MAXCANDDTS], RCopeOrder[MAXCANDDTS];
  uint TCopeland[MAXCANDDTS], TCopeOrder[MAXCANDDTS];
  uint RCopeHist[MAXCANDDTS+1], TCopeHist[MAXCANDDTS+1];
  uint RSmitHist[MAXCANDDTS+1], TSmitHist[MAXCANDDTS+1];
  uint RhardHist[MAXCANDDTS+1], ThardHist[MAXCANDDTS+1];

  printf("TRIES=%d  RANDSEED=%u\n", TRIES, RANDSEED);

  InitRandom(RANDSEED);
  for(j=0; j<4; j++){
    printf("\nrandom uniforms: ");
    for(i=0; i<5; i++){ printf("%f ", RandUniform()); }
  }
  printf("\n\nrandom normals: ");
  for(i=0; i<5; i++){ printf("%f ", RandNormal()); }
  printf("\nrandom exponentials: ");
  for(i=0; i<5; i++){ printf("%f ", RandExponential()); }
  printf("\nrandom gamma(randexp): ");
  for(i=0; i<5; i++){ printf("%f ", RandGamma(RandExponential())); }
  printf("\nrandom beta(randexp, randexp): ");
  for(i=0; i<5; i++){ printf("%f ", RandBeta(RandExponential(), RandExponential())); }
  printf("\n"); 
  fflush(stdout); 

  for(i=0; i<30; i++){
    nc = ElSize[i];
    if(nc >= MAXCANDDTS) goto ALLDONE;
    Tcondwin = Rcondwin = Rbcd = Tbcd = TcoNonUq = RcoNonUq = TSmSum = RSmSum = 0;
    Rbld = Rbsd = Rbnd = Rlsd = Rlnd = Rsnd = 0;
    Tbld = Tbsd = Tbnd = Tlsd = Tlnd = Tsnd = 0;
    for(j=0; j<=MAXCANDDTS; j++){
      RCopeHist[j] = TCopeHist[j] = RSmitHist[j] = TSmitHist[j] = 0;
      ThardHist[j] = RhardHist[j] = 0;
    }
    for(j=0; j<TRIES; j++){
      GenerateRandElectionPairwiseMargins(nc, RandElnMar);
      GenerateTidemanPairwiseMargins(nc, TidemanMar);
      Rblk=RBo=CopelandAndBordaScores(nc, RandElnMar, RBorda, RCopeland, RCopeBord);
      RCW = SortByCopeland(nc, RCopeland, RCopeBord, RCopeOrder);
      RCopeHist[RCW]++;
      RSm = SmithSetSize(nc, RandElnMar, RCopeOrder );
      RTH = TidemanHardnessHack(nc, RSm, RandElnMar, RCopeBord, RCBOrder);
      RhardHist[RTH]++;
      RSmitHist[RSm]++;
      RSmSum += RSm;
      if(RCW==1){ // unique Copeland winner exists
	if( RCopeland[RCopeOrder[0]]==nc+nc-1 ){ // condorcet winner exists
	  Rblk = RCopeOrder[0]; // Black winner
	  Rcondwin++;
	  if( RBo != RCopeOrder[0] ) Rbcd++;  // ...and disagrees with Borda winner
	}
      }else{ // non-unique Copeland winner
	RcoNonUq++;
      }
      Rlrc = BasicCondorcet(nc, RandElnMar);
      Rskc = SimpsonKramerCondorcet(nc, RandElnMar);
      Rnbc = NansonBaldwinCondorcet(nc, RandElnMar);
      if( Rblk != Rlrc ) Rbld++;
      if( Rblk != Rskc ) Rbsd++;
      if( Rblk != Rnbc ) Rbnd++;
      if( Rlrc != Rskc ) Rlsd++;
      if( Rlrc != Rnbc ) Rlnd++;
      if( Rskc != Rnbc ) Rsnd++;
      Tblk=TBo=CopelandAndBordaScores(nc, TidemanMar, TBorda, TCopeland, TCopeBord);
      TCW = SortByCopeland(nc, TCopeland, TCopeBord, TCopeOrder);
      TCopeHist[TCW]++;
      TSm =   SmithSetSize(nc, TidemanMar, TCopeOrder);
      TTH = TidemanHardnessHack(nc, TSm, TidemanMar, TCopeBord, TCBOrder);
      TSmitHist[TSm]++;
      ThardHist[TTH]++;
      TSmSum += TSm;
      if(TCW==1){
	if( TCopeland[TCopeOrder[0]]==nc+nc-1 ){ // Condorcet winner exists
	  Tblk = TCopeOrder[0]; // Black winner
	  Tcondwin++;
	  if( TBo != TCopeOrder[0] ) Tbcd++;  // ...and disagrees with Borda winner
	}
      }else{ // non-unique Copeland winner
	TcoNonUq++;
      }
      Tlrc = BasicCondorcet(nc, TidemanMar);
      Tskc = SimpsonKramerCondorcet(nc, TidemanMar);
      Tnbc = NansonBaldwinCondorcet(nc, TidemanMar);
      // disagreements among different Condorcet cycle-breakers:
      if( Tblk != Tlrc ) Tbld++; 
      if( Tblk != Tskc ) Tbsd++;
      if( Tblk != Tnbc ) Tbnd++;
      if( Tlrc != Tskc ) Tlsd++;
      if( Tlrc != Tnbc ) Tlnd++;
      if( Tskc != Tnbc ) Tsnd++;
    }
    printf("\n%d candidates:\n", nc);
    printf("#TiCWs=%d ", Tcondwin);
    printf("#TiCopNonUq=%d ", TcoNonUq);
    printf("#TiBorConDis=%d ", Tbcd);
    printf("#TiSmiAvg=%.2f\n", (TSmSum+0.000001)/TRIES);
    printf("Tbld=%d, Tbsd=%d, Tbnd=%d, Tlsd=%d, Tlnd=%d, Tsnd=%d\n",
	   Tbld, Tbsd, Tbnd, Tlsd, Tlnd, Tsnd );

    printf("#ReCWs=%d ", Rcondwin);
    printf("#ReCopNonUq=%d ", RcoNonUq);
    printf("#ReBorConDis=%d ", Rbcd);
    printf("#ReSmiAvg=%.2f ", (RSmSum+0.000001)/TRIES);
    printf("\n");
    printf("Rbld=%d, Rbsd=%d, Rbnd=%d, Rlsd=%d, Rlnd=%d, Rsnd=%d\n",
	   Rbld, Rbsd, Rbnd, Rlsd, Rlnd, Rsnd );
    printf("TiSmitHist(%d)=",nc);    PrintHisto(nc,TSmitHist);
    printf("ReSmitHist(%d)=",nc);    PrintHisto(nc,RSmitHist);
    printf("TiHardHist(%d)=",nc);    PrintHisto(nc,ThardHist);
    printf("ReHardHist(%d)=",nc);    PrintHisto(nc,RhardHist);
    printf("TiCopeHist(%d)=",nc);    PrintHisto(nc,TCopeHist);
    printf("ReCopeHist(%d)=",nc);    PrintHisto(nc,RCopeHist);
    fflush(stdout);
  }
 ALLDONE:
  printf("all done.\n");
  fflush(stdout);
}
