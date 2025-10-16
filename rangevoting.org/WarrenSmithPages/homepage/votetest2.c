/*******************
Voting system testing program  votetest.c
  version 2 - incorporates Bucklin system and optional voter ignorance
  version 3 - incorporates Meek voting system(s?) and Carey's IFPP3 system
  (c)  Warren D. Smith  NEC Research Institute NJ USA  October 2000
on LINUX systems, compile with
  gcc votetest.c -lm -lc -O9 -o votetest
on other operating systems, or other C compilers than gcc, may
require some fiddling to get it to compile.
******************************************************
Compares various voting systems
experimentally via monte-carlo... Here is how...
1. generate (for V voters, c candidates)
random in [0,1] utilities for each voter-candidate pair.
(Other utility probability distributions also could be considered.
See below.)

2. since we know the true voter utilities we know the true best
candidate and his utility (summed over all voters).

3. run various voting schemes to elect the candidates various ways.

4. compute the "regret" for that voting system, which is the utility
difference between the true-best candidate and whoever is elected.

Presumably the best voting system is the one with minimal expected
regret.  

With 100,000 simulated elections with [0,1] uniform utilities, the
regret values should be accurate (90% conf error bds) to better than
1%, bootstrapping experiments show.

--------
Following written critique by Prof. Steven J. Brams (NYU dept. of 
Politics) recvd March 2001, I now have modified this program to allow 
"ignorant voters." Brams was worried that range voting might not
work well with ignorant voters and claimed (!?!) that half the
electorate does not even know who the vice-president is.
So, in the simulation, each voter has a true utility for
each election winner, just as before, and these true utility
values are used to assess the Bayesian regrets just as before.
BUT now, each voter does not KNOW his own candidate-election
utilities, instead what he knows is, those utilities, POLLUTED
by the addition of IGNORANCE, i.e., added noise. Specifically
we add a Gaussian random deviate with mean 0 and std. deviation Q.
(If Q=0 this reduces to the old version of the program.)
We may now re-run all the simulations with some Q>0, and see what 
happens.

--------

All this is easy to do if have "honest" voters.  

"Honest approval" voting is somewhat hard to define since there
are c-1 approval votes consistent with any c-candidate ordering. Which
one do we pick? I will, somewhat arbitrarily, assume the honest voters
approve of above-average candidates, in ApHa.

Another kind of honest approval voter chooses his threshhold to
maximize his expected utility assuming all other voters are
completely random uniform, namely he maximizes the
sum of utility differences between his top and bottom 2 groups.
This is ApHu. Theorem (confirmed computationally): This 
actually is the same as ApHa!!

To do the same experiment with "rational" voters is tougher.
We will construct various strategies which seem to approximate
rationality.

PlS (strategic plurality) is: you vote for the best of the 2 frontrunners.

BuS (strategic bullet) is: you vote against the worst of the 2 frontrunners.

RaS (strategic range) is, you use as your threshhold, a utility midway
between the 2 frontrunners and you give +1 to candidates above
threshhold, -1 to candidates below. This seems to be the strategy
which maximizes expected utility in 3-candidate elections, based on
the reasoning that the 2 frontrunners have probability enormously
close to 1 of being elected so you have to give them +1 and -1 to
maximize your expected impact; then in the enormously unlikely event
the 3rd candidate can contend, it is best to then act as though it
were a 3-way tie in the polls so equally likely for you to break a tie
between 3rd candidate and 1st, as between 3rd and 2nd, and given that
you've already committed to +1 and -1 for the 2 frontrunners, the best
choice of +1 and -1 for the 3rd candidate is the one maximizing the
expected utility difference, i.e. use the avg of U1 and U2 as your
threshhold.

For elections with c>=4 candidates, though, I think the RaS strategy
should be modified.  Namely: consider the candidates in decreasing
order of their election chances based on the pre-election polls. Vote
+1,-1 or -1,+1 for the first 2 candidates.  For candidates
i=3,4,5,... select +1 or -1 as your vote, whichever maximizes
   sum                    sum               (Uj - Uk)
  j=cands with +1    k=cands with -1
  j<=i               k<=i
This is the same thing as: pick your vote[i] to be
  sign[       sum   (Ui - Uk)  -      sum    (Uj - Ui)  ]
         k=cands with -1           j=cands with +1
         k<i                       j<i
This in turn is the same thing as: pick your vote[i] to be
  sign[ Ui - average(utility among previous candidates) ].
This acts, when selecting the vote for candidate i, as though it were
an i-way tie in the polls among the i frontrunners, and you want to
maximize expected utility in that situation.  When i=3 this is the
same strategy as the preceding paragraph but it may differ for larger
i. Call this RaS2.

My experiments suggest that in practice, the difference between RaS
and RaS2 is extremely small.

BoS (strategic Borda) is, you give the best of the two frontrunners
the max vote and the worst the min vote. Then you proceed recursively
on the remaining C-2 candidates.

LRS (strategic condorcet) you give the best frontrunner top-vote, other
frontrunner bottom vote, then be honest on other candidates.
Probably this is not as damaging to society as the true
best-strategy condorcet vote, so get a lower bound on regret this way.

STVS: Strategic STV: same strategy as LRS.

Strategic Copeland: same strategy as LRS.

We assume the "frontrunners" are candidates 0,1,2,... in decreasing order
of apparent pre-election likelihood of winning.

Random Dictator and Random Pair Majority are immune to strategy - for
them the best strategy is honesty.

Coombs STV (most least-liked candid is eliminated each round)

Brief results:
If voters are honest, range voting is the best by far.  If strategic
voters, then range voting is also the best.
The regret ratios can be as large as 30:1 for honest range vs
{honest or strategic} bullet... typically range is ahead
of competing systems by factor of 2... strategic range is usually
a little better than honest plurality.

*******************************
Different utility distributions:
Let us assume there are I "issues". The candidates & voters choose
their stances on each of these issues by picking a random number in
[-1,1].  Then the utility of a candidate for a voter is the dot product of
the candidate and voter issue-stance vectors, plus a random noise
term (random uniform in [-1,1]). Note, the all-noise method
above is equivalent to just the special case I=0.

In all cases the utilities are post-affined to make them lie in [0,1],
a normalization that attempts to allow cross-distributional comparisons.
*********Warren D. Smith 2000********************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include<limits.h>

#define uint unsigned
#define uchar unsigned char
#define uint64 unsigned long long
#define bool uint
#define FALSE (1==0)
#define TRUE (!(FALSE))

#define MaxVoters 200
#define MaxCandids 10
#define MaxSystems 50
int NumExperiments = 2000;
#define VERBOSE 0
#define CHEESYCHECK 0
#define HONESTYSTATS 1
#define CONDHARESTATS 0

#define VERBOSECAREY 1

#define RANDOM_UTILITIES 0
#define NORMAL_UTILITIES 5
#define MaxIssues 10

 double sqrt(double x);
 double log(double x); /*base e*/

/*voting systems:*/
#define RaH 0
#define BoH 1
#define LRH 2
#define CoombsH 3
#define STVH 4
#define Hcope 5
#define Hdaba 6
#define BlackH 7
#define Hbuck 8
#define PlRunH 9
#define PlH 10
#define BuH 11
#define RPM 12
#define RandomDict 13
#define RandomWin 14
#define WorstWin 15

#define ApHa 16
/*#define ApHu 11*/

#define RaS 17
#define RaS2 18
#define PlS 19
#define BoS 20
#define BuS 21
#define LRS 22
#define STVS 23
#define BoS2 24
#define CoombsS 25
#define BoS3 26
#define DaS 27
#define CopeS 28
#define BlackS 29

#define HIFPP3 30
#define SIFPP3 31
#define HMeek123 32
#define SMeek123 33
#define Htide 34
#define Stide 35

uint64 RaS2MisOrder[MaxCandids];
uint64 RaS2RaSDiffer[MaxCandids];
uint64 FunnyStatCt[MaxCandids];
int path[MaxCandids][MaxCandids];

double fabs (double x);

/**** UNIX man page description of drand48: 
     ...sequence of 48-bit integer values,
     X[i], according to the linear congruential formula

              X[n+1] = (a X[n] + c)mod m        n>0.

                      48
     The parameter m=2  ; hence 48-bit integer arithmetic is performed.
     Unless lcong48 has been invoked, the multiplier value a and the addend
     value c are given by
           a = 5DEECE66D_16           c = B_16
Therefore Period=2^48, but low-period least signif bits.
*********************************************************/
double drand48(void);
void srand48(long int seedval);

/*********************************************
long int random(void);
void srandom(unsigned int seed);
char *initstate(unsigned int seed, char *state, int n);

UNIX man page falsely claims random() nonlinear. Actually, examining
the source code shows it is a linear subtractive 
lagged fibonacci generator, has very large period if
you go for max-size 256.
***********************************************/
char randstate[256];

int RandSeed = 3456781;
int Voters = 10;
double IGNORANCEQ = 1.0;
double CandIgnorance[MaxCandids];     /* candidate-dependent ignorance values */
double utils[MaxVoters][MaxCandids];  /* utils[i][j] true utility of cand j for voter i */
double igutl[MaxVoters][MaxCandids];  /* utils[][] but with ignorance added */
double voterstance[MaxVoters][MaxIssues];
double candstance[MaxCandids][MaxIssues];
uint hrank[MaxVoters][MaxCandids];
uint drank[MaxVoters][MaxCandids];
double hplurregret;
double hplurRUNregret;
double hbulletregret;
double hbordaregret;
double htidemanregret;
double stidemanregret;
double hIFPP3regret;
double sIFPP3regret;
double hMeek123regret;
double sMeek123regret;
double hblackregret;
double hrangeregret;
double hSTVregret; 
double hLRregret;
double hAPPAregret;
double hAPPUregret;
double hRdictregret;
double hrandregret;
double hworstregret;
double hRPMregret;
double hBucklinregret;
double hCOPEregret;
double hDABAregret;

double srangeregret;
double splurregret;
double sblackregret;
double sbordaregret;
double sbulletregret;
double srangeregret2;
double sLRregret;
double sSTVregret; 
double sBoS2regret; 
double sBoS3regret; 
double hCoombsregret; 
double sCoombsregret; 
double sDABAregret;
double sCOPEregret;

double regrets[MaxCandids][MaxSystems];

uint sclrw[MaxCandids], shstvw[MaxCandids]; 
/* used for collecting stats on strategic condorcet least rev, 
 * strategic hare stv winners */


/*************************************
I.D.Hill, B.A.Wichmann, D.R.Woodall:
Algorithm 123, single transferable vote by Meek's method,
Computer Journal 30,3 (1987) 277-281
available online at
       http://www.bcs.org.uk/election/meek/meekm.htm

Sometimes very high precision floating point #s are required, 
s later demonstrated by Wichmann. For instance, with 69 candidates 
and <1000 votes, one can produce an example requiring
127 decimal places! But in practice high precision
is usually not a problem. I am using double precision and ignoring
the issue.
  Brian L. Meek died in 1997.
His voting method was adopted by the Electoral Reform Society,
the Royal Statistical Society, and the London Mathematical Society.

This implementation has not been tested for NUMSEATS>1.

This may or may not be the same as:
 Brian L. Meek:
A transferable voting system including intensity of preference,
Math\'ematiques et Sciences Humaines
13, 50 (Summer 1975) 23-29.
 This Meek paper was reprinted in the
ERS publication ``Voting Matters,''
Issue 6, May 1996.  
*************************************/

#define NUMSEATS 1 /* more than 1 if want more-than-1-winner elections */

int MeekAlg123( rankings, randperm, C )
uint rankings[MaxVoters][MaxCandids];
uint randperm[MaxCandids]; /* used for tie-breaking */
uint C; /* number of candidates */
/* rankings[3][5] is the 4th voter's 6th ranked candidate (note
 *  counts start from 0). It is left unaltered by this routine. */
{
#define EXCLUDED 0
#define HOPEFUL 1
#define ELECTED 2
#define MaxMeekIters 333
   int i,j,h;
   uint leastliked,iterct,numelected,numexcluded;
   bool Converged,ElectedSomebody;
   double minvots,adjustfactor,totalvotes,quota,excess,w1,w2;
   double CandWeight[MaxCandids];
   double CandV[MaxCandids];
   uchar CandStatus[MaxCandids];

   assert(1<=NUMSEATS);
   assert(NUMSEATS<=C);

   /* initial set up: */
   for(i=0; i<C; i++){ 
      CandStatus[i] = HOPEFUL; 
      CandWeight[i] = 1.0; 
   }
   numelected=0;
   numexcluded=0;

   do{

      /* Meek's nonlinear iteration: */
      iterct=0;
      do{
         iterct++;
	 for(j=0; j<C; j++){
	    CandV[j] = 0.0; 
            if( CandStatus[j] == HOPEFUL ){ CandWeight[j] = 1.0; }
	    else if( CandStatus[j] == EXCLUDED ){ CandWeight[j] = 0.0; }
	 }
	 /*** A vote (pref ordering)  A>B>C  is processed as follows:
          *  A receives from that voter   v[A] += w[A]
          *  B receives from that voter   v[B] += (1-w[A])*w[B]
          *  C receives from that voter   v[C] += (1-w[A])*(1-w[B])*w[C]
          *  the remaining votes go to    excess += (1-w[A])*(1-w[B])*(1-w[C])
          * Then we let  totalvotes = v[0]+v[1]+...+v[C-1]; 
	  ****************************************************/
	 excess = 0.0;
	 for(i=0; i<Voters; i++){
            h = rankings[i][0];
	    w1 = CandWeight[h];
	    CandV[h] += w1;
	    w1 = 1.0-w1;
	    for(j=1; j<C; j++){
	       h = rankings[i][j];
	       w2 = CandWeight[h];
	       CandV[h] += w2*w1;
	       w1 *= (1.0-w2);
	    }
	    excess += w1;
	 }
	 totalvotes = 0.0;
	 for(i=0; i<C; i++){ totalvotes += CandV[i]; }
	 quota = (totalvotes - excess)/(NUMSEATS+1.0);
	 if(quota<0.0001) quota = 0.0001;
	 if(numelected==0) break;

	 /* readjust elected candidate's weights so their votes are equal
          * to the quota: */
	 Converged=TRUE;
	 for(j=0; j<C; j++){  
	    if( CandStatus[j] == ELECTED ){
	       adjustfactor = quota/CandV[j];
	       if(adjustfactor>1.0000001 || adjustfactor<0.9999999){
                  Converged=FALSE;
	       }
	       CandWeight[j] *= adjustfactor;
	       if(CandWeight[j]>1.0) CandWeight[j]=1.0;
	    }
	 }
      }while( iterct<MaxMeekIters && !Converged );

      if(!Converged){
	 printf("Meek123 convergence failure MaxMeekITers=%u C=%u V=%u\n",
            MaxMeekIters,C,Voters);
      }

      ElectedSomebody = FALSE;
      for(j=0; j<C; j++){  
	 if( CandStatus[j] == HOPEFUL && CandV[j] >= quota ){
	    CandStatus[j] = ELECTED;
	    /* CandWeight[j] = quota/CandV[j]; ???*/
	    ElectedSomebody = TRUE;
	    numelected++;
	 }
      }

      leastliked=C;
      if( ElectedSomebody==FALSE ){
         minvots = HUGE;
	 for(j=C-1; j>=0; j--){
            h = randperm[j];
	    if( CandV[h] < minvots ){ 
	       minvots = CandV[h]; 
	       leastliked = h; 
	    }
	 }
	 assert(leastliked<C);
	 CandStatus[leastliked] = EXCLUDED;
	 numexcluded++;
      }

      if(numexcluded + NUMSEATS >= C){
         for(i=0; i<C; i++){
            if( CandStatus[i] == HOPEFUL ){
	       CandStatus[i] = ELECTED;
               numelected++;
            }
         }         
      }

   }while( numelected < NUMSEATS );

   for(i=C-1; i>=0; i--){
      if(numelected<=NUMSEATS) break;
      h = randperm[i];
      if( CandStatus[h] == ELECTED ){
         numelected--;
         CandStatus[h] == EXCLUDED;
      }
   }

   for(i=0; i<C; i++){
      h = randperm[i];
      if( CandStatus[h] == ELECTED ){
         /*printf("MeekAlg123 C=%u randperm=%u %u  winner=%d\n", 
	  *  C, randperm[0],randperm[C-1], h ); */
         assert(0<=h);
         assert(h<C);
         return(h); /* return index 0..C-1 of winning candidate */
      }
   }

   assert( 0==1 ); /* should never reach this point */
#undef EXCLUDED
#undef HOPEFUL
#undef ELECTED
#undef MaxMeekIters
#undef NUMSEATS
}


int ParkMillerX = 34561; /* any nonzero value will do */
/* Iterates X <-- X*48271 mod 2147483647.   X=1..2147483646.
 * Here 48271 is primitive and 2^31-1=2147483647 is prime: */
int ParkMiller(){
    ParkMillerX = 48271*(ParkMillerX%44488) - 3399*(ParkMillerX/44488);
    if(ParkMillerX<0) ParkMillerX = ParkMillerX+2147483647;
    return ParkMillerX;
}

uint64  CoveyouX = 14542; /* anything which is 2 mod 4 will do */
/* iterates X <-- X*(X+1) mod 2^w, with period 2^(w-2), w=wordsize=32.
* The least signif bits are a lot less random 
* [2^(<=k) bits have period max(2^(k-1),1)]. */
uint Coveyou(){
   assert((CoveyouX&3)==2);
   CoveyouX *=   CoveyouX+1;  /* mod 2^64, implied */
   return ((uint)(CoveyouX>>32));
}
     

double myrand(){ /* combines 2 UNIX randgens, and my 2 randgens above,
                  * since don't trust any of them */
   double x,z;
   z = drand48();
   assert(0.<=z);
   assert(z<=1.);
   x = random() / (1. + RAND_MAX);
   assert(0.<=x);
   assert(x<=1.);
   z -= x;
   if(z<0.0) z += 1.0;
   assert(0.<=z);
   assert(z<=1.);
   x = (Coveyou()^((uint)ParkMiller())) * (1/4294967295.7);
   assert(0.<=x);
   assert(x<=1.);
   z -= x;
   if(z<0.0) z += 1.0;
   assert(0.<=z);
   assert(z<=1.);
   return z;
}

/* uses myrand() to find a unit-variance normal deviate
 * Polar method Knuth SNA 3.4.1C */
double mynormalrand(){
   double v1,s;
   static double v2;
   static int normct=0;
   if(normct){ normct=0; return(v2); }
   normct=1;
   do{
     v1 = myrand();
     v2 = myrand();
     v1 = 2.*v1-1.;
     v2 = 2.*v2-1.;
     s = v1*v1+v2*v2;
   }while( s >= 1.0 );
   s = sqrt(-2.*log(s)/s);
   v1 *= s;
   v2 *= s;
   return(v1);
}


/*************************************
G.A.Craig Carey's "IFPP" method
   Carey email: research@ijs.co.nz 
("Improved First Past the Post") for conducting
3-candidate elections, with the 3 candidates named A,B,C.
You can cast exactly one of 9 possible types of votes:

  The 9 kinds        Name of variable that
of allowed votes     counts that kind of vote
----------------     -------------------------
    A                   a0
   A>B                  ab
   A>C                  ac
    B                   b0
   B>C                  bc
   B>A                  ba
    C                   c0
   C>A                  ca
   C>B                  cb
(It seems that Carey intends "A" to mean "A>B and A>C" and "B>C" to mean
"B>C>A"; the candidates un-named in your preference ordering are regarded
as worse, in your view, than the named ones.) Let
  a = a0+ab+ac;   b = b0+ba+bc;   c = c0+ca+cb;
The below C code by WDS is based on Carey's RELOG code.
IFPP is only available for 1,2, and 3 candidate elections,
with either 1 or 2 winners allowed in the 3-candidate case.
******************************************
IFPP3(a,ab,ac, b,bc,ba, c,ca,cb)
int a,ab,ac, b,bc,ba, c,ca,cb;
{
   bool Wa,Wb,Wc, La,Lb,Lc;

   Wa = (b+c<2*a) && ((b+cb<a+ca) || (2*b<a+c)) && ((c+bc<a+ba) || (2*c<a+b));
   Wb = (c+a<2*b) && ((c+ac<b+ab) || (2*c<b+a)) && ((a+ca<b+cb) || (2*a<b+c));
   Wc = (a+b<2*c) && ((a+ba<c+bc) || (2*a<c+b)) && ((b+ab<c+ac) || (2*b<c+a));
   * Wc=TRUE means C is a winner if 1-winner election. *

   La = (b+c>2*a) && ((b+cb>a+ca) || (2*b>a+c)) && ((c+bc>a+ba) || (2*c>a+b));
   Lb = (c+a>2*b) && ((c+ac>b+ab) || (2*c>b+a)) && ((a+ca>b+cb) || (2*a>b+c));
   Lc = (a+b>2*c) && ((a+ba>c+bc) || (2*a>c+b)) && ((b+ab>c+ac) || (2*b>c+a));
    * Lc=TRUE means A and B are both winners if 2-winner election; one could
    * then set    Wa = Lc  ||  Lb;   Wb = Lc  ||  La;   Wc = Lb  ||  La;
    * if desired. **
   return stuff;
}
***********************************
If all voters express full preference orderings 
(i.e. only 6 of the 9 kinds of votes are now allowed)
then IFPP simplifies to the following:
*******************************************/

uint SimpleIFPP3(ab,ac, bc,ba, ca,cb)
uint ab,ac, bc,ba, ca,cb;
{
   uint a,b,c;
   uint Wa,Wb,Wc; /* actually bool */ /*,La,Lb,Lc;*/
   a=ab+ac; b=bc+ba; c=ca+cb;

   Wa = (b+c<2*a) && ((b+cb<a+ca) || (2*b<a+c)) && ((c+bc<a+ba) || (2*c<a+b));
   Wb = (c+a<2*b) && ((c+ac<b+ab) || (2*c<b+a)) && ((a+ca<b+cb) || (2*a<b+c));
   Wc = (a+b<2*c) && ((a+ba<c+bc) || (2*a<c+b)) && ((b+ab<c+ac) || (2*b<c+a));
   /* Wc=TRUE means C is a winner if 1-winner election. */

   /********
   La = (b+c>2*a) && ((b+cb>a+ca) || (2*b>a+c)) && ((c+bc>a+ba) || (2*c>a+b));
   Lb = (c+a>2*b) && ((c+ac>b+ab) || (2*c>b+a)) && ((a+ca>b+cb) || (2*a>b+c));
   Lc = (a+b>2*c) && ((a+ba>c+bc) || (2*a>c+b)) && ((b+ab>c+ac) || (2*b>c+a));
    * Lc=TRUE means A and B are both winners if 2-winner election; one could
    * then set    Wa = Lc  ||  Lb;   Wb = Lc  ||  La;   Wc = Lb  ||  La;
    * if desired. **************/

   if(Wa) return(0);
   if(Wb) return(1);
   if(!Wc){ 
     /*no winner!  
      *this pathology apparently cannot happen if prime # of voters */
#if VERBOSECAREY
     printf("ab=%d ac=%d bc=%d ba=%d ca=%d cb=%d\n", ab,ac, bc,ba, ca,cb);
#endif
     /* not sure below is what Carey really wanted: */
     return( (uint)(myrand()*3.0) );
   }
   assert(Wc); return(2);
}

uint MySimpleIFPP3(ab,ac, bc,ba, ca,cb)
double ab,ac, bc,ba, ca,cb;
{
   uint a,b,c;
   uint Wa,Wb,Wc; /* actually bool */ /*,La,Lb,Lc;*/
   a=ab+ac; b=bc+ba; c=ca+cb;

   Wa = (b+c<2*a) && ((b+cb<a+ca) || (2*b<a+c)) && ((c+bc<a+ba) || (2*c<a+b));
   Wb = (c+a<2*b) && ((c+ac<b+ab) || (2*c<b+a)) && ((a+ca<b+cb) || (2*a<b+c));
   Wc = (a+b<2*c) && ((a+ba<c+bc) || (2*a<c+b)) && ((b+ab<c+ac) || (2*b<c+a));
   /* Wc=TRUE means C is a winner if 1-winner election. */

   if(Wa) return(0);
   if(Wb) return(1);
   assert(Wc); return(2);
}




void genutils(int kind, int CandDepIgnorance, int Cand){
  int i,j,k;
  double dotprod;
  assert(RANDOM_UTILITIES==0);
  if( CandDepIgnorance ){
      for(j=0; j<Cand; j++){
	CandIgnorance[j] = myrand() * IGNORANCEQ;
      }
  }else{
      for(j=0; j<Cand; j++){
	CandIgnorance[j] = IGNORANCEQ;
      }      
  }
  if(kind == RANDOM_UTILITIES){
    for(i=0; i<Voters; i++){
      for(j=0; j<Cand; j++){
	utils[i][j] = myrand();
	igutl[i][j] = utils[i][j] + mynormalrand()*CandIgnorance[j];
      }
    }
  }else if(kind==NORMAL_UTILITIES){
    for(i=0; i<Voters; i++){
      for(j=0; j<Cand; j++){
	utils[i][j] = mynormalrand();
	igutl[i][j] = utils[i][j] + mynormalrand()*CandIgnorance[j];
      }
    }
  }else if(kind > 0){ /* ISSUE_BASED_UTILS; kind=#issues */
    assert(kind<=MaxIssues);
    for(i=0; i<Voters; i++) for(j=0; j<kind; j++) voterstance[i][j] = myrand()*2.-1.;
    for(i=0; i<Cand; i++) for(j=0; j<kind; j++) candstance[i][j] = myrand()*2.-1.;
    for(i=0; i<Voters; i++){
      for(j=0; j<Cand; j++){
	dotprod = kind; /* Offset to prevent negative utilities */
        for(k=0; k<kind; k++) dotprod += voterstance[i][k]*candstance[j][k];
        /* Add random noise term to the dotprod, same scale as 1 issue: */
        dotprod += myrand()*2;
        /* This rescale: causes all utilities to lie in [0,1]: */
	dotprod /= (2+2*kind);
        utils[i][j] = dotprod;
	igutl[i][j] = utils[i][j] + mynormalrand()*CandIgnorance[j];
        assert(0.<=dotprod);
        assert(dotprod<=1.);
      }
    }
  }
}

void dovotes(int C){
  double sumutils[MaxCandids];
  double hrange[MaxCandids];
  double extra[MaxCandids];
  uint elimcand[MaxCandids];
  uint elimd[MaxVoters];
  int elimind[MaxVoters];
  int pairels[MaxCandids][MaxCandids];
  int sumvictories[MaxCandids];
  int summargin[MaxCandids];
  double minutil, maxutil, worstutil, bestutil, utildiff, reciputildiff;
  double uthresh, bestudiff, umsum, udiff, bestD, avg, udiffsum, myutil;
  double topmin, botmax;
  int th, thw, rnd;
  int argmax, argmin, sclrwinner, shstvwinner;
  int Scondwinner, Sbordwinner, Hcondwinner, Hbordwinner;
  int i,j,k,t,c1,c2,maxBvote,minBvote, minDvote, maxDvote;
  int loser, winner, secwinr, best, secbest;
  int round, winner1, winner2, dictator;
  int hplurvotes[MaxCandids], hbullets[MaxCandids];
  int hborda[MaxCandids], hdabagh[MaxCandids];
  int myvotes[MaxCandids], happrov[MaxCandids];
  int splurvotes[MaxCandids], sbullets[MaxCandids];
  int sbordavotes[MaxCandids], srangevotes[MaxCandids];
  int sdabaghvotes[MaxCandids];
  uint randperm[MaxCandids];
  uint tmp;

  assert(1<C);
  assert(C<MaxCandids);

  /* find random perm, used for tiebreaking later */
  for(i=0; i<C; i++) randperm[i]=i;
  for(i=0; i<C; i++){
     j = myrand() * C;
     assert(0<=j);
     assert(j<C);
     tmp = randperm[i];
     randperm[i] = randperm[j];
     randperm[j]=tmp;
  }

  /* find summed utils then find best and worst candidates */
  for(i=0; i<C; i++) sumutils[i]=0.0;
  for(i=0; i<Voters; i++){     for(j=0; j<C; j++){
    sumutils[j] += utils[i][j];
  }}
  bestutil = -99.;
  for(i=0; i<C; i++) if(bestutil < sumutils[i]) bestutil = sumutils[i];
  worstutil = Voters*99.;
  for(i=0; i<C; i++) if(worstutil > sumutils[i]) worstutil = sumutils[i];

  /* honest ranks hrank[0]=best, hrank[1]=2ndbest ... hrank[C-1]=worst: */
  for(i=0; i<Voters; i++){
    for(j=0; j<C; j++){ hrank[i][j] = j; }
    /* cheesy bubble sort: */
    for(j=0; j<C; j++){    for(k=j+1; k<C; k++){
      if( igutl[i][hrank[i][j]] < igutl[i][hrank[i][k]] ){
	 t = hrank[i][j];
	 hrank[i][j] = hrank[i][k];
	 hrank[i][k] = t;
      }
    }}
    /* test sorted in decreasing util order: */
    for(j=1; j<C; j++){ assert( igutl[i][hrank[i][j-1]] >= igutl[i][hrank[i][j]] ); }
  }
  /* dishonest ranks drank[0]=best, drank[1]=2ndbest etc: just like honest ranks except
   * best frontrunner is top ranked and worst is bottom ranked. */
  for(j=0; j<C; j++) extra[j] = 0.;
  for(i=0; i<Voters; i++){
    for(j=0; j<C; j++){ drank[i][j] = j; }
    if( igutl[i][0] < igutl[i][1] ){ extra[0] = -99.; extra[1] =  99.; }
    else{                            extra[0] =  99.; extra[1] = -99.; }
    /* cheesy bubble sort: */
    for(j=0; j<C; j++){    for(k=j+1; k<C; k++){
      if( igutl[i][drank[i][j]]+extra[drank[i][j]] < igutl[i][drank[i][k]]+extra[drank[i][k]] ){
	 t = drank[i][j];
	 drank[i][j] = drank[i][k];
	 drank[i][k] = t;
      }
    }}
  }

  /* honest plurality: */
  for(i=0; i<C; i++) hplurvotes[i] = 0;
  for(i=0; i<Voters; i++){
    maxutil = -99.;
    argmax = -1;
    for(j=0; j<C; j++){
      if(maxutil < igutl[i][j]){ maxutil = igutl[i][j]; argmax = j; }
    }
    assert( argmax >= 0 );
    assert( argmax == hrank[i][0] );
    hplurvotes[argmax] ++;
  }
  best = -99;
  winner = -1;
  for(i=0; i<C; i++) if(best < hplurvotes[randperm[i]]){ best = hplurvotes[randperm[i]]; winner = randperm[i]; }
  assert(winner >= 0);
  hplurregret = bestutil - sumutils[winner];

  /* honest IFPP3 (use honest plurality if c!=3): */
  if(C!=3){
    hIFPP3regret = bestutil - sumutils[winner];
  }

  /* honest plurality with runoff: */
  if(2*hplurvotes[winner] < Voters){ /*only do runoff if no >=50% winner*/
    best = -99;
    secbest = -99;
    for(i=0; i<C; i++){
      if(best < hplurvotes[randperm[i]]){ 
	secbest = best;
	best = hplurvotes[randperm[i]]; 
	secwinr = winner;
	winner = randperm[i]; 
      }else if(secbest < hplurvotes[randperm[i]]){
	best = hplurvotes[randperm[i]]; 
	secwinr = randperm[i];
      }
    }
    assert(winner >= 0);
    assert(secwinr >= 0);
    /* now do runoff between winner, secwinr: */
    c1=c2=0;
    for(i=0; i<Voters; i++){
      if( igutl[i][winner] >= igutl[i][secwinr] ) c1++; else c2++;    
    }
    if(c2>c1) winner = secwinr;
  }
  hplurRUNregret = bestutil - sumutils[winner];

  /* honest IFPP3 (if c==3) */
  if(C==3){
    double ab,ac,bc,ba,ca,cb;
    ab=ac=cb=ca=bc=ba=0.0;
    for(i=0; i<Voters; i++){ 
      if( hrank[i][0]==0 ){
	if( hrank[i][1]==1 ){ ab+=10.0; }
	else if( hrank[i][1]==2 ){ ac+=10.0; }
      }else if( hrank[i][0]==1 ){
	if( hrank[i][1]==0 ){ ba+=10.0; }
	else if( hrank[i][1]==2 ){ bc+=10.0; }
      }else if( hrank[i][0]==2 ){
	if( hrank[i][1]==0 ){ ca+=10.0; }
	else if( hrank[i][1]==1 ){ cb+=10.0; }
      }
    }

    ab -=  myrand();
    ac -=  myrand();
    bc -=  myrand();
    ba -=  myrand();
    ca -=  myrand();
    cb -=  myrand();
    
    winner = MySimpleIFPP3(ab,ac, bc,ba, ca,cb);
    hIFPP3regret = bestutil - sumutils[winner];
  }

  /* honest Meek123: */
  winner = MeekAlg123( hrank, randperm, C );
  hMeek123regret = bestutil - sumutils[winner];

  /* honest bullet: */
  for(i=0; i<C; i++) hbullets[i] = 0;
  for(i=0; i<Voters; i++){
    minutil = 99.;
    argmin = -1;
    for(j=0; j<C; j++){
      if(minutil > igutl[i][j]){ minutil = igutl[i][j]; argmin = j; }
    }
    assert( argmin >= 0 );
    assert( fabs( igutl[i][argmin] - igutl[i][ hrank[i][C-1] ] ) < 0.000001 );
    hbullets[argmin] ++;
  }
  best = Voters*99;
  winner = -1;
  for(i=0; i<C; i++) if(best > hbullets[randperm[i]]){ best = hbullets[randperm[i]]; winner = randperm[i]; }
  assert(winner >= 0);
  hbulletregret = bestutil - sumutils[winner];


  /* honest range (i.e. scaled utility voting): */
  for(i=0; i<C; i++) hrange[i] = 0.0;
  for(i=0; i<Voters; i++){
    maxutil = igutl[i][ hrank[i][0  ] ];
    minutil = igutl[i][ hrank[i][C-1] ];
    utildiff = maxutil-minutil;
    assert(utildiff>0.0);
    reciputildiff = 1.0/utildiff;
    for(j=0; j<C; j++){
      hrange[j] += (igutl[i][j] - minutil) * reciputildiff;
  }}  
  bestD = -99.;
  for(i=0; i<C; i++) if(bestD < hrange[randperm[i]]){ bestD = hrange[randperm[i]]; winner = randperm[i]; }
  hrangeregret = bestutil - sumutils[winner];

  /* honest Hare STV (eliminate candid with fewest best-place votes): */
  for(i=0; i<C; i++) elimcand[i] = 0; /*will change to 1 when eliminate*/
  for(i=0; i<Voters; i++) elimd[i] = 0;
  for(round=1; round<C; round++){
    for(i=0; i<C; i++) hplurvotes[i] = 0;
    for(i=0; i<Voters; i++){
      argmax = hrank[i][elimd[i]];
      hplurvotes[argmax] ++;
    }
    loser = -1;
    best = 99*Voters;
    for(i=C-1; i>=0; i--) if(elimcand[randperm[i]]==0){
       if(best > hplurvotes[randperm[i]]){ best = hplurvotes[randperm[i]]; loser = randperm[i]; }
    }
    assert(loser>=0);
    assert(loser<C);
    assert(elimcand[loser]==0);
    elimcand[loser]=1; /*eliminate loser*/
    for(i=0; i<Voters; i++){ 
       while( elimcand[ hrank[i][elimd[i]] ] != 0 ){ 
	 elimd[i]++; 
	 assert(elimd[i]<C);
       }
    }
  }
  winner = -1;
  for(i=0; i<C; i++) if(elimcand[i]==0){ assert(winner<0); winner=i; }
  assert(winner>=0);
  assert(winner<=C);
  hSTVregret = bestutil - sumutils[winner];

  /* strategic Hare STV - vote one of the 2 frontrunners last place, other top, rest honest. */
  for(i=0; i<C; i++) elimcand[i] = 0; /*will change to 1 when eliminate*/
  for(i=0; i<Voters; i++) elimd[i] = 0;
  for(round=1; round<C; round++){
    for(i=0; i<C; i++) hplurvotes[i] = 0;
    for(i=0; i<Voters; i++){
      argmax = drank[i][elimd[i]];
      hplurvotes[argmax] ++;
    }
    loser = -1;
    best = 99*Voters;
    for(i=C-1; i>=0; i--) if(elimcand[randperm[i]]==0){
       if(best > hplurvotes[randperm[i]]){ best = hplurvotes[randperm[i]]; loser = randperm[i]; }
    }
    assert(loser>=0);
    assert(loser<C);
    assert(elimcand[loser]==0);
    elimcand[loser]=1; /*eliminate loser*/
    for(i=0; i<Voters; i++){ 
       while( elimcand[ drank[i][elimd[i]] ] != 0 ){ 
	 elimd[i]++; 
	 assert(elimd[i]<C);
       }
    }
  }
  winner = -1;
  for(i=0; i<C; i++) if(elimcand[i]==0){ assert(winner<0); winner=i; }
  assert(winner>=0);
  assert(winner<=C);
  sSTVregret = bestutil - sumutils[winner];
  shstvwinner = winner;
  shstvw[shstvwinner]++;


  /* honest Coombs STV (eliminate candid with most worst-place votes): */
  for(i=0; i<C; i++) elimcand[i] = 0; /*will change to 1 when eliminate*/
  for(i=0; i<Voters; i++){ elimind[i] = C-1; }
  for(rnd=1; rnd<C; rnd++){
    for(j=0; j<C; j++){ hplurvotes[j] = 0; }
    for(i=0; i<Voters; i++){
      for(j=0; j<C; j++){ 
         argmax = hrank[i][elimind[i]];
         hplurvotes[argmax] ++;
      }
    }
    loser = -1; best = -1;
    for(j=C-1; j>=0; j--){
       if(best<hplurvotes[randperm[j]]){ best = hplurvotes[randperm[j]]; loser=randperm[j]; }
    }
    assert(loser>=0);
    assert(elimcand[loser]==0);
    elimcand[loser]=1; /*eliminate loser*/
    for(i=0; i<Voters; i++){
       while( elimcand[ hrank[i][elimind[i]] ] != 0 ){ 
	 elimind[i]--; /*eliminate*/
	 assert(elimind[i] >= 0);
       }
       assert( elimcand[hrank[i][elimind[i]]] == 0 );
    }
  }
  winner = -1;
  for(i=0; i<C; i++) if(elimcand[i]==0){ assert(winner<0); winner=i; }
  assert(winner>=0);
  assert(winner<=C);
  hCoombsregret = bestutil - sumutils[winner];


  /* Strategic Coombs STV (eliminate candid with most worst-place votes): */
  for(i=0; i<C; i++) elimcand[i] = 0; /*will change to 1 when eliminate*/
  for(i=0; i<Voters; i++){ elimind[i] = C-1; }
  for(rnd=1; rnd<C; rnd++){
    for(j=0; j<C; j++){ hplurvotes[j] = 0; }
    for(i=0; i<Voters; i++){
      for(j=0; j<C; j++){ 
         argmax = drank[i][elimind[i]];
         hplurvotes[argmax] ++;
      }
    }
    loser = -1; best = -1;
    for(j=C-1; j>=0; j--){
       if(best<hplurvotes[randperm[j]]){ best = hplurvotes[randperm[j]]; loser=randperm[j]; }
    }
    assert(loser>=0);
    assert(elimcand[loser]==0);
    elimcand[loser]=1; /*eliminate loser*/
    for(i=0; i<Voters; i++){
       while( elimcand[ drank[i][elimind[i]] ] != 0 ){ 
	 elimind[i]--; /*eliminate*/
	 assert(elimind[i] >= 0);
       }
       assert( elimcand[drank[i][elimind[i]]] == 0 );
    }
  }
  winner = -1;
  for(i=0; i<C; i++) if(elimcand[i]==0){ assert(winner<0); winner=i; }
  assert(winner>=0);
  assert(winner<=C);
  sCoombsregret = bestutil - sumutils[winner];


  /* strategic IFPP3 (if c==3) */
  if(C==3){
    double ab,ac,bc,ba,ca,cb;
    ab=ac=cb=ca=bc=ba=0.0;
    for(i=0; i<Voters; i++){ 
      if( drank[i][0]==0 ){
	if( drank[i][1]==1 ){ ab+=10.0; }
	else if( drank[i][1]==2 ){ ac+=10.0; }
      }else if( drank[i][0]==1 ){
	if( drank[i][1]==0 ){ ba+=10.0; }
	else if( drank[i][1]==2 ){ bc+=10.0; }
      }else if( drank[i][0]==2 ){
	if( drank[i][1]==0 ){ ca+=10.0; }
	else if( drank[i][1]==1 ){ cb+=10.0; }
      }
    }

    ab -=  myrand();
    ac -=  myrand();
    bc -=  myrand();
    ba -=  myrand();
    ca -=  myrand();
    cb -=  myrand();

    winner = MySimpleIFPP3(ab,ac, bc,ba, ca,cb);
    sIFPP3regret = bestutil - sumutils[winner];
  }

  /* honest Condorcet Least-Reversal: */
  for(i=0; i<C; i++) for(j=0; j<C; j++) pairels[i][j] = 0;
  for(i=0; i<Voters; i++){
    for(j=0; j<C; j++){
      for(k=0; k<C; k++){
        pairels[j][k] += ((igutl[i][j] > igutl[i][k]) ? 1 : -1);
  }}}
  for(i=0; i<C; i++) pairels[i][i] = 0;
  /* pairels[j][k] = margin of victory of j over k */
  for(j=0; j<C; j++){ summargin[j] = 0; }
  for(j=0; j<C; j++){
    for(k=0; k<C; k++){
      summargin[j] += ((pairels[j][k] < 0) ? -pairels[j][k] : 0);
  }}
  best = 99*Voters*C;
  for(i=0; i<C; i++) if(best > summargin[randperm[i]]){ best = summargin[randperm[i]]; winner = randperm[i]; }
  hLRregret = bestutil - sumutils[winner];
  Hcondwinner = winner;
  for( k=0; k<C; k++){ /* Is there a Condorcet winner? */
    if( pairels[winner][k] < 0 ){ Hcondwinner = -1; break; }
  }

  /* Honest Copeland. Ties count 1 point, victory 2, defeat 0 */
  for(j=0; j<C; j++){ sumvictories[j] = 0; }
  for(j=0; j<C; j++){
    for(k=0; k<C; k++){
      if(pairels[j][k] > 0) sumvictories[j] += 2;
      if(pairels[j][k] == 0) sumvictories[j] ++;
  }}
  best = -1;
  for(i=0; i<C; i++) if(best < sumvictories[randperm[i]]){ best = sumvictories[randperm[i]]; winner = randperm[i]; }
  hCOPEregret = bestutil - sumutils[winner];


  /* Tideman ranked pairs with Honest voters.
   * Tideman = Condorcet variant in which you
   * pick the A>B comparison with the largest margin and "lock it in".
   * Then you  pick the next largest one available
   * ("available" means: not already used and not creating a cycle), 
   * and continue on.  This creates an ordering of the candidates. The
   * topmost in the ordering wins.  It is a bit tricky to spot
   * the cycles as we go...
   * (This code based on email from Blake Cretney bcretney@postmark.net):
   ******************************************/

  /*path[][] is used as a changeable copy of pairels.*/
  for(i=0;i<C;i++){
    for(j=0;j<C;j++) path[i][j]=pairels[i][j];
    path[i][i]=INT_MAX;
  }
  /* Whenever a victory
   * is locked in, the appropriate cell is set to INT_MAX.
   * pi,pj are used with randperm to give precedence to victories higher
   * in the random permutation (when tie-breaking).
   * This loop finds the next pair (i,j) to lock in: ********/
  for(;;){
    int maxp,oi,oj,pi,pj;
    maxp=INT_MIN;
    for(pi=0;pi<C;pi++) for(pj=pi+1;pj<C;pj++){
      oi=randperm[pi]; oj=randperm[pj];
      if(   path[oi][oj]!=INT_MAX
	 && path[oj][oi]!=INT_MAX){/*not locked-out*/
	if(path[oi][oj]>maxp){ maxp=path[oi][oj]; i=oi; j=oj; }
	if(path[oj][oi]>maxp){ maxp=path[oj][oi]; i=oj; j=oi; }
      }
    }
    if(maxp==INT_MIN) break;
    /********* print the pair and its margin:
    printf("(%d %d) %d\n",i,j,maxp);
    ***********************/
    /*lock in the pair and clobber future no-good pairs:*/
    for(oi=0;oi<C;oi++) for(oj=0;oj<C;oj++){
      if(path[oi][i]==INT_MAX && path[j][oj]==INT_MAX){
	path[oi][oj]=INT_MAX;
      }
    }
  }
  /* The above code assumes that pairels has been set properly.
   * path[][] ends up with the winning row having all cells
   * set to INT_MAX.  In fact, a complete ranking is given,
   * where path[i][j]==INT_MAX means that i is
   * ranked over j (where i!=j).   So to find the winner: ****/
  winner = -99;
  for(i=0;i<C;i++){
    for(j=0;j<C;j++){
      if(path[i][j] != INT_MAX) break;
    }
    if(j>=C){ winner=i; break; }
  }
  assert(winner >= 0);
  htidemanregret = bestutil - sumutils[winner];  
  
  /* Remarks:
   * B.Cretney: The Ranked Pair algorithm above takes O(C^4) time.
   * An O(C^3) time algorithm is possible, but it is more complex.
   * WD.Smith: An O(C^4) algorithm is perfectly acceptable if C<20,
   * but anyway... for people interested in the true asymptotic
   * complexity for large C...  I suspect an O(C^{2.001}) algorithm is
   * possible, or more precisely O(C^2 * logterms). This would be
   * best possible (except for the logterms).
   * It would work as follows. Regard the already-locked-in ranked pairs
   * as a DAG (directed acyclic graph; the "acyclic" means there
   * are no logical inconsistencies in the ranking). We maintain
   * a spanning forest which sort of "over-represents" the "skeleton"
   * of this DAG. (By "skeleton" I mean the key comparisons from which
   * all the others merely are consequences by transitivity. By
   *"over-represented" I mean: For a DAG   C<--B2<--A-->B1-->C,
   * the  node C is represented twice in the spanning tree shown.)
   * We now go thru the candidate edges in increasing cost order.
   8 (Can use a "heap," also called "priority queue,"  data structure
   * to allow doing that in O(logC) time per candidate edge.)
   * Now at any moment we want to find a new lowest-cost directed-edge
   * to adjoin to the DAG, which does not cause a directed cycle.
   * With the aid of fast "lowest common ancestor" data structures
   * associated with each tree in our forest (and associated
   * algorithm tricks), we may quickly determine for any candidate
   * directed-edge whether it would cause a cycle in the DAG.
   * (Namely: If the candidate edge joins two forest trees, it cannot
   * cause a cycle. If it joins nodes from a single tree, then
   * it causes a cycle if and only if the lowest common ancestor of the
   * two nodes is, in fact, one of the nodes themselves, and
   * the edge is directed toward it.) If it would, we mark that edge
   * as "dead and never to be looked at again ever", i.e. we delete
   * it from the "heap". If it would not cyclize, we adjoin it
   * to the DAG and update all our tree and LCA data structures
   * and lock it in. The total number of such data structure
   * updates is O(C).
   *
   * The above sketch is still a sketch with missing details and
   * it could easily be there is a big problem with it. I.e. it
   * is not to be regarded as a proof such a runtime complexity
   * is achievable, but it is pointing the way to try to
   * devise such a proof. -WDS.
   *****************************************************/
  

  /* honest Borda: */
  for(i=0; i<C; i++) hborda[i] = 0;
  for(i=0; i<Voters; i++){     for(j=0; j<C; j++){
    hborda[ hrank[i][j] ] += C-j-1;
  }}  
  best = -99;
  for(i=0; i<C; i++) if(best < hborda[randperm[i]]){ best = hborda[randperm[i]]; winner = randperm[i]; }
  hbordaregret = bestutil - sumutils[winner];
  Hbordwinner = winner;

  /* Honest Black: if no Condorcet winner then use Borda: */
  if(Hcondwinner >= 0) winner = Hcondwinner;
  else winner = Hbordwinner;
  hblackregret = bestutil - sumutils[winner];

  /* honest Dabagh: */
  for(i=0; i<C; i++) hdabagh[i] = 0;
  for(i=0; i<Voters; i++){    
    hdabagh[ hrank[i][0] ] += 2;
    hdabagh[ hrank[i][1] ] += 1;
  }
  best = -99;
  for(i=0; i<C; i++) if(best < hdabagh[randperm[i]]){ best = hdabagh[randperm[i]]; winner = randperm[i]; }
  hDABAregret = bestutil - sumutils[winner];



  /* strategic Condorcet Least-Reversal: */
  for(i=0; i<C; i++) for(j=0; j<C; j++) pairels[i][j] = 0;
  for(i=0; i<C; i++) extra[i] = 0.;
  for(i=0; i<Voters; i++){
    if( igutl[i][0] < igutl[i][1] ){ extra[0] = -99.; extra[1] =  99.; }
    else{                            extra[0] =  99.; extra[1] = -99.; }
    for(j=0; j<C; j++){
      for(k=0; k<C; k++){
        pairels[j][k] += ((igutl[i][j]+extra[j] > igutl[i][k]+extra[k]) ? 1 : -1);
  }}}
  for(i=0; i<C; i++) pairels[i][i] = 0;
  /* pairels[j][k] = margin of victory of j over k */
  for(j=0; j<C; j++){ summargin[j] = 0; }
  for(j=0; j<C; j++){
    for(k=0; k<C; k++){
      summargin[j] += ((pairels[j][k] < 0) ? -pairels[j][k] : 0);
  }}
  best = 99*Voters*C;
  for(i=0; i<C; i++) if(best > summargin[randperm[i]]){ best = summargin[randperm[i]]; winner = randperm[i]; }
  sLRregret = bestutil - sumutils[winner];
  sclrwinner = winner;
  sclrw[sclrwinner]++;
#if CONDHARESTATS
  if( sclrwinner>1 || shstvwinner>1 || sclrwinner != shstvwinner ){
    printf("sclrwinner=%d != shstvwinner=%d >1\n", sclrwinner, shstvwinner);
  }
#endif
  Scondwinner = winner;
  for( k=0; k<C; k++){ /* Is there a Condorcet winner? */
    if( pairels[winner][k] < 0 ){ Scondwinner = -1; break; }
  }


  /* Strategic Copeland */
  for(j=0; j<C; j++){ sumvictories[j] = 0; }
  for(j=0; j<C; j++){
    for(k=0; k<C; k++){
      if(pairels[j][k] > 0) sumvictories[j] += 2;
      if(pairels[j][k] == 0) sumvictories[j] ++;
  }}
  best = -99;
  for(i=0; i<C; i++) if(best < sumvictories[randperm[i]]){ best = sumvictories[randperm[i]]; winner = randperm[i]; }
  sCOPEregret = bestutil - sumutils[winner];


    /* Tideman ranked pairs with Strategic voters.
   * Tideman = Condorcet variant in which you
   * pick the A>B comparison with the largest margin and "lock it in".
   * Then you  pick the next largest one available
   * ("available" means: not already used and not creating a cycle), 
   * and continue on.  This creates an ordering of the candiates. The
   * topmost in the ordering wins.  It is a bit tricky to spot
   * the cycles as we go...
   * (This code based on email from Blake Cretney bcretney@postmark.net):
   ******************************************/

  /*path[][] is used as a changeable copy of pairels.*/
  for(i=0;i<C;i++){
    for(j=0;j<C;j++) path[i][j]=pairels[i][j];
    path[i][i]=INT_MAX;
  }
  /* Whenever a victory
   * is locked in, the appropriate cell is set to INT_MAX.
   * pi,pj are used with randperm to give precedence to victories higher
   * in the random permutation (when tie-breaking).
   * This loop finds the next pair (i,j) to lock in: ********/
  for(;;){
    int maxp,oi,oj,pi,pj;
    maxp=INT_MIN;
    for(pi=0;pi<C;pi++) for(pj=pi+1;pj<C;pj++){
      oi=randperm[pi]; oj=randperm[pj];
      if(   path[oi][oj]!=INT_MAX
	 && path[oj][oi]!=INT_MAX){/*not locked-out*/
	if(path[oi][oj]>maxp){ maxp=path[oi][oj]; i=oi; j=oj; }
	if(path[oj][oi]>maxp){ maxp=path[oj][oi]; i=oj; j=oi; }
      }
    }
    if(maxp==INT_MIN) break;
    /********* print the pair and its margin:
    printf("(%d %d) %d\n",i,j,maxp);
    ***********************/
    /*lock in the pair and clobber future no-good pairs:*/
    for(oi=0;oi<C;oi++) for(oj=0;oj<C;oj++){
      if(path[oi][i]==INT_MAX && path[j][oj]==INT_MAX){
	path[oi][oj]=INT_MAX;
      }
    }
  }
  /* The above code assumes that pairels has been set properly.
   * path[][] ends up with the winning row having all cells
   * set to INT_MAX.  In fact, a complete ranking is given,
   * where path[i][j]==INT_MAX means that i is
   * ranked over j (where i!=j).   So to find the winner: ****/
  winner = -99;
  for(i=0;i<C;i++){
    for(j=0;j<C;j++){
      if(path[i][j] != INT_MAX) break;
    }
    if(j>=C){ winner=i; break; }
  }
  assert(winner >= 0);
  stidemanregret = bestutil - sumutils[winner];  

  
  /* random dictator; */
  dictator = myrand() * Voters;
  assert(dictator < Voters);
  bestD = -99.;
  for(i=0; i<C; i++) if(bestD < igutl[dictator][randperm[i]]){ bestD = igutl[dictator][randperm[i]]; winner = randperm[i]; }
  hRdictregret = bestutil - sumutils[winner];

  /* random winner: */
  winner = myrand() * C;
  assert(winner<C);
  assert(0<=winner);
  hrandregret = bestutil - sumutils[winner];

  /* random-pair majority winner: */
  winner1 = myrand() * C;
  assert(winner1<C);
  assert(0<=winner1);
  do{
    winner2 = myrand() * C;
  }while(winner2 == winner1);
  assert(winner2<C);
  assert(0<=winner2);
  if(winner2<winner1){ t=winner1; winner1=winner2; winner2=t; }
  assert(winner2<C);
  assert(0<=winner1);
  assert(winner1<winner2);

  c1=c2=0;
  for(i=0; i<Voters; i++){
    if( igutl[i][winner1] < igutl[i][winner2] ) c2++; else c1++;
  }
  if(c2>=c1) winner=winner2; else winner=winner1;
  hRPMregret = bestutil - sumutils[winner];

  /* worst winner (pessimal): */
  hworstregret = bestutil - worstutil;

  /* honest approval: uses avg as threshhold for approval */
  for(i=0; i<C; i++) happrov[i] = 0;
  for(i=0; i<Voters; i++){
    avg = 0.;
    for(j=0; j<C; j++){
       avg += igutl[i][j];
    }
    avg /= C;
    for(j=0; j<C; j++){
       if( igutl[i][j] > avg ) happrov[j]++;
    }
  }
  best = -99;
  for(i=0; i<C; i++) if(best < happrov[randperm[i]]){ best = happrov[randperm[i]]; winner=randperm[i]; }
  hAPPAregret = bestutil - sumutils[winner];

  /***********************************************
   * honest approval: uses max-expected utility to choose threshhold *
   ***** HAS NOW BEEN ELIMINATED SINCE RECOGNIZE IS SAME THING AS ApHa ****
  for(i=0; i<C; i++) happrov[i] = 0;
  for(i=0; i<Voters; i++){
    bestudiff = -9999.9;
    for(th=1; th<C; th++){
       udiff = 0.;
       for(j=0; j<th; j++){
          for(k=th; k<C; k++){
             udiff += igutl[i][ hrank[i][j] ] - igutl[i][ hrank[i][k] ];
       }}
       if(udiff > bestudiff){ bestudiff = udiff; thw = th; }
    }
    for(j=0; j<thw; j++){   happrov[ hrank[i][j] ]++; }
  }
  best = -99;
  for(i=0; i<C; i++) if(best < happrov[randperm[i]]){ best = happrov[randperm[i]]; winner=randperm[i]; }
  hAPPUregret = bestutil - sumutils[winner];
  ******************************************************/

  /*Hbuck (honest Bucklin) ??? */
  for(i=0; i<C; i++) hplurvotes[i] = 0;
  for(rnd=0; rnd<C; rnd++){
    for(i=0; i<Voters; i++){
       hplurvotes[hrank[i][rnd]]++;
    }
    best = -99;
    for(i=0; i<C; i++) if(best < hplurvotes[randperm[i]]){ 
       best = hplurvotes[randperm[i]]; winner = randperm[i]; 
    }
    if(best*2 > Voters || rnd==C-1){
      hBucklinregret = bestutil - sumutils[winner];
      break;
    }
  }

/***********************now for strategic voting********************/

/*PlS (strategic plurality) vote for the best of the 2 frontrunners.*/
  for(i=0; i<C; i++) splurvotes[i] = 0;
  for(i=0; i<Voters; i++){
    if( igutl[i][0] < igutl[i][1] ) splurvotes[1]++; else splurvotes[0]++;
  }
  best = -99;
  for(i=0; i<C; i++) if(best < splurvotes[randperm[i]]){ best = splurvotes[randperm[i]]; winner=randperm[i]; }
  splurregret = bestutil - sumutils[winner];

  /* strategic IFPP3 (use rational plurality if c!=3): */
  if(C!=3){
    sIFPP3regret = splurregret;
  }

/*BuS (strategic bullet): vote against the worst of the 2 frontrunners.*/
  for(i=0; i<C; i++) sbullets[i] = 0;
  for(i=0; i<Voters; i++){
    if( igutl[i][0] < igutl[i][1] ) sbullets[0]++; else sbullets[1]++;
  }
  best = Voters*99;
  for(i=0; i<C; i++) if(best > sbullets[randperm[i]]){ best = sbullets[randperm[i]]; winner=randperm[i]; }
  sbulletregret = bestutil - sumutils[winner];

/*RaS (strategic range): use as your threshhold, a utility midway
* between the 2 frontrunners and you give +1 to candidates above
* threshhold, -1 to candidates below: ****/
  for(i=0; i<C; i++) srangevotes[i] = 0;
  for(i=0; i<Voters; i++){
    uthresh = (igutl[i][0] + igutl[i][1])*0.5;
    for(j=0; j<C; j++) 
       if(igutl[i][j] < uthresh){ srangevotes[j]--;  }
       else{ srangevotes[j]++; }
  }
  best = -Voters*99;
  for(i=0; i<C; i++) if(best < srangevotes[randperm[i]]){ best = srangevotes[randperm[i]]; winner=randperm[i]; }
  srangeregret = bestutil - sumutils[winner];

/*RaS2 (fancier strategic range): Uses moving avg threshhold. */
/*Can be dishonest, note. */
  for(i=0; i<C; i++) srangevotes[i] = 0;
  for(i=0; i<Voters; i++){
    for(j=0; j<C; j++) myvotes[j] = 0;
    if( igutl[i][0] < igutl[i][1] ) myvotes[1]++; else myvotes[0]++;
    umsum = igutl[i][0]+igutl[i][1];
    for(j=2; j<C; j++){
       myutil = igutl[i][j];
       if(myutil*j > umsum) myvotes[j]++;
       umsum += myutil;
    }
    for(j=0; j<C; j++) srangevotes[j] += myvotes[j];
    if(HONESTYSTATS){
      FunnyStatCt[C]++;
      /* check to see if RaS2 ever misorders votes vs utilities: */
      topmin = 999.; botmax = -999.;
      uthresh = (igutl[i][0] + igutl[i][1])*0.5;
      for(j=0; j<C; j++){
	if(myvotes[j]==1 && topmin>igutl[i][j]) topmin = igutl[i][j];
	if(myvotes[j]==0 && botmax<igutl[i][j]) botmax = igutl[i][j]; 
	if((myvotes[j]==1 && igutl[i][j]<uthresh) ||
	   (myvotes[j]==0 && igutl[i][j]>uthresh) ){
	  /********************
	    printf("RaS2 and RaS differ %f\n", fabs(igutl[i][j]-uthresh) );
	    for(k=0; k<C; k++) printf("%.3f ", igutl[i][k]);
	    printf("\n");
	    ****************/
	  RaS2RaSDiffer[C]++;
	}
      }
      if( topmin < botmax ){
	/*****************
        printf("RaS2 can mis-order %f\n", botmax-topmin);
	for(k=0; k<C; k++) printf("%.3f ", igutl[i][k]);
	printf("\n");
	**********************/
        RaS2MisOrder[C]++;
      }    /* End of checks. */
    }
  }
  best = -Voters*99;
  for(i=0; i<C; i++) if(best < srangevotes[randperm[i]]){ best = srangevotes[randperm[i]]; winner=randperm[i]; }
  srangeregret2 = bestutil - sumutils[winner];

/*BoS (strategic Borda) is, you give the best of the two frontrunners
* the max vote and the worst the min vote. Then you proceed recursively
* on the remaining C-2 candidates: ****/
  for(i=0; i<C; i++) sbordavotes[i] = 0;
  for(i=0; i<Voters; i++){
     maxBvote = C-1; minBvote = 0;
     for(j=0; j+1<C; j+=2){
        if( igutl[i][j] < igutl[i][j+1]){
           sbordavotes[j]   += minBvote;
           sbordavotes[j+1] += maxBvote;
        }else{
           sbordavotes[j]   += maxBvote;
           sbordavotes[j+1] += minBvote;
        }
        maxBvote--;
        minBvote++;
     }
     if(C%2 != 0){ 
        assert(minBvote==maxBvote); assert(minBvote*2==C-1); 
        sbordavotes[C-1] = minBvote;
     }
  }
  best = -Voters*C*99;
  for(i=0; i<C; i++) if(best < sbordavotes[randperm[i]]){ best = sbordavotes[randperm[i]]; winner=randperm[i]; }
  sbordaregret = bestutil - sumutils[winner];

/*BoS2 (rational Borda) is, you give the best of the two frontrunners
* the max vote and the worst the min vote. Then on the remaining C-2 
* candidates you vote each the top or bottom possible vote if above or
* below avg of previous candidate utilities (this is general generic rational
* strategy for any COAF system): ****/
  for(i=0; i<C; i++) sbordavotes[i] = 0;
  for(i=0; i<Voters; i++){
     maxBvote = C-1; minBvote = 0;
     if( igutl[i][0] < igutl[i][1]){
       sbordavotes[0] += minBvote;   minBvote++;
       sbordavotes[1] += maxBvote;   maxBvote--;
     }else{
       sbordavotes[0] += maxBvote;   maxBvote--;
       sbordavotes[1] += minBvote;   minBvote++;
     }
     umsum = igutl[i][0]+igutl[i][1];
     for(j=2; j<C; j++){
       myutil = igutl[i][j];
       if(myutil*j > umsum){
	 sbordavotes[j] += maxBvote;   maxBvote--;
       }else{
	 sbordavotes[j] += minBvote;   minBvote++;
       }
       umsum += myutil;
    }
    assert(maxBvote==minBvote-1);
  }
  best = -Voters*C*99;
  for(i=0; i<C; i++) if(best < sbordavotes[randperm[i]]){ best = sbordavotes[randperm[i]]; winner=randperm[i]; }
  sBoS2regret = bestutil - sumutils[winner];

/*DaS (rational Dabagh) is, you give the best of the two frontrunners
* the max vote and the worst the min vote. Then on the remaining C-2 
* candidates you vote each the top or bottom possible vote if above or
* below avg of previous candidate utilities (this is general generic rational
* strategy for any COAF system): ****/
  for(i=0; i<C; i++) sdabaghvotes[i] = 0;
  for(i=0; i<Voters; i++){
     if( igutl[i][0] < igutl[i][1]){ /*assume no utility ties*/
       sdabaghvotes[1] += 2;
     }else{
       sdabaghvotes[0] += 2;
     }
     umsum = igutl[i][0]+igutl[i][1];
     for(j=2; j<C; j++){
       myutil = igutl[i][j];
       if(myutil*j > umsum || j==C-1){
	 sdabaghvotes[j] += 1; break;
       }
       umsum += myutil;
     }
   }
  best = -99;
  for(i=0; i<C; i++) if(best < sdabaghvotes[randperm[i]]){ best = sdabaghvotes[randperm[i]]; winner=randperm[i]; }
  sDABAregret = bestutil - sumutils[winner];

/*BoS3 (strategic Borda) is, you give the best of the two frontrunners
* the max vote and the worst the min vote. Then on the remaining C-2 
* candidates you vote honestly ***/
  for(i=0; i<C; i++) hborda[i] = 0;
  for(i=0; i<Voters; i++){     for(j=0; j<C; j++){
    /* hborda[ drank[i][j] ] += C-j-1; this code should be equiv to below line: */
    hborda[j] += C-1-drank[i][j]; 
  }}  
  best = -99;
  for(i=0; i<C; i++) if(best < hborda[randperm[i]]){ best = hborda[randperm[i]]; winner = randperm[i]; }
  sBoS3regret = bestutil - sumutils[winner];
  Sbordwinner = winner;

  /* Strategic Black: if no Condorcet winner then use Borda:
  * Strat: you give the best of the two frontrunners
  * the max vote and the worst the min vote. Then on the remaining C-2
  * candidates you vote honestly ***/
  
  if(Scondwinner >= 0) winner = Scondwinner;
  else winner = Sbordwinner;
  sblackregret = bestutil - sumutils[winner];

  /* Strategic Meek123: */
  winner = MeekAlg123( drank, randperm, C );
  sMeek123regret = bestutil - sumutils[winner];
}


CheeseTest()
{
  int i,j,k,best,winner;
  double uthresh,umsum,myutil,topmin,botmax,bestutil=100.;
  int C=4;
  int srangevotes[MaxCandids];
  int myvotes[MaxCandids];
  double foou[] = {0.,9.,100.,11.};

/*RaS (strategic range): use as your threshhold, a utility midway
* between the 2 frontrunners and you give +1 to candidates above
* threshhold, -1 to candidates below: ****/
  for(i=0; i<C; i++) srangevotes[i] = 0;
  for(i=0; i<1; i++){
    uthresh = (foou[0] + foou[1])*0.5;
    for(j=0; j<C; j++) 
       if(foou[j] < uthresh){ srangevotes[j]--;  }
       else{ srangevotes[j]++; }
  }
  best = -1*99;
  for(i=0; i<C; i++) if(best < srangevotes[i]){ best = srangevotes[i]; winner=i; }
  srangeregret = bestutil - foou[winner];
  printf("RaS:  winner=%d srangeregret=%f\n", winner,srangeregret );
  printf("srangevotes=%d %d %d %d\n", 
	 srangevotes[0],srangevotes[1],srangevotes[2],srangevotes[3]);

/*RaS2 (fancier strategic range): Uses moving avg threshhold. */
/*Can be dishonest, note. */
  for(i=0; i<C; i++) srangevotes[i] = 0;
  for(i=0; i<1; i++){
    for(j=0; j<C; j++) myvotes[j] = 0;
    if( foou[0] < foou[1] ) myvotes[1]++; else myvotes[0]++;
    umsum = foou[0]+foou[1];
    for(j=2; j<C; j++){
       myutil = foou[j];
       if(myutil*j > umsum) myvotes[j]++;
       umsum += myutil;
    }
    for(j=0; j<C; j++) srangevotes[j] += myvotes[j];
    /* check to see if RaS2 ever misorders votes vs utilities: */
    topmin = 999.; botmax = -999.;
    uthresh = (foou[0] + foou[1])*0.5;
    for(j=0; j<C; j++){
       if(myvotes[j]==1 && topmin>foou[j]) topmin = foou[j];
       if(myvotes[j]==0 && botmax<foou[j]) botmax = foou[j]; 
       if((myvotes[j]==1 && foou[j]<uthresh) ||
          (myvotes[j]==0 && foou[j]>uthresh) ){
             printf("RaS2 and RaS differ %f\n", fabs(foou[j]-uthresh) );
             for(k=0; k<C; k++) printf("%.3f ", foou[k]);
	     printf("\n");
	   }
    }
    if( topmin < botmax ){
      printf("RaS2 can mis-order %f\n", botmax-topmin);
      for(k=0; k<C; k++) printf("%.3f ", foou[k]);
      printf("\n");
    }
    printf("topmin=%f  botmax=%f\n",topmin,botmax);
    /* end of checks. VERY peculiar no check is activated exptly!! */
  }
  best = -1*99;
  for(i=0; i<C; i++) if(best < srangevotes[i]){ best = srangevotes[i]; winner=i; }
  srangeregret2 = bestutil - foou[winner];
  printf("RaS2:  winner=%d srangeregret2=%f\n", winner,srangeregret2 );
  printf("srangevotes=%d %d %d %d\n", 
	 srangevotes[0],srangevotes[1],srangevotes[2],srangevotes[3]);
}


int main(int argc, char *argv[]){
  int i,j,z,C,kind_of_utils;
  int cand_ignorance = 0; /* ignorance not candidate-dependent by default */
  if(argc != 6){ 
     fprintf(stderr, "Need 5 args\nvotetest V E s Q CI for \nV voters, E experiments/datapt, seed s, ignorance Q, candidate-dept ignorance(bool)\n");
     exit(1);
  }
  Voters = atoi(argv[1]);
  if(Voters<0) Voters = -Voters;
  NumExperiments = atoi(argv[2]);
  if(NumExperiments<0) NumExperiments = -NumExperiments;
  RandSeed = atoi(argv[3]);
  IGNORANCEQ = atof(argv[4]);
  cand_ignorance = atoi(argv[5]); /* nonzero means ignorance will be candidate-dependent */

  if(IGNORANCEQ<0.0){
     fprintf(stderr, "Making IGNORANCEQ=0 since you specified negative value, moron.\n");
     fprintf(stderr, "0 means no voter ignorance; 1 is a fairly high ignorance setting.\n");
     IGNORANCEQ=0.0;
  }
  if(IGNORANCEQ>9.0){
     fprintf(stderr, "Making IGNORANCEQ=9 since you specified >9, moron.\n");
     fprintf(stderr, "0 means no voter ignorance; 1 is a fairly high ignorance setting.\n");
     IGNORANCEQ=9.0;
  }
  printf("Seeding with %u\n", RandSeed);
  srand48(RandSeed);
  initstate(RandSeed, randstate, 256);
  CoveyouX ^= (((uint64)RandSeed)<<2);
  ParkMillerX ^= RandSeed;
  ParkMillerX |= 1;
  ParkMillerX &= 2147483647;

  if(NumExperiments<1){
     fprintf(stderr, "At least 1 experiment needed, you said %d\n", NumExperiments);
     exit(1);
  }
  if(Voters>MaxVoters){
     fprintf(stderr, "No more than %d voters allowed, you used %d\n", MaxVoters,Voters);
     exit(1);
  }
  if(Voters<1){
     fprintf(stderr, "Need >=1 voter, fool: used %d\n", Voters);
     exit(1);
  }

  printf("drand48 %lf %lf %lf %lf %lf\n", drand48(),drand48(),drand48(),drand48(),drand48());
  printf("drand48 %f %f %f %f %f\n", drand48(),drand48(),drand48(),drand48(),drand48());

  printf("random %8x %8x %8x %8x %8x\n", random(),random(),random(),random(),random());
  printf("random %8x %8x %8x %8x %8x\n", random(),random(),random(),random(),random());
  printf("random %8x %8x %8x %8x %8x\n", random(),random(),random(),random(),random());
  printf("random %8x %8x %8x %8x %8x\n", random(),random(),random(),random(),random());
  printf("RAND_MAX %8x\n", RAND_MAX);

  printf("Coveyou %8x %8x %8x %8x %8x\n", Coveyou(),Coveyou(),Coveyou(),Coveyou(),Coveyou());
  printf("Coveyou %8x %8x %8x %8x %8x\n", Coveyou(),Coveyou(),Coveyou(),Coveyou(),Coveyou());
  printf("Coveyou %8x %8x %8x %8x %8x\n", Coveyou(),Coveyou(),Coveyou(),Coveyou(),Coveyou());
  printf("Coveyou %8x %8x %8x %8x %8x\n", Coveyou(),Coveyou(),Coveyou(),Coveyou(),Coveyou());
  printf("Coveyou %8x %8x %8x %8x %8x\n", Coveyou(),Coveyou(),Coveyou(),Coveyou(),Coveyou());
  printf("Coveyou %8x %8x %8x %8x %8x\n", Coveyou(),Coveyou(),Coveyou(),Coveyou(),Coveyou());

  printf("ParkMiller %8x %8x %8x %8x %8x\n", ParkMiller(),ParkMiller(),ParkMiller(),ParkMiller(),ParkMiller());
  printf("ParkMiller %8x %8x %8x %8x %8x\n", ParkMiller(),ParkMiller(),ParkMiller(),ParkMiller(),ParkMiller());

  printf("myrand %f %f %f %f %f\n", myrand(),myrand(),myrand(),myrand(),myrand());
  printf("myrand %f %f %f %f %f\n", myrand(),myrand(),myrand(),myrand(),myrand());
  printf("myrand %f %f %f %f %f\n", myrand(),myrand(),myrand(),myrand(),myrand());
  printf("myrand %f %f %f %f %f\n", myrand(),myrand(),myrand(),myrand(),myrand());

  if(CHEESYCHECK){
    printf("CheeseTest:\n");
    CheeseTest();
    printf("end CheeseTest.\n");
  }

  printf("\nsystems\n");
  printf(  "-------\n");
  for(i= -1; i<MaxSystems; i++){
     if(i== -1)  printf("(Baseline=Best-summed-utility winner; regret=0)\n");
     if(i==RaH)  printf("%2d. Honest range voting (scaled utility vote)\n", i);
     if(i==BoH)  printf("%2d. Honest Borda\n", i);
     if(i==LRH)  printf("%2d. Honest Condorcet Least-Reversal (LR)\n", i);
     if(i==STVH) printf("%2d. Honest Hare Single Transferable Vote STV (least most-liked canddt eliminated)\n", i);
     if(i==PlH)  printf("%2d. Honest plurality (1 vote for max-utility canddt)\n", i);
     if(i==BuH)  printf("%2d. Honest bullet (1 vote against min-utility cand)\n", i);
     if(i==RPM)  printf("%2d. Majority vote on random candidate pair\n", i);
     if(i==RandomDict) printf("%2d. Random dictator (chosen from voters) dictates winner\n", i);
     if(i==RandomWin) printf("%2d. Random winner\n", i);
     if(i==WorstWin) printf("%2d. Worst-summed-utility winner\n", i);
     if(i==ApHa) printf("%2d. Honest approval (using threshhold=average candidate utility)\n", i);
     if(i==RaS)  printf("%2d. Strategic range/approval (average of 2 frontrunner utils as thresh)\n", i);
     if(i==RaS2) printf("%2d. Rational range/approval (threshhold=moving average)\n", i);
     if(i==PlS)  printf("%2d. Rational plurality (vote for 1 of 2 frontrunners)\n", i);
     if(i==BoS)  printf("%2d. Strategic Borda I (1 frontrunner top, 1 bottom, rest recursively)\n", i);
     if(i==BuS)  printf("%2d. Rational bullet (vote against 1 of 2 frontrunners)\n", i);
     if(i==LRS)  printf("%2d. Strategic Condorcet LR (strat same as %d)\n",i, BoS3);
     if(i==STVS)  printf("%2d. Strategic Hare STV (strat same as %d)\n",i, BoS3);
     if(i==BoS2)  printf("%2d. Rational Borda (1 frontrunner max, 1 min, rest using moving avg to decide if max or min vote)\n",i);
     if(i==CoombsH)  printf("%2d. Honest Coombs STV (most least-liked candid eliminated each round)\n",i);
     if(i==CoombsS)  printf("%2d. Strategic Coombs STV (strat same as %d)\n",i,BoS3);
     if(i==BoS3)  printf("%2d. Strategic Borda II (1 frontrunner max, 1 min vote, rest honest)\n",i);

     if(i==Hbuck)  printf("%2d. Honest Bucklin \n",i);
     if(i==Hcope)  printf("%2d. Honest Copeland (winner of most pairwise elections)\n",i);
     if(i==Hdaba)  printf("%2d. Honest Dabagh point-and-a-half \n",i);
     if(i==DaS)  printf("%2d. Rational Dabagh point-and-a-half (moving avg)\n",i);
     if(i==CopeS)  printf("%2d. Strategic Copeland (strat same as %d)\n",i,BoS3);
     if(i==PlRunH)  printf("%2d. Honest plurality + runoff for 2 top finishers\n", i);
     if(i==BlackH)  printf("%2d. Honest Black (if no Condorcet winner use Borda)\n", i);
     if(i==BlackS)  printf("%2d. Strategic Black (same strat as %d)\n", i,BoS3);

     if(i==HIFPP3)  printf("%2d. Carey-IFPP3, Honest voters (3 cands only, else plurality)\n", i);
     if(i==SIFPP3)  printf("%2d. Carey-IFPP3, Strategic voters (3 cands only, else use rational plurality)\n", i);

     if(i==HMeek123)  printf("%2d. MeekAlg123, Honest voters\n", i);
     if(i==SMeek123)  printf("%2d. MeekAlg123, Strategic voters (same strat as %d)\n", i,BoS3);
     if(i==Htide)  printf("%2d. Tideman, Honest voters\n", i);
     if(i==Stide)  printf("%2d. Tideman, Stategic voters (same strat as %d)\n", i,BoS3);
   }

  for(kind_of_utils=0; kind_of_utils<=5; kind_of_utils++){
    /* Note it IS possible to have nonzero regret in 2-candidate elections */
    if(kind_of_utils==RANDOM_UTILITIES ){ printf("\nRandom 0-1 Utilities. (0 issues.) "); }
    else if(kind_of_utils==NORMAL_UTILITIES ){ 
      printf("\nRandom-Normal(0,1) Utilities (infinite issue limit). ");
    }
    else{ 
      printf("\nIssue Based Utilities (%d Issue", kind_of_utils); 
      if( kind_of_utils > 1 ) printf("s). "); else printf("). ");
    }
    printf("IgnoranceQ=%.2f. ", IGNORANCEQ);
    if(cand_ignorance){
       printf("Canddt-Dependent. ");
    }else{
       printf("(Identical.) ");
    }

    printf("%d voters.\n", Voters);
    if(kind_of_utils!=NORMAL_UTILITIES ){
       printf("Each candidate utility (for each voter)\n");
       printf("normalized to lie somewhere in [0,1].\n");
    }
    printf("Each regret datapoint averages %u expts.\n", NumExperiments);
      
    for(i=0; i<MaxCandids; i++){ /*HONESTYSTATS*/
       RaS2MisOrder[i]=0; RaS2RaSDiffer[i]=0; FunnyStatCt[i]=0;
    }
    for(C=2; C<=5; C++){
      for(j=0; j<MaxSystems; j++){ regrets[C][j]=0.0; }
      for(z=0; z<C; z++){
	  shstvw[z]=0; sclrw[z]=0;
      }
      for(i=0; i<NumExperiments; i++){
        genutils(kind_of_utils,cand_ignorance,C);
        dovotes(C);
        regrets[C][PlH] += hplurregret; 
        regrets[C][BuH] += hbulletregret; 
        regrets[C][BoH] += hbordaregret; 
        regrets[C][RaH] += hrangeregret; 
        regrets[C][STVH] += hSTVregret; 
        regrets[C][LRH] += hLRregret; 
        regrets[C][RandomDict] += hRdictregret; 
        regrets[C][RandomWin] += hrandregret; 
        regrets[C][WorstWin] += hworstregret; 
        regrets[C][RPM] += hRPMregret; 
        regrets[C][ApHa] += hAPPAregret; 
        regrets[C][CoombsH] += hCoombsregret; 

        regrets[C][Hbuck] += hBucklinregret; 
        regrets[C][Hcope] += hCOPEregret; 
        regrets[C][Hdaba] += hDABAregret; 
        /*regrets[C][ApHu] += hAPPUregret; */

        regrets[C][RaS2] += srangeregret2;
        regrets[C][RaS] += srangeregret;
        regrets[C][PlS] += splurregret;
        regrets[C][BoS] += sbordaregret;
        regrets[C][BuS] += sbulletregret;
        regrets[C][LRS] += sLRregret;
        regrets[C][STVS] += sSTVregret; 
        regrets[C][BoS2] += sBoS2regret; 
        regrets[C][BoS3] += sBoS3regret; 
        regrets[C][CoombsS] += sCoombsregret; 

        regrets[C][DaS] += sDABAregret;
        regrets[C][CopeS] += sCOPEregret;
        regrets[C][PlRunH] += hplurRUNregret; 
        regrets[C][BlackH] += hblackregret; 
        regrets[C][BlackS] += sblackregret; 

        regrets[C][HIFPP3] += hIFPP3regret; 
        regrets[C][SIFPP3] += sIFPP3regret; 

        regrets[C][HMeek123] += hMeek123regret; 
        regrets[C][SMeek123] += sMeek123regret; 
        regrets[C][Htide] += htidemanregret;
        regrets[C][Stide] += stidemanregret;
      }
#if CONDHARESTATS
      printf("#cand. Hare, Condorcet\n");
      for(z=0; z<C; z++){
	printf("%3d.  %9d %9d\n", z, shstvw[z], sclrw[z]);
      }
#endif        
    }
    
    if(HONESTYSTATS){
      printf("HONESTYSTATS (re strategic range voting):\n");
      printf("RaS2 misorder: ");
      for(i=2; i<=5; i++){
	printf("%7.5f ", (double)RaS2MisOrder[i]/(double)FunnyStatCt[i]);
      }
      printf("\nRaS2/1 differ: ");
      for(i=2; i<=5; i++){
	printf("%7.5f ",  (double)RaS2RaSDiffer[i]/(double)FunnyStatCt[i]);
      }
      printf("\nCount:         ");
      for(i=2; i<=5; i++){
	printf("%7g ",  (double)FunnyStatCt[i]);
      }
      printf("\n");
    }

    printf("system|2 canddts 3 canddts 4 canddts 5 canddts\n");
    printf("------+--------- --------- --------- ---------\n");
    for(j=0; j<36; j++){   /* key: # of systems is upper loop limit here */
      printf("   %2d |", j);
      for(C=2; C<=5; C++){
	regrets[C][j] /= NumExperiments; 
	printf("%9.5f ", regrets[C][j] );
      }
      printf("\n");
    }
    if(kind_of_utils==RANDOM_UTILITIES){
       printf("\n%3d RandUtil[4cand] %7.5f %7.5f %5.3f %7.5f %7.5f %5.3f\n",
          Voters,
          regrets[4][RaH],
          regrets[4][PlH],
          regrets[4][PlH]/regrets[4][RaH],
          regrets[4][RaS2],
          regrets[4][PlS],
          regrets[4][PlS]/regrets[4][RaS2] );
    }
    printf(
     "Bootstrap additive error estimates for WorstWin from RandomWin*2\n");
    printf("(from above data with 2-5 cands):\n");
    printf("bterr |");
    for(C=2; C<=5; C++){
       printf("%9.5f ", fabs(regrets[C][RandomWin]*2 - regrets[C][WorstWin]));
    }
    printf("\n");
    fflush(stdout);
  }
}
