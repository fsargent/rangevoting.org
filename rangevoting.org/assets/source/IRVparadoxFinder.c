// gcc IRVparadoxFinder.c -lm -lc -Wall -O6 -o IRVparadoxFinder

//possible other paradoxes to add:
//incentive to vote dishonestly?
//coalitional manipulation?

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>

#define until(x) while(!(x))
#define and &&
#define or ||
#define uint unsigned int
#define TRUE (0<1)
#define FALSE (1<0)

#define uint unsigned int
#define bool char
#define uint32 unsigned int
#define uint64 unsigned long long
#define uchar unsigned char
#define schar signed char
#define real double

#define PI 3.1415926535897932385
#define ROOTPI 1.7724538509055160273
#define ROOT2 1.4142135623730950488
#define ROOT3 1.7320508075688772935

#define MARSGALIA TRUE  //use Marsaglia or me?  Both same speed
/**********Marsgalia rand gen: **********/
#define MAXUINT64  ((uint64)((255ULL<<56)||(255ULL<<48)|(255ULL<<40)|(255ULL<<32)|(255ULL<<24)|(255ULL<<16)|(255ULL<<8)|(255ULL)))

uint64 LFt[256];
uchar LFc;
#define LFIB4 (LFc++, LFt[LFc]=LFt[LFc]+LFt[UC(LFc+58)]+LFt[UC(LFc+119)]+LFt[UC(LFc+178)])
#define UC    (unsigned char)  /*a cast operation*/

/****************** My rand gen: ***********************/
/****** convenient constants: *******/
#define BIGINT 0x7FFFFFFF
#define MAXUINT ((uint)((255<<48)|(255<<40)|(255<<32)|(255<<24)|(255<<16)|(255<<8)|(255)))
/* defn works on 8,16,32, and 64-bit machines */

uint32 BLC32x[60];  /* 32*60=1920 bits of state. Must be nonzero mod P. */
int BLC32NumLeft;

/********************************************************
Warren D. Smith 2001
**********************************************************
Linear congruential pseudo-random number generator mod P,
where P is the enormous prime (578 base-10 digits long;
60 words long if each word is 32 bits)
  P = [2^(48*32) - 1] * 2^(12*32) + 1.
This prime can yield PRNGs suitable for use on
computers with w-bit words, w=8,16,32,48,64,96,128.
The following code is intended for w=32.
The fact that 2^(60*32) = 2^(12*32) - 1 (mod P)
makes modular arithmetic mod P easy, and is the
reason this particular P was chosen.
The period of our generator is P-1.
***********************************************************
Although it is usually easy to detect the nonrandomness of
linear congruential generators because they generate d-tuples
lying on a lattice, in the present case the "grain size" of the
lattice is invisibly small (far less than a least significant bit),
for 32-bit words, if 1<=d<=180. If you want 64-bit words, then need
to concatenate two 32-bit words, and then grain size invisibly small
if 1<=d<=90. These bounds are conservative; in fact I suspect
one is ok, even for 64-bit words, even in up to 1000 dimensions.
***********************************************************
Some even more enormous primes which we have not used are:
[2^{59*32} - 1] * 2^{8 *32} + 1,
[2^{63*32} - 1] * 2^{24*32} + 1,
[2^{69*32} - 1] * 2^{14*32} + 1,
[2^{95*32} - 1] * 2^{67*32} + 1,
[2^{99*32} - 1] * 2^{35*32} + 1;
these would also be suitable for (8,16,32)-bit computers,
and the second of them would also be good for (48,96)-bit computers.
Unfortunately the factorization of P-1 is not known for the last 3
I've listed here, preventing you from being certain you've found a
primitive root mod that P. A still more enormous prime is
  [2^4253 - 1] * 2^4580 + 1    [2659 digits long!]
(where note 2^4253-1 is also prime so that factorizing P-1 is
trivial) but doing arithmetic mod this P is (although still fairly
easy) less pleasant because bit-shifting is required.
*************************************************************/
  uint32 BigLinCong32(){
  uint32 y[120];
  int i;
  uint64 u;

  if(BLC32NumLeft==0){
    /* Need to refill BLC32x[0..59] with 60 new random numbers: */

    /****************************************************************
     * If BLC32x[0..59] is the digits, LS..MS, of a number in base 2^w,
     * then the following code fragment puts A times that number
     * in y[0..119].  Here
     *  A = 1284507170 * 2^(w*3) + 847441413 * 2^(w*44) + 650134147 * 2^(w*59)
     * is a "good" primitive root mod P, if w=32.
     *****************************************************************/
#define lohalf(x) (uint32)(x)
#define A1 (uint64)1284507170
#define A2 (uint64)847441413
#define A3 (uint64)650134147
    for(i=0; i<3; i++){
      y[i] = 0;
    }
    u=0;
    for(/*i=3*/; i<44; i++){
      u += A1 * BLC32x[i-3];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=44*/; i<59; i++){
      u += A1 * BLC32x[i-3];
      u += A2 * BLC32x[i-44];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=59*/; i<60+3; i++){
      u += A1 * BLC32x[i-3];
      u += A2 * BLC32x[i-44];
      u += A3 * BLC32x[i-59];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=60+3*/; i<60+44; i++){
      u += A2 * BLC32x[i-44];
      u += A3 * BLC32x[i-59];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=60+44*/; i<60+59; i++){
      u += A3 * BLC32x[i-59];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=60+44*/; i<60+59; i++){
      u += A3 * BLC32x[i-59];
      y[i] = lohalf(u);
      u = u>>32;
    }
    /*i=60+59=119*/
    y[i] = lohalf(u);
#undef A1
#undef A2
#undef A3
    /*************************************************************
     * If y[0..119] is the digits, LS..MS, of a number in base 2^w,
     * then the following code fragment replaces that number with
     * its remainder mod P in y[0..59]  (conceivably the result will
     * be >P, but this does not matter; it will never be >=2^(w*60)).
     **************************************************************/
    u=1; /*borrow*/
#define AllF 0xffffffff
    /* Step 1: y[0..72] = y[0..59] + y[60..119]shift12 - y[60..119]: */
    for(i=0; i<12; i++){
      u += y[i];
      u += (uint64)~y[60+i];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=12*/; i<60; i++){
      u += y[i];
      u += y[48+i];
      u += (uint64)~y[60+i];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=60*/; i<72; i++){
      u += AllF;
      u += y[48+i];
      y[i] = lohalf(u);
      u = u>>32;
    }
    assert(u>0);
    y[72] = (uint32)(u-1); /*unborrow*/

    /*  Step 2: y[0..60] = y[0..59] + y[60..72]shift12  - y[60..72]: */
    u=1; /*borrow*/
    for(i=0; i<12; i++){
      u += y[i];
      u += (uint64)~y[60+i];
      y[i] = lohalf(u);
      u = u>>32;
    }
    /*i=12*/
    u += y[i] + y[48+i];
    u += (uint64)~y[60+i];
    y[i] = lohalf(u);
    u = u>>32;
    i++;
    for(/*i=13*/; i<25; i++){
      u += AllF;
      u += y[i];
      u += y[48+i];
      y[i] = lohalf(u);
      u = u>>32;
    }
    for(/*i=25*/; i<60; i++){
      u += AllF;
      u += y[i];
      y[i] = lohalf(u);
      u = u>>32;
    }
    /*i=60*/
    assert(u>0);
    y[i] = (uint32)(u-1); /*unborrow*/

    /*It is rare that any iterations of this loop are needed:*/
    while(y[60]!=0){
      printf("rare loop\n");
      /*Step 3+:  y[0..60] = y[0..59] + y[60]shift12 - y[60]:*/
      u=1; /*borrow*/
      u += y[0];
      u += (uint64)~y[60];
      y[0] = lohalf(u);
      u = u>>32;
      for(i=1; i<12; i++){
        u += AllF;
        u += y[i];
        y[i] = lohalf(u);
        u = u>>32;
      }
      /*i=12*/
      u += AllF;
      u += y[i];
      u += y[60];
      y[i] = lohalf(u);
      u = u>>32;
      i++;
      for(/*i=13*/; i<60; i++){
        u += AllF;
        u += y[i];
        y[i] = lohalf(u);
        u = u>>32;
      }
      /*i=60*/
      assert(u>0);
      y[i] = (uint32)(u-1); /*unborrow*/
    }
#undef AllF
#undef lohalf

    /* Copy y[0..59] into BLC32x[0..59]: */
    for(i=0; i<60; i++){
      BLC32x[i] = y[i];
    }
    /*printf("%u\n", BLC32x[44]);*/
    BLC32NumLeft=60;
  }
/* (Else) We have random numbers left, so return one: */
BLC32NumLeft--;
return BLC32x[BLC32NumLeft];
}

#if MARSAGLIA
real Rand01(){ /* returns random uniform in [0,1] */
  return LFIB4 / MAXUINT64;
}
#else
real Rand01(){ /* returns random uniform in [0,1] */
return ((BigLinCong32()+0.5)/(1.0 + MAXUINT) + BigLinCong32())/(1.0 + MAXUINT);
}
#endif


bool RandBool(){ /* returns random boolean */
if( Rand01() > 0.5 ) return TRUE;
return FALSE;
}

real RandExpl(){ /* returns standard exponential (density=exp(-x) for x>0) random deviate */
real x;
do{
x = Rand01();
}while( x==0.0 );
return -log(x);
}


real RandNormal(){ /* returns standard Normal (gaussian variance=1 mean=0) deviate */
real w, x1;
static real x2;
static bool ready = FALSE;
if(ready){
ready = FALSE;
return x2;
}
do{
x1 = 2*Rand01() - 1.0;
x2 = 2*Rand01() - 1.0;
w = x1*x1 + x2*x2;
}while ( w > 1.0 || w==0.0 );
w = sqrt( (-2.0*log(w)) / w );
x1 *= w;
x2 *= w;  /* Now x1 and x2 are two indep normals (Box-Muller polar method) */
ready = TRUE;
return x1;
}


void InitRand(uint seed){ /* initializes the randgen */
int i;
int seed_sec=0, processId=0;
uint seed_usec=0;
#if     MSWINDOWS
tm* locTime;
_timeb currTime;
time_t now;
#else
struct timeval tv;
#endif
if(seed==0){
printf("using time of day and PID to generate seed");
#if MSWINDOWS
now = time(NULL);
locTime = localtime(&now);
_ftime(&currTime);
seed_usec = currTime.millitm*1001;
seed_sec = locTime->tm_sec + 60*(locTime->tm_min + 60*locTime->tm_hour);
processId = _getpid();
#else
gettimeofday(&tv,0);
seed_sec = tv.tv_sec;
seed_usec = tv.tv_usec;
processId = getpid();
#endif
seed = 1000000*(uint)seed_sec + (uint)seed_usec +
(((uint)processId)<<20) + (((uint)processId)<<10);
printf("=%u\n", seed);
}
for(i=0; i<256; i++) LFt[i] = i*i*seed + i*i*i*i*i;
 for(i=0; i<99999; i++) LFIB4;
for(i=0; i<60; i++){ BLC32x[i]=0; }
BLC32x[0] = seed;
BLC32NumLeft = 0;
for(i=0; i<599; i++){ BigLinCong32(); }
printf("Random generator initialized with seed=%u:\n", seed);
for(i=0; i<7; i++){
printf("%.6f ", Rand01());
}
printf("\n");
for(i=0; i<7; i++){
printf("%.6f ", RandNormal());
}
printf("\n");
}

/*---------------------*/

void PrintBinary(uint x){ //10 bits
  uint j;
  for(j=512; j>0; j/=2){
    putchar('0'+((x&j)!=0));
  }
}

#define Classify()     \
                idx=0; dm=IPD;                                  \
		if(PPP){ idx += PPPi; PPPct++; PPPctD+=dm; }		\
		if(NAP){ idx += NAPi; NAPct++; NAPctD+=dm; }		\
		if(NPP){ idx += NPPi; NPPct++; NPPctD+=dm; }		\
		if(PAP){ idx += PAPi; PAPct++; PAPctD+=dm; }			\
		if(PPP or NAP or NPP or PAP){ PFct++; PFctD+=dm; }		\
		if(PPP or NPP){ GMct++; GMctD+=dm; }				\
		if(PAP or NAP){ ABct++; ABctD+=dm; }				\
		if((PPP or NPP) and (PAP or NAP)){ BPAct++; BPActD+=dm; }	\
		if(MLP){ idx += MLPi; MLPct++; MLPctD+=dm; }			\
		/*if(MLP != MLP2){ printf("yikes!\n"); }*/		        \
		if(LMP){ idx += LMPi; LMPct++; LMPctD+=dm; }			\
		if(LMP or MLP){ NMct++; NMctD+=dm; }             		\
		if(LMP and MLP){ BNMct++; BNMctD+=dm; }            		\
		if(CWE){ idx += CWEi; CWEct++; CWEctD+=dm; }			\
		if(REV){ idx += REVi; REVct++; REVctD+=dm; }			\
		if(IPD){ idx += IPDi; PDct++; PDctD+=dm; }			\
		if(CYC){ idx += CYCi; CYct++; CYctD+=dm; }			\
		if(ALL){ idx += ALLi; ALLct++; ALLctD+=dm; }			\
		if(LMP or FBC){ FB2ct++; FB2ctD+=dm; }			        \
		if(LDO){ idx += LDOi; LDOct++; LDOctD+=dm; }			\
		if(ALL or REV or CWE or LMP or MLP or PAP or PPP or NPP or NAP or LDO){ \
                       ANYct++; ANYctD+=dm; }					\
		qu = 6*V + nz;						

//"Standard situation": A is IRVwinner and B the second-placer iff
#define STD (CAB+CBA < ABC+ACB and CAB+CBA < BAC+BCA and ABC+ACB+CAB > BAC+BCA+CBA)

//LDO: loser drop-out paradox.  B can remove self from race, flipping winner from A to C.
#define LDO  (BCA+CAB+CBA > tol+ABC+ACB+BAC)

//FBC: favorite betrayal.  It pays for some BCA voters to betray B, causing C (or B
//in the case LMP, computed separately) to win.   In the above LDO scenario,
//if instead of B dropping out, all the BCA voters simply betray B (now CBA), that has the same 
//effect as B dropping out (unless the remaining BAC voters are so many as to outnumber the
//CAB+CBA+BCA voters in which case C is eliminated as usual and hence A wins as usual,
//but then there would not have been an LDO paradox).  So we conclude that this FBC
//paradox occurs exactly when LDO occurs and we do not need a new test.
//See also http://www.rangevoting.org/IRVStratPf.html
#define FBC LDO

//PPP: possible to add new voters ranking current winner top, such that he now loses.
//Lepelley-Merlin prop1 p57-58, assumes Awins, B2nd. Assumes s!=0 (not IRV).
//Can add ACB voters to make C win and A lose.
#define PPP  (ABC+t*ACB+s*BAC+CBA>tol+s*BCA+t*CAB \
    and s*ABC+BAC+t*BCA>tol+s*ACB+CAB+t*CBA \
    and ABC+ACB+CAB>tol+BAC+BCA+CBA \
    and s*BCA+(1+s)*CAB+CBA>tol+2*s*ABC+(1+s)*BAC+t*BCA)

//NAP: possible to delete voters ranking current winner bottom, such that he now loses.
//Lepelley-Merlin prop2 p58, assumes A=winner, B=2ndplacer.  (Careful: reversed ineq in B5.)
//Can delete some BCA voters causing A to lose and C to win.  Assumes s!=1 (not Coombs).
#define NAP (ABC+t*ACB+s*BAC>tol+s*BCA+t*CAB+CBA \
    and s*ABC+BAC+t*BCA>tol+s*ACB+CAB+t*CBA \
    and ABC+ACB+CAB>tol+BAC+BCA+CBA \
    and s*ACB+(1+t)*CAB+2*t*CBA>tol+ABC+t*ACB+(1+t)*BAC	\
    and s*ACB+CAB+t*CBA>tol+s*ABC+BAC)

//NPP: possible to add new voters ranking a current loser bottom, such that he now wins.
//Lepelley-Merlin prop3 p59, assume IRV A=win B=2nd. Assumes s!=1 (not Coombs).
//Can add CAB voters to make A lose and B win.
#define NPP (ABC+t*ACB+s*BAC>tol+s*BCA+t*CAB+CBA  \
    and t*BAC+(t+s*s)*BCA+(2-s)*s*CBA>tol+(s*s+t)*ABC+t*(s+1)*ACB+s*BAC \
    and ABC+ACB+CAB>tol+BAC+BCA+CBA \
    and t*BAC+BCA+s*CBA>tol+s*ABC+2*t*ACB+s*BAC)

//PAP: possible to delete voters ranking a current loser top, such that he now wins.
//Lepelley-Merlin prop4 p59, assumes IRV A=win, B=2nd. (Careful: reversed ineq in D5.)
//Can delete some BAC voters causing A to lose and B to win. Assumes s!=0 (not IRV).
#define PAP (ABC+t*ACB+s*BAC>tol+s*BCA+t*CAB+CBA \
     and (2-s)*s*BCA+t*CAB+CBA>tol+(s+1)*t*ABC+(s*s+t)*ACB+s*CAB+t*s*CBA \
     and ABC+ACB+CAB>tol+BAC+BCA+CBA \
     and 2*s*BCA+t*CAB+t*CBA>tol+t*ABC+ACB+s*CAB+t*CBA \
     and s*BCA+t*CAB+CBA>tol+ABC+t*ACB)

//MLP if some voters raise their ranking of the current winner, that makes him lose.
//prop1 p136 Lepelley-Chantreuil-Berg1996 assuming IRVwinner=A.
//Proven in Berg&Lepelley 1993.
//Can raise A to top in some BAC and BCA votes to make C win.
#define MLP ( (BAC+BCA+CBA>tol+ABC+ACB+CAB and 4*(BAC+BCA)>V+tol) or (CAB+CBA+BCA>tol+ABC+ACB+BAC and 4*(CAB+CBA)>V+tol) )

//More concise equivalent form of condition by Warren D. Smith.
#define MLP2   ( min(2*(CBA+CAB), ABC+ACB+BAC+BCA, 2*BCA) + CAB+CBA >tol+ ABC+ACB+BAC+BCA )

//Quas says this is equivalent to MLP for Quas model
#define MLPq ( (CBA+CAB)*4>V+tol and CAB+CBA+BCA>tol+ACB+ABC+BAC and CBA+CAB+ACB>tol+BCA+BAC+ABC )

//LMP  if some voters lower their ranking of a current loser, that makes him win.
//prop3 p139 Lepelley-Chantreuil-Berg1996 assuming IRVwinner=A.
//Can lower B to cause B to win.
#define LMP ( (3*(ABC+ACB)<V-tol) and ((BAC+BCA>tol+ABC+ACB and BCA+CAB+CBA>tol+ABC+ACB and 2*ABC+4*ACB<V-tol) \
		   or (CAB+CBA>tol+ABC+ACB and CBA+BAC+BCA>tol+ABC+ACB and 2*ACB+4*ABC<V-tol)) )

// CWE: C (eliminated by IRV in std form) is Condorcet Winner
#define CWE ( CAB+CBA+BCA>tol+ACB+ABC+BAC and CBA+CAB+ACB>tol+BCA+BAC+ABC )

// REV: A is the 'winner' with reversed ballots
#define REV ( (BAC+ABC +tol< CBA+BCA and BAC+ABC +tol< CAB+ACB and CBA+BCA+BAC >tol+ CAB+ACB+ABC) or \
              (CAB+ACB +tol< CBA+BCA and CAB+ACB +tol< BAC+ABC and CBA+BCA+CAB >tol+ ABC+BAC+ACB) )

// IPD: IRV & Plurality winners disagree.  (Assumes IRVwinner=A, IRV2nd=B; demands PlurWinner=B.)
#define IPD (BAC+BCA>tol+ABC+ACB)

// CYC: Condorcet cycle.
#define CYC ((ABC+ACB+CAB>tol+BAC+BCA+CBA and CBA+BCA+CAB>tol+ABC+ACB+BAC and BAC+ABC+BCA>tol+CBA+CAB+ACB) \
  or (ABC+ACB+CAB+tol<BAC+BCA+CBA and CBA+BCA+CAB+tol<ABC+ACB+BAC and BAC+ABC+BCA+tol<CBA+CAB+ACB))

// ALL: All scoring systems agree winner is B.
#define ALL (BAC+BCA>tol+ABC+ACB and CAB+ACB+tol<ABC+BAC and CAB+ACB+tol<CBA+BCA)


//   ABC  ACB  BAC  BCA  CAB  CBA  with sum V = total#voters.
//    1    2    3    4    5    6
//S. Berg & D. Lepelley:
//Note sur le calcul de la probabilite des paradoxes du vote, 
//Mathematiques, Informatique et Sciences Humaines 121 (1993) 33-48.
//William V. Gehrlein & Dominique Lepelley:
//The probability that all weighted scoring rules elect the same winner,
//Economics Letters 66 (2000) 191–197.
//V.Merlin, M.Tataru, F.Valognes: On the probability that all the rules select the same winner,
//1997.
//Dominique Lepelley &amp; Vincent Merlin: Scoring run-off paradoxes for variable electorates,
//Economic Theory 17 (2001) 53-80.
//Dominique Lepelley, Frederic Chantreuil, Sven Berg:
//The likelihood of monotonicity paradoxes in run-off elections,
//Mathematical Social Sciences 31 (1996) 133-146.

#define PPPi 2048
#define PAPi 1024
#define LDOi 512
#define ALLi 256
#define CYCi 128
#define IPDi 64
#define MLPi 32 //works modulo ties
#define LMPi 16 //works modulo ties
#define NAPi 8 //works
#define NPPi 4 //works
#define CWEi 2 //
#define REVi 1 //works

uint qrec[4096];
uint64 ctcn[4096];
uint64 ctNorm[4096];
uint64 ctDirich[4096];
uint64 ctQuas1D[4096];
uint vrec[4096][8];

int min(int x, int y, int z){ if(x<y && x<z) return(x); if(y<z) return(y); return(z); }

#define Sort2(a,b) if(y[a]>y[b]){ t=y[a]; y[a]=y[b]; y[b]=t; }
void Gimme5sortedRands(real y[]){
  real t;
  y[0] = Rand01();
  y[1] = Rand01();
  y[2] = Rand01();
  y[3] = Rand01();
  y[4] = Rand01();
 //damn idiot in Wikipedia had supplied a wrong Sorting network. But this one is ok...
  Sort2(1, 2); 
  Sort2(0, 2);
  Sort2(0, 1);
  Sort2(4, 5);
  Sort2(3, 5);
  Sort2(3, 4);
  Sort2(0, 3);
  Sort2(1, 4);
  Sort2(2, 5);
  Sort2(2, 4);
  Sort2(1, 3);
  Sort2(2, 3);
  //...as is verified by these asserts:
  //assert(y[0]<=y[1]); assert(y[1]<=y[2]); assert(y[3]<=y[3]); 
  //assert(y[3]<=y[4]); assert(y[4]<=y[5]);
}

void Gimme3sortedRands(real y[]){
  real t;
  y[0] = Rand01();
  y[1] = Rand01();
  y[2] = Rand01();
  Sort2(0,1); Sort2(1,2); Sort2(0,1);
}
#undef Sort2

#define UpdateRecs() \
		ctcn[idx]++;						\
		if(qu < qrec[idx]){					\
		  qrec[idx] = qu;					\
		  vrec[idx][0] = ABC;					\
		  vrec[idx][1] = ACB;					\
		  vrec[idx][2] = BAC;					\
		  vrec[idx][3] = BCA;					\
		  vrec[idx][4] = CAB;					\
		  vrec[idx][5] = CBA;					\
		}

real kk[32][32];
real dd[32][32];

#define PrintCounts(j)						   \
  printf(" PFct=%.4f%%\n",  kk[j][0]=  PFct*100.0/(ct+0.0000001)); \
  printf(" NMct=%.4f%%\n",  kk[j][1]=  NMct*100.0/(ct+0.0000001)); \
  printf("LMPct=%.4f%%\n",  kk[j][2]= LMPct*100.0/(ct+0.0000001)); \
  printf("MLPct=%.4f%%\n",  kk[j][3]= MLPct*100.0/(ct+0.0000001)); \
  printf("CWEct=%.4f%%\n",  kk[j][4]= CWEct*100.0/(ct+0.0000001)); \
  printf("REVct=%.4f%%\n",  kk[j][5]= REVct*100.0/(ct+0.0000001)); \
  printf("ALLct=%.4f%%\n",  kk[j][6]= ALLct*100.0/(ct+0.0000001)); \
  printf("ABct =%.4f%%\n",  kk[j][7]=  ABct*100.0/(ct+0.0000001)); \
  printf("GMct =%.4f%%\n",  kk[j][8]=  GMct*100.0/(ct+0.0000001)); \
  printf("PDct =%.4f%%\n",  kk[j][9]=  PDct*100.0/(ct+0.0000001)); \
  printf("CYct =%.4f%%\n",  kk[j][10]= CYct*100.0/(ct+0.0000001)); \
  printf("ANYct=%.4f%%\n",  kk[j][11]=ANYct*100.0/(ct+0.0000001)); \
  printf("BPAct=%.4f%%\n",  kk[j][12]=BPAct*100.0/(ct+0.0000001)); \
  printf("BNMct=%.4f%%\n",  kk[j][13]=BNMct*100.0/(ct+0.0000001)); \
  printf("FB2ct=%.4f%%\n",  kk[j][14]=FB2ct*100.0/(ct+0.0000001)); \
  printf("LDOct=%.4f%%\n",  kk[j][15]=LDOct*100.0/(ct+0.0000001)); 

#define PrintCountsD(j)						   \
  printf(" PFct=%.4f%%\n",  dd[j][0]=  PFctD*100.0/(PDct+0.0000001)); \
  printf(" NMct=%.4f%%\n",  dd[j][1]=  NMctD*100.0/(PDct+0.0000001)); \
  printf("LMPct=%.4f%%\n",  dd[j][2]= LMPctD*100.0/(PDct+0.0000001)); \
  printf("MLPct=%.4f%%\n",  dd[j][3]= MLPctD*100.0/(PDct+0.0000001)); \
  printf("CWEct=%.4f%%\n",  dd[j][4]= CWEctD*100.0/(PDct+0.0000001)); \
  printf("REVct=%.4f%%\n",  dd[j][5]= REVctD*100.0/(PDct+0.0000001)); \
  printf("ALLct=%.4f%%\n",  dd[j][6]= ALLctD*100.0/(PDct+0.0000001)); \
  printf("ABct =%.4f%%\n",  dd[j][7]=  ABctD*100.0/(PDct+0.0000001)); \
  printf("GMct =%.4f%%\n",  dd[j][8]=  GMctD*100.0/(PDct+0.0000001)); \
  printf("PDct =%.4f%%\n",  dd[j][9]=  PDctD*100.0/(PDct+0.0000001)); \
  printf("CYct =%.4f%%\n",  dd[j][10]= CYctD*100.0/(PDct+0.0000001)); \
  printf("ANYct=%.4f%%\n",  dd[j][11]=ANYctD*100.0/(PDct+0.0000001)); \
  printf("BPAct=%.4f%%\n",  dd[j][12]=BPActD*100.0/(PDct+0.0000001)); \
  printf("BNMct=%.4f%%\n",  dd[j][13]=BNMctD*100.0/(PDct+0.0000001)); \
  printf("FB2ct=%.4f%%\n",  dd[j][14]=FB2ctD*100.0/(PDct+0.0000001)); \
  printf("LDOct=%.4f%%\n",  dd[j][15]=LDOctD*100.0/(PDct+0.0000001)); 

FindGoodDenom(real xx[32][32], int j){
  int i,q,r;
  real ns,s,f,recerr,besterr,bestns;
  besterr = 99.9;
  bestns = 99999.9;
  for(q=3; q<990000; q++){
    recerr = 0.0;
    ns = 0.0;
    for(i=0; i<16; i++){
      s = xx[j][i]*q*0.01;
      r = round(s);
      f = s-r;
      ns += f*f;
      if(f<0.0) f = -f;
      assert(f>=0.0);
      assert(f<0.5);
      if(f>recerr){
	recerr=f;
      }
    }
    if(recerr<besterr){
      besterr=recerr; 
      printf("denom=%d  besterr=%f\n", q, besterr);
    }
    if(ns<bestns){
      bestns=ns;
      printf("denom=%d  bestns=%f\n", q, bestns);
    }
  }
}

char * str[] = {
  "<a href=\"PuzzIrvPartic.html\">Participation</a> failure W&cup;X",
  "<a href=\"Monotone.html\">Nonmonotonicity</a> U&cup;V",
  "V: (\"less is more\" nonmonotonicity)",
  "U: (\"more is less\" nonmonotonicity)",
  "Y: Condorcet winner eliminated (\"thwarted majority\")",
  "Z: <a href=\"IrvRFfreq.html\">Reversal</a> \"winner=loser\" failure",
  "R: All scoring rules agree B wins, but IRV says A wins (failure of \"sniff test\")",
  "W: Abstention failure: deleting A-bottom voters <i>stops</i> A from winning",
  "X: Would be strategic mistake for more voters of some single type to come",
  "T: Plurality and IRV winners differ",
  "S: Condorcet cycle",
  "Q&cup;R&cup;U&cup;V&cup;W&cup;X&cup;Y&cup;Z (<b>\"total paradox probability\"</b>)",
  "Both kinds of participation failure simultaneously W&cap;X",
  "Both kinds of nonmonotonicity simultaneously U&cap;V",
  "Q&cup;V: Betraying B makes either B or C win (where either way the betrayers prefer that to A winning)",
  "Q: Loser drop-out paradox: If B drops out, that switches the winner from A to C. Also (which happens in exactly the same set of elections) \"Favorite betrayal\"; voters with favorite B, by betraying B, make C win (whom they prefer as the \"lesser evil\" over current winner A)"
 };

#define PrintHTMLCounts(xx, color)				\
  printf("<table bgcolor=\"%s\" cellspacing=\"2\">\n", color);		\
printf("<tr bgcolor=\"pink\"><th>Phenomenon</th><th>REM</th><th>Dirichlet</th><th>Quas 1D</th></tr>\n"); \
for(j=0; j<16; j++){ \
  printf( \
"<tr><td>%s</td><td align=\"right\">%.4f%%</td><td align=\"right\">%.4f%%</td><td align=\"right\">%.4f%%</td></tr>\n", \
  str[j], xx[1][j], xx[2][j], xx[3][j] );					\
} printf("</table>\n"); \

int IntegerFind(int MX, int BS){
  int   ABC,  ACB,  BAC,  BCA,  CAB,  CBA,  V, tol;
  uint64 ct=0, ALLct=0, ANYct=0, NMct=0, PFct=0, REVct=0, PPPct=0, NPPct=0, LMPct=0, MLPct=0;
  uint64 PAPct=0, NAPct=0, CWEct=0, ABct=0, GMct=0, PDct=0, CYct=0, BPAct=0, BNMct=0;
  uint64 FBCct=0, FB2ct=0, LDOct=0, dm=0;
  uint64 ctD=0, ALLctD=0, ANYctD=0, NMctD=0, PFctD=0, REVctD=0, PPPctD=0, NPPctD=0, LMPctD=0;
  uint64 MLPctD=0, PAPctD=0, NAPctD=0, CWEctD=0, ABctD=0, GMctD=0, PDctD=0, CYctD=0;
  uint64 BPActD=0, BNMctD=0, FBCctD=0, FB2ctD=0, LDOctD=0;
  int oldidx, idx, nz, qu, j;
  //Let s+t=1 and let s=0 for IRV, s=1/2 for Nanson-Baldwin, and s=1 for Coombs.
#define s 0
#define t 1
  ct=0;
  for(j=0; j<4096; j++){ qrec[j] = 99999999; ctcn[j]=(uint64)0; }
  printf("s=%g, t=%g\n", (real)s,(real)t);
  printf("Examining all 3-candidate elections with <=%d votes of each type and <=%d voters.\n", 
	 MX,BS);
  for(ABC=0; ABC<MX; ABC++){
    for(ACB=0; ACB<MX; ACB++){
      for(BAC=0; BAC<MX; BAC++){
	for(BCA=0; BCA<MX; BCA++){
	  for(CAB=0; CAB<MX; CAB++){
	    for(CBA=0; CBA<MX; CBA++){
	      V = ABC+ACB+BAC+BCA+CAB+CBA;
	      if(V<=BS && STD){
		nz = (ABC!=0)+(ACB!=0)+(BAC!=0)+(BCA!=0)+(CAB!=0)+(CBA!=0);
		ct++;
		tol=4;
		Classify();
		oldidx=idx;
		tol= -1;
		Classify();
		if(oldidx==idx){
		  UpdateRecs();
		  if(idx==509 && nz==3){
		    if(BCA-ABC-1+CAB>ABC && BCA-CAB+1+ABC<CAB+CAB-1){
		      printf("V=%d: %d %d %d %d %d %d\n", V,ABC,ACB,BAC, BCA, CAB, CBA);
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
#undef t
#undef s
  PrintCounts(0);
  printf("%lld elections in standard form found.\n", ct);
  for(j=0; j<4096; j++){
    if(ctcn[j]>0){
      printf(
"idx=%3d; #voters=%2d, #types=%d; ABC=%2d, ACB=%2d, BAC=%2d, BCA=%2d, CAB=%2d, CBA=%2d; ct=%2.3f%%=%lld\n", 
	     j, qrec[j]/6, qrec[j]%6, 
	     vrec[j][0],vrec[j][1],vrec[j][2],vrec[j][3],vrec[j][4],vrec[j][5], 
	     (ctcn[j]*100.0)/(ct+0.000000001), ctcn[j] );
    }
  }
  fflush(stdout);
  return(ct);
}

void PrintHTMLtable(uint64 ct){
  uint j;
  printf("<table cellspacing=\"2\" bgcolor=\"pink\">\n");
 printf("<tr bgcolor=\"aqua\"><td><tt>QRSTUVWXYZ</tt></td><th>Election Example</th><th>REM Prob.</th><th>Dirichlet Prob.</th><th>Quas-1D Prob.</th></tr>\n");
  for(j=0; j<4096; j++){
    if(ctcn[j]>0  or ctNorm[j]>0 or ctDirich[j]>0 or ctQuas1D[j]>0){
      printf("<tr><td>");
      printf("<tt>");
      PrintBinary(j);
      printf("</tt>");
      printf("</td><td align=\"center\">ABC=%2d, ACB=%2d, BAC=%2d, BCA=%2d, CAB=%2d, CBA=%2d",
	     vrec[j][0],vrec[j][1],vrec[j][2],vrec[j][3],vrec[j][4],vrec[j][5] );
      printf("</td>");
      printf("<td align=\"right\">%.4f%%</td>", ctNorm[j]*100.0/(ct+0.00000001));
      printf("<td align=\"right\">%.4f%%</td>", ctDirich[j]*100.0/(ct+0.00000001));
      printf("<td align=\"right\">%.4f%%</td>", ctQuas1D[j]*100.0/(ct+0.00000001));
      printf("</tr>\n");
    }
  }
  printf("</table>");
}

int RandElect(uint64 tries){
  real   ABC,  ACB,  BAC,  BCA,  CAB,  CBA,  V, tol, zz;
  uint64 ct=0, ALLct=0, ANYct=0, NMct=0, PFct=0, REVct=0, PPPct=0, NPPct=0, LMPct=0, MLPct=0;
  uint64 PAPct=0, NAPct=0, CWEct=0, ABct=0, GMct=0, PDct=0, CYct=0,  BPAct=0, BNMct=0;
  uint64 FBCct=0, FB2ct=0, LDOct=0, dm=0;
  uint64 ctD=0, ALLctD=0, ANYctD=0, NMctD=0, PFctD=0, REVctD=0, PPPctD=0, NPPctD=0, LMPctD=0;
  uint64 MLPctD=0, PAPctD=0, NAPctD=0, CWEctD=0, ABctD=0, GMctD=0, PDctD=0, CYctD=0;
  uint64 BPActD=0, BNMctD=0, FBCctD=0, FB2ctD=0, LDOctD=0;
  int idx, nz, qu, j;
  //Let s+t=1 and let s=0 for IRV, s=1/2 for Nanson-Baldwin, and s=1 for Coombs.
#define s 0
#define t 1
  ct=0;
  for(j=0; j<4096; j++){ ctNorm[j]=0; }
  printf("s=%g, t=%g\n", (real)s,(real)t);
  printf("trying %llu Random Normal Elections\n", tries);
  ct=0;
  do{
    do{
      ABC = 24.0+RandNormal();
      ACB = 24.0+RandNormal();
      BAC = 24.0+RandNormal();
      BCA = 24.0+RandNormal();
      CAB = 24.0+RandNormal();
      CBA = 24.0+RandNormal();
      if(ABC+ACB<CAB+CBA and ABC+ACB<BAC+BCA){ //swap A,C; this is unnec but here to speed things up
	zz = ABC; ABC = CBA; CBA = zz;
	zz = ACB; ACB = CAB; CAB = zz;
      }
      else if(BAC+BCA<CAB+CBA and BAC+BCA<ABC+ACB){ //swap B,C; this is unnec but here to speed things up
	zz = BAC; BAC = CBA; CBA = zz;
	zz = BCA; BCA = CAB; CAB = zz;
      }
    }until(STD);
    V = ABC+ACB+BAC+BCA+CAB+CBA;
    ct++;
    tol=0.0;
    Classify();
    ctNorm[idx]++;
  }until(ct>=tries);
#undef t
#undef s
  printf("%llu random-normal elections in standard form found:\n", ct);
  PrintCounts(1);
  PrintCountsD(1);
  fflush(stdout);
  return(ct);
}

int DirichletElect(uint64 tries){
  real   ABC,  ACB,  BAC,  BCA,  CAB,  CBA,  V, tol, zz, y[5];
  uint64 ct=0, ALLct=0, ANYct=0, NMct=0, PFct=0, REVct=0, PPPct=0, NPPct=0, LMPct=0, MLPct=0;
  uint64 PAPct=0, NAPct=0, CWEct=0, ABct=0, GMct=0, PDct=0, CYct=0,  BPAct=0, BNMct=0;
  uint64 FBCct=0, FB2ct=0, LDOct=0, dm=0;
  uint64 ctD=0, ALLctD=0, ANYctD=0, NMctD=0, PFctD=0, REVctD=0, PPPctD=0, NPPctD=0, LMPctD=0;
  uint64 MLPctD=0, PAPctD=0, NAPctD=0, CWEctD=0, ABctD=0, GMctD=0, PDctD=0, CYctD=0;
  uint64 BPActD=0, BNMctD=0, FBCctD=0, FB2ctD=0, LDOctD=0;
  int idx, nz, qu, j;
  //Let s+t=1 and let s=0 for IRV, s=1/2 for Nanson-Baldwin, and s=1 for Coombs.
#define s 0
#define t 1
  ct=0;
  for(j=0; j<4096; j++){ ctDirich[j]=0; }
  printf("s=%g, t=%g\n", (real)s,(real)t);
  printf("trying %llu Random Dirichlet Elections\n", tries);
  ct=0;
  V = 1.0;
  do{
    do{
      Gimme5sortedRands(y);
      ABC = y[0];
      ACB = y[1]-y[0];
      BAC = y[2]-y[1];
      BCA = y[3]-y[2];
      CAB = y[4]-y[3];
      CBA = 1.0-y[4];
      if(ABC+ACB<CAB+CBA and ABC+ACB<BAC+BCA){ //swap A,C; this is unnec but here to speed things up
	zz = ABC; ABC = CBA; CBA = zz;
	zz = ACB; ACB = CAB; CAB = zz;
      }
      else if(BAC+BCA<CAB+CBA and BAC+BCA<ABC+ACB){ //swap B,C; this is unnec but here to speed things up
	zz = BAC; BAC = CBA; CBA = zz;
	zz = BCA; BCA = CAB; CAB = zz;
      }
    }until(STD);
    ct++;
    tol=0.0;
    Classify();
    ctDirich[idx]++;
  }until(ct>=tries);
#undef t
#undef s
  printf("%llu random-Dirichlet elections in standard form found:\n", ct);
  PrintCounts(2);
  PrintCountsD(2);
  fflush(stdout);
  return(ct);
}


int Quas1DElect(uint64 tries){
  real   ABC,  ACB,  BAC,  BCA,  CAB,  CBA,  V, tol, y[3], midex, leftmid, rightmid, zz;
  real   ABCt, ACBt, BACt, BCAt, CABt, CBAt;
  uint64 ct=0, ALLct=0, ANYct=0, NMct=0, PFct=0, REVct=0, PPPct=0, NPPct=0, LMPct=0, MLPct=0;
  uint64 PAPct=0, NAPct=0, CWEct=0, ABct=0, GMct=0, PDct=0, CYct=0,  BPAct=0, BNMct=0;
  uint64 FBCct=0, FB2ct=0, LDOct=0, dm=0;
  uint64 ctD=0, ALLctD=0, ANYctD=0, NMctD=0, PFctD=0, REVctD=0, PPPctD=0, NPPctD=0, LMPctD=0;
  uint64 MLPctD=0, PAPctD=0, NAPctD=0, CWEctD=0, ABctD=0, GMctD=0, PDctD=0, CYctD=0;
  uint64 BPActD=0, BNMctD=0, FBCctD=0, FB2ctD=0, LDOctD=0;
  int idx, nz, qu, j;
  //Let s+t=1 and let s=0 for IRV, s=1/2 for Nanson-Baldwin, and s=1 for Coombs.
#define s 0
#define t 1
  ct=0;
  for(j=0; j<4096; j++){ ctQuas1D[j]=0; }
  printf("s=%g, t=%g\n", (real)s,(real)t);
  printf("trying %llu Random Quas1D Elections\n", tries);
  ct=0;
  V = 1.0;
  do{
      Gimme3sortedRands(y);
      // reflect if nec to make the y[0] candidate have more plurality votes than y[2]:
      if( y[0]+2*y[1]+y[2]<2.0 ){  zz = y[2]; y[2] = 1.0-y[0]; y[0] = 1.0-zz; y[1] = 1.0-y[1];  }
      leftmid = 0.5*(y[0]+y[1]);
      ABC = leftmid;
      assert(ABC >= 0.0);
      rightmid = 0.5*(y[2]+y[1]);
      CBA = 1.0-rightmid;
      assert(CBA >= 0.0);
      midex = 0.5*(y[0]+y[2]); //midpoint of extremes
      assert( leftmid <= midex );
      assert( midex <= rightmid );
      BAC = midex - leftmid;
      BCA = rightmid - midex; 
      ACB = 0.0;
      CAB = 0.0;
      assert(BAC >= 0.0 );
      assert(BCA >= 0.0 );
      assert(ABC+ACB >= CBA+CAB);
      if( CBA+CAB < BAC+BCA ){ //C eliminated first
	if( ABC+ACB+CAB > BAC+BCA+CBA ){ // A wins, finish order ABC
	  ABCt=ABC; ACBt=ACB; BACt=BAC; BCAt=BCA; CABt=CAB; CBAt=CBA;
	}else{ // B wins, finish order BAC
	  BACt=ABC; BCAt=ACB; ABCt=BAC; ACBt=BCA; CBAt=CAB; CABt=CBA;
	}
      }else{ //B eliminated first
	if( ABC+ACB+BAC > CBA+CAB+BCA ){ // A wins, finish order ACB
	  ACBt=ABC; ABCt=ACB; CABt=BAC; CBAt=BCA; BACt=CAB; BCAt=CBA;
	}else{ // C wins, finish order CAB
	  BCAt=ABC; BACt=ACB; CBAt=BAC; CABt=BCA; ABCt=CAB; ACBt=CBA;
	}
      }
      /*******old standardization code: buggy?
    //find IRV winner & 2nd
    if( ABC+ACB < BAC+BCA && ABC+ACB<CAB+CBA ){ // A elimd
      if( ABC+BAC+BCA > ACB+CAB+CBA ){ // B wins, BCA (so repl b-->A, c-->B, a-->C):
	CABt=ABC; CBAt=ACB; ACBt=BAC; ABCt=BCA; BCAt=CAB; BACt=CBA;
      }else{ // C wins, CBA
	CBAt=ABC; CABt=ACB; BCAt=BAC; BACt=BCA; ACBt=CAB; ABCt=CBA;
      }
    }
    else if( BAC+BCA < ABC+ACB && BAC+BCA<CAB+CBA ){ // B elimd
      if( ABC+ACB+BAC > CBA+CAB+BCA ){ // A wins, ACB
	ACBt=ABC; ABCt=ACB; CABt=BAC; CBAt=BCA; BACt=CAB; BCAt=CBA;
      }else{ // C wins, CAB
        BCAt=ABC; BACt=ACB; CBAt=BAC; CABt=BCA; ABCt=CAB; ACBt=CBA;
      }
    }
    else if( CAB+CBA < ABC+ACB && CAB+CBA<BCA+BCA ){ // C elimd
      if( ABC+ACB+CAB > BAC+BCA+CBA ){ // A wins, ABC
	ABCt=ABC; ACBt=ACB; BACt=BAC; BCAt=BCA; CABt=CAB; CBAt=CBA;
      }else{ // B wins, BAC
	BACt=ABC; BCAt=ACB; ABCt=BAC; ACBt=BCA; CBAt=CAB; CABt=CBA;
      }
    }
    *****************/
    ABC=ABCt; ACB=ACBt; BAC=BACt; BCA=BCAt; CAB=CABt; CBA=CBAt;
    assert( V <  0.0001 + ABC+ACB+BAC+BCA+CAB+CBA ); 
    assert( V > -0.0001 + ABC+ACB+BAC+BCA+CAB+CBA ); 
    assert(STD);
    assert( ((ABC!=0)+(ACB!=0)+(BAC!=0)+(BCA!=0)+(CAB!=0)+(CBA!=0)) <= 4 );
    if(MLP != MLPq){ printf("Qyikes!\n"); }			
    ct++;
    tol=0.0;
    Classify();
    ctQuas1D[idx]++;
  }until(ct>=tries);
#undef t
#undef s
  printf("%llu random-Quas1D elections in standard form found:\n", ct);
  PrintCounts(3);
  PrintCountsD(3);
  fflush(stdout);
  return(ct);
}

main(){
  int i,j;
  for(i=0;i<32;i++){  for(j=0;j<32;j++){ kk[i][j]=dd[i][j]=0; }}
  InitRand(0);
  system("date"); fflush(stdout);
  IntegerFind(25,45);
IntegerFind(45,117); //under 4 minutes
  system("date");
  //  #define HOWMANY 4000000000  // 4*10^9 ==> under 20000 sec ??
#define HOWMANY 40000000000ULL  // 4*10^10 ==> ??
  //#define HOWMANY 10000000  // 4*10^9 ==> under 20000 sec
  //#define HOWMANY 100000000
  //#define HOWMANY 10000000000ULL
  Quas1DElect(   HOWMANY);
  system("date"); fflush(stdout);
  RandElect(     HOWMANY);
  system("date"); fflush(stdout);
  DirichletElect(HOWMANY);
  system("date"); fflush(stdout);
  PrintHTMLtable(HOWMANY);
#undef HOWMANY
  printf("<br><p>Below is a smaller <b>summary table</b> of some of the\n"
         "most-requested pathology-probability\n"
	 "information  (All of the numbers in the below tables are derivable by adding up\n"
         "appropriate sets of numbers from the above master table):</p><br>\n");
  PrintHTMLCounts(kk, "aqua");
  printf("<br><p>And below is the same table, but <i>restricted</i> to elections in which the IRV process\n"
	 "<i>matters</i>, i.e. in which the IRV and plain-plurality winners differ.\n"
	 "(Warning: The error bars are approximately twice as wide as in the tables above.)\n"
"This almost always makes pathologies substantially more likely:</p><br>\n");
  PrintHTMLCounts(dd, "yellow");
  printf("Good denoms for kk blue Dirichlet table:\n");
  FindGoodDenom(kk, 2);
  printf("Good denoms for kk blue Quas table:\n");
  FindGoodDenom(kk, 3);
  printf("Good denoms for dd yellow Dirichlet table:\n");
  FindGoodDenom(dd, 2);
  printf("Good denoms for dd yellow Quas table:\n");
  FindGoodDenom(dd, 3);
  printf("\nall done.\n");
}

