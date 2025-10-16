#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>

/* 
OK, here is Warren D Smith's simulation hack for the Shentrup-Raphfrk voting experiment.
C, runs under LINUX.  May need a few modifs to run under other OSs.
Set MSWINDOWS to TRUE and that'll go some way toward making it run under MSWINDOWS
but not the whole way.

Compile with
   gcc -Wall -O6 shenhack.c  -DNDEBUG -o shenhack 
Run with
   shenhack > outputfilename

On my machine, with NUMEL 150000 it takes about 8 hours to run.
I recommend 1500000 - 9999999  which would take several days to a few weeks.
For testing purposes to get very high speed and just see if it works, use NUMEL 100

Three modes:
Raphfrk's experiment: 
RAPHFRKMODE TRUE
GENUINEUTIL FALSE

Shentrup's experiment: 
GENUINEUTIL FALSE
RAPHFRKMODE FALSE

Just comparing societal util of Approval vs Range:
GENUINEUTIL TRUE
RAPHFRKMODE FALSE
******************************/

#define TRUE (1==1)
#define FALSE (0==1)

#define MSWINDOWS FALSE
#define NUMEL 150000  /*large for low noise, small for high speed*/
#define STEPME 5
#define MINNUMCANDS 2
#define MAXNUMCANDS 10
#define RAPHFRKMODE TRUE
#define GENUINEUTIL FALSE
#define uint32 uint
#define bool uint
#define uint64 unsigned long long
#define real double

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

void testbiglincong(){
   int i;
   /* lexically minimal permissible seed: */
   for(i=0; i<60; i++){ BLC32x[i]=0; }
   BLC32x[0] = 1; 
   BLC32NumLeft = 0;

   for(i=0; i<599; i++){ BigLinCong32(); }
   printf("%u %u %u %u %u\n", 
       BigLinCong32(),BigLinCong32(),
       BigLinCong32(),BigLinCong32(),BigLinCong32());

   for(i=0; i<12; i++){ 
      printf("%8x %8x %8x %8x %8x\n",
       BigLinCong32(),BigLinCong32(),
       BigLinCong32(),BigLinCong32(),BigLinCong32());
   }
}

/********************************
MAPLE script to check it works:

w := 32;
A := 1284507170 * 2^(w*3) + 847441413 * 2^(w*44) + 650134147 * 2^(w*59);
P := (2^(48*32) - 1) * 2^(12*32) + 1;
for i from 1 to 11 do
   print( floor((A &^ i mod P) / 2^(44*w)) mod (2^w) );
od;
ap := A &^ 11 mod P;
ap mod (2^w);
quit;
********************
MAPLE output:
847441413
4038410930
102374915
470100141
3896743552
243412576
1911259311
1640083353
4014446395
2679952782
4087228849
and
2475219032

C output:
847441413
4038410930
102374915
470100141
3896743552
243412576
oops
2990118053
2614294516
3539391572
1589778147
1758817216
2847725135 1364008005 3563722108 2020382641 1091616930
*************************/

uint32 RandUint(){ /* returns random uint32 */
  return BigLinCong32();
}

real Rand01(){ /* returns random uniform in [0,1] */
  return ((BigLinCong32()+0.5)/(1.0 + MAXUINT) + BigLinCong32())/(1.0 + MAXUINT);
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
   for(i=0; i<60; i++){ BLC32x[i]=0; }
   BLC32x[0] = seed; 
   BLC32NumLeft = 0;
   for(i=0; i<599; i++){ BigLinCong32(); }
   printf("Random generator initialized with seed=%u:\n", seed);
   for(i=0; i<7; i++){ 
     printf("%.6f ", Rand01());
   }
   printf("\n");
}

real SquareReal(real x){ return x*x; }

void TestRand01(){ 
  int i,y, ct[10];
   real x,s,mx,mn,v;
   s=0.0; v=0.0;
   mn = HUGE;
   mx = -HUGE;
   for(i=0; i<10; i++) ct[i]=0; 
   printf("Performing 100000 randgen calls to test that randgen[0,1] behaving ok:\n");
   for(i=0; i<100000; i++){ 
     x = Rand01();
     s += x;
     if(mx<x) mx=x;
     if(x<mn) mn=x;
     v += SquareReal(x-0.5);
     y = (int)(x*10.0);
     if(x>=0 && y<10) ct[y]++;
   }
   printf("mean=%g(should be 1/2)  min=%f  max=%g   variance=%g(should be 1/12=%g)\n", 
      s/100000.0, mn, mx, v/100000.0, 1/12.0);
   printf("counts in 10 bins 0-0.1, 0.1-0.2, etc: ");
   for(i=0; i<10; i++) printf(" %d", ct[i]);
   printf("\n");
   fflush(stdout);
}

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

void TestRandExpl(){ 
  int i, y, ct[10];
  real x,s,mx,mn,v;
  s=0.0; v=0.0;
  mn = HUGE;
  mx = -HUGE;
  for(i=0; i<10; i++) ct[i]=0; 
  printf("Performing 100000 randgen calls to test that expl randgen behaving ok:\n");
  for(i=0; i<100000; i++){ 
     x = RandExpl();
     s += x;
     if(mx<x) mx=x;
     if(x<mn) mn=x;
     v += x*x;
     y = (int)(x*10);
     if(x>=0 && y<10) ct[y]++;
  }
  printf("mean=%g(should be 1)  min=%f  max=%g   variance=%g(should be 2?)\n", 
	 s/100000.0, mn, mx, v/100000.0);
  printf("counts in 10 bins 0-0.1, 0.1-0.2, etc: ");
  for(i=0; i<10; i++) printf(" %d", ct[i]);
  printf("\n");
  fflush(stdout);
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

void TestRandNormal(){ 
  int i, y, ct[10];
  real x,s,mx,mn,v;
  s=0.0; v=0.0;
  mn = HUGE;
  mx = -HUGE;
  for(i=0; i<10; i++) ct[i]=0; 
  printf("Performing 100000 randgen calls to test that normal randgen behaving ok:\n");
  for(i=0; i<100000; i++){ 
     x = RandNormal();
     s += x;
     if(mx<x) mx=x;
     if(x<mn) mn=x;
     v += x*x;
     y = (int)(x*10);
     if(x>=0 && y<10) ct[y]++;
  }
  printf("mean=%g(should be 0)  min=%f  max=%g   variance=%g(should be 1)\n", 
	 s/100000.0, mn, mx, v/100000.0);
  printf("counts in 10 bins 0-0.1, 0.1-0.2, etc: ");
  for(i=0; i<10; i++) printf(" %d", ct[i]);
  printf("\n");
  fflush(stdout);
}


main(){
  int t,i,j,k,ALLAPP,winner,NUMCANDS;
  real mx, mn, av, sum, UtilForApprovalers, UtilForHonests, x, diff, old;
  real Tot[MAXNUMCANDS+1], U[101][MAXNUMCANDS+1], V[101][MAXNUMCANDS+1];

  printf("Approval vs Range Random Electns Model\n");

  InitRand(0);
  TestRandNormal();

#if RAPHFRKMODE
  printf("<p>\nBelow table shows difference in average utility (AvgUtil for strategists minus AvgUtil for honest voters)\n");
  printf("when the first <i>t</i> voters strategize (remaining <i>100-t</i> voters honest),\n");
#else
  printf("<p>\nBelow table shows summed utility difference (util for approval minus util for range)\n");
#if GENUINEUTIL 
  printf("summed for <i>all</i> the voters when the first <i>t</i> voters strategize (rest honest),\n");
#else
  printf("for the <i>100-t</i> honest normalized-range-voters <i>only</i> (the other <i>t</i> are strategic),\n");
#endif
#endif
  printf("for 100-voter elections with %d to %d candidates.\n", MINNUMCANDS, MAXNUMCANDS);
  printf("Each data point based on %d elections.\n", NUMEL);
  printf("The \"strategic\" voters are using mean-candidate-utility as an approval-theshold.\n");
  printf("</p>\n");
  printf("<table cellspacing=\"2\" cellpadding=\"2\">\n");
  printf("<tr>");
  printf("<th>%%strat</th>");
  for(k=MINNUMCANDS; k<=MAXNUMCANDS; k++){ printf("<th>%d</th>", k); }
  printf("</tr>");
  printf("\n");


  for(t=0; t<=100; t+=STEPME){/*t = # of approval-type (i.e. strategic) voters, out of 100 in all*/
    for(NUMCANDS=MINNUMCANDS; NUMCANDS<=MAXNUMCANDS; NUMCANDS++){
#if RAPHFRKMODE
    ALLAPP=0;
#else
    for(ALLAPP=0; ALLAPP<=1; ALLAPP++){ /*1 if ALL voters approval-type; 0 if only the strat ones*/
#endif
      UtilForApprovalers = 0.0;
      UtilForHonests = 0.0;
      for(i=0; i<NUMEL; i++){  /*averaging the results of NUMEL elections*/
	for(k=0; k<NUMCANDS; k++){ Tot[k] = 0.0; } /*start election with zero vote-totals*/
	for(j=0; j<100; j++){  /* there are 100 voters in all */
	  mx = -9999999.9;
	  mn =  9999999.9;
	  av = 0.0;
	  for(k=0; k<NUMCANDS; k++){ /* Generate utility of canddt k for voter j */
	    x = RandNormal();
	    assert(0 <= j);	    assert(j < 100); 
	    assert(0 <= k);	    assert(k < NUMCANDS); 
	    U[j][k] = x;
	    if(x > mx)  mx=x;
	    if(x < mn)  mn=x;
	    av += x;
	  }
	  av /= NUMCANDS;
	  assert(mn<av);	  assert(av<mx);
	  diff = mx-mn;
	  for(k=0; k<NUMCANDS; k++){  /*cast vote j; voters j with j<t are strategic*/
	    assert(0 <= j);	    assert(j < 100); 
	    assert(0 <= k);	    assert(k < NUMCANDS); 
	    x = U[j][k];
	    if(j<t || ALLAPP==1) V[j][k] = ((x<av)?0.:1.);   /*approval vote*/
	    else                 V[j][k] = (x-mn)/diff;  /*normalized range vote*/
	    Tot[k] += V[j][k];
	  }
	}/*end for(j) */
	mx =  -9999999.9;
	winner = -1;
	for(k=0; k<NUMCANDS; k++){  /*find winner*/
	  assert(0 <= k);	    assert(k < NUMCANDS); 
	  if(Tot[k] > mx){ mx = Tot[k]; winner=k; }
	}
	assert(winner >= 0);	assert(winner < NUMCANDS);
	sum = 0.0;
#if GENUINEUTIL 
	for(j=0; j<100; j++){
#else
	for(j=t; j<100; j++){  /* util-sum for honest voters only */
#endif
	  assert(0 <= j);	    assert(j < 100); 
	  sum += U[j][winner]; 
	}
	UtilForHonests += sum;
#if RAPHFRKMODE
	sum = 0.0;
	for(j=0; j<t; j++){
	  assert(0 <= j);	    assert(j < 100); 
	  sum += U[j][winner]; 
	}
	UtilForApprovalers += sum;
#endif
      }/* end for(i) */
#if RAPHFRKMODE
      UtilForHonests /= NUMEL*(100.001-t);  
      UtilForApprovalers /= NUMEL*(t+0.001);  
      if(1){
	if(NUMCANDS==MINNUMCANDS) printf("<tr><td> %d%% </td>", t);
	printf("<td> %.3f </td>\n", UtilForApprovalers-UtilForHonests);
#else
      UtilForHonests /= NUMEL;
      /*      printf("Strat%%=%d    UtilForHonests=%.3f   ALLAPP=%d\n",
      **      t, UtilForHonests, ALLAPP); */
      if(ALLAPP==0)  old = UtilForHonests;
      else{
	if(NUMCANDS==MINNUMCANDS) printf("<tr><td> %d%% </td>", t);
        printf("<td> %.3f </td>\n", UtilForHonests-old);
#endif
        if(NUMCANDS==MAXNUMCANDS){
	  printf("</tr>\n" );
	  fflush(stdout); 
	}
      }
#if !RAPHFRKMODE
    }/* end for(ALLAPP) */
#endif
    }/* end for(NUMCANDS) */
  }/* end for(t) */

  printf("</table>\n");
}



