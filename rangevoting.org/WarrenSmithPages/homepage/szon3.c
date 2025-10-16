
/****************************************************************
* szon3.c By Warren D. Smith July 2002.
* Contains C routines for
* systematic investigation of complexes, quaternions, octonions,
* seizeons/sedenions/hexons, etc.  Compile with
*    gcc -O2 -malign-double -Wall -lm -lc szon3.c -o szon3
* or
*    gcc -Wall -lm -lc szon3.c -o szon3
* or
*    icc -lm -lc -O2 szon3.c -o szon3
********************************************************
* I am making this code publically available.  It may be used
* by anybody for any purpose provided acknowledgment
* is made to me as the original author and modifications you make
* are noted and any errors are pointed out to me.
********************************************************
* No claim is made that code is very efficient.
* C is not a good language for this task.
* LISP would be superior. ML & Ada95 too, maybe Prolog.
* MAPLE is quite nice but much too slow.
********************************************************
* the 2^n-ons ("eon"s) are represented as
* MAXSIZE-element arrays x[0..MAXSIZE-1]
* of 64-bit-signed-integers or 64-bit-reals. The real part is x[0].
* The first half is the 2^{n-1}-on.
* The "level" (value of n=0,1,2,3,4,5,6) is a uint8.
*****************************************************
* To do: 
* other idents: det2, det3, trace2, trace3,
* try finite fields and finite loops
* program division
* seek nonunique division
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
int rand(void);
void srand(unsigned int seed);
double sqrt(double);
double sin(double);
double cos(double);
double log(double);

double myrand(){ return ((double)rand()) / (double)RAND_MAX; }
double normrand(){
  double u1,u2,v1,s;
  static double v2;
  static uint ct;
  if(ct){ ct=0; return v2; }
  do{
    u1 = myrand();  u2 = myrand();  
    v1 = 2*u1-1;    v2 = 2*u2-1;
    s = v1*v1+v2*v2;
  }while( s>=1 );
  if(s > 0) s = sqrt(-2*log(s)/s);
  v1 *= s; v2 *= s;
  ct=1;
  return v1;
}

#define FALSE (1==0)
#define TRUE  (!FALSE)

#define PI     3.14159265358979323844
#define RT3BY2 0.86602540378443864676

typedef long long int64;
typedef unsigned long long uint64;
typedef double real64;
typedef unsigned char uint8;
/*typedef unsigned int uint;*/
typedef uint bool;

#define LGMAXSIZE 6
#define MAXSIZE (1<<LGMAXSIZE)

#define CONTIN TRUE
#define NONZONLY TRUE
#define HEAVYDATA16 TRUE
#define MODULUS 0
#define SLOWSTUFF FALSE
#define TROTTERTEST FALSE

uint RECCUTOFF = 3;
uint MULCHOICE;
uint WDSMULCHOICE;
uint COMBOMULCHOICE;
uint SWIND = 0;
uint TWOFER = 0;

#define CDcode    (0)
#define WDScode   (1)
#define CON1code  (2)
#define CON2code  (3)
#define CON3code  (13*2)
#define CON4code  (13*4)
#define CON5code  (13*13*1)
#define CON6code  (13*13*4)
#define CON7code  (13*13*13*1)
#define CON8code  (13*13*13*3)
#define WEIRDcode (13*13*13*13*5 + 13*13*13*4 + 13*13*0 + 13*7 + 0)
#define PF1code   1000001
#define PF2code   1000002
#define PF3code   1000003
#define JDHScode  1000005
#define ZASScode  1000007
#define WDSeCode  1000008

#define NUMCODES 17
uint MagicCode[NUMCODES] = {
  CDcode, WDScode,
  CON1code, CON2code, CON3code, CON4code,
  CON5code, CON6code, CON7code, CON8code, WEIRDcode,
  PF1code, PF2code, PF3code, JDHScode, ZASScode,
  WDSeCode
};

char MagicName[][30] = {
  "Cayley-Dickson", "WDSmul",
  "Conway0002", "Conway0003", "Conway0020", "Conway0040",
  "Conway0100", "Conway0400", "Conway1000", "Conway3000",
  "WeirdMul", "Pfcomplex", "Pfquat", "PfnearOct",
  "JDHSmith", "Zassenhaus-Eichhorn"
};

char EonName[][30] = {
  "Real-Numbers", "Complex-Numbers", "Quaternions",
  "Octonions", "Hexons/Sedenions", "32ons"
};

#if CONTIN
#define basetype  real64
#define ubasetype real64
#else
#define basetype   int64
#define ubasetype uint64
#endif

typedef struct dummy {
  uint8 lev;                    /* 0,1,2,3,4,5 */
  basetype coord[MAXSIZE];
} eon;

#define CONJSWAP(a,b) tmp=Conj(a); a=Conj(b); b=tmp

#define TESTBOOL(n,ZOG) printf(#ZOG "(%d)=",(n)); printbool(ZOG(n)); printf("\n");

#define TESTREAL1(n,ZOG) printf(#ZOG "(a)=%f,%f,%f\n",\
  ZOG(EErand((n),pureimag,-9)),\
  ZOG(EErand((n),pureimag,-9)),\
  ZOG(EErand((n),pureimag,-9)) )
#define TESTREAL2(n,ZOG) printf(#ZOG "(a,b)=%f,%f,%f\n",\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9) ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9) ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9) ) )
#define TESTREAL3(n,ZOG) printf(#ZOG "(a,b,c)=%f,%f,%f\n",\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9), EErand((n),pureimag,-9)  ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9), EErand((n),pureimag,-9)  ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9), EErand((n),pureimag,-9)  ))
#define TESTREAL4(n,ZOG) printf(#ZOG "(a,b,c,d)=%f,%f,%f\n",\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9),\
     EErand((n),pureimag,-9), EErand((n),pureimag,-9) ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9),\
     EErand((n),pureimag,-9), EErand((n),pureimag,-9) ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9),\
     EErand((n),pureimag,-9), EErand((n),pureimag,-9) ))
#define TESTREAL5(ZOG) printf(#ZOG "(a,b,c,d,e)=%f,%f,%f\n",\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9),\
     EErand((n),pureimag,-9), EErand((n),pureimag,-9), EErand((n),pureimag,-9) ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9),\
     EErand((n),pureimag,-9), EErand((n),pureimag,-9), EErand((n),pureimag,-9) ),\
ZOG(EErand((n),pureimag,-9), EErand((n),pureimag,-9),\
     EErand((n),pureimag,-9), EErand((n),pureimag,-9), EErand((n),pureimag,-9) ))

#define TESTIDENT1(n,pureimag,ZOG){\
  a = EErand((n),(pureimag),-9); \
  t1 = ZOG(a);    \
  a = EErand((n),(pureimag),-9); \
  t2 = ZOG(a);    \
  t1 = mxv(t1,t2);  \
  a = EErand((n),(pureimag),-9); \
  t2 = ZOG(a);    \
  t1 = mxv(t1,t2);  \
   printf("%d&",n);\
  printf(#ZOG "(a)="); if(TWOFER==1) printf("T&"); if(TWOFER==2) printf("Z&");\
  printsgn(t1, 0.005, 0.0001);         \
}

#define TESTIDENT2(n,pureimag,ZOG){\
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  t1 = ZOG(a,b);    \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b);    \
  t1 = mxv(t1,t2);  \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b);    \
  t1 = mxv(t1,t2);  \
   printf("%d&",n);\
  printf(#ZOG "(a,b)="); if(TWOFER==1) printf("T&"); if(TWOFER==2) printf("Z&");\
  printsgn(t1, 0.005, 0.0001);         \
}

#define TESTIDENT3(n,pureimag,ZOG){\
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  t1 = ZOG(a,b,c);    \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c);    \
  t1 = mxv(t1,t2);  \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c);    \
  t1 = mxv(t1,t2);  \
   printf("%d&",n);\
  printf(#ZOG "(a,b,c)="); if(TWOFER==1) printf("T&"); if(TWOFER==2) printf("Z&");\
  printsgn(t1, 0.005, 0.0001);   \
}

#define TESTIDENT4(n,pureimag,ZOG){\
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  t1 = ZOG(a,b,c,d);    \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c,d);    \
  t1 = mxv(t1,t2);  \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c,d);    \
  t1 = mxv(t1,t2);  \
   printf("%d&",n);\
  printf(#ZOG "(a,b,c,d)="); if(TWOFER==1) printf("T&"); if(TWOFER==2) printf("Z&");\
  printsgn(t1, 0.005, 0.0001);         \
}

#define TESTIDENT5(n,pureimag,ZOG){\
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  e = EErand((n),(pureimag),-9); \
  t1 = ZOG(a,b,c,d,e);    \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  e = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c,d,e);    \
  t1 = mxv(t1,t2);  \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  e = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c,d,e);    \
  t1 = mxv(t1,t2);  \
   printf("%d&",n);\
  printf(#ZOG "(a,b,c,d,e)="); if(TWOFER==1) printf("T&"); if(TWOFER==2) printf("Z&");\
  printsgn(t1, 0.005, 0.0001);         \
}

#define TESTIDENT6(n,pureimag,ZOG){\
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  e = EErand((n),(pureimag),-9); \
  f = EErand((n),(pureimag),-9); \
  t1 = ZOG(a,b,c,d,e,f);    \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  e = EErand((n),(pureimag),-9); \
  f = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c,d,e,f);    \
  t1 = mxv(t1,t2);  \
  a = EErand((n),(pureimag),-9); \
  b = EErand((n),(pureimag),-9); \
  c = EErand((n),(pureimag),-9); \
  d = EErand((n),(pureimag),-9); \
  e = EErand((n),(pureimag),-9); \
  f = EErand((n),(pureimag),-9); \
  t2 = ZOG(a,b,c,d,e,f);    \
  t1 = mxv(t1,t2);  \
   printf("%d&",n);\
  printf(#ZOG "(a,b,c,d,e,f)="); if(TWOFER==1) printf("T&"); if(TWOFER==2) printf("Z&");\
  printsgn(t1, 0.005, 0.0001);   \
}

basetype BaseDivide(basetype x, basetype q)
{
#if CONTIN
  return (x / q);
#else
#if MODULUS==0
  return (x / q);
#elsif MODULUS==2
  /*1/0=inf 1/1=1 0/1=0 0/0=1[or0] 1/inf=0 0/inf=0 inf/1=inf inf/0=inf */
  assert(q != 0);
  return x;
#elsif MODULUS==3
  assert(q != 0);
  if(q == 2)
    return 3 - x;
  return x;
#else
#endif
#endif
}


void printmodes()
{
  if(CONTIN)
    printf("CONTIN");
  else
    printf("INT-ONLY");
  printf("\n");
  if(NONZONLY)
    printf("NONZONLY");
  else
    printf("FULLPRINT");
  printf("\n");
  printf("MODULUS=%d\n", MODULUS);
  printf
      ("RECCUTOFF=%d (3 means oct->hex with new formula, 0 means all the way)\n",
       RECCUTOFF);
  printf("rand=%d, %d\n", rand(), rand());
  printf("normrand=%f, %f\n", normrand(), normrand());
  printf("-------------------------------\n");
}

inline uint pow2(uint n)
{
  assert(n < MAXSIZE);
  return (1 << n);              /* 2^n */
}


void print5(uint Q)
{
  uint Q1, Q2, Q3, Q4;
  Q4 = Q % 5;
  Q /= 5;
  Q3 = Q % 5;
  Q /= 5;
  Q2 = Q % 5;
  Q /= 5;
  Q1 = Q % 5;
  Q /= 5;
  printf("%d%d%d%d", Q1, Q2, Q3, Q4);
}

void print9(uint Q)
{
  uint Q1, Q2, Q3, Q4;
  Q4 = Q % 9;
  Q /= 9;
  Q3 = Q % 9;
  Q /= 9;
  Q2 = Q % 9;
  Q /= 9;
  Q1 = Q % 9;
  Q /= 9;
  printf("%d%d%d%d", Q1, Q2, Q3, Q4);
}

void print13(uint Q)
{
  uint Q1, Q2, Q3, Q4;
  Q4 = Q % 13;
  Q /= 13;
  Q3 = Q % 13;
  Q /= 13;
  Q2 = Q % 13;
  Q /= 13;
  Q1 = Q % 13;
  Q /= 13;
  printf("%x%x%x%x", Q1, Q2, Q3, Q4);
}

uint biggest9dig(uint Q)
{
  uint x, y;
  y = 0;
  x = Q % 9;
  Q /= 9;
  if(x > y)
    y = x;
  x = Q % 9;
  Q /= 9;
  if(x > y)
    y = x;
  x = Q % 9;
  Q /= 9;
  if(x > y)
    y = x;
  x = Q % 9;
  Q /= 9;
  if(x > y)
    y = x;
  return y;
}

uint count9dig(uint Q)
{
  uint x, ct;
  ct = 0;
  x = Q % 9;
  Q /= 9;
  if(x > 0)
    ct++;
  x = Q % 9;
  Q /= 9;
  if(x > 0)
    ct++;
  x = Q % 9;
  Q /= 9;
  if(x > 0)
    ct++;
  x = Q % 9;
  Q /= 9;
  if(x > 0)
    ct++;
  return ct;
}

void printnum(basetype y)
{
  double x;
  x = y;
  printf("%.1f", (double) x);
}

void printbool(bool x)
{
  if(x)
    printf("TRUE");
  else
    printf("FALSE");
}

void printeon(eon x)
{
  uint i, j, n;
  n = x.lev;
  j = pow2(n);
  printf("(");
  for (i = 0; i < j; i++){
    printnum(x.coord[i]);
    if(i + 1 < j)
      printf(", ");
  }
  printf(")\n");
}

void printzeon(eon x, double Tlarge, double Tsmall)
{
  uint i, j, n;
  n = x.lev;
  j = pow2(n);
  printf("(");
#if NONZONLY
  for (i = 0; i < j; i++){
    if(fabs((double) (x.coord[i])) > Tlarge)
      printf("N");
    else if(fabs((double) (x.coord[i])) > Tsmall)
      printf("n");
    else
      printf("0");
  }
#else
  for (i = 0; i < j; i++){
    printnum(x.coord[i]);
    if(i + 1 < j)
      printf(", ");
  }
#endif
  printf(")\n");
}

void printTzeon(eon x)
{
  printzeon(x, 0.00001, 0.0000003);
}

/*buggy? last0 seems to foul up*/
void printsgn(eon x, double T, double Tsmall)
{
  uint i, j, n;
  int t;
  bool virgin;
  virgin = TRUE;
  t = -1;
  n = x.lev;
  if(n == 0){
    printf("%.2f\n", x.coord[0]);
    return;
  }
  j = pow2(n);
  printf("{");
  for (i = 0; i < j; i++){
    if(i + 1 == j){
      if(t >= 0){               /*t is index of first zero in string */
        if(t == 0)
          printf("all");
        else
          printf("%d-%u", t, i);
        t = -1;
        virgin = FALSE;
      } else {
        if(fabs((double) (x.coord[i])) < T){
          if(!virgin)
            printf(",");
          printf("%u", i);
          t = -1;
          virgin = FALSE;
        }
      }
    } else {
      if(fabs((double) (x.coord[i])) > T){
        if(t >= 0){
          if(t + 1 == i)
            printf("%d", t);
          else
            printf("%d-%u", t, i - 1);
          t = -1;
        }
      } else {
        if(t < 0){
          t = i;
          if(!virgin)
            printf(",");
          virgin = FALSE;
        }
      }
    }
  }
  if(virgin)
    printf("none");
  printf("}=");
  printzeon(x, T, Tsmall);
}



eon Erand(uint8 n, int overtop)
{
  eon y;
  uint i, j;
  double s;
  basetype t, u;                /*make ints if want ints */
  assert(n <= LGMAXSIZE);
  y.lev = n;
  j = pow2(n);
  if(overtop > 0){              /* non-negative rands only */
    s = overtop;
    for (i = 0; i < j; i++){
      u = (((double) rand()) * s) / ((double) (RAND_MAX));
      y.coord[i] = u;
    }
  } else {
    s = -overtop;
    t = s - 1;
    s = s + s - 1;
    for (i = 0; i < j; i++){
      u = (((double) rand()) * s) / ((double) (RAND_MAX));
      y.coord[i] = u - t;
    }
  }
  for (i = j; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  return y;
}

eon EErand(uint8 n, bool pureimag, int overtop)
{
  eon y;
  uint i, j;
  double s;
  basetype t, u;                /*make ints if want ints */
  assert(n <= LGMAXSIZE);
  y.lev = n;
  j = pow2(n);
  if(overtop > 0){              /* non-negative rands only */
    s = overtop;
    for (i = 0; i < j; i++){
      u = (((double) rand()) * s) / ((double) (RAND_MAX));
      y.coord[i] = u;
    }
  } else {
    s = -overtop;
    t = s - 1;
    s = s + s - 1;
    for (i = 0; i < j; i++){
      u = (((double) rand()) * s) / ((double) (RAND_MAX));
      y.coord[i] = u - t;
    }
  }
  if(pureimag)
    y.coord[0] = 0;
  for (i = j; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  return y;
}

bool NearEq(eon a, eon b){
  uint j,k;
  assert(a.lev==b.lev);
  k = pow2(a.lev);
  assert(k<=MAXSIZE);
  for(j=0; j<k; j++){
    if( fabs(a.coord[j]-b.coord[j]) > 0.001 ) return FALSE;
  }
  return TRUE;
}

basetype maxE(eon x){
  uint j,k;
  basetype m,t;
  m = 0;
  k = pow2(x.lev);
  assert(k<=MAXSIZE);
  for(j=0; j<k; j++){
    t = x.coord[j];
    if(t<0) t = -t;
    if(t>m) m=t;
  }
  return m;
}

/* Cayley-Dickson conjugate:
   * conjugate first half, negate 2nd half;
   * equivalently negate all but real coord */
eon Conj(eon x)
{
  eon z;
  uint i, j, n;
  n = x.lev;
  z.lev = n;
  j = pow2(n);
  z.coord[0] = x.coord[0];
  for (i = 1; i < j; i++)
    z.coord[i] = -x.coord[i];
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

basetype re(eon x)
{                               /* Real part: equiv (x+Conj(x))/2 */
  return (x.coord[0]);
}

eon im(eon x)
{                               /* Imag part: equiv (x-Conj(x))/2 */
  eon z;
  uint i, j, n;
  n = x.lev;
  z.lev = n;
  j = pow2(n);
  z.coord[0] = 0;
  for (i = 1; i < j; i++)
    z.coord[i] = x.coord[i];
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

/* Sum of squares; also sum of norms of two halves */
ubasetype ss(eon x)
{
  basetype t;
  ubasetype s = 0;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    t = x.coord[i];
    s += t * t;
  }
  return s;
}

bool iszero(eon x)
{
  basetype t;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    t = x.coord[i];
    if(t != 0)
      return FALSE;
  }
  return TRUE;
}

bool isnearzero(eon x)
{
  basetype t;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    t = x.coord[i];
    if(t<0) t = -t;
    if(t > 0.01)
      return FALSE;
  }
  return TRUE;
}

bool isreal(eon x)
{
  basetype t;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 1; i < j; i++){
    t = x.coord[i];
    if(t != 0)
      return FALSE;
  }
  return TRUE;
}

bool isimag(eon x)
{
  return (x.coord[0] == 0);
}

eon scal(basetype s, eon x)
{                               /* s*x */
  eon y;
  uint i, j, n;
  n = x.lev;
  y.lev = n;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    y.coord[i] = x.coord[i] * s;
  }
  for (i = j; i < MAXSIZE; i++){
    assert(x.coord[i]==0);
    y.coord[i] = 0;
  }
  return y;
}

eon doub(eon x)
{
  return scal(2, x);
}

eon neg(eon x)
{                               /* -x */
  eon y;
  uint i, j, n;
  n = x.lev;
  y.lev = n;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    y.coord[i] = -x.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  return y;
}

eon divscal(basetype s, eon x)
{                               /* x/s */
  eon y;
  uint i, j, n;
  n = x.lev;
  y.lev = n;
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    y.coord[i] = BaseDivide(x.coord[i], s);
  }
  for (i = j; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  return y;
}

eon normalize(eon x)
{                               /* x/|x| */
  eon y;
  uint i, j, n;
  basetype s,t;
  n = x.lev;
  y.lev = n;
  assert(n <= LGMAXSIZE);
  j = pow2(n);

  t = 0;
  for (i = 0; i < j; i++){
    s = x.coord[i];
    y.coord[i] = s; 
    t += s*s;
  }
  for (i = j; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  s = 1.0/sqrt(t);  
  for (i = 0; i < j; i++){
    y.coord[i] *= s; 
  }
  return y;
}

eon UnitRand(uint8 n, bool pureimag){
  eon y;
  uint i, j;
  assert(n <= LGMAXSIZE);
  y.lev = n;
  j=pow2(n);
  for (i = 0; i < j; i++){
    y.coord[i] = normrand();
  }
  if(pureimag)
    y.coord[0] = 0;
  for (i = j; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  return normalize(y);
}

eon zero(uint8 n)
{
  eon y;
  uint i;
  assert(n <= LGMAXSIZE);
  y.lev = n;
  for (i = 0; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  return y;
}

eon one(uint8 n)
{
  eon y;
  uint i;
  assert(n <= LGMAXSIZE);
  y.lev = n;
  for (i = 1; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  y.coord[0] = 1;
  return y;
}

eon scalar(uint8 n, basetype s)
{
  eon y;
  uint i;
  assert(n <= LGMAXSIZE);
  y.lev = n;
  for (i = 1; i < MAXSIZE; i++){
    y.coord[i] = 0;
  }
  y.coord[0] = s;
  return y;
}

eon add(eon x, eon y)
{
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  for (i = 0; i < j; i++){
    z.coord[i] = x.coord[i] + y.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

basetype ip(eon x, eon y)
{
  basetype s;
  uint i, j, n;
  s = 0;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = 0; i < j; i++){
    s += x.coord[i] * y.coord[i];
  }
  return s;
}

eon basis(uint8 lev, uint8 n)
{
  eon z;
  assert(n < pow2(lev));
  z = zero(lev);
  z.coord[n] = 1;
  return z;
}

eon iunit(uint8 lev)
{
  eon z;
  uint n;
  assert(lev > 0);
  n = pow2(lev - 1);
  z = zero(lev);
  z.coord[n] = 1;
  return z;
}

/* coord[i] = max of absolute values of x.coord[i], y.coord[i] */
eon mxv(eon x, eon y)
{
  eon z;
  uint i, j, n;
  basetype s, t;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  for (i = 0; i < j; i++){
    s = x.coord[i];
    t = y.coord[i];
    if(s < 0)
      s = -s;
    if(t < 0)
      t = -t;
    if(s > t)
      t = s;
    z.coord[i] = t;
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon add3(eon x, eon y, eon q)
{
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n == y.lev);
  assert(n == q.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  for (i = 0; i < j; i++){
    z.coord[i] = x.coord[i] + y.coord[i] + q.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon add4(eon x, eon y, eon q, eon r)
{
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n == y.lev);
  assert(n == q.lev);
  assert(n == r.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  for (i = 0; i < j; i++){
    z.coord[i] = x.coord[i] + y.coord[i] + q.coord[i] + r.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon add5(eon x, eon y, eon q, eon r, eon s)
{
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n == y.lev);
  assert(n == q.lev);
  assert(n == r.lev);
  assert(n == s.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  for (i = 0; i < j; i++){
    z.coord[i] = x.coord[i] + y.coord[i]
        + q.coord[i] + r.coord[i] + s.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon sub(eon x, eon y)
{
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  for (i = 0; i < j; i++){
    z.coord[i] = x.coord[i] - y.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon inv(eon x)
{                               /* reciprocal. If integer-based eons, may not exist.
                                   * Hence recommend avoid use of inv(), prefer Conj() */
  eon y, z;
  ubasetype s;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  if(n == 0){
    y = zero(n);
    y.coord[0] = BaseDivide(1, x.coord[0]);
    return y;
  }
  assert(n <= LGMAXSIZE);
  j = pow2(n);
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  y = Conj(x);
  s = ss(x);
  z = divscal(s, y);
  return z;
}

eon reh(eon x)
{                               /* real half */
  eon z;
  uint i, j, jh, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  assert(n > 0);
  z.lev = n - 1;
  j = pow2(n);
  jh = j / 2;
  for (i = 0; i < jh; i++){
    z.coord[i] = x.coord[i];
  }
  for (i = jh; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon imh(eon x)
{                               /* imaginary half */
  eon z;
  uint i, j, jh, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  assert(n > 0);
  z.lev = n - 1;
  j = pow2(n);
  jh = j / 2;
  for (i = 0; i < jh; i++){
    z.coord[i] = x.coord[i + jh];
  }
  for (i = jh; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon rehalfie(eon x)
{                               /* real half, same level (not decremented) */
  eon z;
  uint i, j, jh, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  assert(n > 0);
  z.lev = n;
  j = pow2(n);
  jh = j / 2;
  for (i = 0; i < jh; i++){
    z.coord[i] = x.coord[i];
  }
  for (i = jh; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon octonly(eon x)
{                          /* octonion only. Same level (not decremented) */
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n > 3)
    n = 3;
  j = pow2(n);
  for (i = 0; i < j; i++){
    z.coord[i] = x.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon twoferize(eon x) 
{
    eon z;
    basetype s1, s2, fac;
    uint i,j,k;
  j = x.lev; 
  k = 1<<j;
  for(i=0; i<k; i++){ z.coord[i] = x.coord[i]; }
  z.lev = j;
  if(j <= 3) return z;
  j = k; 
  k /= 2;  
  s1 = 0;
  for(i=k; i<j; i++){ s1 += x.coord[i]; }
  s2 = 0;
  for(i=0; i<k; i++){ s2 += x.coord[i]; }
  if(s2<0.0) s2 = -s2;
  if(s2<0.01) s2 += 0.01;
  fac = s1/s2;  /* may not work well if not CONTIN */
  x.lev--;
  z=twoferize(x);
  x.lev++;
  z.lev = x.lev;
  z.coord[k] = x.coord[k];
  for(i=k+1; i<j; i++){ z.coord[i] = x.coord[i-k]*fac; }
  return z;
}

eon zeroize(eon x) 
{
    eon z;
    basetype s1, s2, fac;
    uint i,j,k,p;
  j = x.lev; 
  z.lev = j;
  k = 1<<j;
  for(i=0; i<k; i++){ z.coord[i] = x.coord[i]; }
  p=16;
  for(i=9; i<k; i++){
      if( i!=p ){ z.coord[i]=0; }else{ p += p; }
  }
  return z;
}

/* converts to niner, but if TWOFER then converts to twofer. */
eon niner(eon x)
{                         /* niner only. Same level (not decremented) */
  eon z;
  uint i, j, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n > 3)
    n = 3;
  j = pow2(n);
  if(z.lev > 3)
    j++;
  assert(j<=9);
  if(TWOFER==1){
      z=twoferize(x);
  }else if(TWOFER==2){
      z=zeroize(x);
  }else{
      for (i = 0; i < j; i++){
	  z.coord[i] = x.coord[i];
      }
      for (i = j; i < MAXSIZE; i++){
	  z.coord[i] = 0;
      }
  }
  return z;
}

/* imaginary half, same level (not decremented, not moved to re) */
eon imhalfie(eon x)
{
  eon z;
  uint i, j, jh, n;
  n = x.lev;
  assert(n <= LGMAXSIZE);
  assert(n > 0);
  z.lev = n;
  j = pow2(n);
  jh = j / 2;
  for (i = 0; i < jh; i++){
    z.coord[i] = 0;
  }
  for (i = jh; i < j; i++){
    z.coord[i] = x.coord[i];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

eon CombineHalves(eon x, eon y)
{
  eon z;
  uint i, j, jh, m, n;
  m = x.lev;
  assert(m == y.lev);
  n = m + 1;
  assert(n <= LGMAXSIZE);
  z.lev = n;
  j = pow2(n);
  jh = j / 2;
  for (i = 0; i < jh; i++){
    z.coord[i] = x.coord[i];
  }
  for (i = jh; i < j; i++){
    z.coord[i] = y.coord[i - jh];
  }
  for (i = j; i < MAXSIZE; i++){
    z.coord[i] = 0;
  }
  return z;
}

/* Cayley-Dickson multiplication law,
* [J.Baez: Bull. Amer. Math. Soc. 39 (2002) 145-205]
*    (a,b)(c,d)=(a.c-d.Conj(b), c.b+Conj(a).d)
* DIFFERS from  (a.c-Conj(d).b, bConj(c)+da  <-- [Ebbinghaus p257]
*************************************************/
eon CDmul(eon x, eon y)
{
  eon a, b, x1, x2, y1, y2, z;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  a = sub(CDmul(x1, y1), CDmul(y2, Conj(x2)));
  b = add(CDmul(y1, x2), CDmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* (a,b)(c,d) = (a.c-Conj(d).b, bConj(c)+da  <-- [Ebbinghaus p257] */
eon Ebbmul(eon x, eon y)
{
  eon a, b, x1, x2, y1, y2, z;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  a = sub(CDmul(x1, y1), CDmul(Conj(y2), x2));
  b = add(CDmul(x2, Conj(y1)), CDmul(y2, x1));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* Conway-modification of Cayley-Dickson multiplication law, code=1000
* [In Conway-Smith book] (a,b)(c,d)=(ab.inv(b)c-d.conj(b), c.b+conj(a).d)
*/
eon conmul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  s = ss(x2);
  if(isreal(x2) || s == 0)
    c = conmul(x1, y1);
  else
    c = divscal(s, conmul(conmul(x1, x2), conmul(Conj(x2), y1)));
  assert(c.lev == m);
  a = sub(c, conmul(y2, Conj(x2)));
  b = add(conmul(y1, x2), conmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* Fancier modification of Cayley-Dickson multiplication law,
* (a,b)(c,d)=(a(dConj(b)d).inv(dConj(b)d)c-d.conj(b), c.b+conj(a).d)
*/
eon conDBDmul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z, dbd;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  dbd = conDBDmul(y2, conDBDmul(Conj(x2), y2));
  s = ss(dbd);
  if(isreal(dbd) || s == 0)
    c = conDBDmul(x1, y1);
  else
    c = divscal(s,
                conDBDmul(conDBDmul(x1, dbd), conDBDmul(Conj(dbd), y1)));
  assert(c.lev == m);
  a = sub(c, conDBDmul(y2, Conj(x2)));
  b = add(conDBDmul(y1, x2), conDBDmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/***   (a, b) (c, d)  =  ( ac - bd,   bc + b [conj(a) . b d] / |b|^2 ),
     conj(a, b) = (conj(a),  -conj(b)).
Unlike Pfister original, this gives the complexes starting from the reals.
It does not give the quaternions.
***********************************/
eon pf1mul(eon x, eon y)
{
  eon a, b, c, d, t1, t2, t3, t4, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  if(n <= RECCUTOFF){
    return CDmul(x, y);
  }                             /*RECCUTOFF=3 for oct->hex via Pf */
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  a = reh(x);
  b = imh(x);
  c = reh(y);
  d = imh(y);
  s = ss(b);
  if(isreal(b) || s == 0)
    t4 = pf1mul(Conj(a), d);
  else
    t4 = divscal(s, pf1mul(b, pf1mul(Conj(a), pf1mul(b, d))));
  t1 = pf1mul(a, c);
  t2 = pf1mul(b, d);
  t3 = pf1mul(b, c);
  assert(t1.lev == m);
  assert(t2.lev == m);
  assert(t3.lev == m);
  assert(t4.lev == m);
  z = CombineHalves(sub(t1, t2), add(t3, t4));
  assert(z.lev == n);
  return z;
}

/***  
    M(x) = ( M(a)     -M(b)^T                   )    
           ( M(b)     M(b)^{-T} M(a)^T M(b)^T   ).
  (a, b) (c, d)  =  ( ac - conj(b)d,   bc + b [conj(a) . b^{-1} d] ).
     conj(a, b) = (conj(a),  -b)
Unlike Pfister original, this gives the complexes & quaternions
starting from the reals.  It does not give the octonions.
***********************************/
eon pf2mul(eon x, eon y)
{
  eon a, b, c, d, t1, t2, t3, t4, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  if(n <= RECCUTOFF){
    return CDmul(x, y);
  }                             /*RECCUTOFF=3 for oct->hex via Pf */
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  a = reh(x);
  b = imh(x);
  c = reh(y);
  d = imh(y);
  s = ss(b);
  if(isreal(b) || s == 0)
    t4 = pf2mul(Conj(a), d);
  else
    t4 = divscal(s, pf2mul(b, pf2mul(Conj(a), pf2mul(Conj(b), d))));
  t1 = pf2mul(a, c);
  t2 = pf2mul(Conj(b), d);
  t3 = pf2mul(b, c);
  assert(t1.lev == m);
  assert(t2.lev == m);
  assert(t3.lev == m);
  assert(t4.lev == m);
  z = CombineHalves(sub(t1, t2), add(t3, t4));
  assert(z.lev == n);
  return z;
}

/***    Finally, here is a further modification, which also works,
     M(x) = ( M(a)             -K M(b) K                      )
            ( K M(b)^T K       K M(b)^{-1} K M(a)^T K M(b) K  ).
and which gives you something still more closely resembling
the Cayley-Dickson (and Conway) formula:
    (a, b) (c, d)  =  ( ac - Conj(bConj(d)),
           Conj(Conj(b)Conj(c)) + Conj(b^{-1}Conj(Conj(a)Conj(bConj(d)))) ),
     conj(a, b) = (conj(a),  -b).
Here K is a matrix which obeys K^2 = I and K^T = K^{-1} = K
and which (by acting on a vector) performs conjugation.
This last formula becomes the same as the Cayley-Dickson
formula if one magically gets rid of the M(b) and M(b)^{-1} in
the 4-term-product term.  But it is not exactly the same, so this,
starting from the reals, gives you the complexes and quaternions
but does NOT give you the octonions.
***********************************/
eon pf3mul(eon x, eon y)
{
  eon a, b, c, d, t1, t2, t3, t4, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  if(n <= RECCUTOFF){
    return CDmul(x, y);
  }                             /*RECCUTOFF=3 for oct->hex via Pf */
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  a = reh(x);
  b = imh(x);
  c = reh(y);
  d = imh(y);
  s = ss(b);
  if(isreal(b) || s == 0)
    t4 = pf3mul(Conj(a), d);
  else
    t4 = divscal(s, Conj(pf3mul(Conj(b),
                                Conj(pf3mul
                                     (Conj(a),
                                      Conj(pf3mul(b, Conj(d))))))));
  t1 = pf3mul(a, c);
  t2 = Conj(pf3mul(b, Conj(d)));
  t3 = Conj(pf3mul(Conj(b), Conj(c)));
  assert(t1.lev == m);
  assert(t2.lev == m);
  assert(t3.lev == m);
  assert(t4.lev == m);
  z = CombineHalves(sub(t1, t2), add(t3, t4));
  assert(z.lev == n);
  return z;
}

/***     The RIGHT formula, modestly named after myself.
Consider the following 2x2 block matrix:
       [    A           -K B K             ]
   M = [ K B^T K    K B^T K A^T K B^{-T} K ]
Here K is an nXn matrix  K=diag(1,-1,-1,...,-1,-1);  K^T = K^{-1} = K.
(K's effect on a vector is to perform "conjugation.")
giving
    (a, b) (c, d)  =  ( ac-Conj(bConj(d)),
Conj(Conj(b)Conj(c))+Conj(Conj(b)Conj(Conj(a)Conj(Conj(b^{-1})Conj(d)))) )
     conj(a, b) = (conj(a),  -b).
***********************************/
eon WDSmul(eon x, eon y)
{
  eon a, b, c, d, t1, t2, t3, t4, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  if(n <= RECCUTOFF){
    z = CDmul(x, y);
    assert(z.lev == n);
    return z;
  }                             /*RECCUTOFF=3 for oct->hex via Pf */
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  a = reh(x);
  b = imh(x);
  c = reh(y);
  d = imh(y);
  s = ss(b);
  if(isreal(b) || s == 0)
    t4 = WDSmul(Conj(a), d);
  else
    t4 = divscal(s,
                 Conj(WDSmul(Conj(b),
                             Conj(WDSmul
                                  (Conj(a), Conj(WDSmul(b, Conj(d))))))));
  t1 = WDSmul(a, c);
  t2 = Conj(WDSmul(b, Conj(d)));
  t3 = Conj(WDSmul(Conj(b), Conj(c)));
  assert(t1.lev == m);
  assert(t2.lev == m);
  assert(t3.lev == m);
  assert(t4.lev == m);
  z = CombineHalves(sub(t1, t2), add(t3, t4));
  assert(z.lev == n);
  return z;
}

/***    The WDSmul formula
    (a, b) (c, d)  =  ( ac-Conj(bConj(d)),
Conj(Conj(b)Conj(c))+Conj(Conj(b)Conj(Conj(a)Conj(Conj(b^{-1})Conj(d)))) )
     conj(a, b) = (conj(a),  -b)
* can be simplified.
* Try #1:  Consider replacing every element (x,y) of
* the algebra by (x, Conj(y)).  Then WDSmul becomes
* [upon replacing b by Conj(b), d by Conj(d) and conjugating the
* 2nd output octonion]
    (a, b) (c, d)  =  ( ac-Conj(Conj(b)d),
                        bConj(c)+bConj(Conj(a)Conj(b^{-1}d)) )
* Try #2:  Consider replacing every element (x,y) of
* the algebra by (Conj(x), y).  Then WDSmul becomes
    (a, b) (c, d)  =  ( Conj(Conj(a)Conj(c))-bConj(d),
Conj(Conj(b)c)+Conj(Conj(b)Conj(aConj(Conj(b^{-1})Conj(d)))) )
* Try #3:  Consider replacing every element (x,y) of
* the algebra by (Conj(x), Conj(y)).  Then WDSmul becomes
    (a, b) (c, d)  =  ( Conj(Conj(a)Conj(c))-Conj(b)d,
                        bc+bConj(aConj(b^{-1}d)) )
* Of these 4, try #1 and try #3 are tied for having the fewest Conj()'s
* with 6 each. I somehow prefer try #1 aesthetically.
************************************************************
* A: Now in try #2 consider replacing every element (x,y) of
* the algebra by (x, -y).  Result is the same as doing nothing:
    (a, b) (c, d)  =  ( ac-Conj(Conj(b)d),
                        +bConj(c)+bConj(Conj(a)Conj(b^{-1}d)) )
* B: Now in try #2 consider replacing every element (x,y) of
* the algebra by (-x, y).  Result
    (a, b) (c, d)  =  ( -ac+Conj(Conj(b)d),
                        -bConj(c)-bConj(Conj(a)Conj(b^{-1}d)) )
* C: Now in try #2 consider replacing every element (x,y) of
* the algebra by (-x, -y).  Result is same as B.
* Of these A = doing nothing have fewest - signs with 1 each.
************************************************************
* We could also consider redefining (a,b)(c,d) to be (c,d)(a,b).
* The result (acting on try #1) is
    (a, b) (c, d)  =  ( CA-Conj(Conj(D)B),
                        DConj(A)+DConj(Conj(C)Conj(D^{-1}B))  )
* The result (acting on try #3) is
    (a, b) (c, d)  =  ( Conj(Conj(C)Conj(A))-Conj(D)B,
                        DA+DConj(CConj(D^{-1}B))  )
* but I don't see why either should be considered better.
CONCLUSION: The simplest/nicest is
    (a, b) (c, d)  =  ( ac-Conj(Conj(b)d),
                        bConj(c)+bConj(Conj(a)Conj(b^{-1}d)) )
     conj(a, b) = (conj(a),  -b).
Now it is possible to make 192 variations of this formula,
all of which are equivalent for octonion->hexon purposes but
inequivalent at 16->32.  We use (or not)  x*y = Conj(Conj(y)Conj(x))
on each of the 6 muls (2^6=64) and use Moufang on the 4th term
as desired. The simplest among these 192 variations, which
I shall use as a "base", are [these 3 are all Moufang-equivalent]
    (a, b) (c, d)  =  ( ac-Conj(d)b, bConj(c)+b(b^{-1}d . a) )
    (a, b) (c, d)  =  ( ac-Conj(d)b, bConj(c)+d b . b^{-1} a )   NICEST
    (a, b) (c, d)  =  ( ac-Conj(d)b, bConj(c)+(d . ab)b^{-1} )
************************************************************
Results: after exploring all 192, exactly two produced 32-ons with
multiplicative norm:
[3]   (a, b) (c, d)  =  ( ac - Conj(d)b,
fails at 64               bConj(c) + bConj(Conj(a) . Conj(d)Conj(b^{-1}) ) )

[17]  (a, b) (c, d)  =  ( ac - Conj(Conj(b)d),
works forever             bConj(c) + bConj(Conj(a) Conj(b^{-1}d) ) )
***********************************/
eon WDSmulE(eon x, eon y);

eon WDSmulEopt(bool flip, eon x, eon y)
{
  if(flip) return Conj( WDSmulE( Conj(y), Conj(x) ) );
  return WDSmulE(x,y);
}

eon WDSmulE(eon x, eon y)
{
  eon a, b, c, d, t1, t2, t3, t4, z;
  basetype s;
  bool f1,f2,f3,f4,f5,f6;
  uint n, m, Q;

  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);

  Q = WDSMULCHOICE;
  assert(Q<192);
  f1 = Q%2;
  Q /= 2;
  f2 = Q%2;
  Q /= 2;
  f3 = Q%2;
  Q /= 2;
  f4 = Q%2;
  Q /= 2;
  f5 = Q%2;
  Q /= 2;
  f6 = Q%2;
  Q /= 2;
  assert(Q<=2);
  /*if(n>=5) printf("%d%d%d%d%d%d;%d\n",f1,f2,f3,f4,f5,f6,Q);*/
  /* The two that work at 32 are  100010;0=17,
   * 110000;0=3 [but fails at 64 and not Lalt at 32 and not R-linear at 32] */

  if(n <= RECCUTOFF){
    z = CDmul(x, y);
    assert(z.lev == n);
    return z;
  }                             /*RECCUTOFF=3 for oct->hex via Pf */
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  a = reh(x);
  b = imh(x);
  c = reh(y);
  d = imh(y);
  s = ss(b);
  if(isreal(b) || s == 0)
    t4 = WDSmulEopt(f1, d, Conj(a));
  else{
    if(Q==0)
       t4 = divscal(s, WDSmulEopt(f3, b, WDSmulEopt(f1, WDSmulEopt(f2, Conj(b),d), a)) );
    else{ if(Q==1)
       t4 = divscal(s, WDSmulEopt(f1, WDSmulEopt(f2, d,b), WDSmulEopt(f3, Conj(b),a)) );
          else /*if(Q==2)*/
            t4 = divscal(s, WDSmulEopt(f2,WDSmulEopt(f1,d,WDSmulEopt(f3,a,b)),Conj(b)) );
        }
  }
  t1 = WDSmulEopt(f4, a, c);
  t2 = WDSmulEopt(f5, Conj(d), b);
  t3 = WDSmulEopt(f6, b, Conj(c));
  assert(t1.lev == m);
  assert(t2.lev == m);
  assert(t3.lev == m);
  assert(t4.lev == m);
  z = CombineHalves(sub(t1, t2), add(t3, t4));
  assert(z.lev == n);
  return z;
}

eon NaiveDiv(eon cd, eon XY)
{
  eon ab;
  assert(cd.lev == XY.lev);
  ab = divscal(ss(cd), WDSmul(XY, Conj(cd)));
  return ab;
}

/* finds (a,b) so (a,b)(c,d)=(X,Y).  Uses

                     |c|^2 Y.Conj(d) - |d|^2 c.Conj(X)
         b = Conj(c) --------------------------------- d
                       |c|^2 ( |c|^2 + |d|^2 ) |d|^2

         a = (X+dConj(b))c^{-1}
THIS IS WRONG. MY DERIVATION WAS WRONG. SHOULD BE OK IN QUATS/OCTS
JUST NOT HEXONS.
**************************************/
eon WDSdiv(eon cd, eon XY)
{
  eon a, b, c, d, X, Y, z, ab;
  uint n;
  basetype ssc, ssd;
  n = cd.lev;
  assert(n == XY.lev);
  if(n <= 1)
    /*RECCUTOFF*/ return NaiveDiv(cd, XY);
  c = reh(cd);
  d = imh(cd);
  X = reh(XY);
  Y = imh(XY);
  assert(1 + c.lev == n);
  assert(1 + d.lev == n);
  assert(1 + X.lev == n);
  assert(1 + Y.lev == n);
  ssc = ss(c);
  ssd = ss(d);
  z = sub(scal(ssc, WDSmul(Y, Conj(d))), scal(ssd, WDSmul(c, Conj(X))));
  assert(1 + z.lev == n);
  b = divscal(ssc * (ssc + ssd) * ssd, WDSmul(WDSmul(Conj(c), z), d));
  assert(1 + b.lev == n);
  a = divscal(ssc, WDSmul(add(X, WDSmul(d, Conj(b))), Conj(c)));
  assert(1 + a.lev == n);
  ab = CombineHalves(a, b);
  assert(ab.lev == n);
  return ab;
}

eon divcheck1(eon cd, eon XY)
{                              
  eon ab, z;
  ab = WDSdiv(cd, XY);
  z = WDSmul(ab, cd);
  return sub(XY, z);
}

eon divcheck2(eon cd, eon XY)
{
  eon ab, z;
  ab = WDSdiv(cd, XY);
  z = NaiveDiv(cd, XY);
  return sub(ab, z);
}

eon selfback(eon x)
{
  return x;
}

/***  Hans Zassenhaus & Wolfgang Eichhorn:
Herleitung von Acht- und Sechzehn-Quadrate-Identit\"aten mit
Hilfe von Eigenschaften der verallgemeinerten Quaternionen und der
Cayley-Dicksonschen Zahlen,
Archiv der Mathematik  17 (1966) 492-496
    Z&E follow Taussky's lead and for the 2x2 octonion matrix
    a b
    c d
define
    det = a d - (a c) (a^{-1} b)
and find
    (det M) conj(det M) = det(M M^H)
[noting as an aside, following Taussky,
that  conj(det M)  is generally not equal det(M^H)]
which upon rearranging terms becomes the 16-square identity
   ( |a|^2 + |b|^2 ) ( |c|^2 + |d|^2 ) =
      | a conj(c) + b conj(d) |^2 +
      | a d - (a c)(a^{-1} b) |^2
where a,b,c,d are octonions.  This suggests
(a,b)(c,d) = ( a Conj(c) + b Conj(d) , a d - ac.a^{-1}b )
***********************************************************/
eon zassmul(eon x, eon y)
{
  eon a, b, c, d, t1, t2, t3, t4, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  if(n <= RECCUTOFF){
    return CDmul(x, y);
  }                             /*RECCUTOFF=3 for oct->hex via zass */
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  a = reh(x);
  b = imh(x);
  c = reh(y);
  d = imh(y);
  s = ss(a);
  if(isreal(a) || s == 0)
    t4 = zassmul(c, b);
  else
    t4 = divscal(s, zassmul(zassmul(a, c), zassmul(Conj(a), b)));
  t1 = zassmul(a, Conj(c));
  t2 = zassmul(b, Conj(d));
  t3 = zassmul(a, d);
  assert(t1.lev == m);
  assert(t2.lev == m);
  assert(t3.lev == m);
  assert(t4.lev == m);
  z = CombineHalves(add(t1, t2), sub(t3, t4));
  assert(z.lev == n);
  return z;
}

/****************************************/
eon g1cmul(eon x, eon y);

eon wedgein(eon a, eon b, eon c)
{
  basetype s;
  s = ss(b);
  if(isreal(b) || s == 0)
    return g1cmul(a, c);
  return divscal(s, g1cmul(g1cmul(a, b), g1cmul(Conj(b), c)));
}

eon opt(uint Q, eon a, eon b, eon c, eon d)
{
  if(Q == 0)
    return g1cmul(a, d);
  if(Q == 1)
    return wedgein(a, b, d);
  if(Q == 2)
    return wedgein(a, Conj(b), d);
  if(Q == 3)
    return wedgein(a, c, d);
  if(Q == 4)
    return wedgein(a, Conj(c), d);
  if(Q == 5)
    return wedgein(a, a, d);
  if(Q == 6)
    return wedgein(a, Conj(a), d);
  if(Q == 7)
    return wedgein(a, d, d);
  if(Q == 8)
    return wedgein(a, Conj(d), d);
  fprintf(stderr, "opt given bad choice parameter Q=%d\n", Q);
  exit(1);
}

void printopt(uint Q)
{
  /* (a,b)(c,d) = */
  uint Q1, Q2, Q3, Q4;
  Q4 = Q % 9;
  Q /= 9;
  Q3 = Q % 9;
  Q /= 9;
  Q2 = Q % 9;
  Q /= 9;
  Q1 = Q % 9;
  Q /= 9;
  printf("( a ");
  if(Q1 == 0){;
  }
  if(Q1 == 1)
    printf("b . b^{-1} ");
  if(Q1 == 2)
    printf("b^{-1} . b ");
  if(Q1 == 3)
    printf("d . d^{-1} ");
  if(Q1 == 4)
    printf("d^{-1} . d ");
  if(Q1 == 5)
    printf("a . a^{-1} ");
  if(Q1 == 6)
    printf("a^{-1} . a ");
  if(Q1 == 7)
    printf("c . c^{-1} ");
  if(Q1 == 8)
    printf("c^{-1} . c ");
  printf("c - d ");
  if(Q2 == 0){;
  }
  if(Q2 == 1)
    printf("a . a^{-1} ");
  if(Q2 == 2)
    printf("a^{-1} . a ");
  if(Q2 == 3)
    printf("c . c^{-1} ");
  if(Q2 == 4)
    printf("c^{-1} . c ");
  if(Q2 == 5)
    printf("b . b^{-1} ");
  if(Q2 == 6)
    printf("b^{-1} . b ");
  if(Q2 == 7)
    printf("d . d^{-1} ");
  if(Q2 == 8)
    printf("d^{-1} . d ");
  printf("\\conj{b} , \\, c ");
  if(Q3 == 0){;
  }
  if(Q3 == 1)
    printf("a . a^{-1} ");
  if(Q3 == 2)
    printf("a^{-1} . a ");
  if(Q3 == 3)
    printf("d . d^{-1} ");
  if(Q3 == 4)
    printf("d^{-1} . d ");
  if(Q3 == 5)
    printf("b . b^{-1} ");
  if(Q3 == 6)
    printf("b^{-1} . b ");
  if(Q3 == 7)
    printf("c . c^{-1} ");
  if(Q3 == 8)
    printf("c^{-1} . c ");
  printf("b + \\conj{a} ");
  if(Q4 == 0){;
  }
  if(Q4 == 1)
    printf("b . b^{-1} ");
  if(Q4 == 2)
    printf("b^{-1} . b ");
  if(Q4 == 3)
    printf("c . c^{-1} ");
  if(Q4 == 4)
    printf("c^{-1} . c ");
  if(Q4 == 5)
    printf("a . a^{-1} ");
  if(Q4 == 6)
    printf("a^{-1} . a ");
  if(Q4 == 7)
    printf("d . d^{-1} ");
  if(Q4 == 8)
    printf("d^{-1} . d ");
  printf("d )");
}

eon g1cmul(eon x, eon y)
{
  eon a, b, a1, a2, b1, b2, x1, x2, y1, y2, z;
  uint Q1, Q2, Q3, Q4, Q;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  Q = MULCHOICE;
  Q4 = Q % 9;
  Q /= 9;
  Q3 = Q % 9;
  Q /= 9;
  Q2 = Q % 9;
  Q /= 9;
  Q1 = Q % 9;
  Q /= 9;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  assert(x1.lev == m);
  a1 = opt(Q1, x1, x2, y2, y1);
  assert(a1.lev == m);
  a2 = opt(Q2, y2, x1, y1, Conj(x2));
  a = sub(a1, a2);
  b1 = opt(Q3, y1, x1, y2, x2);
  b2 = opt(Q4, Conj(x1), x2, y1, y2);
  b = add(b1, b2);
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/****************************************/

eon g2cmul(eon x, eon y);
eon grecmul(eon x, eon y);

eon wedgein2(eon a, eon b, eon c)
{
  basetype s;
  s = ss(b);
  if(isreal(b) || s == 0)
    return grecmul(a, c);
  return divscal(s, grecmul(grecmul(a, b), grecmul(Conj(b), c)));
}

eon wedgeinR(eon a, eon b, eon c)
{
  basetype s;
  s = ss(b);
  if(isreal(b) || s == 0)
    return grecmul(a, c);
  return divscal(s, grecmul(grecmul(a, b), grecmul(c, Conj(b))));
}

eon wedgeinL(eon a, eon b, eon c)
{
  basetype s;
  s = ss(b);
  if(isreal(b) || s == 0)
    return grecmul(a, c);
  return divscal(s, grecmul(grecmul(b, a), grecmul(Conj(b), c)));
}

eon opt2(uint Q, eon a, eon b, eon c, eon d)
{
  if(Q == 0)
    return grecmul(a, d);
  if(Q == 1)
    return wedgein2(a, b, d);
  if(Q == 2)
    return wedgein2(a, Conj(b), d);
  if(Q == 3)
    return wedgein2(a, c, d);
  if(Q == 4)
    return wedgein2(a, Conj(c), d);

  if(Q == 5)
    return wedgeinL(a, b, d);
  if(Q == 6)
    return wedgeinL(a, Conj(b), d);
  if(Q == 7)
    return wedgeinL(a, c, d);
  if(Q == 8)
    return wedgeinL(a, Conj(c), d);

  if(Q == 9)
    return wedgeinR(a, b, d);
  if(Q == 10)
    return wedgeinR(a, Conj(b), d);
  if(Q == 11)
    return wedgeinR(a, c, d);
  if(Q == 12)
    return wedgeinR(a, Conj(c), d);

  fprintf(stderr, "opt2 given bad choice parameter Q=%d\n", Q);
  exit(1);
}

/* This recurses to grecmul - sometimes want to recurse to g2cmul,
* other times want CDmul instead. Change grecmul to change that.
******/
eon g2cmul(eon x, eon y)
{
  eon a, b, a1, a2, b1, b2, x1, x2, y1, y2, z, tmp;

  uint MySWIND, Q1, Q2, Q3, Q4, Q;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  Q = COMBOMULCHOICE;
  Q4 = Q % 13;
  Q /= 13;
  Q3 = Q % 13;
  Q /= 13;
  Q2 = Q % 13;
  Q /= 13;
  Q1 = Q % 13;
  Q /= 13;
  MySWIND = Q;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  assert(x1.lev == m);
  if(MySWIND & 8){
    CONJSWAP(x1, y1);
  }
  a1 = opt2(Q1, x1, x2, y2, y1);
  if(MySWIND & 8){
    CONJSWAP(x1, y1);
  }
  assert(a1.lev == m);

  if(MySWIND & 4){
    CONJSWAP(y2, x2);
  }
  a2 = opt2(Q2, y2, x1, y1, Conj(x2));
  if(MySWIND & 4){
    CONJSWAP(y2, x2);
  }
  a = sub(a1, a2);

  if(MySWIND & 2){
    CONJSWAP(y1, x2);
  }
  b1 = opt2(Q3, y1, x1, y2, x2);
  if(MySWIND & 2){
    CONJSWAP(y1, x2);
  }

  if(MySWIND & 1){
    CONJSWAP(x1, y2);
  }
  b2 = opt2(Q4, Conj(x1), x2, y1, y2);
  if(MySWIND & 1){
    CONJSWAP(x1, y2);
  }
  b = add(b1, b2);
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

eon grecmul(eon x, eon y)
{
  if(x.lev <= RECCUTOFF)
    return CDmul(x, y);
  return g2cmul(x, y);
}

/*********************************************************/

/* Try to swap binv and b in conway mul  - NOT norm-multip at 16ons */
eon conRmul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  s = ss(x2);
  if(isreal(x2) || s == 0)
    c = conRmul(x1, y1);
  else
    c = divscal(s, conRmul(conRmul(x1, Conj(x2)), conRmul(x2, y1)));
  assert(c.lev == m);
  a = sub(c, conRmul(y2, Conj(x2)));
  b = add(conRmul(y1, x2), conRmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* replace b.binv by b^2.binv^2 in conway formula - NOT norm-multip at 16ons*/
eon con2mul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z, x2s;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  x2s = con2mul(x2, x2);
  s = ss(x2s);
  if(isreal(x2) || s == 0)
    c = con2mul(x1, y1);
  else
    c = divscal(s, con2mul(con2mul(x1, x2s), con2mul(Conj(x2s), y1)));
  assert(c.lev == m);
  a = sub(c, con2mul(y2, Conj(x2)));
  b = add(con2mul(y1, x2), con2mul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* replace b by b+17 in conway formula - NOT norm-mult at 16ons */
eon con17mul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z, x2k;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  x2k = add(x2, scalar(x2.lev, 17));
  s = ss(x2k);
  if(isreal(x2) || s == 0)
    c = con17mul(x1, y1);
  else
    c = divscal(s, con17mul(con17mul(x1, x2k), con17mul(Conj(x2k), y1)));
  assert(c.lev == m);
  a = sub(c, con17mul(y2, Conj(x2)));
  b = add(con17mul(y1, x2), con17mul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* replace b.binv by (b+d).(b+d)inv in conway formula.
* Still norm-multip at 16 but not alternative nor distributive
* on left or right or middle. Still, flexible in last 8 coords. */
eon conAmul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z, x2s;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  x2s = add(x2, y2);
  s = ss(x2s);
  if(isreal(x2) || s == 0)
    c = conAmul(x1, y1);
  else
    c = divscal(s, conAmul(conAmul(x1, x2s), conAmul(Conj(x2s), y1)));
  assert(c.lev == m);
  a = sub(c, conAmul(y2, Conj(x2)));
  b = add(conAmul(y1, x2), conAmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* replace b.binv by (3b+5d).(3b+5d)inv in conway formula.
* Still norm-multip at 16 but not alternative nor distributive
* on left or right or middle. Not flexible in last 8 coords either. */
eon conWmul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z, x2s;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  x2s = add(scal(3, x2), scal(5, y2));
  s = ss(x2s);
  if(isreal(x2) || s == 0)
    c = conWmul(x1, y1);
  else
    c = divscal(s, conWmul(conWmul(x1, x2s), conWmul(Conj(x2s), y1)));
  assert(c.lev == m);
  a = sub(c, conWmul(y2, Conj(x2)));
  b = add(conWmul(y1, x2), conWmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* replace b.binv by (b-d).(b-d)inv in conway formula.
* Seems to have about same properties as conAmul().
* Still norm-multip at 16 but not alternative nor distributive
* on left or right or middle. Still, flexible in last 8 coords. */
eon conSmul(eon x, eon y)
{
  eon a, b, c, x1, x2, y1, y2, z, x2s;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  x2s = sub(x2, y2);
  s = ss(x2s);
  if(isreal(x2) || s == 0)
    c = conSmul(x1, y1);
  else
    c = divscal(s, conSmul(conSmul(x1, x2s), conSmul(Conj(x2s), y1)));
  assert(c.lev == m);
  a = sub(c, conSmul(y2, Conj(x2)));
  b = add(conSmul(y1, x2), conSmul(Conj(x1), y2));
  assert(a.lev == m);
  assert(b.lev == m);
  z = CombineHalves(a, b);
  assert(z.lev == n);
  return z;
}

/* JDH Smith's modification of Cayley-Dickson multiplication law
* [J.Algeb 176 (1995) 128-138, eq 4.2 p132]   code=900a - NOT!
*  (a,b)(c,d)=(ab.cInv(b)-b.Conj(d), b.Conj(c)+dInv(b).ab)
* notice this DIFFERS from CDmul's  [Baez, ZSSS p28, Schafer p45]
*             (a.c       -d.Conj(b), c.b      +Conj(a).d)
* & DIFFERS   (a.c       -Conj(d).b, b.Conj(c)+da  <-- [Ebbinghaus p257]
*  JDHS is not a modified cayley-dickson formula at all!
*****************************************************************
* Trouble is, JDHS does not want to use his mul law recursively;
* he wants to use Cayley-Dickson unmodified law at every recursive
* level below the top. I do so here as he desires. */
eon jdhsmul(eon x, eon y)
{
  eon a, b, c, d, x1, x2, y1, y2, z;
  basetype s;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  s = ss(x2);
  if(isreal(x2) || s == 0){
    c = CDmul(x1, y1);
    d = CDmul(y2, x1);
  } else {
    c = divscal(s, CDmul(CDmul(x1, x2), CDmul(y1, Conj(x2))));
    d = divscal(s, CDmul(CDmul(y2, Conj(x2)), CDmul(x1, x2)));
  }
  a = sub(c, CDmul(x2, Conj(y2)));
  b = add(CDmul(x2, Conj(y1)), d);
  z = CombineHalves(a, b);
  return z;
}

/* Expts indicate this is NOT a valid CD formula - does not yield
* multiplicative norm at the 8ons  !! */
eon jdhsBaseMul(eon x, eon y)
{
  eon a, b, c, d, x1, x2, y1, y2, z;
  uint n, m;
  n = x.lev;
  assert(n == y.lev);
  assert(n <= LGMAXSIZE);
  z.lev = n;
  if(n == 0){
    z = zero(n);
    z.coord[0] = x.coord[0] * y.coord[0];
    return z;
  }
  assert(n > 0);
  m = n - 1;
  x1 = reh(x);
  x2 = imh(x);
  y1 = reh(y);
  y2 = imh(y);
  c = CDmul(x1, y1);
  d = CDmul(y2, x1);
  a = sub(c, CDmul(x2, Conj(y2)));
  b = add(CDmul(x2, Conj(y1)), d);
  z = CombineHalves(a, b);
  return z;
}

/* The big question is: Which mul shall we use? */
eon mul(eon x, eon y)
{
  switch (COMBOMULCHOICE){
  case WDScode:
    return WDSmul(x, y);
  case WDSeCode:
    return WDSmulE(x, y);
  case CDcode:
  case CON1code:
  case CON2code:
  case CON3code:
  case CON4code:
  case CON5code:
  case CON6code:
  case CON7code:
  case CON8code:
  case WEIRDcode:
    return g2cmul(x, y);

  case PF1code:
    return pf1mul(x, y);
  case PF2code:
    return pf2mul(x, y);
  case PF3code:
    return pf3mul(x, y);
  case JDHScode:
    return jdhsmul(x, y);
  case ZASScode:
    return zassmul(x, y);
  }
  /* default */ return CDmul(x, y);
}

/* multiply then zero the real part */
eon imul(eon x, eon y)
{
  eon z;
  z = mul(x, y);
  z.coord[0] = 0;
  return z;
}

eon m3(eon a, eon b, eon c)
{
  return mul(a, mul(b, c));
}

eon m4(eon a, eon b, eon c, eon d)
{
  return m3(a, b, mul(c, d));
}

/***********************************************************/

/* x*x */
eon sq(eon x)
{
  return mul(x, x);
}

eon cube(eon x)
{
  return m3(x, x, x);
}

basetype reab(eon a, eon b)
{
    return ( re(mul(a,b)) - ip(a,Conj(b)) );
}

eon bimul(eon x, eon y)
{
    return mul(mul(x,y),x);
}

eon comul(eon x, eon y)
{
    return mul(mul(x,y),Conj(x));
}

eon linbimul(eon x, eon y, eon z)
{
    return sub( add(mul(mul(x,y),x), mul(mul(x,z),x)),
		mul(mul(x,add(y,z)),x) );
}

eon lincomul(eon x, eon y, eon z)
{
    return sub( add(mul(mul(x,y),Conj(x)), mul(mul(x,z),Conj(x))),
		mul(mul(x,add(y,z)),Conj(x)) );
}

eon ninlinbimul(eon x, eon y, eon z)
{
    eon q; q = niner(x);
    return sub( add(mul(mul(q,y),q), mul(mul(q,z),q)),
		mul(mul(q,add(y,z)),q) );
}

eon ninlincomul(eon x, eon y, eon z)
{
    eon q; q = niner(x);
    return sub( add(mul(mul(q,y),Conj(q)), mul(mul(q,z),Conj(q))),
		mul(mul(q,add(y,z)),Conj(q)) );
}

eon bi2mul(eon x, eon y)
{
    return mul(x,mul(y,x));
}

eon co2mul(eon x, eon y)
{
    return mul( x,mul(y,Conj(x)) );
}

eon linbi2mul(eon x, eon y, eon z)
{
    return sub( add(mul(x,mul(y,x)), mul(x,mul(z,x))),
		mul(x,mul(add(y,z),x)) );
}

eon linco2mul(eon x, eon y, eon z)
{
    return sub( add(mul(x,mul(y,Conj(x))), mul(x,mul(z,Conj(x)))),
		mul(x,mul(add(y,z),Conj(x))) );
}

eon ninlinbi2mul(eon x, eon y, eon z)
{
    eon q; q = niner(x);
    return sub( add(mul(q,mul(y,q)), mul(q,mul(z,q))),
		mul(q,mul(add(y,z),q)) );
}

eon ninlinco2mul(eon x, eon y, eon z)
{
    eon q; q = niner(x);
    return sub( add(mul(q,mul(y,Conj(q))), mul(q,mul(z,Conj(q)))),
		mul(q,mul(add(y,z),Conj(q))) );
}

eon biCDmul(eon x, eon y)
{
    return CDmul(CDmul(x,y),x);
}

eon coCDmul(eon x, eon y)
{
    return CDmul(CDmul(x,y),Conj(x));
}

eon CDvsme(eon x, eon y)  /* compares Cayley-Dickson CDmul vs the real mul */
{
  return sub( CDmul(x,y), mul(x,y) );
}

eon CDsqsame(eon x)  /* compares Cayley-Dickson square vs real */
{
  return sub( CDmul(x,x), mul(x,x) );
}

eon bcx;
eon basCDvsme(eon y)
{
  return sub( CDmul(y,bcx), mul(y,bcx) );
}

eon linbasCDvsme(eon y, eon z)
{
  return sub( mul(add(y,z),bcx), add(mul(y,bcx), mul(z,bcx)) );
}

eon biCDvsme(eon x, eon y)  /* compares Cayley-Dickson biCDmul vs the real bimul */
{
  return sub( biCDmul(x,y), bimul(x,y) );
}

eon coCDvsme(eon x, eon y)  /* compares Cayley-Dickson biCDmul vs the real bimul */
{
  return sub( coCDmul(x,y), comul(x,y) );
}

eon bi2CDvsme(eon x, eon y)  /* compares Cayley-Dickson biCDmul vs the real bimul */
{
  return sub( biCDmul(x,y), bi2mul(x,y) );
}

eon co2CDvsme(eon x, eon y)  /* compares Cayley-Dickson biCDmul vs the real bimul */
{
  return sub( coCDmul(x,y), co2mul(x,y) );
}

eon ninCDvsme(eon x, eon y)  /* compares Cayley-Dickson CDmul vs the real mul */
{
    eon z;  z = niner(y);
  return sub( CDmul(x,z), mul(x,z) );
}

eon nin1CDvsme(eon x, eon y)  /* compares Cayley-Dickson CDmul vs the real mul */
{
    eon z;  z = niner(x);
  return sub( CDmul(z,x), mul(z,x) );
}

eon ninninCDvsme(eon x, eon y)  /* compares Cayley-Dickson CDmul vs the real mul */
{
    eon q,z;  z = niner(y);  q = niner(x);
  return sub( CDmul(q,z), mul(q,z) );
}

eon ninbiCDvsme(eon x, eon y)  /* compares Cayley-Dickson biCDmul vs the real bimul */
{
    eon z;  z = niner(y);
  return sub( biCDmul(x,z), bimul(x,z) );
}

eon nincoCDvsme(eon x, eon y)  /* compares Cayley-Dickson coCDmul vs the real comul */
{
    eon z;  z = niner(y);
  return sub( coCDmul(x,z), comul(x,z) );
}

eon ninbi2CDvsme(eon x, eon y)  /* compares Cayley-Dickson biCDmul vs the real bi2mul */
{
    eon z;  z = niner(y);
  return sub( biCDmul(x,z), bi2mul(x,z) );
}

eon ninco2CDvsme(eon x, eon y)  /* compares Cayley-Dickson coCDmul vs the real co2mul */
{
    eon z;  z = niner(y);
  return sub( coCDmul(x,z), co2mul(x,z) );
}

/* (xy+yx)   unscaled Jordan multiplication */
eon ujor(eon x, eon y)
{
  return add(mul(x, y), mul(y, x));
}

/* (xy+yx)/2   Jordan multiplication, field with char not 2 */
eon jor(eon x, eon y)
{
  return divscal(2, ujor(x, y));
}

/* (xy-yx)   Lie multiplication */
eon lie(eon x, eon y)
{
  return sub(mul(x, y), mul(y, x));
}

eon sqlie(eon x, eon y)
{
  return sq(lie(x, y));
}

eon sqjor(eon x, eon y)
{
  return sq(jor(x, y));
}

eon squjor(eon x, eon y)
{
  return sq(ujor(x, y));
}

eon pow4lie(eon x, eon y)
{
  return sq(sqlie(x, y));
}

eon lieC(eon x)
{
  return lie(x, Conj(x));
}

eon sqlieC(eon x)
{
  return sqlie(x, Conj(x));
}

eon pow4lieC(eon x)
{
  return pow4lie(x, Conj(x));
}

/* (xy-yx)/2  half Lie multiplication = cross product for pure imag */
eon cross(eon x, eon y)
{
  return divscal(2, sub(mul(x, y), mul(y, x)));
}

/* xy.z - x.yz */
eon assoc(eon x, eon y, eon z)
{
  eon a, b, c;
  a = mul(mul(x, y), z);
  b = mul(x, mul(y, z));
  c = sub(a, b);
  assert(a.lev == x.lev);
  assert(b.lev == x.lev);
  assert(c.lev == x.lev);
  return c;
}

eon sqassoc(eon x, eon y, eon z)
{
  return sq(assoc(x, y, z));
}

/* xy.z - x.yz using imul */
eon iassoc(eon x, eon y, eon z)
{
  eon a, b, c;
  a = imul(imul(x, y), z);
  b = imul(x, imul(y, z));
  c = sub(a, b);
  assert(a.lev == x.lev);
  assert(b.lev == x.lev);
  assert(c.lev == x.lev);
  return c;
}

/* xy.z - x.yz using cross */
eon Cassoc(eon x, eon y, eon z)
{
  eon a, b, c;
  a = cross(cross(x, y), z);
  b = cross(x, cross(y, z));
  c = sub(a, b);
  assert(a.lev == x.lev);
  assert(b.lev == x.lev);
  assert(c.lev == x.lev);
  return c;
}

/* xy.z - x.yz using ujor */
eon Jassoc(eon x, eon y, eon z)
{
  eon a, b, c;
  a = ujor(ujor(x, y), z);
  b = ujor(x, ujor(y, z));
  c = sub(a, b);
  assert(a.lev == x.lev);
  assert(b.lev == x.lev);
  assert(c.lev == x.lev);
  return c;
}

/* Khalil-Yiu Lem 1.4.1 attributed Adem1975. 
 * Holds for CayleyDickson algebras.*/
eon ky141(eon x, eon y, eon z, eon w)
{
    return sub(
	add3(mul(x,assoc(y,z,w)), mul(assoc(x,y,z),w), assoc(x,mul(y,z),w)),
	add(assoc(mul(x,y),z,w), assoc(x,y,mul(z,w)) ) );
}

/* Khalil-Yiu Lem 1.4.2 */
eon ky142a(eon x, eon y)
{
    return sub(
	add( jor(x,y), scal(ip(x,y), one(x.lev)) ),
	add( scal(re(x),y), scal(re(y),x) ) );
}


/* yxy = 2<Conj(x),y>y - <y,y>Conj(x)   triple product for quats & octs */
eon trip(eon x, eon y)
{
  return sub(add(mul(y, mul(x, y)), mul(mul(y, x), y)),
             sub(scal(4 * ip(Conj(x), y), y), scal(2 * ss(y), Conj(x))));
}

eon sqtrip(eon x, eon y)
{
  return sq(trip(x, y));
}

/* flexible xy.x-x.yx */
eon flex(eon x, eon y)
{
  return assoc(x, y, x);
}

eon flex2(eon x, eon y)
{
  return sub(mul(mul(x, y), x), mul(x, mul(y, x)));
}

eon sqflex(eon x, eon y)
{
  return sq(flex(x, y));
}

eon pow4flex(eon x, eon y)
{
  return sq(sqflex(x, y));
}

eon octflex1(eon x, eon y)
{
  return flex(niner(x), y);
}

eon octflex2(eon x, eon y)
{
  return flex(x, niner(y));
}

eon ninflex1(eon x, eon y)
{
  return flex(niner(x), y);
}

eon ninflex2(eon x, eon y)
{
  return flex(x, niner(y));
}

/* left alternative */
eon Lalt(eon x, eon y)
{
  return assoc(x, x, y);
}

eon octLalt2(eon x, eon y)
{
  return Lalt(x, niner(y));
}

eon octLalt1(eon x, eon y)
{
  return Lalt(niner(x), y);
}

eon ninLalt2(eon x, eon y)
{
  return Lalt(x, niner(y));
}

eon ninLalt1(eon x, eon y)
{
  return Lalt(niner(x), y);
}

eon sqLalt(eon x, eon y)
{
  return sq(Lalt(x, y));
}

eon pow4Lalt(eon x, eon y)
{
  return sq(sqLalt(x, y));
}                               /*Mikheev */

/* right alternative */
eon Ralt(eon x, eon y)
{
  return assoc(x, y, y);
}

eon octRalt1(eon x, eon y)
{
  return Ralt(x, niner(y));
}

eon octRalt2(eon x, eon y)
{
  return Ralt(niner(x), y);
}

eon ninninRalt(eon x, eon y)
{
  return Ralt(niner(x), niner(y));
}

eon ninRalt1(eon x, eon y)
{
  return Ralt(x, niner(y));
}

eon ninRalt2(eon x, eon y)
{
  return Ralt(niner(x), y);
}

eon sqRalt(eon x, eon y)
{
  return sq(Ralt(x, y));
}

eon pow4Ralt(eon x, eon y)
{
  return sq(sqRalt(x, y));
}                               /*Mikheev */

eon eyeshiftA(eon a, eon b)
{
  uint n;
  eon e;
  n = a.lev;
  assert(n == b.lev);
  assert(n>=2);
  assert(n<=LGMAXSIZE);
  e = basis(n,pow2(n-1));
  return sub(mul(b, mul(e,a)), mul(e, mul(Conj(b),a)));
}

eon eyeshiftB(eon a, eon b)
{
  uint n;
  eon e;
  n = a.lev;
  assert(n == b.lev);
  assert(n>=2);
  assert(n<=LGMAXSIZE);
  e = basis(n,pow2(n-1));
  return sub(mul(mul(e,b),a), mul(e, mul(a,b)));
}

eon eyeshiftC(eon a, eon b)
{
  uint n;
  eon e;
  n = a.lev;
  assert(n == b.lev);
  assert(n>=2);
  assert(n<=LGMAXSIZE);
  e = basis(n,pow2(n-1));
  return add( mul(mul(e,b),mul(e,a)),  mul(a,Conj(b)) );
}

eon octeyeshiftA(eon a, eon b){  return eyeshiftA(octonly(a), octonly(b)); }
eon octeyeshiftB(eon a, eon b){  return eyeshiftB(octonly(a), octonly(b)); }
eon octeyeshiftC(eon a, eon b){  return eyeshiftC(octonly(a), octonly(b)); }

eon nineyeshiftA(eon a, eon b){  return eyeshiftA(niner(a), niner(b)); }
eon nineyeshiftB(eon a, eon b){  return eyeshiftB(niner(a), niner(b)); }
eon nineyeshiftC(eon a, eon b){  return eyeshiftC(niner(a), niner(b)); }

eon ReyeshiftA(eon a, eon b)
{                            /* ae.b = aConj(b).e */
  uint n;
  eon e;
  n = a.lev;
  assert(n == b.lev);
  assert(n>=2);
  assert(n<=LGMAXSIZE);
  e = basis(n,pow2(n-1));
  return sub(mul(mul(a,e),b), mul(mul(a,Conj(b)),e));
}

eon ReyeshiftB(eon a, eon b)
{                           /* a.be = ba.e */
  uint n;
  eon e;
  n = a.lev;
  assert(n == b.lev);
  assert(n>=2);
  assert(n<=LGMAXSIZE);
  e = basis(n,pow2(n-1));
  return sub(mul(a,mul(b,e)), mul(mul(b,a),e));

}

eon ReyeshiftC(eon a, eon b)
{                      /*   ae.be = -Conj(b)a  */
  uint n;
  eon e;
  n = a.lev;
  assert(n == b.lev);
  assert(n>=2);
  assert(n<=LGMAXSIZE);
  e = basis(n,pow2(n-1));
  return add( mul(mul(a,e),mul(b,e)),  mul(Conj(b),a) );
}

eon octReyeshiftA(eon a, eon b){  return ReyeshiftA(octonly(a), octonly(b)); }
eon octReyeshiftB(eon a, eon b){  return ReyeshiftB(octonly(a), octonly(b)); }
eon octReyeshiftC(eon a, eon b){  return ReyeshiftC(octonly(a), octonly(b)); }

eon ninReyeshiftA(eon a, eon b){  return ReyeshiftA(niner(a), niner(b)); }
eon ninReyeshiftB(eon a, eon b){  return ReyeshiftB(niner(a), niner(b)); }
eon ninReyeshiftC(eon a, eon b){  return ReyeshiftC(niner(a), niner(b)); }

eon dashL(eon x, eon y)
{
  eon z;
  z = mul(x, y);
  return add3(assoc(x, y, z), assoc(y, z, y), assoc(z, x, y));
}

eon dashR(eon x, eon y)
{
  eon z;
  z = mul(y, x);
  return add3(assoc(z, y, x), assoc(y, z, y), assoc(y, x, z));
}

eon assoclie1R(eon a, eon b)
{
  return assoc(lie(a, b), a, b);
}

eon assoclie2R(eon a, eon b)
{
  return assoc(lie(a, b), a, a);
}

eon assoclie1L(eon a, eon b)
{
  return assoc(b, a, lie(b, a));
}

eon assoclie2L(eon a, eon b)
{
  return assoc(a, a, lie(b, a));
}

eon s3(eon a, eon b, eon c)
{
  return add3(assoc(a, b, c), assoc(b, c, a), assoc(c, a, b));
}

eon sident1R(eon x, eon y, eon z)
{
  return sub(mul(s3(z, y, x), x), s3(z, mul(x, y), x));
}

eon sident2R(eon x, eon y, eon z)
{
  return sub(s3(z, mul(y, x), x), mul(x, s3(z, y, x)));
}

eon sident1L(eon x, eon y, eon z)
{
  return sub(mul(x, s3(x, y, z)), s3(x, mul(y, x), z));
}

eon sident2L(eon x, eon y, eon z)
{
  return sub(s3(x, mul(x, y), z), mul(s3(x, y, z), x));
}

eon gaR(eon a, eon b, eon c, eon d)
{
  return
      sub(add
          (assoc(mul(a, b), c, d), assoc(a, b, lie(c, d))),
          add(mul(a, assoc(b, c, d)), mul(assoc(a, c, d), b)));
}

eon gaL(eon a, eon b, eon c, eon d)
{
  return
      sub(add
          (assoc(d, c, mul(b, a)), assoc(lie(d, c), b, a)),
          add(mul(assoc(d, c, b), a), mul(b, assoc(d, c, a))));
}

eon lieassoc(eon a, eon b, eon c, eon x)
{
  return lie(assoc(a, b, c), x);
}

eon lieflex2(eon a, eon b)
{
  return lie(flex(a, b), a);
}

eon lieflex3(eon a, eon b, eon c)
{
  return lie(flex(a, b), c);
}

eon HPQ9L(eon a, eon b)
{
  return sub(assoc(sq(a), a, b), assoc(a, sq(a), b));
}

eon HPQ9R(eon a, eon b)
{
  return sub(assoc(b, sq(a), a), assoc(b, a, sq(a)));
}

eon HPQ11(eon a, eon b)
{                               /*central for quadratic algebra */
  return
      sub(add4
          (scal(2, mul(a, assoc(a, b, b))),
           scal(2, mul(b, assoc(b, a, a))), ujor(a,
                                                 flex(b, a)),
           ujor(b, flex(a, b))), scal(2,
                                      add(assoc
                                          (a, b, mul(a, b)),
                                          assoc(b, a, mul(b, a)))));
}

eon HPQ12(eon a, eon b)
{                               /*central for quadratic algebra */
  return
      sub(scal
          (2,
           add3(ujor(sq(a), sq(b)), mul(a, Jassoc(a, b, b)),
                mul(b, Jassoc(b, a, a)))), sq(ujor(a, b)));
}

eon Rasmyslov(eon a, eon b)
{  return lie(a, sqlie(a, b));
}

eon rr(eon x, eon a, eon b, eon c)
{
  return ujor(a, ujor(b, ujor(c, x)));

}

eon p3(eon a, eon b, eon c, eon x)
{
  return
      sub(add3
          (rr(x, a, b, c), rr(x, b, c, a), rr(x, c, a, b)),
          add3(rr(x, c, b, a), rr(x, a, c, b), rr(x, b, a, c)));
}

eon bigRacine(eon a, eon b, eon c, eon x)
{
  return sub(p3(a, b, c, sq(x)), ujor(p3(a, b, c, x), x));
}

eon HP23(eon a, eon b, eon c, eon d)
{
  return
      sub(add
          (scal
           (2,
            add(assoc(flex(a, b), c, d),
                assoc(c, flex(a, b), d))), assoc(lie(a, b),
                                                 ujor(a, c),
                                                 d)),
          add(assoc(ujor(lie(a, b), a), c, d),
              assoc(ujor(lie(a, b), c), a, d)));
}

eon hp26term(eon a, eon d, eon b, eon c, eon e)
{
  return add(assoc(a, Jassoc(c, d, e), b), assoc(a, b, Jassoc(c, d, e)));
}

eon HP26(eon a, eon b, eon c, eon d, eon e)
{
  return
      sub(add3
          (hp26term(a, d, e, b, c), hp26term(a, d, b, c, e),
           hp26term(a, d, c, e, b)), add3(hp26term(a, d, c,
                                                   b, e),
                                          hp26term(a, d, e,
                                                   c, b),
                                          hp26term(a, d, b, e, c)));
}


/* power-associativity for x^3. Hentzel&Peresi: J.Algebra 206(1998)1-16
* Thm1 claim all degree-3 idents for quadratic algebras come from this=0. */
eon powass3(eon x)
{
  return assoc(x, x, x);
}

/* power-associativity for x^4 */
eon powass4(eon x)
{
  return sub(mul(mul(x, mul(x, x)), x), mul(x, mul(x, mul(x, x))));
}

eon sqpowass3(eon x)
{
  return sq(powass3(x));
}

eon sqpowass4(eon x)
{
  return sq(powass4(x));
}

eon jorlies(eon a, eon b, eon c, eon d)
{
  return ujor(lie(a, b), lie(c, d));
}

eon expt31(eon a, eon b, eon c)
{
  return
      add3(add(m3(Conj(a), b, Conj(c)), m3(Conj(b), a, c)),
           add(m3(Conj(b), c, Conj(a)), m3(Conj(a), c, b)),
           add(m3(Conj(c), a, Conj(b)), m3(Conj(c), b, a)));
}

eon expt32(eon a, eon b, eon c)
{
  return
      add3(add(m3(Conj(c), b, Conj(a)), m3(Conj(c), a, b)),
           add(m3(Conj(a), c, Conj(b)), m3(Conj(b), c, a)),
           add(m3(Conj(b), a, Conj(c)), m3(Conj(a), b, c)));
}

eon expt3(eon a, eon b, eon c)
{                               /*0 for octonion a,b,c ???! */
  return add(expt31(a, b, c), expt32(a, b, c));
}

/* finds determinant of nXn matrix m, which is destroyed */
real64 det( real64 m[MAXSIZE][MAXSIZE], int n ){
  int p,pp,c,r,k,nm1, nm2;
  real64 t,mp,det,mult,piv;
  
  assert(n<=MAXSIZE);
  assert(n>0);
  if(n==1) return m[0][0];
  if(n==2) return ( m[0][0]*m[1][1] - m[0][1]*m[1][0] );
  det = 1;
  nm1 = n-1;
  nm2 = n-2;
  for(c=0; c<nm2; c++){ /* reduce column c to upper tri form */
    /* find pivot */
    pp = c;
    mp = m[c][c];
    if(mp<0) mp = -mp;    
    for(p=c+1; p<n; p++){ /* find partial pivot */
      t = m[p][c];
      if(t<0) t = -t;
      if(mp < t){
        pp = p;
        mp = t;
      }
    }
    if(mp==0) return(0);
    piv = m[pp][c];
    if(pp>c){
      det *= -piv;
      /* swap rows c, pp */
      for(k=c; k<n; k++){
        t = m[c][k];
        m[c][k] = m[pp][k];
        m[pp][k] = t;
      }
    }
    else{
      det *=  piv;
    }
    
    for(r=c+1; r<n; r++){    /* subtract mult of the pivot row from row r */
      mult = m[r][c]/piv;
      for(k=c+1; k<n; k++){
        m[r][k] -= m[c][k]*mult;
      }    
    }
  }
  return ( det * ( m[nm1][nm1] * m[nm2][nm2] - m[nm1][nm2] * m[nm2][nm1] ) );
}

basetype jacdet(eon a, eon b){
  eon y,z,q;
  uint n,j,i;
  real64 mat[MAXSIZE][MAXSIZE];
  real64 tm;
  assert(a.lev == b.lev);
  n = pow2(a.lev);
  assert(n>0);
  for(j=0; j<n; j++){
    tm = a.coord[j];
    a.coord[j] = tm + 0.0000001;
    z = mul(a, b);
    a.coord[j] = tm - 0.0000001;
    y = mul(a, b);
    a.coord[j] = tm;
    q = divscal(  0.0000002, sub(z,y)  );
    for(i=0; i<n; i++){
       mat[j][i] = q.coord[i];
    }
  }
  return det( mat, n );
}

/**************  
     [d1      a     b ]
     [                ]
     [cja    d2     c ]
     [                ]
     [cjb    cjc    d3]
**************/
basetype hermdet3(basetype d1, basetype d2, basetype d3,
                  eon a, eon b, eon c)
{
  return d1 * d2 * d3 - d1 * ss(c) - d2 * ss(b) - d3 * ss(a)
      + re(m3(Conj(a), b, Conj(c))) + re(m3(Conj(b), a, c));
}

/* det = a d - (a c) (a^{-1} b);  if a=0 then  -cb.  */
eon zassdet2(eon a, eon b, eon c, eon d)
{
  basetype s;
  s = ss(a);
  if(fabs(s) < 0.00000000000001)
    return scal(-1, mul(c, b));
  return sub(mul(a, d), divscal(s, mul(mul(a, c), mul(Conj(a), b))));
}

/* zassident ok in octonions but fails for conway hexons in all coords. */
eon zassident(eon a, eon b, eon c, eon d)
{
  eon A, B, C, D;
  A = add(mul(a, Conj(a)), mul(b, Conj(b)));
  B = add(mul(a, Conj(c)), mul(b, Conj(d)));
  C = add(mul(c, Conj(a)), mul(d, Conj(b)));
  D = add(mul(c, Conj(c)), mul(d, Conj(d)));
  /* A,B, C,D is M M^H where M=(a,b, c,d) */
  return sub(zassdet2(A, B, C, D),
             mul(zassdet2(a, b, c, d), Conj(zassdet2(a, b, c, d))));
}

/* zassext works up thru quaternions but fails for 8-ons in all coords. */
eon zassext(eon a, eon b, eon c, eon d, eon p, eon q)
{
  eon A, B, C, D, pc;
  eon A2, B2, C2, D2;
  basetype d1,d2;
  d1 = q.coord[2];
  d2 = q.coord[1];
  pc = Conj(p);
  A = add( scal(d1, a), mul(b, pc)  );
  B = add( mul(a, p),   scal(d2, b) );
  C = add( scal(d1,c),  mul(d, pc)  );
  D = add( mul(c, p),   scal(d2, d) );

  A2 = add(mul(A, Conj(a)), mul(B, Conj(b)));
  B2 = add(mul(A, Conj(c)), mul(B, Conj(d)));
  C2 = add(mul(C, Conj(a)), mul(D, Conj(b)));
  D2 = add(mul(C, Conj(c)), mul(D, Conj(d)));

  /* A2,B2, C2,D2 is M H M^H where M=(a,b, c,d) and H=(d1,p, pc,d2) */
  return sub(zassdet2(A2, B2, C2, D2),
             scal(d1*d2-ss(p), 
                  mul(zassdet2(a, b, c, d), Conj(zassdet2(a, b, c, d)))));
}

/* Works up thru octonions but fails for conway hexons. */
ubasetype sszassident(eon a, eon b, eon c, eon d)
{
  eon A, B, C, D;
  A = add(mul(a, Conj(a)), mul(b, Conj(b)));
  B = add(mul(a, Conj(c)), mul(b, Conj(d)));
  C = add(mul(c, Conj(a)), mul(d, Conj(b)));
  D = add(mul(c, Conj(c)), mul(d, Conj(d)));
  /* A,B, C,D is M M^H where M=(a,b, c,d) */
  return ss(zassdet2(A, B, C, D)) -
      ss(mul(zassdet2(a, b, c, d), Conj(zassdet2(a, b, c, d))));
}


/***********************************
Determinants of 3x3 octonion matrices.
Hans Freudenthal: Zur ebenen Oktavengeometrie,
Proc. Konink. Nederl. Akad. Wetensch. A56 (1953) 195-200;
Oktaven, Ausnahmegruppen, und Oktavengeometrie,
Geometriae Dedicata 19 (1985) 1-63 [written in 1951].
jor(A,B) = (AB+BA)/2
freud(A,B) = jor(A,B) - (Atr(B)+Btr(A))/2 + (tr(A)tr(A)-tr(jor(A,B)))
det(A) = tr(jor(freud(A,A),A))/3
det(A)I = Ijor(freud(A,A),A)
sigma(A) = tr(freud(A,A))
A^3 - (trA)A^2 + sigma(A)A - (detA)I = 0
In the 2x2 case, what are these?
p a Cb
Ca m c
b Cc n
has
det = pmn-n|a|^2-m|b|^2-p|c|^2+b.ac+Conj(b.ac)

The matrix product A B A^H is then involving
all terms being of the form  x y Conj(x)
which do not care about parenthesization?
It is claimed that 2x2 octonion matrices are
an alternative algebra, I presume under Jordan product,
regular product is not even power-associative.
No sorry, there also are cross terms x y z.
And if B is hermitian so is this.
(In the hexons this "reflex" property is no longer true, sadly.)
***********************************************/

basetype freud(eon a, eon b, eon c, eon pmn){
  basetype p,m,n;
  eon z;
  p = pmn.coord[1];
  m = pmn.coord[2];
  n = pmn.coord[3];
  z = add( mul(b,mul(a,c)), Conj(mul(b,mul(a,c))) );
  return (  z.coord[0] + p*m*n - ( n*ss(a)+m*ss(b)+p*ss(c) )  );
}

eon jordiJ(eon x, eon y)
{                               /* xy.xx = x(y.xx) */
  eon sx;
  sx = sq(x);
  assert(x.lev == y.lev);
  return sub(jor(jor(x, y), sx), jor(x, jor(y, sx)));
}

eon jordi1(eon x)
{
  return jordiJ(x, x);
}

eon jordiMR(eon x, eon y)
{                               /* xy.xx = x(y.xx)  if x pureimag */
  eon sx;
  sx = sq(x);
  assert(x.lev == y.lev);
  return assoc(x,y,sx);
  /*  return sub(mul(mul(x, y), sx), mul(x, mul(y, sx)));  */
}

eon jordiMRN(eon x, eon y)
{                              
  eon sx;
  sx = scalar(x.lev, ss(x));
  assert(x.lev == y.lev);
  return assoc(x,y,sx);
}

eon jordiMRR(eon x, eon y)
{                               /*TRUE for hexons if x pureimag*/
  eon sx;
  basetype rex;
  rex = re(x);
  sx = scal(2*rex, x);
  assert(x.lev == y.lev);
  return assoc(x,y,sx);
}

eon jordiMLR(eon x, eon y) /* other diriction jordiMRR */
{                      
  eon sx;
  basetype rex;
  rex = re(x);
  sx = scal(2*rex, x);
  assert(x.lev == y.lev);
  return assoc(sx,y,x);
}

eon jordiMRRcons2(eon x, eon y)
{                               /*FALSE for hexons (equiv. flex)*/
  eon sx;
  basetype rex;
  rex = 2;
  sx = scal(rex,x);
  assert(x.lev == y.lev);
  return assoc(x,y,sx);
}

eon jordiMRRcons3(eon x, eon y)
{                       /*TRUE for hexons if x pureimag*/
  eon sx;
  basetype rex;
  rex = re(x);
  sx = scal(rex, x);
  assert(x.lev == y.lev);
  assert(sx.lev == y.lev);
  return assoc(x,y,sx);
}

eon jordiML(eon x, eon y) /*Schafer p.141*/
{                               /* xx.yx = (xx.y)x  */
  eon sx;
  sx = sq(x);
  assert(x.lev == y.lev);
  return assoc(sx,y,x);
  /*  return sub(mul(sx, mul(y, x)), mul(mul(sx, y), x)); */
}

eon jordiMRlinA(eon w, eon x, eon y, eon z) /* Schafer p140, p91 */
{               /*true for octonions. Conseq of linearizing jordiMR.*/  
  eon sx;
  sx = sq(x);
  assert(x.lev == y.lev);
  return add( assoc(x,y,ujor(z,x)), assoc(z,y,sx) );
}

eon jordiMRlin(eon w, eon x, eon y, eon z)
{     /* incorrectly stated Schafer p140.
       * As corrected here becomes true for octonions.
       * w=1 in a linear algebra in which this holds ==> flexibility.  */
  assert(x.lev == y.lev);
  return
    add3( assoc(x,y,ujor(w,z)), assoc(w,y,ujor(z,x)), assoc(z,y,ujor(w,x)) );
/* Schafer p140 wrongly wanted
* add3( assoc(x,y,mul(w,z)), assoc(w,y,mul(z,x)), assoc(z,y,mul(x,w)) ); */
}

eon jordiMLlin(eon w, eon x, eon y, eon z)
{
  assert(x.lev == y.lev);
  return
    add3( assoc(ujor(z,w),y,x), assoc(ujor(x,z),y,w), assoc(ujor(w,x),y,z) );
}

eon jordiMR1(eon x, eon y) /*Schafer p.141*/
{                               /* xx.xy = x(xx.y) */
  eon sx;
  sx = sq(x);
  assert(x.lev == y.lev);
  return sub(mul(sx, mul(x, y)), mul(x, mul(sx,y)));
}

eon jordiML1(eon x, eon y) /*Schafer p.141*/
{                               /* yx.xx = (y.xx)x */
  eon sx;
  sx = sq(x);
  assert(x.lev == y.lev);
  return sub(mul(mul(y, x), sx), mul(mul(y,sx),x));
}

/* Schafer p140 shows JordanIdent (if unit 1 exists) ==> Flexible by
1. linearize,
      (x,y,wz) + (w,y,zx) + (z,y,xw) = 0  --incorrect, but ok if prod is ujor
2. use fact we have unit 1 make w=1, 3. done.
This required distributivity.
Hexons do not do this.  
Right-linearity implies
xy.(a+b)=xy.a+xy.b
x.y(a+b)=x.(ya+yb)=x.ya+x.yb
But quadratic ident implies x^2 = scalar+linear = norm(x) + 2*re(x)*x
so
xy.x^2 = xy.scalar + xy.linear
x.yx^2 = x.yscalar + x.ylinear
so JordiMR, quadratic ident, and right-linearity implies flexibility???!!!

He then says flexible + any of the following ==> NonCommut Jordan:
1.  x^2 y.x = x^2 . yx,  2.  x^2 . xy = x.x^2 y,  3. yx.x^2 = y x^2.x.
He then says  tr( xy.z ) = tr( x.yz )  and <xy,z> = <x,yz>.
Hexons do not do this.                  Hexons do not do this.
*********************************************/

eon JacLie(eon x, eon y, eon z)
{
  return add3(lie(lie(x, y), z), lie(lie(z, x), y), lie(lie(y, z), x));
}

eon sqJacLie(eon x, eon y, eon z)
{
  return sq(JacLie(x, y, z));
}

eon JacLie21(eon x, eon z)
{
  return JacLie(x, x, z);
}

eon JacMul(eon x, eon y, eon z)
{
  return add3(mul(mul(x, y), z), mul(mul(z, x), y), mul(mul(y, z), x));
}

eon JacI(eon x, eon y, eon z)
{
  return add3(imul(imul(x, y), z), imul(imul(z, x), y),
              imul(imul(y, z), x));
}

eon MalcevLie(eon x, eon y, eon z)
{
  return sub(lie(lie(x, y), lie(x, z)),
             add3(lie(lie(lie(x, y), z), x),
                  lie(lie(lie(y, z), x), x), lie(lie(lie(z, x), x), y)));
}

eon MalcevMul(eon x, eon y, eon z)
{
  return sub(mul(mul(x, y), mul(x, z)),
             add3(mul(mul(mul(x, y), z), x),
                  mul(mul(mul(y, z), x), x), mul(mul(mul(z, x), x), y)));
}

eon MalcevI(eon x, eon y, eon z)
{
  return sub(imul(imul(x, y), imul(x, z)),
             add3(imul(imul(imul(x, y), z), x),
                  imul(imul(imul(y, z), x), x),
                  imul(imul(imul(z, x), x), y)));
}

eon cycassI(eon x, eon y, eon z)
{
  return add3(imul(imul(x, y), z),
              imul(imul(y, z), x),
              imul(imul(z, x), y));
}

eon cycassL(eon x, eon y, eon z)
{
  return add3(lie(lie(x, y), z),
              lie(lie(y, z), x),
              lie(lie(z, x), y));
}

eon varMal1I(eon x, eon y, eon z)
{
  return sub(cycassI(x, y, imul(x, z)),
             imul(cycassI(x, y, z), x) );
}

eon varMal2I(eon x, eon y, eon z)
{
  return sub(cycassI(x, imul(x, y), z), imul(cycassI(x, y, z), x));
}

eon varMal3I(eon x, eon y, eon z)
{
  return sub(cycassI(x, imul(x, y), z), cycassI(x, y, imul(x, z)));
}

eon cycassLie(eon x, eon y, eon z)
{
  return add3(lie(lie(x, y), z), lie(lie(y, z), x), lie(lie(z, x), y));
}

eon varMal1Lie(eon x, eon y, eon z)
{
  return sub(cycassL(x, y, lie(x, z)), lie(cycassL(x, y, z), x));
}

eon varMal2Lie(eon x, eon y, eon z)
{
  return sub(cycassL(x, lie(x, y), z), lie(cycassL(x, y, z), x));
}

eon varMal3Lie(eon x, eon y, eon z)
{
  return sub(cycassL(x, lie(x, y), z), cycassL(x, y, lie(x, z)));
}


/*********************
"determinant" (sum of 3! terms) of
a b c
a b c
a b c
*******************/
eon sym3ident(eon a, eon b, eon c)
{
  return sub(add3(m3(a, b, c), m3(b, c, a), m3(c, a, b)),
             add3(m3(c, b, a), m3(a, c, b), m3(b, a, c)));
}

/******************
"determinant" (sum of 4! terms) of
a b c d
a b c d
a b c d
a b c d
is 0 for quaternion a,b,c,d.
****************************/
eon sym4ident(eon a, eon b, eon c, eon d)
{
  return
      add5(add5
           (m4(a, b, c, d), neg(m4(a, b, d, c)),
            neg(m4(a, c, b, d)), m4(a, c, d, b), m4(a, d, b,
                                                    c)),
           add5(neg(m4(a, d, c, b)), neg(m4(b, a, c, d)),
                m4(b, a, d, c), m4(b, c, a, d),
                neg(m4(b, c, d, a))),
           add5(neg(m4(b, d, a, c)), m4(b, d, c, a),
                m4(c, a, b, d), neg(m4(c, a, d, b)),
                neg(m4(c, b, a, d))), add5(m4(c, b, d, a),
                                           m4(c, d, a, b),
                                           neg(m4
                                               (c, d, b, a)),
                                           neg(m4
                                               (d, a, b, c)),
                                           m4(d, a, c, b)),
           add4(m4(d, b, a, c), neg(m4(d, b, c, a)),
                neg(m4(d, c, a, b)), m4(d, c, b, a)));
}

eon crosstest(eon x, eon y)
{
  return sub(imul(x, y), cross(x, y));
}

/* uX(vXw) = <u,w>v - <u,v>w  Ebbingp199 */
eon Grassman1(eon x, eon y, eon z)
{
  return
      sub(add(imul(x, imul(y, z)), scal(ip(x, y), z)), scal(ip(x, z), y));
}

/* uX(vXw) = <u,w>v - <u,v>w   allows us to deduce
* (uXv)Xw = - wX(uXv) = - <w,v>u + <w,u>v
* (uXv)Xw - (uXv)Xw = - <w,v>u + <w,u>v - <u,w>v + <u,v>w
* (uXv)Xw - (uXv)Xw = <u,v>w - <w,v>u
******************************************/
eon Grassman2(eon x, eon y, eon z)
{
  return sub(add(iassoc(x, y, z), scal(ip(z, y), x)), scal(ip(x, y), z));
}

eon Grassman1C(eon x, eon y, eon z)
{
  return
      sub(add(cross(x, cross(y, z)), scal(ip(x, y), z)),
          scal(ip(x, z), y));
}

eon Grassman2C(eon x, eon y, eon z)
{
  return sub(add(Cassoc(x, y, z), scal(ip(z, y), x)), scal(ip(x, y), z));
}

eon crossanti(eon x, eon y)
{
  return add(cross(x, y), cross(y, x));
}

eon imulanti(eon x, eon y)
{
  return add(imul(x, y), imul(y, x));
}

basetype crossxchg(eon w, eon x, eon y)
{
  return ip(w, cross(x, y)) - ip(x, cross(y, w));
}

basetype imulxchg(eon w, eon x, eon y)
{
  return ip(w, imul(x, y)) - ip(x, imul(y, w));
}

/* Moufang identity   xy.ax = (x.ya)x */
eon MoufmidR(eon x, eon y, eon a)
{
  return sub(mul(mul(x, y), mul(a, x)), mul(mul(x, mul(y, a)), x));
}

/* Moufang identity   xy.ax = x(ya.x) */
eon MoufmidL(eon x, eon y, eon a)
{
  return sub(mul(mul(x, y), mul(a, x)), mul(x, mul(mul(y, a), x)));
}

/* Middle Moufang identity   (xy.a)x = x(ya.x) */
eon MoufmidM(eon x, eon y, eon a)
{
  return sub(mul(mul(x, mul(y, a)), x), mul(x, mul(mul(y, a), x)));
}

eon pflM2(eon x, eon y, eon a)
{
    eon z;
    z = sq(x);
  return sub(mul(mul(x, mul(y, a)), z), mul(x, mul(mul(y, a), z)));
}

eon pflM3(eon x, eon y, eon a)
{
    eon z;
    z = cube(x);
  return sub(mul(mul(x, mul(y, a)), z), mul(x, mul(mul(y, a), z)));
}

eon pflM4(eon x, eon y, eon a)
{
    eon z;
    z = sq(sq(x));
  return sub(mul(mul(x, mul(y, a)), z), mul(x, mul(mul(y, a), z)));
}

eon pflM5(eon x, eon y, eon a)
{
    eon z;
    z = mul(x,sq(sq(x)));
  return sub(mul(mul(x, mul(y, a)), z), mul(x, mul(mul(y, a), z)));
}

eon pflM6(eon x, eon y, eon a)
{
    eon z;
    z = cube(sq(x));
  return sub(mul(mul(x, mul(y, a)), z), mul(x, mul(mul(y, a), z)));
}

eon pflM7(eon x, eon y, eon a)
{
    eon z;
    z = mul(x,cube(sq(x)));
  return sub(mul(mul(x, mul(y, a)), z), mul(x, mul(mul(y, a), z)));
}

/* Fenyves Extra identity   x(y.zx) = ((xy.z)x 
 * fails in all coords in 2^n-ons with n>=4.
 * Works if n<=2 in all coords.
 * If n=3 then it works in the real coord only (all coords if pure imag octons)
 *******************************************/
eon FenyvesExtra(eon x, eon y, eon z)
{
    return sub( mul(x, mul(y, mul(z, x))), mul(mul(mul(x,y), z), x) );
}

/* All 3 nucsq identities hold for pure-imaginary 2^n-ons for all n
 * but for fully general ones hold only for quaternions & below, (and
 * octonions also if only require validity in real part only)  */
/* Left nuclear square identity   xx.yz = (xx.y)z */
eon Lnucsq(eon x, eon y, eon z)
{
   return sub( mul(mul(x,x), mul(y,z)), mul(mul(mul(x,x),y),z) );
}

/* Right nuclear square identity   x(y.zz) = xy.zz */
eon Rnucsq(eon x, eon y, eon z)
{
   return sub( mul(mul(x,y), mul(z,z)), mul(x,mul(y,mul(z,z))) );
}

/* Middle nuclear square identity   x(yy.z) = (x.yy)z */
eon Mnucsq(eon x, eon y, eon z)
{
   return sub( mul(x,mul(mul(y,y),z)), mul(mul(x,mul(y,y)),z) );
}


eon W1a(eon x, eon y, eon z) /* W1:   $zx \cdot (yxy) = z(xyx) \cdot y$  */
{
    return sub( mul(mul(z,x), mul(mul(y,x), y)) , 
                mul(mul(z, mul(mul(x,y),x)), y) );
}

eon W1b(eon x, eon y, eon z) /* W1:   $zx \cdot (yxy) = z(xyx) \cdot y$  */
{
    return sub( mul(mul(z,x), mul(y,mul(x,y)) ) , 
                mul(mul(z, mul(x,mul(y,x))), y) );
}

eon W1c(eon x, eon y, eon z) /* W1:   $zx \cdot (yxy) = z(xyx) \cdot y$  */
{
    return sub( mul(mul(z,x), mul(y,mul(x,y)) ) , 
                mul(mul(z, mul(mul(x,y),x)), y) );
}

eon W1d(eon x, eon y, eon z) /* W1:   $zx \cdot (yxy) = z(xyx) \cdot y$  */
{
    return sub( mul(mul(z,x), mul(mul(y,x), y)), 
                mul(mul(z, mul(x,mul(y,x))), y) );
}

eon W2a(eon x, eon y, eon z) /* W2: $(yxy) \cdot xz = y \cdot (xyx) z$ */
{
   return sub( mul( mul(mul(y,x),y), mul(x,z) ) ,
               mul(y, mul(mul(mul(x,y),x),z) ) );
}

eon W2b(eon x, eon y, eon z) /* W2: $(yxy) \cdot xz = y \cdot (xyx) z$ */
{
   return sub( mul( mul(y,mul(x,y)), mul(x,z) ) ,
               mul(y, mul(mul(x,mul(y,x)),z) ) );
}

eon W2c(eon x, eon y, eon z) /* W2: $(yxy) \cdot xz = y \cdot (xyx) z$ */
{
   return sub( mul( mul(y,mul(x,y)), mul(x,z) ) ,
               mul(y, mul(mul(mul(x,y),x),z) ) );
}

eon W2d(eon x, eon y, eon z) /* W2: $(yxy) \cdot xz = y \cdot (xyx) z$ */
{
   return sub( mul( mul(mul(y,x),y), mul(x,z) ) ,
               mul(y, mul(mul(x,mul(y,x)),z) ) );
}

/* non-loop quasigroup versions: ??? */
/* LG1: x(y.zz) = (x.yz)z */
/* LG2: xy.zz   = (x.yz)z */
/* LG3: x(y.zy) = (x.yz)y */
/* LC2: x(x.yz) = (x.xy)z */
/* LC3: x(x.yz) = (xx.y)z */
/* LC4: x(y.yz) = (x.yy)z */

/* RG1: x(xy.z) = (xx.y)z */
/* RG2: x(x.yz) = xx.yz   true by left-alt */
/* RG3: x(yx.z) = (xy.x)z */
/* RC2: x(yz.z) = (xy.z)z */
/* RC3: x(y.zz) = (xy.z)z */
/* RC4: x(yy.z) = (xy.y)z */

/* LC, LC1, LC2, LC3: work in pure-imag 2^n-ons, but only work
   in real part of general octonions, and fail in general 16-ons. */
eon LC2(eon x, eon y, eon z) /* LC2: x(x.yz) = (x.xy)z */
{
    return sub( mul(x, mul(x, mul(y,z))),  mul(mul(x,mul(x,y)),z) );
}

eon LC3(eon x, eon y, eon z) /* LC3: x(x.yz) = (xx.y)z */
{
    return sub( mul(x, mul(x, mul(y,z))),  mul(mul(mul(x,x),y),z) );
}

eon LC4(eon x, eon y, eon z) /* LC4: x(y.yz) = (x.yy)z */ 
{
    return sub( mul(x, mul(y, mul(y,z))),  mul(mul(x,mul(y,y)),z) );
}

/* LG1,LG2,LG3 only ok in real coord in octonions, total falsity if n>=4. */
eon LG1(eon x, eon y, eon z)  /* LG1: x(y.zz) = (x.yz)z */
{
    return sub( mul(x, mul(y, mul(z,z))),  mul(mul(x,mul(y,z)),z) );
}

eon LG2(eon x, eon y, eon z) /* LG2: xy.zz   = (x.yz)z */
{
    return sub( mul(mul(x,y), mul(z,z)),  mul(mul(x,mul(y,z)),z) );
}

eon LG3(eon x, eon y, eon z) /* LG3: x(y.zy) = (x.yz)y */
{
    return sub( mul(x, mul(y, mul(z,y))),  mul(mul(x,mul(y,z)),y) );
}

/* RG1: x(xy.z) = (xx.y)z */
/* RG2: x(x.yz) = xx.yz   */
/* RG3: x(yx.z) = (xy.x)z */
/* RC2: x(yz.z) = (xy.z)z */
/* RC3: x(y.zz) = (xy.z)z */
/* RC4: x(yy.z) = (xy.y)z */

eon RG1(eon x, eon y, eon z) /* RG1:  x(xy.z) = (xx.y)z */
{
   return sub( mul(x, mul(mul(x,y),z)),  mul(mul(mul(x,x),y),z) );
}

eon RG2(eon x, eon y, eon z) /* RG2: x(x.yz) = xx.yz   */ 
{
   return sub( mul(x, mul(x,mul(y,z))),  mul(mul(x,x), mul(y,z)) );
}

eon RG3(eon x, eon y, eon z) /* RG3: x(yx.z) = (xy.x)z */
{
   return sub( mul(x, mul(mul(y,x),z)),  mul(mul(mul(x,y),x),z) );
}

eon RC2(eon x, eon y, eon z) /* RC2: x(yz.z) = (xy.z)z */ 
{
    return sub( mul(x, mul(mul(y,z),z)),  mul(mul(mul(x,y),z),z) );
}

eon RC3(eon x, eon y, eon z) /* RC3: x(y.zz) = (xy.z)z */ 
{
    return sub( mul(x, mul(y,mul(z,z))),  mul(mul(mul(x,y),z),z) );
}

eon RC4(eon x, eon y, eon z) /* RC4: x(yy.z) = (xy.y)z */
{
    return sub( mul(x, mul(mul(y,y),z)),  mul(mul(mul(x,y),y),z) );
}

eon RWIP(eon x, eon y)
{
    return sub( scal(ss(y), Conj(x)),  mul(y, Conj(mul(x,y))) );
}

eon LWIP(eon x, eon y)
{
    return sub( scal(ss(y), Conj(x)),  mul(Conj(mul(y,x)), y) );
}

/* Cloop identity   x(y.yz) = (xy.z)y 
 * This and RCloop totally fail in the 16-ons and beyond (even if pure imag);
 * works for n<=2, and works for octonions if real part only needed
 * or if pure imag octonions */
eon Cloop(eon x, eon y, eon z)
{
    return sub( mul(x,mul(y,mul(y,z))), mul(mul(mul(x,y),y),z) );
}

/* RCloop identity   x(yz.z) = xy.zz */
eon RCloop(eon x, eon y, eon z)
{
    return sub( mul(x,mul(mul(y,z),z)), mul(mul(x,y), mul(z,z)) );
}

/* LCloop identity   xx.yz = (x.xy)z 
 * Holds for pure imaginary 2^n-ons, all n.
 * If fully general 2^n-ons, then fails in the 16-ons and beyond,
 * works for n<=2, and works for octonions if real part only */
eon LCloop(eon x, eon y, eon z)
{
    return sub( mul(mul(x,x), mul(y,z)), mul(mul(x,mul(x,y)),z) );
}

/* left Bol identity   (x.yx)z = x(y.xz) */
eon Lbol(eon x, eon y, eon z)
{
  return sub(mul(mul(x, mul(y, x)), z), mul(x, mul(y, mul(x, z))));
}

/* right Bol identity  z(xy.x)=(zx.y)x */
eon Rbol(eon x, eon y, eon z)
{
  return sub(mul(z, mul(mul(x, y), x)), mul(mul(mul(z, x), y), x));
}

/* left Moufang identity   (xy.x)z = x(y.xz) */
eon Lmouf(eon x, eon y, eon z)
{
  return sub(mul(mul(mul(x, y), x), z), mul(x, mul(y, mul(x, z))));
}

/* right Moufang identity  z(x.yx)=(zx.y)x */
eon Rmouf(eon x, eon y, eon z)
{
  return sub(mul(z, mul(x, mul(y, x))), mul(mul(mul(z, x), y), x));
}

/*  R.Moufang: Math Ann 110 (1934/5) 416-430  * p420eq10 & relatives */
eon mouf5L(eon x, eon y, eon z)
{                               /* (x.yConj(z))z = xConj(z).zy  */
  return sub(mul(mul(x, mul(y, Conj(z))), z),
             mul(mul(x, Conj(z)), mul(z, y)));
}

eon mouf5R(eon x, eon y, eon z)
{
  return sub(mul(z, mul(mul(Conj(z), y), x)),
             mul(mul(y, z), mul(Conj(z), x)));
}

/*  R.Moufang: Math Ann 110 (1934/5) 416-430  * p420eq10 & relatives */
eon mouf6L(eon x, eon y, eon z)
{                               /* Conj(z)(zx.y) = xConj(z).zy  */
  return sub(mul(Conj(z), mul(mul(z, x), y)),
             mul(mul(x, Conj(z)), mul(z, y)));
}

eon mouf6R(eon x, eon y, eon z)
{
  return sub(mul(mul(y, mul(x, z)), Conj(z)),
             mul(mul(y, z), mul(Conj(z), x)));
}

eon mouf7L(eon x, eon y, eon z)
{                               /* (x.yConj(z))z = xConj(z).zy  */
  return sub(mul(mul(x, mul(y, Conj(z))), z),
             mul(Conj(z), mul(mul(z, x), y)));
}

eon mouf7R(eon x, eon y, eon z)
{
  return sub(mul(z, mul(mul(Conj(z), y), x)),
             mul(mul(y, mul(x, z)), Conj(z)));
}

/**************************************???
Variant-Moufang laws:
(x \cdot yx^{-1})z = x(y \cdot x^{-1}z) = (xy \cdot x^{-1})z
     %var left Moufang/Bol

z(xy \cdot x^{-1}) = (zx \cdot y)x^{-1} = z(x \cdot yx^{-1})
   %var right Moufang/Bol

%each of the flexators can be done in 2 ways...

xy \cdot z = x z^{-1} \cdot z y z = x z \cdot z^{-1} y z
= x(y x^{-1} \cdot x z)

z \cdot yx = z y z \cdot z^{-1} x = z y z^{-1} \cdot z x
= (z x \cdot x^{-1} y)x
********************************************/

eon varMouf2R(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(scal(ss(z), mul(mul(x, y), z)),
             mul(mul(x, z), mul(Conj(z), mul(y, z))));
}

eon varMouf2Ra(eon x, eon y, eon z)
{                               /* xy.z = xInv(z).zyz */
  return sub(scal(ss(z), mul(mul(x, y), z)),
             mul(mul(x, Conj(z)), mul(z, mul(y, z))));
}

eon varMouf2Rc(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(scal(ss(z), mul(mul(x, y), z)),
             mul(mul(x, z), mul(mul(Conj(z), y), z)));
}

eon varMouf2Rd(eon x, eon y, eon z)
{                               /* xy.z = xInv(z).zyz */
  return sub(scal(ss(z), mul(mul(x, y), z)),
             mul(mul(x, Conj(z)), mul(mul(z, y), z)));
}

/* works for 2^n-ons! */
eon varMouf2Re(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(mul(mul(x, Conj(z)), mul(z, mul(y, z))),
             mul(mul(x, z), mul(Conj(z), mul(y, z))));
}

eon varMouf2Le(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(mul(mul(mul(z, y), z), mul(Conj(z), x)),
             mul(mul(mul(y, z), Conj(z)), mul(z, x)));
}

/*** do L equivalents??? */

eon varMouf3R(eon a, eon b, eon c)
{                               /* (b.ca)Inv(a) = ba.Inv(a)c */
  return sub(mul(mul(b, mul(c, a)), Conj(a)),
             mul(mul(b, a), mul(Conj(a), c)));
}

eon varMouf2L(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(scal(ss(z), mul(z, mul(y, x))),
             mul(mul(mul(z, y), Conj(z)), mul(z, x)));
}

eon varMouf3L(eon a, eon b, eon c)
{                               /* (b.ca)Inv(a) = ba.Inv(a)c */
  return sub(mul(Conj(a), mul(mul(a, c), b)),
             mul(mul(c, Conj(a)), mul(a, b)));
}

eon varMouf2Rv(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(scal(ss(z), mul(mul(x, y), z)),
             mul(mul(x, z), mul(mul(Conj(z), y), z)));
}

eon varMouf2Lv(eon x, eon y, eon z)
{                               /* xy.z = xz.Inv(z)yz */
  return sub(scal(ss(z), mul(z, mul(y, x))),
             mul(mul(z, mul(y, Conj(z))), mul(z, x)));
}

/* Middle Moufang identity   xy.aX = (x.ya)X */
eon MoufmidRC(eon x, eon y, eon a)
{
  return sub(mul(mul(x, y), mul(a, Conj(x))),
             mul(mul(x, mul(y, a)), Conj(x)));
}

/* Middle Moufang identity   xy.aX = x(ya.X) */
eon MoufmidLC(eon x, eon y, eon a)
{
  return sub(mul(mul(x, y), mul(a, Conj(x))),
             mul(x, mul(mul(y, a), Conj(x))));
}

/* Middle Moufang identity   (xy.a)X = x(ya.X) */
eon MoufmidMC(eon x, eon y, eon a)
{
  return sub(mul(mul(x, mul(y, a)), Conj(x)),
             mul(x, mul(mul(y, a), Conj(x))));
}

/* left Bol identity   (x.yX)z = x(y.Xz) where X=Conj(x) */
eon LbolC(eon x, eon y, eon z)
{
  return sub(mul(mul(x, mul(y, Conj(x))), z),
             mul(x, mul(y, mul(Conj(x), z))));
}

/* right Bol identity  z(xy.X)=(zx.y)X  where X=conj(x)*/
eon RbolC(eon x, eon y, eon z)
{
  return sub(mul(z, mul(mul(x, y), Conj(x))),
             mul(mul(mul(z, x), y), Conj(x)));
}

eon sqLbol(eon x, eon y, eon z)
{
  return sq(Lbol(x, y, z));
}

eon sqRbol(eon x, eon y, eon z)
{
  return sq(Rbol(x, y, z));
}

eon sqMoufmidL(eon x, eon y, eon z)
{
  return sq(MoufmidL(x, y, z));
}

eon sqMoufmidR(eon x, eon y, eon z)
{
  return sq(MoufmidR(x, y, z));
}

eon sqMoufmidM(eon x, eon y, eon z)
{
  return sq(MoufmidM(x, y, z));
}

eon comm3(eon x, eon y, eon z)
{
  return sub(mul(mul(x, y), z), mul(mul(y, x), z));
}

eon comm3A(eon x, eon y, eon z)
{
  return add3(mul(mul(x, y), z), mul(mul(y, x), z), scal(2 * ip(x, y), z));
}

eon comm3B(eon x, eon y, eon z)
{
  return add3(mul(y, mul(x, z)), mul(y, mul(z, x)), scal(2 * ip(x, z), y));
}

basetype a4(eon x, eon y, eon z)
{
  return re(mul(x, mul(y, z))) + re(mul(z, mul(y, x)))
      - re(mul(mul(x, y), z)) - re(mul(mul(z, y), x));
}

/*  Ebb p276 has gram law
|aXb|^2 = |a|^2 |b|^2 - <a,b>^2 = det( gram matrix of a,b inner prods )
************/
basetype gramC(eon a, eon b)
{
  basetype t;
  t = ip(a, b);
  return ss(cross(a, b)) - ss(a) * ss(b) + t * t;
}

/* gramI still works for pure imag a,b even at 32-ons and 64-ons;
 * gramC stops working at 16-ons */
basetype gramI(eon a, eon b)
{
  basetype t;
  t = ip(a, b);
  return ss(imul(a, b)) - ss(a) * ss(b) + t * t;
}

/*  aX(aXb) = <a,b>a - |a|^2 b  */
eon aXaXb(eon a, eon b)
{
  return sub(cross(a, cross(a, b)),
             sub(scal(ip(a, b), a), scal(ss(a), b)));
}

/*  aX(aXb) = <a,b>a - |a|^2 b  */
eon aXaXbI(eon a, eon b)
{
  return sub(imul(a, imul(a, b)),
             sub(scal(ip(a, b), a), scal(ss(a), b)));
}

/* aX(bXa) = <a,a> b - <a,b> a */
eon aXbXaI(eon a, eon b){
    return sub( imul(a, imul(b,a)), 
                sub( scal(ss(a), b), scal(ip(a,b), a) ));
}

/*  (aXb)Xa = <a,a>b - <a,b>a  */
eon aXbXaI2(eon a, eon b){
    return sub( imul(imul(a,b),a),
		sub( scal(ss(a),b), scal(ip(a,b),a) ));
}

/*  (aXb)Xa = <a,b>a - |a|^2 b  */
eon bXaXa(eon a, eon b)
{
  return sub(cross(cross(b, a), a),
             sub(scal(ip(a, b), a), scal(ss(a), b)));
}

/*  (aXb)Xa = <a,b>a - |a|^2 b  */
eon bXaXaI(eon a, eon b)
{
  return sub(imul(imul(b, a), a),
             sub(scal(ip(a, b), a), scal(ss(a), b)));
}

eon cyc3(eon x, eon y, eon z)
{
  return sub(mul(mul(x, y), z), mul(mul(z, x), y));
}

eon sqcyc3(eon x, eon y, eon z)
{
  return sq(cyc3(x, y, z));
}

eon qbarq(eon q)
{
  return mul(Conj(q), q);
}

eon antiaut(eon a, eon b)
{
  return sub(mul(Conj(a), Conj(b)), Conj(mul(b, a)));
}

eon autinv(eon a, eon b)
{
  return sub(mul(Conj(b), Conj(a)), Conj(mul(b, a)));
}

eon octantiaut(eon a, eon b)
{
  return antiaut(a, niner(b));
}

eon ninantiaut(eon a, eon b)
{
  return antiaut(a, niner(b));
}

eon antiaut1(eon a)
{
  return antiaut(a, a);
}

eon Cantiaut(eon a)
{
  return antiaut(a, Conj(a));
}

eon sqantiaut1(eon a)
{
  return sq(antiaut1(a));
}

eon sqantiaut(eon a, eon b)
{
  return sq(antiaut(a, b));
}

eon sqCantiaut(eon a)
{
  return sq(Cantiaut(a));
}

eon Rdistrib(eon a, eon b, eon c)
{                               /* a(b+c)=ab+ac */
  return sub(mul(a, add(b, c)), add(mul(a, b), mul(a, c)));
}

eon Ldistrib(eon a, eon b, eon c)
{                               /* (b+c)a=ba+ca */
  return sub(mul(add(b, c), a), add(mul(b, a), mul(c, a)));
}

eon RdistribPow(eon b)
{
  return Rdistrib(b, b, sq(b));
}

eon LdistribPow(eon b)
{
  return Ldistrib(b, b, sq(b));
}

eon RdistribPow3(eon b)
{
  return Rdistrib(b, b, cube(b));
}

eon LdistribPow3(eon b)
{
  return Ldistrib(b, b, cube(b));
}

eon Rdistrib2(eon a, eon b)
{
  return Rdistrib(a, sq(b), b);
}

eon Ldistrib2(eon a, eon b)
{
  return Ldistrib(a, sq(b), b);
}

eon Rdistrib3(eon a, eon b)
{
  return Rdistrib(a, sq(a), b);
}

eon Ldistrib3(eon a, eon b)
{
  return Ldistrib(a, sq(a), b);
}

eon Rdistrib4(eon a, eon b)
{
  return Rdistrib(sq(a), a, b);
}

eon Ldistrib4(eon a, eon b)
{
  return Ldistrib(sq(a), a, b);
}

eon Rdistrib2c(eon a, eon b)
{
  return Rdistrib(a, cube(b), b);
}

eon Ldistrib2c(eon a, eon b)
{
  return Ldistrib(a, cube(b), b);
}

eon Rdistrib3c(eon a, eon b)
{
  return Rdistrib(a, cube(a), b);
}

eon Ldistrib3c(eon a, eon b)
{
  return Ldistrib(a, cube(a), b);
}

eon Rdistrib4c(eon a, eon b)
{
  return Rdistrib(cube(a), a, b);
}

eon Ldistrib4c(eon a, eon b)
{
  return Ldistrib(cube(a), a, b);
}

eon Rdistrib5(eon a, eon b)
{
  return Rdistrib(mul(a, b), a, b);
}

eon Ldistrib5(eon a, eon b)
{
  return Ldistrib(mul(a, b), a, b);
}

eon Rdistrib6(eon a, eon b)
{
  return Rdistrib(a, mul(b, a), b);
}

eon Ldistrib6(eon a, eon b)
{
  return Ldistrib(a, mul(b, a), b);
}

eon Rdistrib7(eon a, eon b)
{
  return Rdistrib(a, mul(a, b), b);
}

eon Ldistrib7(eon a, eon b)
{
  return Ldistrib(a, mul(a, b), b);
}

eon Rdistribimh(eon z, eon b, eon c)
{
  eon a;
  a = imhalfie(z);
  return sub(mul(a, add(b, c)), add(mul(a, b), mul(a, c)));
}

eon Ldistribimh(eon z, eon b, eon c)
{
  eon a;
  a = imhalfie(z);
  return sub(mul(add(b, c), a), add(mul(b, a), mul(c, a)));
}

eon Rdistribreh(eon z, eon b, eon c)
{
  eon a;
  a = rehalfie(z);
  return sub(mul(a, add(b, c)), add(mul(a, b), mul(a, c)));
}

eon Ldistribreh(eon z, eon b, eon c)
{
  eon a;
  a = rehalfie(z);
  return sub(mul(add(b, c), a), add(mul(b, a), mul(c, a)));
}


/* EQ 49, Corinne+Manogue, Mod Phys Lett A14 (1999) 1243-1256.  */
eon G2aut49(eon a, eon b, eon y){
  eon ab;
  basetype sab;
  ab = mul(a,b);
  sab = ss(ab);
  return divscal(sab*sab,
     m3( Conj(ab), m3(b, m3(a,y,Conj(a)), Conj(b)), ab));
}

eon G2aut50(eon a, eon b, eon y){
  eon ab;
  basetype sab;
  ab = mul(a,b);
  sab = ss(ab);
  return divscal( sab*sab*sab,
     m3( Conj(ab), m3(b, m3(a,y,sq(a)), sq(b)), sq(Conj(ab))) );
}

eon G2aut51(eon a, eon b, eon y){
  eon ab;
  basetype sab;
  ab = mul(a,b);
  sab = ss(ab);
  return divscal( sab*sab*sab,
    m3( sq(ab), m3(sq(Conj(b)), m3(sq(Conj(a)),y,Conj(a)), Conj(b)), ab ) );
}

eon G2aut57(eon a, eon f, eon b, eon y){
  eon c,d,L;
  basetype sn;
  c = octonly(a); d = octonly(b); L = octonly(f);
  c.coord[0]=0;
  c.coord[2]=0;
  c.coord[3]=0;
  c.coord[4]=0;
  c.coord[5]=0;
  c.coord[6]=0;
  c.coord[7]=0;
  d.coord[0]=0;
  d.coord[1]=0;
  d.coord[3]=0;
  d.coord[4]=0;
  d.coord[5]=0;
  d.coord[6]=0;
  d.coord[7]=0;
  L.coord[0]=0;
  L.coord[1]=0;
  L.coord[2]=0;
  L.coord[3]=0;
  sn = ss(c)*ss(d)*ss(L);
  return divscal( sn,
                 m3( mul(d,L), mul(c,L), mul(d,mul(c,y)) ) );
}

/* G2aut49 works in quaternions but only true in real part in octs, fails utterly 16-ons */
eon G2t49(eon a, eon b, eon x, eon y){
  return sub( mul( G2aut49(a,b,x), G2aut49(a,b,y) ), G2aut49(a,b,mul(x,y)) );
}

eon G2t50(eon a, eon b, eon x, eon y){
  return sub( mul( G2aut50(a,b,x), G2aut50(a,b,y) ), G2aut50(a,b,mul(x,y)) );
}

eon G2t51(eon a, eon b, eon x, eon y){
  return sub( mul( G2aut51(a,b,x), G2aut51(a,b,y) ), G2aut51(a,b,mul(x,y)) );
}

/* G2t57 works in all the 2^k-ons */
eon G2t57(eon a, eon b, eon c, eon x, eon y){
  return sub( mul( G2aut57(a,b,c,x), G2aut57(a,b,c,y) ),
                   G2aut57(a,b,c,mul(x,y)) );
}

eon weakRdistrib(eon a, eon b)
{
  return sub(mul(a, add(b, one(a.lev))),
             add(mul(a, b), mul(a, one(a.lev))));
}

eon weakLdistrib(eon a, eon b)
{
  return sub(mul(add(b, one(a.lev)), a),
             add(mul(b, a), mul(one(a.lev), a)));
}

eon octRdistrib(eon a, eon b, eon c)
{
  return sub(mul(a, add(b, niner(c))),
             add(mul(a, b), mul(a, niner(c))));
}

eon octLdistrib(eon a, eon b, eon c)
{                               /*fails for hexons */
  return sub(mul(add(b, niner(c)), a),
             add(mul(b, a), mul(niner(c), a)));
}

eon octRdistrib2(eon a, eon b, eon c)
{
  return sub(mul(niner(a), add(b, c)),
             add(mul(niner(a), b), mul(niner(a), c)));
}

eon octLdistrib2(eon a, eon b, eon c)
{
  return sub(mul(add(b, c), niner(a)),
             add(mul(b, niner(a)), mul(c, niner(a))));
}

eon octRdistrib3(eon a, eon b, eon c)
{
  return sub(mul(a, add(niner(b), niner(c))),
             add(mul(a, niner(b)), mul(a, niner(c))));
}

eon octLdistrib3(eon a, eon b, eon c)
{                               /*fails for hexons */
  return sub(mul(add(niner(b), niner(c)), a),
             add(mul(niner(b), a), mul(niner(c), a)));
}

eon ninRdistrib(eon a, eon b, eon c)
{
  return sub(mul(a, add(b, niner(c))),
             add(mul(a, b), mul(a, niner(c))));
}

eon ninLdistrib(eon a, eon b, eon c)
{                               /*fails for hexons */
  return sub(mul(add(b, niner(c)), a),
             add(mul(b, a), mul(niner(c), a)));
}

eon ninRdistrib2(eon a, eon b, eon c)
{
  return sub(mul(niner(a), add(b, c)),
             add(mul(niner(a), b), mul(niner(a), c)));
}

eon ninLdistrib2(eon a, eon b, eon c)
{
  return sub(mul(add(b, c), niner(a)),
             add(mul(b, niner(a)), mul(c, niner(a))));
}

eon ninRdistrib3(eon a, eon b, eon c)
{
  return sub(mul(a, add(niner(b), niner(c))),
             add(mul(a, niner(b)), mul(a, niner(c))));
}

eon ninLdistrib3(eon a, eon b, eon c)
{                               /*works for niners and zeroers*/
  return sub(mul(add(niner(b), niner(c)), a),
             add(mul(niner(b), a), mul(niner(c), a)));
}

/*Kantor-S:    ae.b = aConj(b).e,    a.be=ba.e,    ae.be=-Conj(b)a */
eon KSe1(eon x, eon y)
{
  eon e, a, b;
  e = iunit(x.lev);
  a = rehalfie(x);
  b = rehalfie(y);
  return sub(mul(mul(a, e), b), mul(mul(a, Conj(b)), e));
}

eon KSe2(eon x, eon y)
{
  eon e, a, b;
  e = iunit(x.lev);
  a = rehalfie(x);
  b = rehalfie(y);
  return sub(mul(a, mul(b, e)), mul(mul(b, a), e));
}

eon KSe3(eon x, eon y)
{
  eon e, a, b;
  e = iunit(x.lev);
  a = rehalfie(x);
  b = rehalfie(y);
  return add(mul(mul(a, e), mul(b, e)), mul(Conj(b), a));
}

/* name "braid" may be bad here */
basetype Rbraid(eon x, eon y, eon z)
{
  return ip(mul(x, y), z) - ip(x, mul(z, Conj(y)));
}

basetype Lbraid(eon x, eon y, eon z)
{
  return ip(mul(Conj(y), x), z) - ip(x, mul(y, z));
}

basetype Lip(eon x, eon y, eon z)
{
  return (ip(mul(x, y), mul(x, z)) - ss(x) * ip(y, z));
}

basetype Rip(eon x, eon y, eon z)
{
  return (ip(mul(y, x), mul(z, x)) - ss(x) * ip(y, z));
}

eon rere(eon x, eon y)
{
  return mul(rehalfie(x), rehalfie(y));
}

eon reim(eon x, eon y)
{
  return mul(rehalfie(x), imhalfie(y));
}

eon imre(eon x, eon y)
{
  return mul(imhalfie(x), rehalfie(y));
}

eon imim(eon x, eon y)
{
  return mul(imhalfie(x), imhalfie(y));
}

/* <ab,cd> + <ad,cb> = 2<a,c><b,d> */
basetype ip4(eon a, eon b, eon c, eon d)
{
  return
      ip(mul(a, b), mul(c, d)) + ip(mul(a, d),
                                    mul(c, b)) - 2 * ip(a, c) * ip(b, d);
}

basetype ip3L(eon a, eon b, eon d)
{
  return ip4(a, b, a, d);
}

basetype ip3R(eon a, eon b, eon c)
{
  return ip4(a, b, c, b);
}


basetype ipcheck1(eon x, eon y)
{
  return ip(x, y) - re(mul(x, Conj(y)));
}

basetype ipcheck2(eon x, eon y)
{
  return ip(x, y) - re(mul(Conj(x), y));
}

basetype ipcheck3(eon x, eon y)
{
  return ip(x, y) + re(mul(x, y));
}

basetype ipcheck4(eon x, eon y)
{
  return 2 * ip(x, y) - re(mul(x, Conj(y))) - re(mul(y, Conj(x)));
}

basetype ipcheck5(eon x, eon y)
{
  return ip(x, y) + re(jor(x, y)) - 2 * re(x) * re(y);
}

/* refl(1) maps  x --> -Conj(x). 
 * refl(t) maps  x -->  x  - 2*ip(x,t)/ss(t) * t.  **/
eon refl(eon x, eon y)
{ /* reflection of x in hyperplane thru 0 normal to y */
    return sub(x, scal(2*ip(x,y)/ss(y), y));
}

/*false*/
eon reflident1(eon x, eon b)  /*Conway-Smith B operator p73*/
{ 
    eon a;   a = normalize(b);
    return  add( Conj(refl(x,a)), mul(mul(a,x),a) );
}

/*false*/
eon reflident2(eon x, eon b)  /*Conway-Smith B operator p73*/
{ 
    eon a;   a = normalize(b);
    return  add( Conj(refl(x,a)), mul(a,mul(x,a)) );
}

eon reflident1R(eon x, eon b)  /*Conway-Smith B operator p73*//*true*/
{
    eon a;   a = normalize(b);
    return  add( refl(Conj(x),a), mul(mul(a,x),a) );
}

eon reflident2R(eon x, eon b)  /*Conway-Smith B operator p73*/ 
{ 
    eon a;   a = normalize(b);
    return  add( refl(Conj(x),a), mul(a,mul(x,a)) );
}


eon reflident1C(eon x, eon b)  /*Conway-Smith B operator p73*/ /*true*/
{ 
    eon a;   a = normalize(b);
    return  add( Conj(refl(x,a)), mul(mul(Conj(a),x),Conj(a)) );
}

eon reflident2C(eon x, eon b)  /*Conway-Smith B operator p73*/ 
{ 
    eon a;   a = normalize(b);
    return  add( Conj(refl(x,a)), mul(Conj(a),mul(x,Conj(a))) );
}

/*false*/
eon reflident1RC(eon x, eon b)  /*Conway-Smith B operator p73*/
{ 
    eon a;   a = normalize(b);
    return  add( refl(Conj(x),a), mul(mul(Conj(a),x),Conj(a)) );
}

/*false*/
eon reflident2RC(eon x, eon b)  /*Conway-Smith B operator p73*/
{ 
    eon a;   a = normalize(b);
    return  add( refl(Conj(x),a), mul(Conj(a),mul(x,Conj(a))) );
}

eon nin1reflident1R(eon x, eon b)  /*Conway-Smith B operator p73*//*true*/
{ 
    return  reflident1R(niner(x), b);
}

eon nin1reflident2R(eon x, eon b)  /*Conway-Smith B operator p73*/ /*true*/
{ 
    return  reflident2R(niner(x), b);
}

eon nin1reflident1C(eon x, eon b)  /*Conway-Smith B operator p73*/ /*true*/
{ 
    return  reflident1C(niner(x), b);
}

eon nin1reflident2C(eon x, eon b)  /*Conway-Smith B operator p73*/ /*true*/
{
    return  reflident2C(niner(x), b); 
}


eon nin2reflident1R(eon x, eon b)  /*Conway-Smith B operator p73*//*true*//*powon 0-16+32*/
{ /*5&nin2reflident1R(a,b)=Z&{0-16}=(00000000000000000NNNNNNNNNNNNNNN)*/
    return  reflident1R(x, niner(b));
}

eon nin2reflident2R(eon x, eon b)  /*Conway-Smith B operator p73*/ /*true*/
{ 
    return  reflident2R(x, niner(b));
}

eon nin2reflident1C(eon x, eon b)  /*Conway-Smith B operator p73*//*true*//*powon 0-16+32*/
{ 
    return  reflident1C(x, niner(b));
}

eon nin2reflident2C(eon x, eon b)  /*Conway-Smith B operator p73*/ /*true*/
{
    return  reflident2C(x, niner(b)); 
}





eon reflex(eon x, eon y)
{
  return sub(mul(x, mul(y, Conj(x))), mul(mul(x, y), Conj(x)));
}

eon sqreflex(eon x, eon y)
{
  return sq(reflex(x, y));
}

eon pow4reflex(eon x, eon y)
{
  return sq(sqreflex(x, y));
}

eon octreflex1(eon x, eon y)
{
  return reflex(niner(x), y);
}

eon octreflex2(eon x, eon y)
{
  return reflex(x, niner(y));
}

eon ninreflex1(eon x, eon y)
{
  return reflex(niner(x), y);
}

eon ninreflex2(eon x, eon y)
{
  return reflex(x, niner(y));
}

eon bruck51ii(eon x, eon y)
{
    eon a,b;
    a = mul(Conj(y), mul(x, mul(y, Conj(x))));
    b = mul(mul(Conj(y), mul(x, y)), Conj(x));
    return(sub(a,b));
}

eon linsq(eon x, eon y)
{
  return sub(sq(add(x, y)), add4(sq(x), sq(y), mul(x, y), mul(y, x)));
}

/* s is scalar, but we use real part of eon to make test macros happy */
eon weaklinsq(eon s, eon y)
{
  eon x;
  x = scalar(y.lev, re(s));
  return sub(sq(add(x, y)), add4(sq(x), sq(y), mul(x, y), mul(y, x)));
}

/* orthogcross1 works if x,y pure imag up thru 64-ons */
basetype orthogcross1(eon x, eon y)
{
  return ip(x, mul(y, x));
}

/* orthogcross2 works if x,y pure imag up thru 64-ons */
basetype orthogcross2(eon x, eon y)
{
  return ip(x, mul(x, y));
}

/* orthogL1 still works for general 64-ons!! */
basetype orthogL1(eon x, eon y)
{
  return ip(x, lie(y, x));
}

/* orthogL2 works for pure-imag 64-ons */
basetype orthogL2(eon x, eon y)
{
  return ip(x, jor(x, y));
}

eon cross3a(eon x, eon y, eon z){
  return sub(
	     add( scal(ip(x,y),z),
		  scal(ip(y,z),x) ),
	     add( scal(ip(x,z),y),
		  mul(x,mul(Conj(y),z)) ) );
}

/* orthog3a1 keeps working thru 64-ons */
basetype orthog3a1(eon x, eon y, eon z){
  return ip( cross3a(x,y,z), x );
}

/* orthog3a2 stops, last valid in 8-ons */
basetype orthog3a2(eon x, eon y, eon z){
  return ip( cross3a(x,y,z), y );
}

/* orthog3a3 stops, last valid in 8-ons */
basetype orthog3a3(eon x, eon y, eon z){
  return ip( cross3a(x,y,z), z );
}

/* gram3a is last valid in 8-ons */
basetype gram3a(eon x, eon y, eon z){
  eon c;
  basetype m[3][3], det;
  c = cross3a(x,y,z);
  m[0][0] = ip(x,x);
  m[0][1] = ip(x,y);
  m[0][2] = ip(x,z);
  m[1][0] = ip(y,x);
  m[1][1] = ip(y,y);
  m[1][2] = ip(y,z);
  m[2][0] = ip(z,x);
  m[2][1] = ip(z,y);
  m[2][2] = ip(z,z);
  det = 
     m[1][1]*m[2][2]*m[0][0]
    -m[1][1]*m[2][0]*m[0][2]
    -m[2][1]*m[1][2]*m[0][0]
    +m[2][1]*m[1][0]*m[0][2]
    +m[0][1]*m[1][2]*m[2][0]
    -m[0][1]*m[1][0]*m[2][2];
  return (ss(c)-det);
}

eon cross3b(eon x, eon y, eon z){
  return sub(
	     add( scal(ip(x,y),z),
		  scal(ip(y,z),x) ),
	     add( scal(ip(x,z),y),
		  mul(mul(x,Conj(y)),z) ) );
}

/* orthog3b1 stops, last valid in 8-ons */
basetype orthog3b1(eon x, eon y, eon z){
  return ip( cross3b(x,y,z), x );
}

/* orthog3b2 stops, last valid in 8-ons */
basetype orthog3b2(eon x, eon y, eon z){
  return ip( cross3b(x,y,z), y );
}

/* orthog3b3 keeps working thru 64-ons */
basetype orthog3b3(eon x, eon y, eon z){
  return ip( cross3b(x,y,z), z );
}

/* gram3b is last valid in 8-ons */
basetype gram3b(eon x, eon y, eon z){
  eon c;
  basetype m[3][3], det;
  c = cross3b(x,y,z);
  m[0][0] = ip(x,x);
  m[0][1] = ip(x,y);
  m[0][2] = ip(x,z);
  m[1][0] = ip(y,x);
  m[1][1] = ip(y,y);
  m[1][2] = ip(y,z);
  m[2][0] = ip(z,x);
  m[2][1] = ip(z,y);
  m[2][2] = ip(z,z);
  det = 
     m[1][1]*m[2][2]*m[0][0]
    -m[1][1]*m[2][0]*m[0][2]
    -m[2][1]*m[1][2]*m[0][0]
    +m[2][1]*m[1][0]*m[0][2]
    +m[0][1]*m[1][2]*m[2][0]
    -m[0][1]*m[1][0]*m[2][2];
  return (ss(c)-det);
}

/* aut3 is last valid in the 8-ons.
 * Thm by Conway & D.Smith (?) is that
 * conjugation by c is an automorphism of the octonions iff
 * c^3=1.  Our expts confirm this condition suffices - but c^3=1 is NOT
 * enough to get an automorphism of the 16-ons. */
eon aut3(eon x, eon y, eon c){
  eon tx, ty, tz, z, c3;
  c.coord[0] = 0;
  c = normalize(c);
  assert(  ss(c) < 1.001 );
  assert(  ss(c) > 0.999 );
  c = scal( RT3BY2, c );
  c.coord[0] = -0.5;
  c3 = m3(c,c,c);
  printeon(c3);
  assert(  NearEq( c3, one(c.lev) )  );
  z = mul(x,y);
  tx = mul(mul(Conj(c), x), c);
  ty = mul(mul(Conj(c), y), c);
  tz = mul(mul(Conj(c), z), c);
  return sub(mul(tx,ty), tz);
}

eon Wilson(eon x, eon y, eon z)
{
    return sub( scal(ss(z), mul(x, Conj(mul(x,y)))), 
	        mul(mul(x,z), Conj(mul(x,mul(y,z)))) );
}

eon LCC2(eon x, eon y, eon z)
{
    return sub( scal(ss(z), mul(mul(x,y), z)),  
                mul(mul(x,z), mul(Conj(z), mul(y,z)))
	);
}


/* LCC identity: (zx)^(-1) * zy.x   is indep of z */
eon LCC(eon x, eon y, eon z, eon w)
{
    return sub( divscal( ss(z), mul( Conj(mul(z,x)), mul(mul(z,y),x) ) ) ,
                divscal( ss(w), mul( Conj(mul(w,x)), mul(mul(w,y),x) ) ) );
}

eon Aloop1(eon x, eon y, eon u, eon w)
{
    return sub(
	mul( mul(mul(mul(u,x),y), Conj(mul(x,y))) ,
	     mul(mul(mul(w,x),y), Conj(mul(x,y))) ), 
        mul(mul(mul(mul(u,w),x),y), Conj(mul(x,y))));
}

eon Aloop2(eon x, eon y, eon u, eon w)
{
    return sub(
	mul(  mul(Conj(mul(y,x)), mul(y,mul(x,u))),
              mul(Conj(mul(y,x)), mul(y,mul(x,w))) ),
        mul(Conj(mul(y,x)), mul(y,mul(x,mul(u,w)))) );
}

eon Aloop3(eon x, eon y, eon u, eon w)
{
    return sub(
	mul(  mul(Conj(x), mul(u,x)), 
              mul(Conj(x), mul(w,x)) ),
        mul(Conj(x), mul(mul(u,w),x)) );
}

/* RIF identity: $xy \cdot (z \cdot xy ) = (x \cdot yz ) x \cdot\ y$ */
eon RIFkkp(eon x, eon y, eon z)
{
    return sub( mul(mul(x,y), mul(z, mul(x,y))),
                mul(mul(mul(x, mul(y,z)), x), y) );
}

eon RIFtest(eon a, eon b)
{
    eon d;
    return sub( mul(b, mul(      a,Conj(b))) , 
           Conj(mul(b, mul(Conj(a),Conj(b)))) );
}

eon RIFtest2(eon a, eon b, eon c)
{
    eon d;
    d = Conj(mul(b,c));
    return sub( mul(mul(b, mul(a,c)), d), 
           Conj(mul(mul(b, mul(Conj(a),c)), d)) );
}

eon RIFtest3(eon a, eon b, eon c)
{
    eon d;
    d = Conj(mul(b,c));
    return sub( mul(d,mul(b, mul(a,c))), 
           Conj(mul(d,mul(b, mul(Conj(a),c)))) );
}

eon RIFtest4(eon a, eon b, eon c)
{
    eon d;
    d = Conj(mul(b,c));
    return sub( mul( mul(mul(a,b),c), d), 
		Conj(mul( mul(mul(Conj(a),b),c), d)) );
}

eon RIFtest5(eon a, eon b, eon c)
{
    eon d;
    d = Conj(mul(c,b));
    return sub( mul(d,mul(c,mul(b,a))),
                Conj(mul(d,mul(c,mul(b,Conj(a))))) );
}

uint BV1, BV2;
eon aut3BV(eon x, eon y){
  return aut3(x,y,  add( basis(x.lev, BV1),  basis(x.lev, BV2) )  );
}

/*  xy.z + x.yz = 2re(yz)x - 2re(xz)y + 2re(xy)z + 2re(xyz)  */
eon Ebbing1p254(eon x, eon y, eon z)
{
  eon s;
  uint8 n;
  n = x.lev;
  assert(n == y.lev);
  assert(n == z.lev);
  s = add(mul(mul(x, y), z), mul(x, mul(y, z)));
  return
      sub(add3
          (s, scalar(n, -re(s)), scal(2 * re(mul(x, z)), y)),
          add(scal(2 * re(mul(y, z)), x), scal(2 * re(mul(x, y)), z)));
}

/******************
*   xy.z + x.yz =  2re(x)yz + 2re(y)xz + 2re(z)xy
*         - 2<y,z>x - 2re(xz)y - 2<x,y>z + 4re(z)<x,y> - 2<x,yz>
*****************************************/
eon Ebbing6p255(eon x, eon y, eon z)
{
  eon s;
  uint8 n;
  n = x.lev;
  assert(n == y.lev);
  assert(n == z.lev);
  s = add(mul(mul(x, y), z), mul(x, mul(y, z)));
  return add(s,
             sub(add5
                 (scal(2 * ip(y, z), x),
                  scal(2 * ip(x, y), z), scalar(n,
                                                2 * ip(x,
                                                       mul(y,
                                                           z))),
                  scalar(n, -4 * re(z) * ip(x, y)),
                  scal(2 * re(mul(x, z)), y)),
                 add3(scal(2 * re(x), mul(y, z)),
                      scal(2 * re(y), mul(x, z)),
                      scal(2 * re(z), mul(x, y)))));
}

eon ZSSSlem5p130(eon x, eon y, eon z)
{
  return sub(scal(2, assoc(jor(mul(z, x), mul(x, z)), y, z)),
             assoc(add(mul(mul(x, z), x), mul(x, mul(z, x))), y, sq(z)));
}

eon ZSSSp145eq15(eon x, eon y, eon z)
{
  return ujor(assoc(x, y, z), lie(x, y));
}

eon ZSSSp145eq13(eon x, eon y, eon z)
{
  return sub(assoc(x, y, assoc(x, y, z)), mul(lie(x, y), assoc(x, y, z)));
}

/* Kleinfeld1 = Kleinfeld2, ZSSSp139.
* In alt algebra, if any 2 arguments equal, Kleinfeld=0. */
eon Kleinfeld1(eon w, eon x, eon y, eon z)
{
  return sub(assoc(mul(w, x), y, z),
             add(mul(x, assoc(w, y, z)), mul(assoc(x, y, z), w)));
}

eon Kleinfeld2(eon w, eon x, eon y, eon z)
{
  return add(assoc(lie(w, x), y, z), assoc(lie(y, z), w, x));
}

eon Kleindiff(eon w, eon x, eon y, eon z)
{
  return sub(Kleinfeld2(w, x, y, z), Kleinfeld1(w, x, y, z));
}

/* I.Kleinfeld: PAMS 4 (1953) 939-944 */
eon KleinIdent(eon w, eon x, eon y, eon z)
{
  return
      sub(add
          (assoc(mul(w, x), y, z), assoc(w, x, lie(y, z))),
          add(mul(w, assoc(x, y, z)), mul(assoc(w, y, z), z)));
}

bool imaganticomm(uint8 n)  /* true only when n<=3. false when n>=4. */
{
    uint x, y, p, hp;
    double theta,s,c;
    basetype z;
  eon a,b;
  p = pow2(n);
  hp = p/2;
  a = zero(n);
  b = zero(n);
  for (x = 1; x < hp; x++){ a.coord[x] = 1.0; }
  for (x = hp; x < p; x++){ b.coord[x] = 1.0; }
  for (x=2; x<p; x++){
      for (y=1; y<x; y++){
	  theta = (2.0*PI)*myrand();
	  s = sin(theta); c = cos(theta);
          z          = c*a.coord[x] + s*a.coord[y];
          a.coord[y] = c*a.coord[y] - s*a.coord[x];
          a.coord[x] = z;
          z          = c*b.coord[x] + s*b.coord[y];
          b.coord[y] = c*b.coord[y] - s*b.coord[x];
          b.coord[x] = z;
      }
  }
  z = ip(a,b); if(z<0.0) z = -z;
  assert( z < 0.0001 );
  if( !isnearzero( add(mul(a,b), mul(b,a)) ) )   return FALSE;
  return TRUE;
}

bool basiscomm(uint8 n)
{
  uint x, y, p;
  p = pow2(n);
  for (x = 1; x < p; x++){
    for (y = 1; y < p; y++)
      if(y != x){
        if(!isnearzero(lie(basis(n, x), basis(n, y))))
          return FALSE;
      }
  }
  return TRUE;
}

bool basisanticomm(uint8 n)
{
  uint x, y, p;
  p = pow2(n);
  for (x = 1; x < p; x++){
    for (y = 1; y < p; y++)
      if(y != x){
        if(!isnearzero(ujor(basis(n, x), basis(n, y))))
          return FALSE;
      }
  }
  return TRUE;
}

bool basis3cyc(uint8 n)
{
  uint x, y, z, p;
  p = pow2(n);
  for (x = 1; x < p; x++){
    for (y = 1; y < p; y++)
      if(y != x){
        for (z = 1; z < p; z++)
          if(z != x && z != y){
            if(!isnearzero
               (sub
                (mul
                 (basis(n, x),
                  mul(basis(n, y), basis(n, z))),
                 mul(basis(n, z), mul(basis(n, x), basis(n, y))))))
              return FALSE;
          }
      }
  }
  return TRUE;
}

bool basisFlex(uint8 n)
{
  uint x, y, p;
  p = pow2(n);
  for (x = 1; x < p; x++){
    for (y = 1; y < p; y++)
      if(y != x){
        if(!isnearzero
           (sub
            (mul
             (basis(n, x),
              mul(basis(n, y), basis(n, x))),
             mul(mul(basis(n, x), basis(n, y)), basis(n, x)) )))
          return FALSE;
      }
  }
  return TRUE;
}

bool basisRalt(uint8 n)
{
  uint x, y, p;
  p = pow2(n);
  for (x = 1; x < p; x++){
    for (y = 1; y < p; y++)
      if(y != x){
        if(!isnearzero
           (sub
            (mul
             (basis(n, x),
              mul(basis(n, y), basis(n, y))),
             mul(mul(basis(n, x), basis(n, y)), basis(n, y)) )))
          return FALSE;
      }
  }
  return TRUE;
}

bool basisMoufmid(uint8 n)
{
    uint x, y, p, z;
  p = pow2(n);
  for (x = 1; x < p; x++){
      for (y = 1; y < p; y++){
	  if(y!=x)
	      for (z = 1; z < p; z++)
		  if(z!=x && z!=y){
		      if(!isnearzero(
    sub(mul(mul(basis(n,x), mul(basis(n,y), basis(n,z))), basis(n,x)), 
        mul(basis(n,x), mul(mul(basis(n,y), basis(n,z)), basis(n,x))))
			     ));
		      return FALSE;
		  }
      }
  }
  return TRUE;
}

/* ???something wrong here */
bool basisantiassoc(uint8 n)
{
  uint x, y, z, p;
  p = pow2(n);
  for (x = 0; x < p; x++){
    for (y = 0; y < p; y++)
      if(y != x){
        for (z = 0; z < p; z++)
          if(z != x && z != y){
            if(!isnearzero
               (add
                (mul
                 (basis(n, x),
                  mul(basis(n, y), basis(n, z))),
                 mul(basis(n, z), mul(basis(n, y), basis(n, x))))))
              return FALSE;
          }
      }
  }
  return TRUE;
}

bool basisassoc(uint8 n)
{
  uint x, y, z, p;
  p = pow2(n);
  for (x = 0; x < p; x++){
    for (y = 0; y < p; y++)
      if(y != x){
        for (z = 0; z < p; z++)
          if(z != x && z != y){
            if(!isnearzero
               (sub
                (mul
                 (basis(n, x),
                  mul(basis(n, y), basis(n, z))),
                 mul(mul(basis(n, x), basis(n, y)), basis(n, z)))))
              return FALSE;
          }
      }
  }
  return TRUE;
}

/* Ebb (3) p277.  0 if x,y pure imag */
eon jorip(eon x, eon y)
{
  return add(ujor(x, y), scalar(x.lev, 2 * ip(x, y)));
}


/* Ebbing 4 p251 */
eon jorcheck1(eon x, eon y)
{
  return sub(ujor(x, y),
             add3(scal(2 * re(x), y), scal(2 * re(y), x),
                  scalar(x.lev, -2 * ip(x, y))));
}

eon alt1(eon x, eon y, eon z)
{
  return sub(assoc(x, y, z), assoc(y, z, x));
}

eon Lalt2(eon x, eon y, eon z)
{
  return add(assoc(x, y, z), assoc(y, x, z));
}

eon Ralt2(eon x, eon y, eon z)
{
  return add(assoc(x, y, z), assoc(x, z, y));
}

/*1st arg of scalcheck1...scalcheck3 is scalar. We use real part of an eon to
make test macros happier */
eon scalcheck1(eon s, eon a, eon b)
{                                 /* s.ab = sa.b */
  basetype res;
  res = re(s);
  return sub(scal(res, mul(a, b)), mul(scal(res, a), b));
}

eon scalcheck2(eon s, eon a, eon b)
{                                /* s.ab = a.sb */
  basetype res;
  res = re(s);
  return sub(scal(res, mul(a, b)), mul(a, scal(res, b)));
}

eon scalcheck3(eon s, eon a, eon b)
{                                /* sa.b = a.sb */
  basetype res;
  res = re(s);
  return sub(mul(scal(res, a), b), mul(a, scal(res, b)));
}

/* Conj(x) = sum_b  b.xb / (2-p)    for the (2^n=p)-ons, b=basis elmtns
* sum is over p terms */
eon Lconjlin(eon x)
{
  eon b, z;
  uint n, j, p;
  basetype s;
  n = x.lev;
  p = pow2(n);
  s = p - 2;
  z = x;
  for (j = 1; j < p; j++){
    b = basis(n, j);
    z = add(z, mul(b, mul(x, b)));
  }
  /* at this point z = -(p-2)*Conj(x) */
  return add(scal(s, Conj(x)), z);
}

/* Conj(x) = sum_b  bx.b / (2-p)    for the (2^n=p)-ons, b=basis elmtns
* sum is over p terms */
eon Rconjlin(eon x)
{
  eon b, z;
  uint n, j, p;
  basetype s;
  n = x.lev;
  p = pow2(n);
  s = p - 2;
  z = x;
  for (j = 1; j < p; j++){
    b = basis(n, j);
    z = add(z, mul(mul(b, x), b));
  }
  /* at this point z = -(p-2)*Conj(x) */
  return add(scal(s, Conj(x)), z);
}

/* 2Conj(x) = sum_b  (b.xb+bx.b) / (2-p)    for the (2^n=p)-ons,
* b=basis elmtns; sum is over p terms */
eon conjlin(eon x)
{
  eon b, z;
  uint n, j, p;
  basetype s;
  n = x.lev;
  p = pow2(n);
  s = p - 2;
  s *= 2;
  z = scal(2, x);
  for (j = 1; j < p; j++){
    b = basis(n, j);
    z = add3(z, mul(b, mul(x, b)), mul(mul(b, x), b));
  }
  /* at this point z = -2(p-2)*Conj(x) */
  return add(scal(s, Conj(x)), z);
}

/* none of the above will work since nonlinear hexon mul.
*   sum_m Conj(Conj(em) Conj(em x)) = -sum_m Conj(em Conj(em x))  
* however seems to do the job right, excpet you get x back not conj(x). */
eon conjnonlin(eon x)
{
  eon b, z;
  uint n, j, p;
  basetype s;
  n = x.lev;
  p = pow2(n);
  s = p - 2;
  z = x;
  for (j = 1; j < p; j++){
    b = basis(n, j);
    z = add(z, Conj(mul(b, Conj(mul(b, x)))));
  }
  /* at this point z = (p-2)*x */
  return sub(scal(s, x), z);
}

eon simil(eon a, eon b){
  return mul(Conj(b), mul(a, b));
}

eon sqsimil(eon a, eon b){
  return sq(simil(a, b));
}

eon sim(eon a, eon b){
  return sub(mul(Conj(b), mul(a, b)), scal(ss(b), a));
}

eon sim2(eon x, eon y, eon b){
  return sub(mul(mul(Conj(b), x), mul(y, b)), scal(ss(b), mul(x, y)));
}

eon sim3(eon x, eon y, eon b){
  return sub(mul(mul(mul(Conj(b), x), y), b), scal(ss(b), mul(x, y)));
}

eon sqsim(eon x, eon y){
  return sq(sim(x, y));
}

eon sqsim2(eon x, eon y, eon z){
  return sq(sim2(x, y, z));
}

eon sqsim3(eon x, eon y, eon z){
  return sq(sim3(x, y, z));
}

eon Lcan(eon a, eon b){
  return sub(mul(mul(a, Conj(a)), b), mul(a, mul(Conj(a), b)));
}

eon sqLcan(eon x, eon y){
  return sq(Lcan(x, y));
}

eon pow4Lcan(eon x, eon y){
  return sq(sqLcan(x, y));
}

eon Rcan(eon a, eon b){
  return sub(mul(mul(b, Conj(a)), a), mul(b, mul(Conj(a), a)));
}

eon sqRcan(eon x, eon y){
  return sq(Rcan(x, y));
}

eon pow4Rcan(eon x, eon y){
  return sq(sqRcan(x, y));
}

eon octLcan2(eon x, eon y){
  return Lcan(x, niner(y));
}

eon octLcan1(eon x, eon y){
  return Lcan(niner(x), y);
}

eon octRcan1(eon x, eon y){
  return Rcan(x, niner(y));
}

eon octRcan2(eon x, eon y){
  return Rcan(niner(x), y);
}

eon ninLcan2(eon x, eon y){
  return Lcan(x, niner(y));
}

eon ninLcan1(eon x, eon y){
  return Lcan(niner(x), y);
}

eon ninRcan1(eon x, eon y){
  return Rcan(x, niner(y));
}

eon ninRcan2(eon x, eon y){
  return Rcan(niner(x), y);
}


eon Lsand(eon a, eon b, eon c){
  return sub(mul(mul(a, b), mul(c, Conj(b))),
             mul(a, mul(mul(b, c), Conj(b))));
}

eon Lsand2(eon a, eon b){
  return Lsand(a, b, a);
}

eon Lsand2C(eon a, eon b){
  return Lsand(a, b, Conj(a));
}

basetype foo(eon a, eon b){
  basetype t;
  t = ip(a, b);
  return ss(imul(a, b)) + t * t - ss(a) * ss(b);
}

basetype Nmultip(eon a, eon b){
  return (ss(a) * ss(b) - ss(mul(a, b)));
}

basetype Nmulsim(eon a, eon b){
  ubasetype x;
  x = ss(b);
  x *= x;
  return (ss(a) * x - ss(mul(mul(Conj(b), a), b)));
}

/* quadratic ident  2*re(x)*x - x*x = norm(x) = x*Conj(x) */
eon quadident(eon x){
  uint8 n;
  n = x.lev;
  return sub(sub(scal(2 * re(x), x), sq(x)), scalar(n, ss(x)));
}

basetype requad(eon x){
  basetype r;
  r = re(x);
  return ss(x) - 2 * r * r + re(sq(x));
}


/***************************************************/

real64 jacsamp(uint n){
  eon a,b;
  real64 jd;
  a=UnitRand(n,FALSE);
  b=UnitRand(n,FALSE);
  assert( ss(a) > 0.999 );
  assert( ss(a) < 1.001 );
  assert( ss(b) > 0.999 );
  assert( ss(b) < 1.001 );
  jd = jacdet(a,b);
  return jd;
}

void testjac(){
    uint n,k,cthi,ctlo,ctneg;
  real64 jd, sj, jx, sa, jn, jn2, jda, jpa, jna;
  printf("JACOBIAN DET TESTER:\n");
  for(n=0; n<=LGMAXSIZE; n++){
    cthi=ctlo=ctneg=0;
    jpa = jna = sa = 0;
    printf("UnitRand Jacobian Samples(%d)=", pow2(n));  
    sj = 0;
    jx=0;
    jn=99999;
    jn2=jn;
#define NUMSAMP 10000
    for(k=0; k<NUMSAMP; k++){
      do{
       jd = jacsamp(n);
       if(jd<0) jda = -jd;
       else     jda = jd;

      }while( jda > 2500 );  /* avoid super-large jd's probly /0 */
      if(NUMSAMP<1001){
	printf("%.3f ", jd);
	fflush(stdout);
      }
      sj += jd;
      if(jd<jn2) jn2=jd;
      if(jd>1.6) cthi++;
      if(jd<0.2) ctlo++;
      if(jd < -0.000001) ctneg++;
      if(jd>jx) jx=jd;      
      if(jd>0) jpa += jd;
      else     jna += jd;
      jd = jda;
      if(jd<jn) jn=jd;
      sa += jda;
    }
    printf("\n");
    printf(
"AvgJacobian(%d)=%f Avg|J|=%f A+=%f A-=%f minJ=%f min|J|=%f maxJ=%f (#>1.6)=%u/%u (#<0.2)=%u/%u (#<0)=%u/%u\n",
        pow2(n), sj/NUMSAMP, sa/NUMSAMP, 
        jpa/(NUMSAMP-ctneg), jna/(ctneg), jn2, jn, jx, 
        cthi, NUMSAMP, ctlo, NUMSAMP, ctneg, NUMSAMP);
    fflush(stdout);    
  }
}
/************************
**************************/


int testidents(uint n, bool pureimag){
  eon a, b, c, d, e, f, t1, t2;
  uint points = 0;
  uint k, p2n, be;
  real64 jd;
  p2n = pow2(n);
  printf
      ("TESTING IDENTS(2^%d=%d) [EXERCISE CARE RE NUMERICS!!) ",
       n, p2n);
  if(pureimag)
    printf("IMAG ONLY\n");
  else
    printf("UNRESTRICTED\n");
  /*printf("Conj(a)-a="); printsgn(sub(Conj(a),a), 0.0001, 0.000003); */
  a = EErand((n), (pureimag), -9);
  b = EErand((n), (pureimag), -9);
  c = EErand((n), (pureimag), -9);
  d = EErand((n), (pureimag), -9);
  e = EErand((n), (pureimag), -9);
  f = EErand((n), (pureimag), -9);
  printf("(%d)reab: %g\n", n, reab(a,b));
  TESTREAL2(n, Nmultip);
  TESTREAL2(n, Nmulsim);
  TESTIDENT1(n, pureimag, selfback);
  TESTIDENT1(n, pureimag, CDsqsame);
  TESTIDENT2(n, pureimag, ky142a);
  TESTIDENT2(n, pureimag, CDvsme);
  TESTIDENT2(n, pureimag, biCDvsme);
  TESTIDENT2(n, pureimag, coCDvsme);
  TESTIDENT2(n, pureimag, bi2CDvsme);
  TESTIDENT2(n, pureimag, co2CDvsme);
  TESTIDENT2(n, pureimag, reflident1);
  TESTIDENT2(n, pureimag, reflident2);
  TESTIDENT2(n, pureimag, reflident1R);
  TESTIDENT2(n, pureimag, reflident2R);
  TESTIDENT2(n, pureimag, reflident1C);
  TESTIDENT2(n, pureimag, reflident2C);
  TESTIDENT2(n, pureimag, reflident1RC);
  TESTIDENT2(n, pureimag, reflident2RC);
  TESTIDENT2(n, pureimag, divcheck1);
  TESTIDENT2(n, pureimag, divcheck2);
  TESTIDENT1(n, pureimag, quadident);
  TESTIDENT1(n, pureimag, qbarq);
  TESTIDENT2(n, pureimag, lie);
  TESTIDENT2(n, pureimag, ujor);
  TESTIDENT2(n, pureimag, squjor);
  TESTIDENT2(n, pureimag, sqlie);
  TESTIDENT2(n, pureimag, pow4lie);
  TESTIDENT1(n, pureimag, lieC);
  TESTIDENT1(n, pureimag, sqlieC);
  TESTIDENT1(n, pureimag, pow4lieC);
  TESTIDENT3(n, pureimag, assoc);
  TESTIDENT3(n, pureimag, sqassoc);
  TESTIDENT2(n, pureimag, flex);
  TESTIDENT2(n, pureimag, flex2);
  TESTIDENT2(n, pureimag, sqflex);
  TESTIDENT2(n, pureimag, pow4flex);
  TESTIDENT2(n, pureimag, octflex1);
  TESTIDENT2(n, pureimag, octflex2);
  TESTIDENT2(n, pureimag, ninninCDvsme);
  TESTIDENT2(n, pureimag, ninCDvsme);
  TESTIDENT2(n, pureimag, nin1CDvsme);
  TESTIDENT2(n, pureimag, ninbiCDvsme);
  TESTIDENT2(n, pureimag, nincoCDvsme);
  TESTIDENT2(n, pureimag, ninbi2CDvsme);
  TESTIDENT2(n, pureimag, ninco2CDvsme);
  TESTIDENT2(n, pureimag, nin1reflident1R);
  TESTIDENT2(n, pureimag, nin1reflident2R);
  TESTIDENT2(n, pureimag, nin1reflident1C);
  TESTIDENT2(n, pureimag, nin1reflident2C);
  TESTIDENT2(n, pureimag, nin2reflident1R);
  TESTIDENT2(n, pureimag, nin2reflident2R);
  TESTIDENT2(n, pureimag, nin2reflident1C);
  TESTIDENT2(n, pureimag, nin2reflident2C);
  TESTIDENT2(n, pureimag, ninflex1);
  TESTIDENT2(n, pureimag, ninflex2);
  TESTIDENT2(n, pureimag, reflex);
  TESTIDENT2(n, pureimag, sqreflex);
  TESTIDENT2(n, pureimag, pow4reflex);
  TESTIDENT2(n, pureimag, octreflex1);
  TESTIDENT2(n, pureimag, octreflex2);
  TESTIDENT2(n, pureimag, ninreflex1);
  TESTIDENT2(n, pureimag, ninreflex2);
  TESTIDENT2(n, pureimag, bruck51ii);
  TESTIDENT4(n, pureimag, LCC);
  TESTIDENT3(n, pureimag, RIFkkp);
  TESTIDENT4(n, pureimag, Aloop1);
  TESTIDENT4(n, pureimag, Aloop2);
  TESTIDENT4(n, pureimag, Aloop3);
  TESTIDENT2(n, pureimag, RIFtest);
  TESTIDENT3(n, pureimag, RIFtest2);
  TESTIDENT3(n, pureimag, RIFtest3);
  TESTIDENT3(n, pureimag, RIFtest4);
  TESTIDENT3(n, pureimag, RIFtest5);
  TESTIDENT2(n, pureimag, Lalt);
  TESTIDENT2(n, pureimag, sqLalt);
  TESTIDENT2(n, pureimag, pow4Lalt);
  TESTIDENT2(n, pureimag, Ralt);
  TESTIDENT2(n, pureimag, sqRalt);
  TESTIDENT2(n, pureimag, pow4Ralt);
  TESTIDENT2(n, pureimag, Lcan);
  TESTIDENT2(n, pureimag, sqLcan);
  TESTIDENT2(n, pureimag, pow4Lcan);
  TESTIDENT2(n, pureimag, Rcan);
  TESTIDENT2(n, pureimag, sqRcan);
  TESTIDENT2(n, pureimag, pow4Rcan);
  TESTIDENT2(n, pureimag, octLcan1);
  TESTIDENT2(n, pureimag, octRcan1);
  TESTIDENT2(n, pureimag, octLcan2);
  TESTIDENT2(n, pureimag, octRcan2);
  TESTIDENT2(n, pureimag, ninLcan1);
  TESTIDENT2(n, pureimag, ninRcan1);
  TESTIDENT2(n, pureimag, ninLcan2);
  TESTIDENT2(n, pureimag, ninRcan2);
  TESTIDENT1(n, pureimag, powass3);
  TESTIDENT1(n, pureimag, sqpowass3);
  TESTIDENT1(n, pureimag, powass4);
  TESTIDENT1(n, pureimag, sqpowass4);
  TESTIDENT2(n, pureimag, antiaut);
  TESTIDENT2(n, pureimag, autinv);
  TESTIDENT2(n, pureimag, sqantiaut);
  TESTIDENT1(n, pureimag, Cantiaut);
  TESTIDENT1(n, pureimag, sqCantiaut);
  TESTIDENT1(n, pureimag, antiaut1);
  TESTIDENT1(n, pureimag, sqantiaut1);
  TESTIDENT2(n, pureimag, octantiaut);
  TESTIDENT2(n, pureimag, ninantiaut);
  TESTIDENT2(n, pureimag, trip);
  TESTIDENT2(n, pureimag, sqtrip);
  TESTIDENT3(n, pureimag, Ldistrib);
  TESTIDENT4(n, pureimag, ky141);
  TESTIDENT4(n, pureimag, G2t49);
  TESTIDENT4(n, pureimag, G2t50);
  TESTIDENT4(n, pureimag, G2t51);
  TESTIDENT5(n, pureimag, G2t57);
  TESTIDENT3(n, pureimag, Rdistrib);
  TESTIDENT1(n, pureimag, LdistribPow);
  TESTIDENT1(n, pureimag, RdistribPow);
  TESTIDENT1(n, pureimag, LdistribPow3);
  TESTIDENT1(n, pureimag, RdistribPow3);
  TESTIDENT2(n, pureimag, Ldistrib2);
  TESTIDENT2(n, pureimag, Rdistrib2);
  TESTIDENT2(n, pureimag, Ldistrib3);
  TESTIDENT2(n, pureimag, Rdistrib3);
  TESTIDENT2(n, pureimag, Ldistrib4);
  TESTIDENT2(n, pureimag, Rdistrib4);
  TESTIDENT2(n, pureimag, Ldistrib2c);
  TESTIDENT2(n, pureimag, Rdistrib2c);
  TESTIDENT2(n, pureimag, Ldistrib3c);
  TESTIDENT2(n, pureimag, Rdistrib3c);
  TESTIDENT2(n, pureimag, Ldistrib4c);
  TESTIDENT2(n, pureimag, Rdistrib4c);
  TESTIDENT2(n, pureimag, Ldistrib5);
  TESTIDENT2(n, pureimag, Rdistrib5);
  TESTIDENT2(n, pureimag, Ldistrib6);
  TESTIDENT2(n, pureimag, Rdistrib6);
  TESTIDENT2(n, pureimag, Ldistrib7);
  TESTIDENT2(n, pureimag, Rdistrib7);
  if(n>=2){
    TESTIDENT2(n, pureimag, eyeshiftA);
    TESTIDENT2(n, pureimag, eyeshiftB);
    TESTIDENT2(n, pureimag, eyeshiftC);
    TESTIDENT2(n, pureimag, octeyeshiftA);
    TESTIDENT2(n, pureimag, octeyeshiftB);
    TESTIDENT2(n, pureimag, octeyeshiftC);
    TESTIDENT2(n, pureimag, nineyeshiftA);
    TESTIDENT2(n, pureimag, nineyeshiftB);
    TESTIDENT2(n, pureimag, nineyeshiftC);

    TESTIDENT2(n, pureimag, ReyeshiftA);
    TESTIDENT2(n, pureimag, ReyeshiftB);
    TESTIDENT2(n, pureimag, ReyeshiftC);
    TESTIDENT2(n, pureimag, octReyeshiftA);
    TESTIDENT2(n, pureimag, octReyeshiftB);
    TESTIDENT2(n, pureimag, octReyeshiftC);
    TESTIDENT2(n, pureimag, ninReyeshiftA);
    TESTIDENT2(n, pureimag, ninReyeshiftB);
    TESTIDENT2(n, pureimag, ninReyeshiftC);
  }
  TESTIDENT3(n, pureimag, linbimul);
  TESTIDENT3(n, pureimag, lincomul);
  TESTIDENT3(n, pureimag, ninlinbimul);
  TESTIDENT3(n, pureimag, ninlincomul);
  TESTIDENT3(n, pureimag, linbi2mul);
  TESTIDENT3(n, pureimag, linco2mul);
  TESTIDENT3(n, pureimag, ninlinbi2mul);
  TESTIDENT3(n, pureimag, ninlinco2mul);

  TESTIDENT3(n, pureimag, octLdistrib);
  TESTIDENT3(n, pureimag, octRdistrib);
  TESTIDENT3(n, pureimag, octLdistrib2);
  TESTIDENT3(n, pureimag, octRdistrib2);
  TESTIDENT3(n, pureimag, octLdistrib3);
  TESTIDENT3(n, pureimag, octRdistrib3);
  TESTIDENT3(n, pureimag, ninLdistrib);
  TESTIDENT3(n, pureimag, ninRdistrib);
  TESTIDENT3(n, pureimag, ninLdistrib2);
  TESTIDENT3(n, pureimag, ninRdistrib2);
  TESTIDENT3(n, pureimag, ninLdistrib3);
  TESTIDENT3(n, pureimag, ninRdistrib3);
  TESTIDENT2(n, pureimag, weakLdistrib);
  TESTIDENT2(n, pureimag, weakRdistrib);
  if(n > 0){
    TESTIDENT3(n, pureimag, Ldistribimh);
    TESTIDENT3(n, pureimag, Rdistribimh);
    TESTIDENT3(n, pureimag, Ldistribreh);
    TESTIDENT3(n, pureimag, Rdistribreh);
  }
  TESTIDENT2(n, pureimag, linsq);
  TESTIDENT2(n, pureimag, weaklinsq);
  TESTIDENT2(n, pureimag, jordiJ);
  TESTIDENT1(n, pureimag, jordi1);
  TESTIDENT2(n, pureimag, jordiMR);
  TESTIDENT2(n, pureimag, jordiML);
  TESTIDENT4(n, pureimag, jordiMRlin);
  TESTIDENT4(n, pureimag, jordiMRlinA);
  TESTIDENT4(n, pureimag, jordiMLlin);
  TESTIDENT2(n, pureimag, jordiMR1);
  TESTIDENT2(n, pureimag, jordiML1);
  TESTIDENT4(n, pureimag, zassident);
  TESTIDENT6(n, pureimag, zassext);
  TESTREAL4(n, sszassident);
  TESTIDENT3(n, pureimag, JacLie);
  TESTIDENT3(n, pureimag, sqJacLie);
  TESTIDENT2(n, pureimag, JacLie21);
  TESTIDENT3(n, pureimag, MalcevLie);
  TESTIDENT3(n, pureimag, MalcevMul);
  TESTIDENT3(n, pureimag, MalcevI);
  TESTIDENT3(n, pureimag, varMal1Lie);
  TESTIDENT3(n, pureimag, varMal2Lie);
  TESTIDENT3(n, pureimag, varMal3Lie);
  TESTIDENT3(n, pureimag, varMal1I);
  TESTIDENT3(n, pureimag, varMal2I);
  TESTIDENT3(n, pureimag, varMal3I);
  TESTIDENT3(n, pureimag, FenyvesExtra);
  TESTIDENT3(n, pureimag, Cloop);
  TESTIDENT3(n, pureimag, LCloop);
  TESTIDENT3(n, pureimag, RCloop);
  TESTIDENT3(n, pureimag, W1a);
  TESTIDENT3(n, pureimag, W1b);
  TESTIDENT3(n, pureimag, W1c);
  TESTIDENT3(n, pureimag, W1d);
  TESTIDENT3(n, pureimag, W2a);
  TESTIDENT3(n, pureimag, W2b);
  TESTIDENT3(n, pureimag, W2c);
  TESTIDENT3(n, pureimag, W2d);
  TESTIDENT3(n, pureimag, LG1);
  TESTIDENT3(n, pureimag, LG2);
  TESTIDENT3(n, pureimag, LG3);
  TESTIDENT3(n, pureimag, LC2);
  TESTIDENT3(n, pureimag, LC3);
  TESTIDENT3(n, pureimag, LC4);
  TESTIDENT3(n, pureimag, RG1);
  TESTIDENT3(n, pureimag, RG2);
  TESTIDENT3(n, pureimag, RG3);
  TESTIDENT3(n, pureimag, RC2);
  TESTIDENT3(n, pureimag, RC3);
  TESTIDENT3(n, pureimag, RC4);
  TESTIDENT3(n, pureimag, Lnucsq);
  TESTIDENT3(n, pureimag, Rnucsq);
  TESTIDENT3(n, pureimag, Mnucsq);
  TESTIDENT3(n, pureimag, pflM2);
  TESTIDENT3(n, pureimag, pflM3);
  TESTIDENT3(n, pureimag, pflM4);
  TESTIDENT3(n, pureimag, pflM5);
  TESTIDENT3(n, pureimag, pflM6);
  TESTIDENT3(n, pureimag, pflM7);
  TESTIDENT3(n, pureimag, MoufmidL);
  TESTIDENT3(n, pureimag, MoufmidM);
  TESTIDENT3(n, pureimag, MoufmidR);
  TESTIDENT3(n, pureimag, MoufmidLC);
  TESTIDENT3(n, pureimag, MoufmidMC);
  TESTIDENT3(n, pureimag, MoufmidRC);
  TESTIDENT3(n, pureimag, varMouf2R);
  TESTIDENT3(n, pureimag, varMouf3R);
  TESTIDENT3(n, pureimag, varMouf2L);
  TESTIDENT3(n, pureimag, varMouf3L);
  TESTIDENT3(n, pureimag, varMouf2Rv);
  TESTIDENT3(n, pureimag, varMouf2Lv);
  TESTIDENT3(n, pureimag, varMouf2Ra);
  TESTIDENT3(n, pureimag, varMouf2Rc);
  TESTIDENT3(n, pureimag, varMouf2Rd);
  TESTIDENT3(n, pureimag, varMouf2Re);
  TESTIDENT3(n, pureimag, varMouf2Le);
  TESTIDENT3(n, pureimag, mouf5L);
  TESTIDENT3(n, pureimag, mouf6L);
  TESTIDENT3(n, pureimag, mouf7L);
  TESTIDENT3(n, pureimag, mouf5R);
  TESTIDENT3(n, pureimag, mouf6R);
  TESTIDENT3(n, pureimag, mouf7R);
  TESTIDENT3(n, pureimag, Lbol);
  TESTIDENT3(n, pureimag, Rbol);
  TESTIDENT3(n, pureimag, Lmouf);
  TESTIDENT3(n, pureimag, Rmouf);
  TESTIDENT3(n, pureimag, LbolC);
  TESTIDENT3(n, pureimag, RbolC);
  TESTIDENT3(n, pureimag, sqMoufmidL);
  TESTIDENT3(n, pureimag, sqMoufmidM);
  TESTIDENT3(n, pureimag, sqMoufmidR);
  TESTIDENT3(n, pureimag, sqLbol);
  TESTIDENT3(n, pureimag, sqRbol);
  TESTIDENT3(n, pureimag, comm3);
  TESTIDENT3(n, pureimag, comm3A);
  TESTIDENT3(n, pureimag, comm3B);
  TESTIDENT4(n, pureimag, lieassoc);
  TESTIDENT3(n, pureimag, lieflex3);
  TESTIDENT2(n, pureimag, lieflex2);
  TESTIDENT2(n, pureimag, HPQ9L);
  TESTIDENT2(n, pureimag, HPQ9R);
  TESTIDENT2(n, pureimag, HPQ11);
  TESTIDENT2(n, pureimag, HPQ12);
  TESTIDENT4(n, pureimag, jorlies);
  TESTIDENT3(n, pureimag, expt31);
  TESTIDENT3(n, pureimag, expt32);
  TESTIDENT3(n, pureimag, expt3);
  TESTIDENT3(n, pureimag, sym3ident);
  TESTIDENT4(n, pureimag, sym4ident);
  TESTIDENT2(n, pureimag, Rasmyslov);
  TESTIDENT4(n, pureimag, bigRacine);
  TESTIDENT4(n, pureimag, HP23);
  TESTIDENT5(n, pureimag, HP26);
  if(n > 0){
    TESTIDENT2(n, pureimag, rere);
    TESTIDENT2(n, pureimag, reim);
    TESTIDENT2(n, pureimag, imre);
    TESTIDENT2(n, pureimag, imim);
  }
  TESTREAL3(n, a4);
  TESTIDENT3(n, pureimag, cyc3);
  TESTIDENT3(n, pureimag, sqcyc3);
  TESTIDENT2(n, pureimag, crosstest);
  TESTREAL2(n, gramC);
  TESTREAL2(n, gramI);
  TESTIDENT2(n, pureimag, jorip);
  TESTIDENT2(n, pureimag, crossanti);
  TESTIDENT2(n, pureimag, imulanti);
  TESTREAL3(n, crossxchg);
  TESTREAL3(n, imulxchg);
  TESTIDENT2(n, pureimag, aXaXb);
  TESTIDENT2(n, pureimag, aXaXbI);
  TESTIDENT2(n, pureimag, aXbXaI);
  TESTIDENT2(n, pureimag, aXbXaI2);
  TESTIDENT2(n, pureimag, bXaXa);
  TESTIDENT2(n, pureimag, bXaXaI);
  TESTIDENT3(n, pureimag, Grassman1);
  TESTIDENT3(n, pureimag, Grassman2);
  TESTIDENT3(n, pureimag, Grassman1C);
  TESTIDENT3(n, pureimag, Grassman2C);
  TESTIDENT3(n, pureimag, cross3a);
  TESTREAL3(n, orthog3a1);
  TESTREAL3(n, orthog3a2);
  TESTREAL3(n, orthog3a3);
  TESTREAL3(n, gram3a);
  TESTIDENT3(n, pureimag, cross3b);
  TESTREAL3(n, orthog3b1);
  TESTREAL3(n, orthog3b2);
  TESTREAL3(n, orthog3b3);
  TESTREAL3(n, gram3b);
  if(n > 0){
    TESTIDENT3(n, pureimag, aut3);
    TESTIDENT2(n, pureimag, KSe1);
    TESTIDENT2(n, pureimag, KSe2);
    TESTIDENT2(n, pureimag, KSe3);
    for(BV1=1; BV1<p2n; BV1++){
      BV2 = BV1;
      printf("%d:", BV1); TESTIDENT2(n, pureimag, aut3BV);
    }
  }
  TESTREAL3(n, Lbraid);
  TESTREAL3(n, Rbraid);
  TESTREAL3(n, Lip);
  TESTREAL3(n, Rip);
  TESTREAL4(n, ip4);
  TESTREAL3(n, ip3L);
  TESTREAL3(n, ip3R);
  TESTREAL2(n, ipcheck1);
  TESTREAL2(n, ipcheck2);
  TESTREAL2(n, ipcheck3);
  TESTREAL2(n, ipcheck4);
  TESTREAL2(n, ipcheck5);
  TESTREAL2(n, orthogcross1);
  TESTREAL2(n, orthogcross2);
  TESTREAL2(n, orthogL1);
  TESTREAL2(n, orthogL2);
  TESTREAL1(n, requad);
  TESTIDENT3(n, pureimag, Ebbing1p254);
  TESTIDENT3(n, pureimag, Ebbing6p255);
  TESTIDENT3(n, pureimag, ZSSSlem5p130);
  TESTIDENT3(n, pureimag, ZSSSp145eq13);
  TESTIDENT3(n, pureimag, ZSSSp145eq15);
  TESTIDENT2(n, pureimag, jorcheck1);
  TESTIDENT3(n, pureimag, Wilson);
  TESTIDENT3(n, pureimag, LCC2);
  TESTIDENT2(n, pureimag, RWIP);
  TESTIDENT2(n, pureimag, LWIP);
  TESTIDENT3(n, pureimag, alt1);
  TESTIDENT3(n, pureimag, Lalt2);
  TESTIDENT3(n, pureimag, Ralt2);
  TESTIDENT2(n, pureimag, octLalt2);
  TESTIDENT2(n, pureimag, octRalt2);
  TESTIDENT2(n, pureimag, octLalt1);
  TESTIDENT2(n, pureimag, octRalt1);
  TESTIDENT2(n, pureimag, ninninRalt);
  TESTIDENT2(n, pureimag, ninLalt2);
  TESTIDENT2(n, pureimag, ninRalt2);
  TESTIDENT2(n, pureimag, ninLalt1);
  TESTIDENT2(n, pureimag, ninRalt1);

  TESTIDENT3(n, pureimag, scalcheck1);
  TESTIDENT3(n, pureimag, scalcheck2);
  TESTIDENT3(n, pureimag, scalcheck3);

  TESTIDENT1(n, pureimag, conjlin);
  TESTIDENT1(n, pureimag, conjnonlin);
  TESTIDENT1(n, pureimag, Lconjlin);
  TESTIDENT1(n, pureimag, Rconjlin);
  TESTREAL2(n, foo);
  TESTIDENT2(n, pureimag, simil);
  TESTIDENT2(n, pureimag, sqsimil);
  TESTIDENT2(n, pureimag, sim);
  TESTIDENT2(n, pureimag, sqsim);
  TESTIDENT3(n, pureimag, sim2);
  TESTIDENT3(n, pureimag, sqsim2);
  TESTIDENT3(n, pureimag, sim3);
  TESTIDENT3(n, pureimag, sqsim3);
  TESTIDENT3(n, pureimag, Lsand);
  TESTIDENT2(n, pureimag, Lsand2);
  TESTIDENT2(n, pureimag, Lsand2C);
  TESTIDENT1(n, pureimag, sq);
  TESTIDENT2(n, pureimag, dashL);
  TESTIDENT2(n, pureimag, dashR);
  TESTIDENT2(n, pureimag, assoclie1R);
  TESTIDENT2(n, pureimag, assoclie1L);
  TESTIDENT2(n, pureimag, assoclie2R);
  TESTIDENT2(n, pureimag, assoclie2L);
  TESTIDENT3(n, pureimag, sident1R);
  TESTIDENT3(n, pureimag, sident2R);
  TESTIDENT3(n, pureimag, sident1L);
  TESTIDENT3(n, pureimag, sident2L);
  TESTIDENT4(n, pureimag, gaR);
  TESTIDENT4(n, pureimag, gaL);
  TESTIDENT4(n, pureimag, Kleindiff);
  TESTIDENT4(n, pureimag, KleinIdent);
  fflush(stdout);
  for(be=0; be<p2n; be++){
      bcx = basis(n, be);
      printf("basis(%d,%d): ",n,be);
      TESTIDENT1(n, pureimag, basCDvsme);
  }
  for(be=0; be<p2n; be++){
      bcx = basis(n, be);
      printf("basis(%d,%d): ",n,be);
      TESTIDENT2(n, pureimag, linbasCDvsme);
  }
  if(!pureimag){
    printf("UnitRand Jacobian Samples(%d)=", pow2(n));
    for(k=0; k<10; k++){
       jd = jacsamp(n);
       printf("%.3f ", jd);
     }
    printf("\n");
    fflush(stdout);    
    TESTBOOL(n, imaganticomm);
    fflush(stdout);
    TESTBOOL(n, basiscomm);
    fflush(stdout);
    TESTBOOL(n, basisFlex);
    fflush(stdout);
    TESTBOOL(n, basisRalt);
    fflush(stdout);
    TESTBOOL(n, basisanticomm);
    fflush(stdout);
    if(SLOWSTUFF){
      TESTBOOL(n, basis3cyc);
      fflush(stdout);
      TESTBOOL(n,  basisMoufmid);
      fflush(stdout);
      TESTBOOL(n, basisassoc);
      fflush(stdout);
      TESTBOOL(n, basisantiassoc);
      fflush(stdout);
    }
  }
  return points;
}

/*******************************************************/

void printMouf(uint M){
  assert(M<512);
  if(M&256){
    if(M&16){ printf("\\Inv{a} "); }
  }else{
    if(M&16){ printf("a "); }
  }
  if(M&(16|32)) printf("( ");

  if(M&128){
    if(M&1){  printf("\\Inv{a} "); }
    printf("x ");
    if(M&2){  printf("a "); }
  }else{
    if(M&1){  printf("a "); }
    printf("x ");
    if(M&2){  printf("\\Inv{a} "); }
  }

  printf("\\cdot ");

  if(M&64){
    if(M&4){  printf("\\Inv{a} "); }
    printf("y ");
    if(M&8){  printf("a "); }
  }else{
    if(M&4){  printf("a "); }
    printf("y ");
    if(M&8){  printf("\\Inv{a} "); }
  }

  if(M&(16|32)) printf(") ");
  if(M&256){
    if(M&32){ printf("a "); }
  }else{
    if(M&32){ printf("\\Inv{a} "); }
  }
}

eon evalMouf(uint M, eon a, eon x, eon y){
  eon x1,x2,y1,y2,z,z1,z2,z3,c,q;
  basetype s, ssp;
  assert(M<512);
  q = zero(a.lev);
  c = Conj(a);
  s = ss(a);
  ssp = 1;
  if(M&128){
    if(!(M&(1|2))) return q;
    if(M&1){  x1 = mul(c,x);   ssp *= s; }else x1=x;
    if(M&2){  x2 = mul(x1,a);            }else x2=x1;
  }else{
    if(M&1){  x1 = mul(a,x);             }else x1=x;
    if(M&2){  x2 = mul(x1,c);  ssp *= s; }else x2=x1;
  }
  if(M&64){
    if(!(M&(4|8))) return q;
    if(M&4){  y1 = mul(c,y);   ssp *= s; }else y1=y;
    if(M&8){  y2 = mul(y1,a);            }else y2=y1;
  }else{
    if(M&4){  y1 = mul(a,y);             }else y1=y;
    if(M&8){  y2 = mul(y1,c);  ssp *= s; }else y2=y1;
  }
  z = mul(x2,y2);
  if(M&256){
    if(!(M&(16|32))) return q;
    if(M&16){ z1 = mul(c,z);   ssp *= s; }else z1=z;
    if(M&32){ z2 = mul(z1,a);            }else z2=z1;
  }else{
    if(M&16){ z1 = mul(a,z);             }else z1=z;
    if(M&32){ z2 = mul(z1,c);  ssp *= s; }else z2=z1;
  }
  z3 = divscal( ssp, z2 );
  return z3;
}

uint FindFirstMouf(uint M, eon a, eon x, eon y){
  uint j;
  eon z;
  assert(a.lev==x.lev);
  assert(a.lev==y.lev);
  z = evalMouf(M,a,x,y);
  for(j=0; j<M; j++){
    if( NearEq(z, evalMouf(j,a,x,y)) ) return j;
  }
  return M;
}

void exploreMouf(uint n){
  eon a, x, y;
  eon a2, x2, y2;
  eon a3, x3, y3;
  uint fm3,fm2,fm,M;
  a = Erand(n, -9);
  x = Erand(n, -9);
  y = Erand(n, -9);
  a2 = Erand(n, -9);
  x2 = Erand(n, -9);
  y2 = Erand(n, -9);
  a3 = Erand(n, -9);
  x3 = Erand(n, -9);
  y3 = Erand(n, -9);
  printf("exploring 512 candidate Moufang-like laws for %d-ons\n", pow2(n));
  for(M=0; M<512; M++){
    fm  = FindFirstMouf(M,a,x,y);
    fm2 = FindFirstMouf(M,a2,x2,y2);
    fm3 = FindFirstMouf(M,a3,x3,y3);
    printf("%3d %3d %3d %4d ", fm, fm2, fm3, M); printMouf(M); printf("\n");
  }
   printf("done with exploreMouf\n");
   fflush(stdout);
}

void exploreg1c(){
  eon a, b, c, x, y;
  double n32, n32min;
  uint8 n;
  uint bestmc = 0;
  double ssab;
  n32min = HUGE;
  n = 4;
  a = Erand(n, -99);
  b = Erand(n, -99);
  c = Erand(n, -99);
  x = Erand(5, -99);
  y = Erand(5, -99);
  ssab = ss(a);
  ssab *= ss(b);
  printf("exploring 6561 candidate mul laws \n");
  for (MULCHOICE = 0; MULCHOICE < 6561; MULCHOICE++){
    /* 6561 does not find a 32-on multip norm formula, just as 625 does not */
    /* use this to prevent huge search: */
    if(biggest9dig(MULCHOICE) < 5 && count9dig(MULCHOICE) == 1){
      if(fabs(ssab - ss(mul(a, b))) < 0.003){
        if(fabs(Nmultip(b, c)) < 0.003){
          print9(MULCHOICE);
          printf(" ");
          printopt(MULCHOICE);
          printf(": ");
          printf("\n");
#if HEAVYDATA16
          testidents(n, FALSE);
          testidents(n, TRUE);
#endif
          printf("Nmultip(32)=%f  Nmulsim=%f\n",
                 n32 = (double) Nmultip(x, y), (double) Nmulsim(x, y));
          n32 = fabs(n32);
          if(n32 < n32min){
            bestmc = MULCHOICE;
            n32min = n32;
          }
          fflush(stdout);
        }
      }
    }
  }
  printf("best was MULCHOICE=%d with n32=%g: i.e. ", bestmc, n32min);
  print9(bestmc);
  printf("\n");
}

void exploreg2c(){
  eon a2, b2, c2, a, b, c, x, y, z;
  double n32, n32min;
  double ssab;
  uint8 n;
  uint bestmc = 0;
  uint ct = 0;
  n32min = HUGE;
  n = 4;
  a = Erand(n, -99);
  b = Erand(n, -99);
  c = Erand(n, -99);
  a2 = Erand(n, -99);
  b2 = Erand(n, -99);
  c2 = Erand(n, -99);
  x = Erand(5, -99);
  y = Erand(5, -99);
  z = Erand(5, -99);
  ssab = ss(a);
  ssab *= ss(b);
  printf("exploring 16*28561 candidate mul laws \n");
  for (SWIND = 0; SWIND < 16; SWIND++){
    for (MULCHOICE = 0; MULCHOICE < 28561; MULCHOICE++){
      COMBOMULCHOICE = MULCHOICE + (13 * 13 * 13 * 13) * SWIND;
      if(fabs(ssab - ss(mul(a, b))) < 0.0005){
        if(fabs(Nmultip(c, c2)) < 0.0005){
          if(fabs(Nmultip(a2, b2)) < 0.0005){
            printf("%x.", SWIND);
            print13(MULCHOICE);
            printf(":\n");
            ct++;
#if HEAVYDATA16
            testidents(n, FALSE);
            testidents(n, TRUE);
#endif
            printf("Nmultip=%f  Nmulsim=%f\n",
                   n32 = (double) Nmultip(x, y), (double) Nmulsim(x, y));
            n32 = fabs(n32);
            if(n32 < n32min){
              bestmc = MULCHOICE;
              n32min = n32;
            }
            fflush(stdout);
          }
        }
      }
    }
  }
  printf
      ("ct=%d found: best was MULCHOICE=%d with n32=%g: i.e. ",
       ct, bestmc, n32min);
  print13(bestmc);
  printf("\n");
}

void exploreWDS192(){
  eon a2, b2, c2, a, b, c, x, y;
  double ssab;
  uint8 n;
  assert(COMBOMULCHOICE==WDSeCode);
  n = 5;
  a = Erand(n, -99);
  b = Erand(n, -99);
  c = Erand(n, -99);
  a2 = Erand(n, -99);
  b2 = Erand(n, -99);
  c2 = Erand(n, -99);
  x = Erand(6, -99);
  y = Erand(6, -99);
  ssab = ss(a);
  ssab *= ss(b);
  printf("exploring 192 candidate mul laws \n");
  for (WDSMULCHOICE = 0; WDSMULCHOICE < 192; WDSMULCHOICE++){
    if(fabs(ssab - ss(mul(a, b))) < 0.0005){
      if(fabs(Nmultip(c, c2)) < 0.0005){
        if(fabs(Nmultip(a2, b2)) < 0.0005){
          printf("%d, ", WDSMULCHOICE);
          if(maxE(Lalt(a,b)) < 0.0005&& maxE(Lalt(a2,b2)) < 0.0005){
            printf("%d is Lalt at 32\n", WDSMULCHOICE);
          }else{
            printf("%d Lalt FAILS at 32\n", WDSMULCHOICE);
          }          
          if(fabs(Nmultip(x,y)) < 0.0005){
            printf("%d works at 64 too\n", WDSMULCHOICE);
          }else{
            printf("%d FAILS at 64\n", WDSMULCHOICE);
          }          
          fflush(stdout);
        }
      }
    }
    printf("\n"); fflush(stdout);    
  }
  printf("done with exploreWDS192\n");
}

void findNasty(){
  uint a,b,c;
  eon ea,eb,ec,z;
  for(a=1; a<8; a++){
    ea = basis(3,a);
    for(b=1; b<8; b++){
      eb = basis(3,b);
      for(c=1; c<8; c++){
        ec = basis(3,c);
        z = sub(mul(mul(ea, Conj(eb)), mul(eb, ec)), mul(ea,ec));
        if( ss(z) > 0.9 ) printf("nasty: %d %d %d\n", a,b,c);
      }
    }
  }
}

#if 0
void exploreidents(uint8 n){
  eon q[103525];
  char *name[103525];
  eon v, w, x, y, z;
  uint i, j, k, p;
  assert(n <= LGMAXSIZE);
  v = Erand(n, -99);
  w = Erand(n, -99);
  x = Erand(n, -99);
  y = Erand(n, -99);
  z = Erand(n, -99);
  p = pow2(n);
  printf("searching for idents(%d)...\n", n);
#include "expout.c"
  printf("done computing quantities(%d). Now seeking idents...\n", n);
  for (i = 0; i < 103525; i++){
    for (j = i + 1; j < 103525; j++){
      for (k = 0; k < p; k++){
        if(fabs(q[i].coord[k] - q[j].coord[k]) > 0.00001)
          goto NODEAL;
      }
      printf("q[%d]=%s = %s=q[%d]\n", i, name[i], name[j], j);
    NODEAL:;
    }
  }
}
#endif

eon tay500(eon z){ /* exp(z), maclaurin series up to z^50 */
  int k;
  eon z2,z3;
  z3 = one(z.lev);
  for(k=500; k>0; k--){
    z2 = mul(z, z3);
    z3 = add( one(z.lev), divscal((basetype)k, z2) );
  }
  return z3;
}

void trotter(n){
  uint p2,i,p,p2p;
  eon ep,tt,tts,x,y,xpy, xs, ys;
  printf("test of trotter formula in 2^%d-ons\n", n);
  p2 = pow2(n);  
  for(i=0; i<p2; i++){
    x.coord[i] = (int)(myrand()*19.999 - 9.999);    
    y.coord[i] = (int)(myrand()*19.999 - 9.999);
  }
  x.coord[0] = 2.0;
  y.coord[0] = 2.0;
  x.lev = y.lev = n;
  printeon(x);
  printeon(y);
  xpy = add(x,y);
  ep = tay500(xpy);
  printf("0: "); printeon(ep);
  for(p=1; p<=14; p++){
    p2p = pow2(p);
    xs = divscal( (basetype)p2p, x );
    ys = divscal( (basetype)p2p, y );
    tt = mul( tay500(xs), tay500(ys) );
    for(i=0; i<p; i++){
      tts = sq(tt);
      tt = tts;
    }
    printf("%d: ", p2p); printeon(tt);
  }
}

main(){
  eon a, b, c, d, e, t1, t2, t3;
  uint8 n;
  uint be;
  bool imagonly;
  uint pts, CodeIndex;

  COMBOMULCHOICE = WDScode;
  srand(time(0));
  printmodes();

printf("sample 5twofer before twoferize:\n");
a = Erand(5, -99); printeon(a); 
printf("after twoferize:\n");
b=twoferize(a); printeon(b); 
printf("after zeroize:\n");
c=zeroize(a); printeon(c); 
printf("end of sample 5twofer.\n");

  /*** findNasty(); ***/

  /*******************
   exploreMouf(1);
   exploreMouf(2);
   exploreMouf(3);
   exploreMouf(4);
  ******************/

   /*********
   printf("exploreg2c:\n");
   printf(
     "exploring 13^4=28561 possible ways to modify Cayley-Dickson doubling\n");
   exploreg2c();
   printf("done exploring 28561.\n");
   fflush(stdout);
   exit(0);
  *************/

   /*******
   printf("exploreg1c:\n");
  exploreg1c();
  exit(0);
  ********/

  /*exploreidents(); */

/*****************
   printf("comparative survey:\n");
   for(n=0; n<=5; n++){
     for(CodeIndex=0; CodeIndex<NUMCODES; CodeIndex++){
       if(CodeIndex==0 || n>RECCUTOFF){
         COMBOMULCHOICE=MagicCode[CodeIndex];
         printf("********************************************************\n");
         printf("TESTING %s %dons(%s), code=",
               MagicName[CodeIndex], pow2(n), EonName[n]);
         SWIND = COMBOMULCHOICE/(13*13*13*13);
         printf("%x.",SWIND); print13( COMBOMULCHOICE%(13*13*13*13) );
         printf(", RECCUTOFF=%d:\n", RECCUTOFF);      
         for(imagonly=0; imagonly<=1; imagonly++){
           testidents(n, imagonly);
           fflush(stdout);
         }
       }
     }
   }
**************/

/********************
  printf("CDmul tests:\n");
  for (n = 0; n <= LGMAXSIZE; n++){
    COMBOMULCHOICE = CDcode;
    printf("********************************************************\n");
    printf("TESTING WDSmul %dons(%s), code=%d", pow2(n),
           EonName[n], WDScode);
    printf(", RECCUTOFF=%d:\n", RECCUTOFF);
    for (imagonly = 0; imagonly <= 1; imagonly++){
      testidents(n, imagonly);
      fflush(stdout);
    }
  }
*********************/

/********************
  COMBOMULCHOICE = WDSeCode;
  MULCHOICE = WDSeCode;
  exploreWDS192();
  ****************/

  printf("WDSmul tests:\n");
  COMBOMULCHOICE = WDScode;

#if TROTTERTEST
  for (n = 0; n <= LGMAXSIZE; n++){
    trotter(n);
    trotter(n);
  }
#endif

  /*************************
  testjac();
  ***************************/

  for (n = 0; n <= LGMAXSIZE; n++){
    COMBOMULCHOICE = WDScode;
    printf("********************************************************\n");
    printf("TESTING WDSmul %dons(%s), code=%d", pow2(n),
           EonName[n], WDScode);
    printf(", RECCUTOFF=%d:\n", RECCUTOFF);
    for (imagonly = 0; imagonly <= 1; imagonly++){
	for(TWOFER=0; TWOFER<=2; TWOFER++){
           testidents(n, imagonly);
	}
      fflush(stdout);
    }
  }

  printf("all done.\n");
}


/****
13^4 search on all-levels-same recursion, 8 norm16mul found:
0002,
0003,
0020,
0040,
0100,
0400,
1000,
3000.
*******/

/********
I did another computer search for modified Cayley-Dickson formulae
which, starting from the octonions, generate hexons with
multiplicative euclidean norm.  456976=16*13^4  candidate
formulae were screened, and of these, 120 had multiplicative norm
for hexons, but NONE (if recursed 1 step further) gave a
multiplicative norm for 32-ons.

The unmodified CD formula,   code=0.0000,  is
(a,b)(c,d) = (ac - dConj(b),  cb+Conj(a)d).

Now for each of the 4 terms xy on the right, I consider
either leaving xy alone (code=0)
or replacing it with Conj(y)Conj(x) (code=1)
with the concatenation of these 4 code-bits being a number, in hex.
Example: 5 means, do this conj-swap operation on the 2nd and 4th terms.

Now, for each term (as modified according to the previous paragraph) pq
I consider (where r and s are the remaining two letters from the set
{a,b,c,d}, in alphabetical order)
code   what to do
0  leave it alone
1  p r . rinv q
2  p rinv . r q
3  p s . sinv q
4  p sinv . s q
5  r p . rinv q
6  rinv p . r q
7  s p . sinv q
8  sinv p . s q
9  p r . q rinv
a  p rinv . q r
b  p s . q sinv
c  p sinv . q s
getting a 4-digit (in base 13) code number.  The two code numbers
are then concatentated with a "." between.
EXAMPLE:  code=5.4070 means
(a,b)(c,d) = (a dinv . d c - Conj(b) d,  d c . dinv b + Conj(d) a).

Here are the 120 code numbers (separated by colons :):
0.0002 : 0.0003 : 0.0020 : 0.0040 : 0.0100 : 0.0400 : 0.074b : 0.084c : 0.0992 : 0.0aa2 :
0.1000 : 0.1990 : 0.1aa0 : 0.1c6b : 0.1cc6 : 0.3000 : 0.358a : 0.35a7 : 0.370b : 0.380c :
0.5025 : 0.5105 : 0.5c4a : 0.5cc2 : 0.6026 : 0.6106 : 0.718a : 0.71a7 : 0.7927 : 0.7963 :
0.946b : 0.94c6 : 0.9827 : 0.9863 : 0.b073 : 0.b470 : 0.b54a : 0.b5c2 : 0.c083 : 0.c480 :
1.005a : 1.00b7 : 1.0500 : 1.0a50 : 1.1a53 : 1.3502 : 1.3503 : 1.3510 : 1.3530 : 1.4540 :
1.5c07 : 1.5c70 : 1.7000 : 1.7102 : 1.7103 : 1.7110 : 1.7130 : 1.7220 : 1.790a : 1.7990 :
1.7aa0 : 1.980a : 1.9890 : 1.b507 : 1.b570 : 1.c0b0 : 1.c4b2 : 1.c580 : 4.000a : 4.00c0 :
4.014a : 4.01c2 : 4.044a : 4.04c2 : 4.056b : 4.05c6 : 4.0c8a : 4.0ca7 : 4.10c1 : 4.204a :
4.20c2 : 4.303a : 4.404a : 4.40c2 : 4.50c5 : 4.606b : 4.60c6 : 4.8070 : 4.8173 : 4.8c00 :
4.a005 : 4.a425 : 4.a500 : 4.b07a : 4.c08a : 4.c0a7 : 5.0037 : 5.0073 : 5.0107 : 5.0170 :
5.0407 : 5.0470 : 5.071c : 5.0c02 : 5.0c03 : 5.0c10 : 5.0c30 : 5.2007 : 5.2070 : 5.270c :
5.4007 : 5.4070 : 5.4c00 : 5.70c2 : 5.71c0 : 5.c002 : 5.c003 : 5.c010 : 5.c030 : 5.c400 :

Among these 120, there were no flexible algebras (xy.x=x.yx).
There were no algebras obeying the antiaut property  conj(ab)=conj(b)conj(a).
There were 2 left-alternative algebras 0.1000 & 0.0002 [Conway] (xx.y=x.xy).
There were 3 right-alternative algebras 0.0040 & 0.3000 [Conway]
plus the following weird one:  5.4070 [used as example above].
(I admit to puzzlement about why 2 and 3 are unequal, probably has to do
with eyes-one-way convention, not the-other-way.)
5.4070 in many ways seems a little less attractive than 0.3000, but in some
ways it seems superior. Specifically, it obeys the right-Bol z(xy.x)=(zx.y)x
(also called right Moufang) identity, while 0.3000 doesn't.
************************/



