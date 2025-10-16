/* compile with
 * gcc rhc.c -lm -lc -Wall -O5 -o rhc 
PROGRAM FOR COMPUTING THE c_n MACLAURIN SERIES COEFFS FOR F(z)
FUNCTION RELATED TO RIEMANN ZETA FUNCTION BY  F(z) = ln( z/(1-z) * Zeta(1/(1-z)) )
 ******* Warren D. Smith Feb 2005 **********/

#define TWOPI    6.2831853071795864769252867665590057683943387987502
#define PI       3.1415926535897932384626433832795028841971693993751
#define ZETA3    1.2020569031595942853997381615114499907649862923405
#define ROOT2    1.4142135623730950488016887242096980785696718753769

#define ALLCOEFFS 0
#define CNSQONLY 1
#define LOGPOW 0.0

/*#define MAXP2  524288*/
/*#define MAXP2   1048576*/
/*#define MAXP2      2097152*/
#define MAXP2       4194304 
/*#define MAXP2        8388608 
#define MAXP2   16777216*****/

/* Complex Numbers, taken from Press, Teukolsky, Vettering and Flannery
 * Numerical Recipes in C, 2nd Edition, Cambridge 1992 with a few 
 * changes & additions by WDS. */
#include <stdio.h>
#include <math.h>

double sq(double x){ return x*x; }
double cu(double x){ return x*x*x; }

typedef struct DCMPLX {double r,i;} dcmplx; 

dcmplx Cadd(dcmplx a, dcmplx b)
{
    dcmplx c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

dcmplx Csub(dcmplx a, dcmplx b)
{
    dcmplx c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


dcmplx Cmul(dcmplx a, dcmplx b)
{
    dcmplx c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

dcmplx Complex(double re, double im)
{
    dcmplx c;
    c.r = re;
    c.i = im;
    return c;
}

dcmplx Conjg(dcmplx z)
{
    dcmplx c;
    c.r = z.r;
    c.i = -z.i;
    return c;
}

dcmplx Cdiv(dcmplx a, dcmplx b)
{
    dcmplx c;
    double r,den;
    if (fabs(b.r) >= fabs(b.i)) {
      r   = b.i/b.r;
      den = b.r+r*b.i;
      c.r = (a.r+r*a.i)/den;
      c.i = (a.i-r*a.r)/den;
    }
    else {
      r   = b.r/b.i;
      den = b.i+r*b.r;
      c.r = (a.r*r+a.i)/den;
      c.i = (a.i*r-a.r)/den;
    }
    return c;
}

dcmplx Crecip(dcmplx b)
{
    dcmplx c;
    double r,den;
    if (fabs(b.r) >= fabs(b.i)) {
      r   = b.i/b.r;
      den = b.r+r*b.i;
      c.r = 1.0/den;
      c.i = -r/den;
    }
    else {
      r   = b.r/b.i;
      den = b.i+r*b.r;
      c.r = r/den;
      c.i = -1.0/den;
    }
    return c;
}

double Cabsq(dcmplx z)
{
  return z.r * z.r + z.i * z.i;
}

double Cabs(dcmplx z)
{
    double x,y,ans;
    double temp;
    x = (double)fabs(z.r);
    y = (double)fabs(z.i);
    if (x == 0.0)
      ans  = y;
    else if (y == 0.0)
      ans  = x;
    else if (x > y) {
      temp = (double)(y/x);
      ans  = x*(double)sqrt(1.0+temp*temp);
    }
    else {
      temp = (double)(x/y);
      ans  = y*(double)sqrt(1.0+temp*temp);
    }
    return ans;
}

dcmplx Csqrt(dcmplx z)
{
    dcmplx c;
    double w;
    double x,y,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
      c.r=0.0;
      c.i=0.0;
      return c;
    }
    else {
      x=fabs(z.r);
      y=fabs(z.i);
      if (x >= y) {
        r   = y/x;
        w   = (double)(sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r))));
      }
      else {
        r   = x/y;
        w   = (double)(sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r))));
      }
      if (z.r >= 0.0) {
        c.r = w;
        c.i = z.i/(2.0*w);
      } else {
        c.i = (z.i >= 0.0) ? w : -w;
        c.r = z.i/(2.0*c.i);
      }
      return c;
    }
}

dcmplx RCmul(double x, dcmplx a)
{
        dcmplx c;
        c.r=x*a.r;
        c.i=x*a.i;
        return c;
}

dcmplx RCadd(double x, dcmplx a)
{
        dcmplx c;
        c.r=x+a.r;
        c.i=  a.i;
        return c;
}

dcmplx Cexp(dcmplx z)
{
        dcmplx c;
	double e;
	e = exp(z.r);
	c.r = cos(z.i)*e;
        c.i = sin(z.i)*e;
        return c;
}

static int argadjust;
dcmplx Clog(dcmplx z, double oldarg)
{
        dcmplx c;
	double numrots;
	int nri;
	c.r = log( Cabs(z) );
	c.i = atan2( z.i, z.r );
	numrots = (oldarg - c.i)/TWOPI;
	if(numrots > 0.51){
	  nri = numrots+0.5;
	  c.i += TWOPI * nri;
	  argadjust++; 
	}else if(numrots < -0.51){
	  nri = 0.5-numrots;
	  c.i -= TWOPI * nri;	
	  argadjust++; 
	}
        return c;
}

dcmplx RCpow(double p, dcmplx z)  /* p^z where p>0. */
{
	return Cexp( RCmul( log(p), z ) );
}

dcmplx RCnegpow(double p, dcmplx z)  /* p^(-z) where p>0. */
{
	return Cexp( RCmul( -log(p), z ) );
}

dcmplx CCircle(double theta)
{
	return Complex(cos(theta), sin(theta));
}

static double emz[] = { (1/12.0),
  -(1/720.0),
  (1/30240.0),
  -(1/1209600.0),
  (1/47900160.0),
  -(691/1307674368000.0),
  (1/74724249600.0),
  -(3617/10670622842880000.0),
  (43867/5109094217170944000.0),
  -(174611/802857662698291200000.0),
  (77683/14101100039391805440000.0),
  -(236364091/1693824136731743669452800000.0),
  (657931/186134520519971831808000000.0),
 -(3392780147/37893265687455865519472640000000.0) };

/***This Zeta function routine based on Euler-Maclaurin is
not used except for testing purposes and has insufficient precision.  Anyhow,
if one were to use Euler-Maclaurin it would be better to use it on the alternating
sum (-1)^n/n^s  by regarding each consecutive pair of terms as "one term".
That idea for some reason was not mentioned in the literature, which is why
I am mentioning it here.
 ***/
dcmplx EulerMacZeta(dcmplx s)  
{
  int k;
  dcmplx t;

  t = Complex(0.0, 0.0);
  for(k=13; k>=1; k--){
    t = RCadd(emz[k], t);
    t = RCmul( 0.01, Cmul( Cadd(Complex(2.0*k-1.0, 0.0), s), Cadd(Complex(2.0*k, 0.0), s) ) );
  }
  t = RCadd(emz[0], t);
  t = Cmul( RCmul(0.1, s), t );
  t = Cadd( Cdiv( Complex(10.0, 0.0),  Csub( s, Complex(1.0, 0.0) ) ), t );
  t = RCadd( 0.5, t );
  t = Cmul( RCnegpow(10, s), t );
  for( k=9; k>=2; k--){
    t = Cadd( RCnegpow((double)k, s),  t );
  }
  t = RCadd( 1.0,  t );
  return t;
}



/*** An efficient algorithm for the Riemann zeta function 
 *** Peter Borwein  http://www.cecm.sfu.ca/preprints/1995pp.html  95:043
 *******/
dcmplx BorweinZeta(dcmplx s)  /* s.r > 0  is assumed. Actually works if s.r > -(n-1). */
{
  int n,n1,n2,j,k,corrneeded;
	static int bwrecn = 0;
	dcmplx term, scalefac, asum, csum;
	double binom, rp2, mysig, bwe;
	if(s.r==0.0 && s.i==0.0) return Complex(-0.5, 0.0);
	/****
	if(s.r < 0.0){
	  printf("ERROR: BorweinZeta called on s=%g,%g with negative real part\n", s.r, s.i);
	  exit(1);
	  }***/
	/* find n needed: */
	corrneeded = 1;
	n = 0.4808983470 * log(2.9 + 2.0 * fabs(s.i));
	if(s.r<0.0) n += 0.666666667 * fabs(s.r);
	n += 21.4 + 0.7553933572 * fabs(s.i);
	/* 21.4 is for 64-bit reals.  For 128-bit reals, double this to 42.8. */
	if(s.r > 7.4){
	  n1 = exp(41.6/s.r);
	  if(n > n1*0.45){ n = n1; corrneeded = 0; }
	}

	if(n > bwrecn){ 
	  bwrecn = n; 
	  printf("BorweinZeta record n=%d s=%g,%g\n", bwrecn, s.r, s.i); 
	}

	csum = Complex(0.0, 0.0);
	if( corrneeded ){
	  rp2 = 1.0;
	  for(j=1; j<n; j += 2){ rp2 *= 2.0; } 
	  if(n&1){ rp2 *= ROOT2; }  /* rp2 = sqrt(2^n). */
	  /* compute correction sum 1..n: */
	  mysig = -1.0;
	  bwe = binom = 1.0/rp2;
	  for(k=n+n; k>n; k--){
	    term = RCmul( mysig*bwe, RCnegpow( (double)k, s ) );
	    binom *= k-n;
	    binom /= n+n+1-k;  /* binom = (n choose k-n) / rp2. */
	    bwe += binom;
	    mysig = -mysig;
	    csum = Cadd( csum, term );
	  }
	  csum = RCmul( 1.0/rp2, csum );
	}

	/* compute alternating sum 1..n: */
	mysig = -1.0;
	if(n&1) mysig = 1.0;  /* sig = (-1)^(n-1). */
	asum = csum;
	for(k=n; k>1; k--){
	  term = RCmul( mysig, RCnegpow( (double)k, s ) );
	  mysig = -mysig;
	  asum = Cadd( asum, term );
	}
	/* compute scalefac = 1 - 2^(1-s): */
	scalefac = Csub( Complex(1.0, 0.0), RCpow( 2.0, Csub( Complex(1.0, 0.0), s ) ) );
	return Cdiv(Cadd(Complex(1.0, 0.0), asum), scalefac);
}

/*** FFT of nn complex data points 
 ** (data[1], data[2]) = first complex (real, imag),
 ** (data[3], data[4]) = second complex, etc.
 ** Use isign = +1 for forward and -1 for inverse FFT (if inverse must divide by nn afterwards).
 *********/
void four1(data, nn, isign)
double data[];
int nn, isign;
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}

/* Calculates FFT of 2n real data values data[1..2n], n=power of 2.  (If isign = +1.)
 * Replaces by positive freq half of complex FT.
 * The real-valued first and last components of the complex transform are returned
 * as data[1] and data[2].  
 * Same routine computes inverse transform (if is of real data); in that case
 * must post-divide result by n and use isign = -1.
 *  Based on Press et al Numerical recipes 1992.
 *****************************/
void realft(data,n,isign)
double data[];
int n,isign;
{
    int i, i1, i2, i3, i4, n2p3;
    double c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    double wr, wi, wpr, wpi, wtemp, theta;
    void four1();

    theta = PI/(double) n;
    if (isign == 1) {
	c2 = -0.5;
	four1(data, n, 1);
    } 
    else {
	c2 = 0.5;
	theta = -theta;
    }
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;
    n2p3 = 2*n+3;
    for (i = 2; i <= n/2; i++) {
	i4 = 1 + (i3 = n2p3 - (i2 = 1 + ( i1 = i + i - 1)));
	h1r =  c1*(data[i1] + data[i3]);
	h1i =  c1*(data[i2] - data[i4]);
	h2r = -c2*(data[i2] + data[i4]);
	h2i =  c2*(data[i1] - data[i3]);
	data[i1] =  h1r + wr*h2r - wi*h2i;
	data[i2] =  h1i + wr*h2i + wi*h2r;
	data[i3] =  h1r - wr*h2r + wi*h2i;
	data[i4] = -h1i + wr*h2i + wi*h2r;
	wr = (wtemp = wr)*wpr - wi*wpi+wr;
	wi = wi*wpr + wtemp*wpi + wi;
    }
    if (isign == 1) {
	data[1] = (h1r = data[1]) + data[2];
	data[2] = h1r - data[2];
    } else {
	data[1] = c1*((h1r = data[1]) + data[2]);
	data[2] = c1*(h1r - data[2]);
	four1(data, n, -1);
    }
}

double cheapo(int k)
{
  return( sqrt(k-1.0) * pow( log((double)k), LOGPOW ) );
}

void DoIt(int p2, double r)
{
  dcmplx z,w,zw,f;
	double rpow, rrecip, recmax, recmin, recsmax, recsmin, recslmax, recslmin;
	double sqsum, t, maxdat, oldarg; 
	int k, p2h, kmin,kmax, ksmin, ksmax, kslmin, kslmax, kmd, maxct, nextp2;
	static double wirec = 0.0;
	static double w1rec = 0.0;
	static double data[MAXP2+4];
	p2h = p2/2;
	maxdat = 0.0;
	oldarg = 0.0;
	argadjust = 0;
	for(k=0; k<=p2h; k++){
	  z = RCmul(r, CCircle(TWOPI * k/p2));
	  w = Crecip( Csub( Complex(1.0, 0.0), z ) );  /* w = 1/(1-z). */
	  if( w.i > wirec ){  wirec = w.i; printf("ComputeFValues record w.i=%g\n", wirec); }
	  if( w.i > w1rec && w.r <= 1.0 ){  
	    w1rec = w.i; printf("ComputeFValues record w1=%g\n", w1rec); }
	  zw = Cmul(z, w);
	  f = Cmul( zw, BorweinZeta(w) );
	  f = Clog(f, oldarg);  /* watch out for branch!! */
	  oldarg = f.i;
	  data[k+k+1] = f.r;
	  data[k+k+2] = f.i;
	  t = Cabsq(f);
	  if(maxdat < t){ maxdat = t; kmd = k; }
	}
	if( fabs( data[2] ) > 0.0001 ){
	  printf("ERROR: first data point not real?  %g\n", data[2]);
	}
	data[2] = f.r;
	if( fabs( f.i ) > 0.0001 ){
	  printf("ERROR: last data point not real?  %g\n", f.i);
	}

	sqsum = 0.0;
        for(k=1; k<=p2; k++){
	  t = sq(data[k]);
	  sqsum += t;
	}
	sqsum /= p2h;

	realft(data, p2h, -1);

	rpow = 1.0;
	rrecip = 1.0/r;
        for(k=1; k<=p2h; k++){
	  data[k] *= rpow/p2h;
	  rpow *= rrecip;
	}
	printf("p2=%d  rpow=%g argadjust=%d\n", p2, rpow, argadjust);

	printf("        n       c_n      c_n*sqrt(n)   c_n*sqrt(n)*log(n+1)^3\n");
#if ALLCOEFFS
        for(k=1; 7*k<=p2h/2; k++){
#if CNSQONLY
	  printf("%20.16f\n", data[k]*cheapo(k) );
#else
	  printf("%5d %20.16f %20.16f  %20.16f\n", k-1, data[k], data[k]*cheapo(k),
		 data[k]*sqrt(k-1.0)*cu(log((double)k))  );
#endif

	}
#else
	recmax  = data[2]; /*0.57721566490153286;*/  kmax=2;
	recsmax = -999.;
	recsmin =  999.; 
	recslmax = -999.;
	recslmin =  999.; 
	recmin  =  999.;
	maxct=0;
	nextp2 = 32;
        for(k=1; k<11 && k<p2h; k++){
	    printf("    %5d  %20.16f  %20.16f  %20.16f\n", k-1, data[k], data[k]*cheapo(k),
		   data[k]*sqrt(k-1.0)*log((double)k)  );
	    if( data[k]*cheapo(k) > recsmax ){ recsmax = data[k]*cheapo(k); ksmax = k; }
	    if( data[k]*sqrt(k-1.0)*cu(log((double)k)) > recslmax ){ recslmax = data[k]*sqrt(k-1.0)*cu(log((double)k)); kslmax = k; }
	}
        for(k=10; 7*k<p2h/2; k++){
	  if(data[k]>data[k-1] && data[k]>data[k+1]){
	    maxct++;
	    if(k>=nextp2){ 
	      nextp2 += nextp2;
	      printf("maxct=%d k=%d period=%.2f fit=%.2f\n", maxct, k-1, (k-2.0)/maxct, 55+6.6*log(k-1.0)); 
	    }
	    printf("%d max %5d %20.16f %20.16f ", maxct, k-1, data[k], data[k]*cheapo(k) );
	    if( data[k] > recmax ){ recmax = data[k]; kmax = k; }
	    if( data[k]*cheapo(k) > recsmax ){ recsmax = data[k]*cheapo(k); ksmax = k; }
	    if( data[k]*sqrt(k-1.0)*cu(log((double)k)) > recslmax ){ 
	      recslmax = data[k]*sqrt(k-1.0)*cu(log((double)k)); kslmax = k; 
	      printf(" recL\n");
	    }else{ printf("\n"); }
	  }
	  if(data[k]<data[k-1] && data[k]<data[k+1]){
	    printf("min %5d %20.16f %20.16f", k-1, data[k], data[k]*cheapo(k) );
	    if( data[k] < recmin ){ recmin = data[k]; kmin = k; }
	    if( data[k]*cheapo(k) < recsmin ){ recsmin = data[k]*cheapo(k); ksmin = k; }
	    if( data[k]*sqrt(k-1.0)*cu(log((double)k)) < recslmin ){ 
	      recslmin = data[k]*sqrt(k-1.0)*cu(log((double)k)); kslmin = k; 
	      printf(" recL\n");
	    }else{ printf("\n"); }
	  }
	}
	printf("sqsum >= %20.16f using %d points at radius=%20.16f \n", sqsum, p2, r);
	printf("maxdat = %20.16f using %d points at radius=%20.16f, k=%d \n", maxdat, p2, r, kmd);
	printf("recmax = %20.16f at %d\n", recmax, kmax-1);
	printf("recmin = %20.16f at %d\n", recmin, kmin-1);
	printf("recsmax = %20.16f at %d\n", recsmax, ksmax-1);
	printf("recsmin = %20.16f at %d\n", recsmin, ksmin-1);
	printf("recslmax = %20.16f at %d\n", recslmax, kslmax-1);
	printf("recslmin = %20.16f at %d\n", recslmin, kslmin-1);
#endif
}

void CheckZeta(){
  dcmplx z,y;
  int j,k;
  double q, qmax;
  printf("performing checks:\n");
  z = BorweinZeta( Complex(0.0, 1.0) );
  printf("zeta(i) = %.18g,%.18g\nZeta(i) = %.18g,%.18g\n", z.r, z.i, 
	 0.0033002236853241028742171142101345659714896472402784,
	 -0.41815544914132167668927423984336106083595018690104 );
  z = BorweinZeta( Complex(0.5, 0.0) );
  printf("zeta(0.5) = %.18g,%.18g\nZeta(0.5) = %.18g\n", z.r, z.i, 
	 -1.4603545088095868128894991525152980124672293310126 );
  z = BorweinZeta( Complex(0.5, 1.0) );
  printf("zeta(0.5+i) = %.18g,%.18g\nZeta(0.5+i) = %.18g,%.18g\n", z.r, z.i, 
	 0.14393642707718906032438966648372157903562010555575, 
        -0.72209974353167308912617513458032492501318439535370 );
  z = BorweinZeta( Complex(1.0, 1.0) );
  printf("zeta(1+i) = %.18g,%.18g\nZeta(1+i) = %.18g,%.18g\n", z.r, z.i, 
     0.58215805975200364819946316791425920187798931682653, 
	 -0.92684856433080707653642431391750077405345489387394 );
  z = BorweinZeta( Complex(1.5, 0.0) );
  printf("zeta(1.5) = %.18g,%g\nZeta(1.5) = %.18g\n", z.r, z.i,  2.6123753486854883433485675679 );
  z = BorweinZeta( Complex(2.0, 0.0) );
  printf("zeta(2) = %.18g,%g\n pi^2/6 = %.18g\n", z.r, z.i, PI*PI/6);
  z = BorweinZeta( Complex(3.0, 0.0) );
  printf("zeta(3) = %.18g,%g\nZeta(3) = %.18g\n", z.r, z.i, ZETA3);
  z = BorweinZeta( Complex(4.0, 0.0) );
  printf("zeta(4) = %.18g,%g\npi^4/90 = %.18g\n", z.r, z.i, PI*PI*PI*PI/90);
  z = BorweinZeta( Complex(0.0, 0.0) );
  printf("zeta(0) = %.18g,%.18g\nZeta(0) = %.18g,%.18g\n", z.r, z.i, -0.5, 0.0 );
  z = BorweinZeta( Complex(-1.0, 0.0) );
  printf("zeta(-1) = %.18g,%.18g\nZeta(-1) = %.18g,%.18g\n", z.r, z.i,  -0.08333333333, 0.0 );
  z = BorweinZeta( Complex(-2.0, 0.0) );
  printf("zeta(-2) = %.18g,%.18g\nZeta(-2) = %.18g,%.18g\n", z.r, z.i,  0.0, 0.0 );
  z = BorweinZeta( Complex(-3.0, 0.0) );
  printf("zeta(-3) = %.18g,%.18g\nZeta(-3) = %.18g,%.18g\n", z.r, z.i,   0.008333333333 , 0.0 );
  z = BorweinZeta( Complex(-4.0, 0.0) );
  printf("zeta(-4) = %.18g,%.18g\nZeta(-4) = %.18g,%.18g\n", z.r, z.i,   0.0 , 0.0 );
  z = BorweinZeta( Complex(-2.5, 0.0) );
  printf("zeta(-2.5) = %.18g,%.18g\nZeta(-2.5) = %.18g,%.18g\n", z.r, z.i, 0.00851692877785, 0.0 );
  z = BorweinZeta( Complex(0.5, 6.0) );
  printf("zeta(0.5+6i) = %.18g,%.18g\nZeta(0.5+6i) = %.18g,%.18g\n", z.r, z.i,
	 0.837223808066879519392085379,    0.34021839694376641529294 );
  z = BorweinZeta( Complex(0.5, 9.0) );
  printf("zeta(0.5+9i) = %.18g,%.18g\nZeta(0.5+9i) = %.18g,%.18g\n", z.r, z.i,
	 1.4476424519337558111, 0.191803012762665405488);
  z = BorweinZeta( Complex(0.5, 9.5) );
  printf("zeta(0.5+9.5i) = %.18g,%.18g\nZeta(0.5+9.5i) = %.18g,%.18g\n", z.r, z.i,
	 1.51716040681491, 0.0534140302686796 );
  z = BorweinZeta( Complex(0.5, 7.5) );
  printf("zeta(0.5+7.5i) = %.18g,%.18g\nZeta(0.5+7.5i) = %.18g,%.18g\n", z.r, z.i,
	 1.1291091838090642 , 0.392406583197541295 );
  z = BorweinZeta( Complex(1.0, 9.0) );
  printf("zeta(1+9i) = %.18g,%.18g\nZeta(1+9i) = %.18g,%.18g\n", z.r, z.i,
	 1.3396141551524297978,     0.1224817717747397 );
  z = BorweinZeta( Complex(1.0, 9.5) );
  printf("zeta(1+9.5i) = %.18g,%.18g\nZeta(1+9.5i) = %.18g,%.18g\n", z.r, z.i,
	 1.381615485927153574641,     0.0153775951333729900 );

  printf("WDS SUMMARY: seems ok now.\n");
  printf("  Perhaps Borwein error estimates were a bit wrong;\n");
  printf("  & I have added slight extra safety factors inside BorweinZeta that fix problems.\n");

  printf("now for 400 checks on a grid:\n");
  qmax = -1.0;
  for(j=1; j<20; j++){
    for(k=0; k<20; k++){
      z = Complex( j*0.5, k*0.5 );
      y = BorweinZeta(z);
      z = EulerMacZeta(z);
      q = Cabs( Csub(y, z) );
      printf(" (%.1f,%.1f), %g -- Bor=%f,%f,  Eul=%f,%f\n", j*0.5, k*0.5, q, y.r, y.i, z.r, z.i );
      if(q > qmax){
	qmax = q;
	printf("Zeta(%.1f,%.1f):  |diff| = %g\n", j*0.5, k*0.5, qmax);
      }
    }
  }
  printf("checks done.\n\n");
  printf("LOGPOW = %g\n\n", LOGPOW);
}

main(){
  int p2;
  double r, offset;

  CheckZeta();
  for(p2 = 8; p2<=MAXP2; p2 += p2){
    offset = sqrt(4.0/p2);
    if(p2 > 512){
      offset = 285./p2;
    }
    if(p2 > 3000000){
      offset = 500./p2;  /* attempt to get the 1st 10^5 coeffs, fails? */
    }
    r = 1.0 - offset;
    printf("running p2=%d  r=%g\n", p2, r);
    DoIt(p2, r);
  }
}

