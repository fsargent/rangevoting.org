#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

     double drand48 (void);  /* random in [0,1) */
     void srand48 (long seedval);

#define real double
#define uint unsigned int

/*********************************************************************
* Implementation of algorithm for sampling from a discrete
* probability N-vector X[1], X[2], ..., X[N].  (N>=1.)
* Runs on O(1) worst case time per sample,
* and uses one integer and one real N-element array for storage.
* Preprocessing consumes O(N) time and temporarily uses one 
* additional integer array (B[0..N+1]) for bookkeeping. 
* X[0] and X[N+1] are also used as sentinels in the Build() algorithm.
*
*  compile with    gcc -O9 -lm -lc  WDSsampler.c -o WDSsampler
****************Warren D. Smith August 2000**************************/

/* returns random sample i=1..N with prob X[i] */
uint WalkerSample(real Y[], uint N, uint A[]){
   uint i; real r;
   /* Let i = random uniform integer from {1,2,...N};  */
   i = 1 + (uint)(N*drand48()); 
   r = drand48();
   if(r > Y[i]){ i=A[i]; }
   return i;
}

/*overwrites X[1..N] by Y[1..N], creates A[1..N]. B[1..N] used for bookkeeping: */
void BuildSampler(real X[], uint N, uint A[], uint B[]){
   uint i,j,k;
   assert(1<=N);
   for(i=1; i<=N; i++){ 
     A[i]=i; B[i]=i; /* initial destins=stay there */
     assert(X[i]>=0.0);
     X[i]=X[i]*N; /* scale probvec */
   }
   B[0]=0; X[0]=0.0; B[N+1]=N+1; X[N+1]=2.0; /* sentinels */
   i=0; j=N+1;
   for(;;){	
      do{ i++; }while( X[B[i]]< 1.0 );  /* find i so X[B[i]] needs more */
      do{ j--; }while( X[B[j]]>=1.0 );  /* find j so X[B[j]] wants less */
      if(i>=j) break;
      k = B[i]; B[i] = B[j]; B[j] = k;   /* swap B[i], B[j] */
   }
   /*** printf("partitioned i=%u j=%u:\n",i,j);
   * assert: X[B[k]]<1.0 for k=1..j, 1.0<=X[B[k]] for k=j+1..N. 
   * for(i=1; i<=N; i++){ 
   *   assert( (i<=j && X[B[i]]<1.0) || (i>j && X[B[i]]>=1.0) );
   *   printf("%3u %u %f\n", i, B[i], X[B[i]]); if(i==j) printf("--------\n");
   * }
   ************************/

   i=j; j++;
   while(i>0){
     while(X[B[j]]<=1.0){ j++; }  /* find j so X[B[j]] needs more */
     assert(X[B[i]]<1.0);         /* meanwhile X[B[i]] wants less */
      if(j>N) break;
      assert(j<=N);
      assert(X[B[j]]>1.0);
      X[B[j]] -= 1.0-X[B[i]];     /* B[i] will donate to B[j] to fix up */
      A[B[i]] = B[j];             
      if( X[B[j]] < 1.0 ){        /* X[B[j]] now wants less so readjust ordering */
         assert(i<j);
         k = B[i]; B[i] = B[j]; B[j] = k;         /* swap B[j], B[i] */
	 j++;
      }else{ i--; }
   }
}


/*******************************************************************
* Test driver, applies above sampling algorithms to Zipfian distributions
*   (kth histogram bar has height proportional to  k^p)
* with parameter p input by user:
********************************************************************/

void WalkerVerify(real Y[], uint N, uint A[], real X[]){
   real Z[10002];
   uint i,j,k;
   real t, sum=0.0;
   for(i=1; i<=N; i++){   Z[i] = Y[i]/N;  }
   for(i=1; i<=N; i++){   Z[A[i]] += (1.0-Y[i])/N;  }
   for(i=1; i<=N; i++){   t = Z[i] - X[i]; if(t<0.0) t = -t; sum += t; }
   printf("Verification of data structure: this value should be 0: %g\n", sum);
}

void normalize(real X[], uint N){
   uint i;
   real sum;
   assert(N>=1);
   sum = 0.0;
   for(i=1; i<=N; i++){ sum += X[i]; }
   assert(sum > 0.0);
   sum = 1.0/sum;
   for(i=1; i<=N; i++){ X[i] *= sum; }
}

void makezipf(real X[], uint N, real p){
   uint i;
   for(i=1; i<=N; i++){
      X[i] = exp( log((double)i) * p );
   }
   normalize(X,N);
}

void randompermute(real X[], uint N){
   uint i,j;
   real t;
   for(i=1; i<N; i++){
      j = i+1+ (uint)( (N-i)*drand48() );
      assert(i<j);
      assert(j<=N);
      t = X[i]; X[i] = X[j]; X[j] = t;      
   }
}
    
void TestSampler(){ /* test driver, accepts seed, N, p as input */
   uint i, j, N, seed;
   real p, t, sum, expval;
   real X[10002];
   real Y[10002];
   uint A[10002];
   uint B[10002];
   uint count[10002];
   printf("please input RandSeed, #points, and Zipf power parameter p:\n");
   scanf("%u",&seed);
   srand48(seed); for(i=0; i<100; i++) drand48();
   printf("seed=%u  rands=%f %f %f\n", seed, drand48(), drand48(), drand48());
   scanf("%u",&N);
   printf("N=%u\n", N);
   assert(N<10002);
   scanf("%lf",&p);
   printf("p=%f\n", p);

   makezipf(X,N,p);
   printf("made zipf:\n");
   for(i=1; i<=N; i++){ printf("%f\n", X[i]); }
   randompermute(X,N);
   printf("randompermuted:\n");
   for(i=1; i<=N; i++){ printf("%3u %f\n", i, X[i]); }

   for(i=1; i<=N; i++){ Y[i]=X[i]; }
   printf("copied.\n");

   BuildSampler(Y,N,A,B);
   printf("built: i Y[i] A[i]:\n");
   for(i=1; i<=N; i++){ printf("%3u %f %u\n", i, Y[i], A[i]); }

   for(i=1; i<=N; i++){ count[i]=0; }
   for(i=100000; i>0; i--){
      j = WalkerSample(Y,N,A);
      count[j]++;
   }

   sum=0.0;
   printf("100000 samples:\n");
   printf("prob   #samples:\n");
   for(i=1; i<=N; i++){ 
      printf("%f %6u\n", X[i], count[i]); 
      expval = X[i]*100000;
      t = expval-count[i];
      sum += t*t/expval;
   }
   sum /= N;
   printf("sumvar=%f (should be about 1)\n", sum);

   WalkerVerify(Y,N,A,X);
   printf("done.\n");
}

main(){ 
   TestSampler();
}
