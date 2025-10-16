#include <math.h>
#include <stdio.h>

// BiasFinder.c  by  Warren D Smith   Nov 2014
// gcc -Wall -O6 BiasFinder.c

#define uint unsigned int
#define real double
#define popcnt(x) __builtin_popcount(x)

int bias(uint x){
  int i,b;
  i=0; b=0;
  while(x){
    b += (x&1)*i;
    i++; x >>= 1;
  }
  return(b);
}

main(){
  int A, N, i, c, x, b, bx, k, wc;
  uint ew, bw, p, j;
  scanf("%d %d", &N, &A);
  x=0; b=0; k=0; wc=0;
  for(k=0; k<N; ){ //read N-bit word
    c = getchar();
    if(c=='0' || c=='1'){ x *= 2;  x += c-'0'; k++; wc += c-'0'; }
  }
  if(wc!=A){ printf("MORON! #1bits=%d which is NOT %d\n", wc,A); }
  if(k!=N){ printf("MORON! #1bits=%d which is NOT %d\n", wc,A); }
  bx = bias(x);
  bw = 0; ew = 0; p = 1; p = p<<N;
  for(j=0; j<p; j++){
    if(popcnt(j)==A){ 
      ew++; if(bias(j)>=bx){ bw++; }
    }
  }
  printf("Among the %u binary %d-bit words with %d bits '1', the number\n", ew, N, A);
  printf("whose bias score is at least as great as %d\n", bx);  
  printf("(which is the bias of ");
  for(i=N-1; i>=0; i--){
    if((x>>i)&1) putchar('1'); else putchar('0');
  }
  printf("), is %d.\nAs a fraction this is %d/%d=%f.\n", bw, bw, ew, (real)bw/ew);
}
