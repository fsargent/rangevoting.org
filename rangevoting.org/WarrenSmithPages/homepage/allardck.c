/*****
It has been pointed out that STV single-winner elections can be nonmonotonic.
One cite giving an example is
G. Doron \& R. Kronick: 
Single Transferable Vote: An Example of a Perverse Social Choice Function,
American J. Political Science 21 (May 1997) 303-311.
Another is:
Steven J. Brams \& Peter C. Fishburn:
Some logical defects of the single transferable vote,
Chapter 14, pp.  147-151, in
Choosing an Electoral System:
Issues and Alternatives
(Arend Lijphart and Bernard Grofman, eds.)
Praeger, New York 1984.

But how common is this?

Crispin Allard:
Estimating the Probability of Monotonicity Failure in a UK General Election,
Voting matters 5 (January 1996)

Let the three candidates be $A$, $B$ and $C$.
Sufficient conditions for a single-winner STV monotonicity failure 
with 3 candidates be $A,B,C$ are:
\Benumerate
\item
$A$ is ahead of $B$ who is ahead of $C$;
\item
When $C$ is eliminated, his transfers put $B$ ahead of $A$ so $B$ is elected;
\item
If some number of voters switch their relevant preference from $A$ to $C$, 
so that both $A$ and $C$ are ahead of B, then when $B$ is eliminated, 
$A$ goes ahead of $C$, so that $A$ is elected.
\Eenumerate

Writing these conditions down in mathematical terms we get:
\Benumerate
\item
$a > b > c$.
\item
$a < b + \alpha c$,
\item
Some $x>0$ exists so that
$a - x > b$, $c + x > b$, and $a > c + 2x + \beta b$
where
$\alpha = T_{CB} - T_{CA}$, $\beta = T_{BC} - T_{BA}$
where $T_{ij}$ is the proportion of $i$'s votes which transfer to $j$ if $i$ is eliminated.
\Eenumerate
Assuming $a>b>c$, by letting $x$ be slightly larger than $b-c$ we find that
these conditions are equivalent to
\BE
a < b + \alpha c ,\;\;\;
%I get a-b > b-c  which is equiv to a+c>2b, OK:
a+c > 2 b ,\;\;\;
%I get by choosing x = b-c
%that  a>c+2b-2c+beta*b = 2b-c+beta*b  so a+c > (2+beta)b, OK.
a+c > (2+\beta) b .
\EE


the same as the probability that 

a > b > c
a < b + alpha*c
a+c > 2*b
a+c > (2+beta)*b

where beta=y-z, alpha=u-w,
and where the 7 real variables a,b,c,u,w,y,z are chosen
uniformly from the following 4-dimensional subset of R7:
a+b+c=1,
y+z=1,
u+w=1,
u,w,y,z,a,b,c>0.

In fact, nonmonotonicity also can happen in other ways, so this is only a lower 
bound on the true probability.  However, this may be the most important
kind of nonmonotonicity, in which case Allard's lower bound would be a good
approximation.  M.A.E.Dummett in his book "Principles of electoral reform" Cambridge 1997
attacked Allard's estimate as likely to be erroneously too small, although
he gave no evidence to back that accusation up.

The first criticism I have of Allard is that we may as well ASSUME a>b>c by
appropriately renaming the three candidates.  This immediately would increase the
Allard's lower bound by about a factor of 6.

The second criticism I have is, Allard got the wrong number.
He made an inexact geometric estimate and got 0.00025 (without the 6 increase).
I instead attempted to find the exact probability (to 1% accuracy or better)
by using Monte Carlo integration.
My program generates random a,b,c, y,z, u,w>0,  satisfying
a+b+c=1, y+z=1, u+w=1,  and a>b>c.  It then computes
beta=y-z  and  alpha=u-w  and  checks whether
a < b + alpha*c  and   a+c > 2*b   and    a+c > (2+beta)*b
(mono failure) hold. 

My C program (on 6 runs with different random seeds each time) gets:
10419 out of 1000000 experiments are nonmonotonic, fraction=0.010419
10485 out of 1000000 experiments are nonmonotonic, fraction=0.010485
10686 out of 1000000 experiments are nonmonotonic, fraction=0.010686
10580 out of 1000000 experiments are nonmonotonic, fraction=0.010580
10559 out of 1000000 experiments are nonmonotonic, fraction=0.010559
10356 out of 1000000 experiments are nonmonotonic, fraction=0.010356
and in some larger independent runs:
52540 out of 5000000 experiments are nonmonotonic, fraction=0.010508
105038 out of 10000000 experiments are nonmonotonic, fraction=0.010504
so that the correct number is  (1.05 +- 0.01)%  of STV elections (at least)
lead to nonmonotonicity.  This is two orders of magnitude higher than
Allard's wrong number.

Compile with   gcc -lm -lc allardck.c -O9 -o allardck
**************Warren D. Smith***Aug 2004**********/

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

unsigned long RandSeed = 2719467481;

#define VERBOSE (0==1)

main(){
  unsigned long ct, NumExp;
  unsigned int i,j;
  double t,a,b,c,s,u,w,x,y,z,alpha,beta;

  srand48(RandSeed);
  printf("some rands %f %f  %f\n", drand48(), drand48(), drand48() );
  NumExp = 1;
  for(j=0; j<20; j++){
    if(j&1) NumExp *= 2;
    else NumExp *= 5;
    ct=0;
    for(i=0; i<NumExp; i++){
      y = drand48(); z = 1.0-y;
      u = drand48(); w = 1.0-u;
      do{
	a = drand48();
	b = drand48();
	c = drand48();
	s = a+b+c;
      }while(s>1.0 || s==0.0);
      if(a<b){ t=a; a=b; b=t; }
      if(b<c){ t=b; b=c; c=t; }
      if(a<b){ t=a; a=b; b=t; }
      /* now a>b>c holds */
      a /= s; b /= s; c /= s;
      /* now a+b+c=1 holds */
      beta = y-z;
      alpha = u-w;
      if(VERBOSE) printf("a=%f b=%f c=%f alpha=%f beta=%f : ",a,b,c,alpha,beta);
      if( a < b + alpha*c ){
	if( a+c > 2*b ){
	  if( a+c > (2+beta)*b ){
	    ct++;
	    if(VERBOSE) printf("success\n");
	  }else{ if(VERBOSE) printf("3fails\n"); }
	}else{ if(VERBOSE) printf("2fails\n"); }
      }else{ if(VERBOSE) printf("1fails\n"); }
    }
    printf("%.0f out of %.0f experiments are nonmonotonic, fraction=%f\n", 
	   (double)ct, (double)NumExp, ((double)ct)/((double)NumExp) );
  }
}






