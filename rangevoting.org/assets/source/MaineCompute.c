#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

// Author: Warren D. Smith, 13 Nov. 2014.  Two BUGS found & corrected 14 Nov.
// Compile:   gcc -Wall -O6 MaineCompute.c
// Run:       ./a.out < MaineBallots.txt
// Run to produce output file:   ./a.out < MaineBallots.txt > outfile
// Cheapo analysis of output file:  this unix 1-liner counts WGT elections of type LCMC:
//    cat outfile | fgrep WGT | fgrep LCMC | wc -l

#define uint uint32_t
#define real double
#define pf printf
#define rand() drand48()

#define MXBALL 999

real fuzz(){ return( (rand()-0.5)*0.0000001 ); }

char who(int i){
  switch(i){
  case 1: return('C');
  case 2: return('L');
  case 3: return('M');
  case 4: return('X');
  }
  return('_');
}

int whichcand(char c){
  if(c=='C') return(1);
  if(c=='L') return(2);
  if(c=='M') return(3);
  if(c=='X') return(4);
  return(0);
}

uint offT, offC, offL, offM, RandSeed, NumBallots, SubsetSz, NumBoots;
real percC, percL, percM;
real WT[6];
uint PV[MXBALL], AV1[MXBALL], AV2[MXBALL], AV3[MXBALL], 
  I1V[MXBALL], I2V[MXBALL], I3V[MXBALL], GEO[MXBALL];
uint ptot[5];

int ReadBallots(){  // assumes MaineBallots.txt fed in to stdin
  int line, field, c, w, n, nc,np,nr,ns, plur, app1, app2, app3, irv1, irv2, irv3, j, b, b2;
  pf("Reading ballots.\n");
  scanf("%u\t%u\t%u\t%u\t%u", &RandSeed, &NumBoots, &offC, &offL, &offM);
  offT= offC+offL+offM;
  percC = (100.0*offC)/offT;
  percL = (100.0*offL)/offT;
  percM = (100.0*offM)/offT;
  pf("RandSeed=%u\tNumBoots=%u\toffC=%u\toffL=%u\toffM=%u\toffT=%u\nOffic.Percentages:\tC=%.2f\tL=%.2f\tM=%.2f\n",
     RandSeed, NumBoots, offC, offL, offM, offT, percC, percL, percM);
  while( (c=getchar())!='*' ){;}
  getchar();
  b=0; b2=0; ptot[0]=0; ptot[1]=0; ptot[2]=0; ptot[3]=0; ptot[4]=0;
  pf(
     "            Plur Appv Ranked GeographicLocation\n"
     "LINE#  #voters: P CLM x>y>z p+c+s+r\n"
     "========================================\n");
  for(line=1; line<=56; line++){
    np=0; ns=0; nc=0; nr=0; irv1=0; irv2=0; irv3=0; app1=0; app2=0; app3=0; plur=0;
    for(field=0;;){
      c = getchar();
      if(c=='\n'){ field=0; break; }
      if(c=='\t'){ field++; }
      else{
	w=whichcand(c);
	switch( field ){
	case 0:
	  plur=w;
	  break;
	case 1:
	  if(w==1) app1=1; else app1=0;
	  break;
	case 2:
	  if(w==2) app2=1; else app2=0;
	  break;
	case 3:
	  if(w==3) app3=1; else app3=0;
	  break;
	case 4:
	  irv1=w;
	  break;
	case 5:
	  irv2=w;
	  break;
	case 6:
	  irv3=w;
	  break;
	case 7:
	  np = np*10 + c-'0';
	  break;
	case 8:
	  nc = nc*10 + c-'0';
	  break;
	case 9:
	  ns = ns*10 + c-'0';
	  break;
	case 10:
	  nr = nr*10 + c-'0';
	  break;
	default:
	  pf("HUH???\n");
	}
      }
    } //end of for(field)
    n = np+nc+ns+nr; b2+=n;
    pf("line(%2d)   %3d: %c %d%d%d %c>%c>%c %d+%d+%d+%d\n", 
       line, n, who(plur), app1, app2, app3, who(irv1), who(irv2), who(irv3), np, nc, ns, nr);
    ptot[plur] += n;
    for(j=0; j<np; j++){
      PV[b]=plur; AV1[b]=app1; AV2[b]=app2; AV3[b]=app3; I1V[b]=irv1; I2V[b]=irv2; I3V[b]=irv3;
      GEO[b]=0; b++;
    }
    for(j=0; j<nc; j++){
      PV[b]=plur; AV1[b]=app1; AV2[b]=app2; AV3[b]=app3; I1V[b]=irv1; I2V[b]=irv2; I3V[b]=irv3;
      GEO[b]=1; b++;
    }
    for(j=0; j<ns; j++){
      PV[b]=plur; AV1[b]=app1; AV2[b]=app2; AV3[b]=app3; I1V[b]=irv1; I2V[b]=irv2; I3V[b]=irv3;
      GEO[b]=2; b++;
    }
    for(j=0; j<nr; j++){
      PV[b]=plur; AV1[b]=app1; AV2[b]=app2; AV3[b]=app3; I1V[b]=irv1; I2V[b]=irv2; I3V[b]=irv3;
      GEO[b]=3; b++;
    }
  } //end of for(line)
  pf("Read %d=%d ballots.\n", b2,b);
  pf("Raw plurality totals: _=%u\tC=%u\tL=%u\tM=%u\tX=%u\n", 
	 ptot[0],ptot[1],ptot[2],ptot[3],ptot[4]);
  WT[1] = percC/(ptot[1]+0.5*ptot[4]+ptot[0]/3.0);
  WT[2] = percL/(ptot[2]+ptot[0]/3.0);
  WT[3] = percM/(ptot[3]+0.5*ptot[4]+ptot[0]/3.0);
  WT[4] = 0.5*(WT[1]+WT[3]);
  WT[0] = (WT[1]+WT[2]+WT[3])/3.0;
  pf("Weights: _=%f\tC=%f\tL=%f\tM=%f\tX=%f\n", WT[0],WT[1],WT[2],WT[3],WT[4]);
  return(b);
}

uint avtot1, avtot2, avtot3, pvtot[5], irvtot1[5];
uint pairw[5][5];
real Wavtot1, Wavtot2, Wavtot3, Wpvtot[5], Wirvtot1[5];
real Wpairw[5][5];
int AVwin, PVwin, IRV1loss, IRVwin;
int WAVwin, WPVwin, WIRV1loss, WIRVwin;

DoElection(int N, int S){
  int Nx,Sx,i,j,k,ct,  inset[MXBALL]={0};
  real r,w;
  avtot1 = 0;  avtot2 = 0;  avtot3 = 0;
  Wavtot1 = 0;  Wavtot2 = 0;  Wavtot3 = 0;
  AVwin=0;PVwin=0;IRV1loss=0;IRVwin=0;
  WAVwin=0;WPVwin=0;WIRV1loss=0;WIRVwin=0;
  for(i=0; i<5; i++){ pvtot[i]=0; irvtot1[i]=0; Wpvtot[i]=0; Wirvtot1[i]=0; }
  for(j=0; j<5; j++){ for(k=0; k<5; k++){ pairw[j][k] = 0; Wpairw[j][k] = 0.0; }}
  Nx=N; Sx=S; ct=0;
  for(i=0; i<N; i++){
    r = rand();
    if(Nx<=0){ break; }
    if(Nx*r>Sx){ 
      inset[i]=0;  Nx--;
    }else{ //include ballot i in random subset of S out of the N total ballots
      inset[i]=1;
      ct++;
      pvtot[PV[i]]++;
      avtot1 += AV1[i];
      avtot2 += AV2[i];
      avtot3 += AV3[i];
      irvtot1[I1V[i]]++;
      for(j=1; j<=3; j++){ for(k=1; k<=3; k++){
	  if( k!=j && (I1V[i]==j || (I2V[i]==j && I1V[i]!=k) ||
			(I3V[i]==j && I2V[i]==0 && I1V[i]==0)) ){
	    pairw[j][k]++;
	  }
      }}
      Nx--; Sx--;
    }
  }
  //pf("ct=%d  S=%d  N=%d\n", ct,S,N);
  assert(ct==S);
  if(avtot1>avtot2 && avtot1>avtot3){ AVwin=1; }
  if(avtot2>avtot1 && avtot2>avtot3){ AVwin=2; }
  if(avtot2==avtot1 && avtot2>avtot3){ AVwin=(rand()>0.5)?2:1; }
  if(avtot3>avtot1 && avtot3>avtot2){ AVwin=3; }
  if(pvtot[1]>pvtot[2] && pvtot[1]>pvtot[3]){ PVwin=1; }
  if(pvtot[2]>pvtot[1] && pvtot[2]>pvtot[3]){ PVwin=2; }
  if(pvtot[3]>pvtot[1] && pvtot[3]>pvtot[2]){ PVwin=3; }
  if(irvtot1[1]<irvtot1[2] && irvtot1[1]<irvtot1[3]){ IRV1loss=1; }
  if(irvtot1[2]<irvtot1[1] && irvtot1[2]<irvtot1[3]){ IRV1loss=2; }
  if(irvtot1[3]<irvtot1[1] && irvtot1[3]<irvtot1[2]){ IRV1loss=3; }
  if(IRV1loss==1 && pairw[2][3]>pairw[3][2]+fuzz() ){ IRVwin=2; }
  if(IRV1loss==1 && pairw[2][3]<pairw[3][2]+fuzz() ){ IRVwin=3; }
  if(IRV1loss==2 && pairw[1][3]>pairw[3][1] ){ IRVwin=1; }
  if(IRV1loss==2 && pairw[1][3]<pairw[3][1] ){ IRVwin=3; }
  if(IRV1loss==2 && pairw[1][3]==pairw[3][1] ){ IRVwin=(rand()>0.5)?3:1; }
  if(IRV1loss==3 && pairw[2][1]>pairw[1][2]+fuzz() ){ IRVwin=2; }
  if(IRV1loss==3 && pairw[2][1]<pairw[1][2]+fuzz() ){ IRVwin=1; }

  WT[1] = percC/(pvtot[1]+0.5*pvtot[4]+pvtot[0]/3.0);
  WT[2] = percL/(pvtot[2]+pvtot[0]/3.0);
  WT[3] = percM/(pvtot[3]+0.5*pvtot[4]+pvtot[0]/3.0);
  WT[4] = 0.5*(WT[1]+WT[3]);
  WT[0] = (WT[1]+WT[2]+WT[3])/3.0;
  for(i=0; i<N; i++){
    if(inset[i]){
      w = WT[PV[i]];
      Wpvtot[PV[i]] += w;
      Wavtot1 += AV1[i]*w;
      Wavtot2 += AV2[i]*w;
      Wavtot3 += AV3[i]*w;
      Wirvtot1[I1V[i]] += w;
      for(j=1; j<=3; j++){ for(k=1; k<=3; k++){
	  if( k!=j && (I1V[i]==j || (I2V[i]==j && I1V[i]!=k) ||
		       (I3V[i]==j && I2V[i]==0 && I1V[i]==0)) ){
	    Wpairw[j][k] += w;
	  }
      }}
    }
  }
  if(Wavtot1>Wavtot2 && Wavtot1>Wavtot3){ WAVwin=1; }
  if(Wavtot2>Wavtot1 && Wavtot2>Wavtot3){ WAVwin=2; }
  if(Wavtot2==Wavtot1 && Wavtot2>Wavtot3){ WAVwin=(rand()>0.5)?2:1; }
  if(Wavtot3>Wavtot1 && Wavtot3>Wavtot2){ WAVwin=3; }
  if(Wpvtot[1]>Wpvtot[2] && Wpvtot[1]>Wpvtot[3]){ WPVwin=1; }
  if(Wpvtot[2]>Wpvtot[1] && Wpvtot[2]>Wpvtot[3]){ WPVwin=2; }
  if(Wpvtot[3]>Wpvtot[1] && Wpvtot[3]>Wpvtot[2]){ WPVwin=3; }
  if(Wirvtot1[1]<Wirvtot1[2] && Wirvtot1[1]<Wirvtot1[3]){ WIRV1loss=1; }
  if(Wirvtot1[2]<Wirvtot1[1] && Wirvtot1[2]<Wirvtot1[3]){ WIRV1loss=2; }
  if(Wirvtot1[2]==Wirvtot1[1] && Wirvtot1[2]<Wirvtot1[3]){ WIRV1loss=(rand()>0.5)?2:1; }
  if(Wirvtot1[3]<Wirvtot1[1] && Wirvtot1[3]<Wirvtot1[2]){ WIRV1loss=3; }
  if(WIRV1loss==1 && Wpairw[2][3]>Wpairw[3][2]+fuzz() ){ WIRVwin=2; }
  if(WIRV1loss==1 && Wpairw[2][3]<Wpairw[3][2]+fuzz() ){ WIRVwin=3; }
  if(WIRV1loss==2 && Wpairw[1][3]>Wpairw[3][1]+fuzz() ){ WIRVwin=1; }
  if(WIRV1loss==2 && Wpairw[1][3]<Wpairw[3][1]+fuzz() ){ WIRVwin=3; }
  if(WIRV1loss==3 && Wpairw[2][1]>Wpairw[1][2]+fuzz() ){ WIRVwin=2; }
  if(WIRV1loss==3 && Wpairw[2][1]<Wpairw[1][2]+fuzz() ){ WIRVwin=1; }
}

void PrintElectionResults(int N, int S){
  int j,k;
  pf("#voters=%d,  Random Subset Size=%d\n", N,S);
  pf("UNWEIGHTED (raw) RESULTS:\n");
  pf("=========================\n");
  pf("PLURALITY: \t _=%d\tC=%d\tL=%d\tM=%d\tX=%d\n", 
     pvtot[0], pvtot[1], pvtot[2], pvtot[3], pvtot[4] );
  pf("APPROVAL: C=%d\tL=%d\tM=%d\n", 
     avtot1,     avtot2,     avtot3 );
  pf("IRV ROUND ONE:  \t_=%d\tC=%d\tL=%d\tM=%d\tX=%d\n", 
     irvtot1[0], irvtot1[1], irvtot1[2], irvtot1[3], irvtot1[4] );
  pf("PAIRWISE:\n");
  for(j=1; j<=3; j++){ for(k=1; k<=3; k++){
      if( k!=j ){
	pf("(%c>%c)=%u,  ", who(j), who(k), pairw[j][k]);
      }
  }}
  pf("\n");
  pf("PVwin=%c\tAVwin=%c\t IRV1loss=%c\tIRVwin=%c\n",
     who(PVwin), who(AVwin), who(IRV1loss), who(IRVwin) );


  pf("WEIGHTED RESULTS:\n");
  pf("=========================\n");
  pf("WEIGHTS:  \t_=%.8f\tC=%.8f\tL=%.8f\tM=%.8f\tX=%.8f\n", 
     WT[0], WT[1], WT[2], WT[3], WT[4] );
  pf("PLURALITY:  \t_=%.3f\tC=%.3f\tL=%.3f\tM=%.3f\tX=%.3f\n", 
     Wpvtot[0], Wpvtot[1], Wpvtot[2], Wpvtot[3], Wpvtot[4] );
  pf("APPROVAL: \tC=%.3f\tL=%.3f\tM=%.3f\n", 
     Wavtot1,     Wavtot2,     Wavtot3 );
  pf("IRV ROUND ONE:  \t_=%.3f\tC=%.3f\tL=%.3f\tM=%.3f\tX=%.3f\n", 
     Wirvtot1[0], Wirvtot1[1], Wirvtot1[2], Wirvtot1[3], Wirvtot1[4] );
  pf("PAIRWISE:\n");
  for(j=1; j<=3; j++){ for(k=1; k<=3; k++){
      if( k!=j ){
	pf("(%c>%c)=%.3f,  ", who(j), who(k), Wpairw[j][k]);
      }
  }}
  pf("\n");
  pf("PVwin=%c\tAVwin=%c\t IRV1loss=%c\tIRVwin=%c\n",
     who(WPVwin), who(WAVwin), who(WIRV1loss), who(WIRVwin) );
}

void OneLineElRes(int N, int S){
  int j,k;
  pf("RAW");
  pf(" %d %d %d %d %d", 
     pvtot[0], pvtot[1], pvtot[2], pvtot[3], pvtot[4] );
  pf(" %d %d %d", 
     avtot1,     avtot2,     avtot3 );
  pf(" %d %d %d %d %d", 
     irvtot1[0], irvtot1[1], irvtot1[2], irvtot1[3], irvtot1[4] );
  for(j=1; j<=3; j++){ for(k=1; k<=3; k++){
      if( k!=j ){
	pf(" %u",  pairw[j][k]);
      }
  }}
  pf(" %c%c%c%c\n",
     who(PVwin), who(AVwin), who(IRV1loss), who(IRVwin) );


  pf("WGT %.8f %.8f %.8f %.8f %.8f", 
     WT[0], WT[1], WT[2], WT[3], WT[4] );
  pf(" %.3f %.3f %.3f %.3f %.3f", 
     Wpvtot[0], Wpvtot[1], Wpvtot[2], Wpvtot[3], Wpvtot[4] );
  pf(" %.3f %.3f %.3f", 
     Wavtot1,     Wavtot2,     Wavtot3 );
  pf(" %.3f %.3f %.3f %.3f %.3f", 
     Wirvtot1[0], Wirvtot1[1], Wirvtot1[2], Wirvtot1[3], Wirvtot1[4] );
  for(j=1; j<=3; j++){ for(k=1; k<=3; k++){
      if( k!=j ){
	pf(" %.3f", Wpairw[j][k]);
      }
  }}
  pf(" %c%c%c%c\n",
     who(WPVwin), who(WAVwin), who(WIRV1loss), who(WIRVwin) );
}

main(){
  int i;
  NumBallots = ReadBallots();  
  srand48(RandSeed); for(i=0; i<20; i++){ rand(); }
  //SubsetSz = 3*NumBallots/4;   WDS 13 Nov 2014: I now think 3/4 was a mistake.
  //Instead I now prefer 2/3 for reasons explained in  http://rangevoting.org/PuzzBootstrap.html .
  SubsetSz = (NumBallots*2)/3;
  DoElection(NumBallots, NumBallots);
  PrintElectionResults(NumBallots, NumBallots);
  pf("Here follow %d Random-Ballot-Subset elections, %d ballots (that's %.0f%% of them) in each subset:\n",
     NumBoots, SubsetSz, (SubsetSz*100.0)/NumBallots);
  for(i=0; i<NumBoots; i++){
    DoElection(NumBallots, SubsetSz);
    OneLineElRes(NumBallots, SubsetSz);
  }
}

/********
Some sample output using 14 Nov 2014 version of program
with 2/3 as subset fractional size and 100K bootstrap elections:

Raw Unweighted:
./a.out < MaineBallots.txt | fgrep RAW | fgrep MMLM | wc -l
   97519
./a.out < MaineBallots.txt | fgrep RAW | fgrep MMLC | wc -l
    2481
97519+2481=100000.

Weighted:
./a.out < MaineBallots.txt | fgrep WGT | fgrep LCMC | wc -l
   91294
./a.out < MaineBallots.txt | fgrep WGT | fgrep LCCM | wc -l
    6238
./a.out < MaineBallots.txt | fgrep WGT | fgrep LCCL | wc -l
    2160
./a.out < MaineBallots.txt | fgrep WGT | fgrep LMMC | wc -l
     150
./a.out < MaineBallots.txt | fgrep WGT | fgrep LMCM | wc -l
      93
./a.out < MaineBallots.txt | fgrep WGT | fgrep LMCL | wc -l
      59
./a.out < MaineBallots.txt | fgrep WGT | fgrep LCLC | wc -l
       3
./a.out < MaineBallots.txt | fgrep WGT | fgrep LLMC | wc -l
       2
./a.out < MaineBallots.txt | fgrep WGT | fgrep LLCL | wc -l
       1
91294+6238+2160+150+93+59+3+2+1=100000.
***********************/
