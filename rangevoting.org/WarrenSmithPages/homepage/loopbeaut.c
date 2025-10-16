/******************
 * Compile with   gcc loopbeaut.c -O -o loopbeaut
 *
 * To do:  add more identities, add parser,  
 * clean up and modularize yukky code, add modes for other things 
 * besides loop, like Hadamard matrices...
 ********/

#include <stdio.h>
#include <ctype.h>
#define until(x) while(!(x))
#define VERSION 1.0
#define MAXSIZE 115

/*
     int isalpha(int c);
     int isupper(int c);
     int islower(int c);
     int isdigit(int c);
     int isxdigit(int c);
     int isalnum(int c);
     int isspace(int c);
     int ispunct(int c);
     int isprint(int c);
     int isgraph(int c);
     int iscntrl(int c);
     int isascii(int c);
*/

char bok[] = "ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuwxyz!@#$^&|:;<,~_-.`%?????";
int Lident=-1, Rident=-2;
int Linv[MAXSIZE], Rinv[MAXSIZE];
int pow[MAXSIZE];
int mark[MAXSIZE];
int L[MAXSIZE][MAXSIZE];
int Z[MAXSIZE][MAXSIZE];
char commhist[4999];
int OSPACE = 0;
int VERBOS = 1;
int record = 0;
int recguy = 0;
int rec2ord = 0;
int recguy1 = 0;
int recguy2 = 0;
int N;

char toalnum(int x){
  if(x<0) return '?';
  if(x<=9) return '0'+x;
  x -= 10;
  if(x>67) return '?';
  return bok[x];
}

xxx(int x, int y){
    int i,j,r2,r1,z;
    printf("x%c%c ", x,y);
      /* printf("relabeling chars %c %c\n", x,y );*/
      for(i=0; i<=N; i++){
	for(j=0; j<=N; j++){
          if(L[i][j]==x) L[i][j]=y;
          else if(L[i][j]==y) L[i][j]=x;
      }}


      r2=0; r1=0;
      for(i=1; i<=N; i++){ 
         if(L[i][0]==x) r1=i;
         if(L[i][0]==y) r2=i;
	 /*printf("%c", L[0][i]);*/
      }
      /*printf("swapping rows %d %d\n", r1,r2 );*/
      if(r1>0 && r2>0) for(j=0; j<=N; j++){
	z = L[r1][j]; L[r1][j] = L[r2][j]; L[r2][j] = z;
      }
      r2=0; r1=0;
      for(i=1; i<=N; i++){ 
         if(L[0][i]==x) r1=i;
         if(L[0][i]==y) r2=i;
	 /*printf("%c", L[i][0]);*/
      }
      /*printf("swapping cols %d %d\n", r1,r2 );*/
      if(r1>0 && r2>0) for(j=0; j<=N; j++){
	z = L[j][r1]; L[j][r1] = L[j][r2]; L[j][r2] = z;
      }
}

main(){
  char inptype, ouptype;
  int i,j,k,m,c,x,y,z,r1,r2,si;
  int hh=0, pct=0;

  printf(" Loopbeaut.c version %.2f by  Warren D. Smith Feb 2004 \n",VERSION);
  printf(" Interactive loop beautification tool\n");
  printf(
"(also applicable to quasigroups, groups, hadamard matrices, graphs, etc)\n");
  printf("\n");
  printf(" Distributed under the GPL (Gnu public license).  \n");
  printf("Please credit Warren D. Smith if use this code & please\n");
  printf("give me any worthwhile code modifs you make WDSmith@fastmail.fm\n");

  for(j=0; j<67; j++) printf("%c", bok[j]);
  printf("\n\n");

  printf("How many elements?\n");
  scanf("%d", &N);
  printf("  You chose N=%d.\n", N);
  if(N>=MAXSIZE) printf("too large\n");
  if(N<0) printf("negative\n");
  printf("Would you like to input the loop using\n");
  printf("  a: 1-digit-wide decimal no spaces\n");
  printf("  b: 2-digit-wide decimal no spaces (mace tabular output format)\n");
  printf("  c: space separated decimal \n");
  printf("  d: alnum \n");
  do{ c=getchar(); }until(c>='a' && c<='d');
  inptype = c;
  printf("  You chose %c.\n", inptype);
  printf(" Input the loop mul table now (1st nonwhite char=* please):\n");
  for(i=0; i<=N; i++){
    if(i==0){
      do{ c=getchar(); }until(c=='*'); L[0][0]='*'; 
    }
    else{
      if(inptype=='a'){
	do{ c=getchar(); }until(isdigit(c)); x = c-'0';
	x = toalnum(x);
      }
      if(inptype=='b'){
	do{ c=getchar(); }until(isdigit(c)); x = c-'0';
	c=getchar(); if(c!=' '){ y=c-'0'; x = 10*x+y; }
	x = toalnum(x);
      }
      if(inptype=='c'){     
	do{ c=getchar(); }until(isdigit(c)); x = c-'0';
	for(;;){ 
	  c=getchar(); if(isdigit(c)){ y=c-'0'; x = 10*x+y; } else break; 
	}
	x = toalnum(x);
     }
     if(inptype=='d'){     
       do{ c=getchar(); }until(isalnum(c)); x = c;
     }
     L[i][0] = x;
    }
    do{ c=getchar(); }until(c=='|');
    for(j=1; j<=N; j++){
      if(inptype=='a'){
	c=getchar(); x = c-'0';
	x = toalnum(x);
      }
      if(inptype=='b'){
	c=getchar(); x=0; if(c!=' ') x = c-'0';
	c=getchar(); x = 10*x + c-'0';
	x = toalnum(x);
      }
      if(inptype=='c'){     
	do{ c=getchar(); }until(isdigit(c)); x = c-'0';
	for(;;){ 
	  c=getchar(); if(isdigit(c)){ y=c-'0'; x = 10*x+y; } else break; 
	}
	x = toalnum(x);
      }
      if(inptype=='d'){     
	do{ c=getchar(); }until(c!=' '); x = c;
      }
      L[i][j] = x;       
    }
    do{ c=getchar(); }until(c=='\n'); 
    printf("**read in line %d:",i);
    for(j=0; j<=N; j++) printf("%c ", L[i][j]);
    printf("**\n");
    if(i==0){ do{ c=getchar(); }until(c=='\n'); }
  }
  L[0][0] = '*';

  /* print the loop: */
  for(;;){  
    for(i=0; i<=N; i++){
      if(OSPACE) printf("\n%c |", L[i][0]);
      else       printf("\n%c|", L[i][0]);
      for(j=1; j<=N; j++){ 
	  if(OSPACE) printf(" %c", L[i][j]); 
	  else       printf("%c", L[i][j]); 
      }
      if(i==0){
	  if(OSPACE){
	      printf("\n--+");
	      for(j=1; j<=N; j++){ printf("--"); }
	  }else{
	      printf("\n-+");
	      for(j=1; j<=N; j++){ printf("-"); }
	  }
      }
    } 

    printf("\n");
    /* compute some loop info: */
    for(i=1; i<=N; i++){
	for(j=1; j<=N; j++){
	    Z[i][j] = 0;
	    for(k=1; k<=N; k++){
	      if( L[i][j]==L[0][k] ){
		  Z[i][j] = k; 
                  break;
	      }
	    }
	    if(Z[i][j]<=0 || Z[i][j]>N) 
		printf("bogus loop %c.%c=%c\n",
				  L[0][i], L[0][j], L[i][j] );
	}
    }
    
    /* read command and do it: */  
    printf("\n\nwhat now (? for help):  ");
    do{ c=getchar(); }until(
   	   c=='s' || c=='x' || c=='r' || c=='?' 
	|| c=='L' || c=='v' || c=='h' || c=='e' 
	|| c=='i' || c=='q' || c=='O' || c=='V'
        || c=='a'
    );
    commhist[hh++] = c;
    if(c=='r' || c=='s' || c=='x'){ /* 2 char args */
       x=getchar(); y=getchar();
       printf("  You chose %c%c%c.\n", c,x,y);
       commhist[hh++] = x;
       commhist[hh++] = y;
    }
    if(c=='h' || c=='v'){ /* 1 char args */
       x=getchar(); 
       printf("  You chose %c%c.\n", c,x);
       commhist[hh++] = x;
    }
    if(c=='e' || c=='i'){ /* balanced paren string as arg */
	do{ x=getchar(); }until(x=='(');
	pct=1; y=0;
        bok[y++]=x;
	commhist[hh++] = x;
        do{  
	    x=getchar(); 
	    bok[y++]=x; 
	    commhist[hh++] = x;
	    if(x=='(') pct++; 
	    else if(x==')') pct--;
	}until(pct==0);
	bok[y]=0;
	printf("  You chose %c%s.\n", c,bok);
    }
    commhist[hh++] = '\n';

    if(c=='?'){
	printf("?:  help (prints this command summary)\n");
	printf("c:  command history so far - recounts it.\n");
	printf("a:  auto-beautify.\n");
	printf(
"rAB: replaces symbol A with symbol B and vice versa everywhere.\n");
	printf(
"sAB: swaps row with header A with row with header B, ditto for columns.\n");
	printf("xAB: same as doing both sAB and rAB.\n");
	printf("q:  prints facts about current table\n");
	printf("O:  toggles space-consuming output mode\n");
	printf("V:  toggles verbosity\n");
	printf("e(1.2(76)):  evaluates expression.\n");
	printf(
"i(a.b(cd)=(ab)c.d):  checks identity, finds violations if any. Lowercase vars.\n");
	printf(
"vA: if vertical col headed by A contains only 2 symbols, it swaps them.\n");
	printf(
	    "hA: if horizontal row headed by A contains only 2 symbols, it swaps them.\n");
	printf("L: latex loop table output\n");
    }

    if(c=='r' || c=='x'){
      printf("relabeling chars %c %c\n", x,y );
      for(i=0; i<=N; i++){
	for(j=0; j<=N; j++){
          if(L[i][j]==x) L[i][j]=y;
          else if(L[i][j]==y) L[i][j]=x;
      }}
    }
    if(c=='s' || c=='x'){
      r2=0; r1=0;
      for(i=1; i<=N; i++){ 
         if(L[i][0]==x) r1=i;
         if(L[i][0]==y) r2=i;
	 printf("%c", L[0][i]);
      }
      printf("swapping rows %d %d\n", r1,r2 );
      if(r1>0 && r2>0) for(j=0; j<=N; j++){
	z = L[r1][j]; L[r1][j] = L[r2][j]; L[r2][j] = z;
      }
      r2=0; r1=0;
      for(i=1; i<=N; i++){ 
         if(L[0][i]==x) r1=i;
         if(L[0][i]==y) r2=i;
	 printf("%c", L[i][0]);
      }
      printf("swapping cols %d %d\n", r1,r2 );
      if(r1>0 && r2>0) for(j=0; j<=N; j++){
	z = L[j][r1]; L[j][r1] = L[j][r2]; L[j][r2] = z;
      }
    }
    if(c=='h'){
	j=0;
	for(i=1; i<=N; i++){ if(L[i][0]==x) j=i; }
        if(j){
	    r1 = L[j][1]; r2=0;
	    for(i=2; i<=N; i++){ 
		if(L[j][i]!=r1){
		    if(r2==0) r2=L[j][i]; 
		    if(L[j][i]!=r2) goto LOSEh;
		}
	    }
	    for(i=1; i<=N; i++){ 
		if(L[j][i]==r1) L[j][i]==r2;
		else    	L[j][i]==r1;
	    }
	}
    LOSEh: ;
    }
    if(c=='v'){
	j=0;
	for(i=1; i<=N; i++){ if(L[0][i]==x) j=i; }
        if(j){
	    r1 = L[1][j];
	    for(i=2; i<=N; i++){ 
		if(L[i][j]!=r1){
		    if(r2==0) r2=L[i][j]; 
		    if(L[i][j]!=r2) goto LOSEv;
		}
	    }
	    for(i=1; i<=N; i++){ 
		if(L[i][j]==r1) L[i][j]==r2;
		else    	L[i][j]==r1;
	    }
	}
    LOSEv: ;
    }
    if(c=='c'){
	for(i=0; i<=hh; i++){
	    putchar(commhist[i]);
	}
    }
    if(c=='q'){
        printf("\nNumber of elements=%d\n",N);
        printf("Same row headers as column headers in same order? ");
        for(i=1; i<=N; i++){
	    if(L[0][i] != L[i][0]){
		printf("NO! %c != %c\n", L[0][i], L[i][0]);
		goto BADHEADS;
	    }
	}
        printf("YES!\n");
    BADHEADS: ;
        printf("Is it a quasigroup (latin square)? ");
        for(i=0; i<=N; i++){
	    for(j=2; j<=N; j++){
		for(k=1; k<j; k++){
		    if(L[i][j]==L[i][k]){
			printf("NO! row %c: entries %c %c agree as %c\n",
			       L[i][0], L[0][k], L[0][j], L[i][j]);
			goto LOSEquasi;
		    }
		}
	    }
	}
        for(i=0; i<=N; i++){
	    for(j=2; j<=N; j++){
		for(k=1; k<j; k++){
		    if(L[j][i]==L[k][i]){
			printf("NO! col %c: entries %c %c agree as %c\n",
			       L[0][i], L[k][0], L[j][0], L[j][i]);
			goto LOSEquasi;
		    }
		}
	    }
	}
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		x = L[i][j];
		for(k=1; k<=N; k++){
		    if(L[0][k]==x) goto FOUNDitq;
		}
		printf("NO! entry %c in row %c col %c is strange\n", 
		       x, L[i][0], L[0][j]);
		goto LOSEquasi;
	    FOUNDitq: ;
	    }
	}
	printf("YES!\n");
    LOSEquasi: ;
        printf("Left-identity(ies):  {");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if(L[i][j]!=L[0][j]) goto NOTLident;
	    }
	    printf("%c ", L[i][0]);
	    Lident = i;
	NOTLident: ;
	}
	printf("}\n");
        printf("Right-identity(ies):  {");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if(L[j][i]!=L[j][0]) goto NOTRident;
	    }
	    printf("%c ", L[0][i]);
	    Rident = i;
	NOTRident: ;
	}
	printf("}\n");
	printf("Left and right idents same? ");
	if(Lident==Rident) printf("YES!\n");
	else printf("NO!\n");
	printf("Right-inverses:\n");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if(Z[i][j]==Lident){
                   printf("%c%c=%c ", L[i][0], L[j][0], L[Lident][0] );
		   Rinv[i] = j;
		   break;
                }
	    }
        }
        printf("\n");
	printf("Left-inverses:\n");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if(Z[j][i]==Rident){
                   printf("%c%c=%c ", L[j][0], L[i][0], L[Rident][0] );
		   Linv[i] = j;
		   break;
                }
	    }
        }
        printf("\n");
	printf("commutative? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if(L[j][i]!=L[i][j]){
		    printf("NO! %c.%c = %c != %c = %c.%c\n",
			   L[i][0], L[0][j], L[i][j],
			   L[j][i], L[j][0], L[0][i] );
		    goto NOTcommut;
		}
	    }
	}
	printf("YES!\n");
	NOTcommut: ;
	printf("associative? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		for(k=1; k<=N; k++){
		    if( Z[i][Z[j][k]] != Z[Z[i][j]][k] ){
			printf("NO! %c.%c%c = %c != %c = %c%c.%c\n",
			   L[i][0], L[j][0], L[k][0],
                           L[ Z[i][Z[j][k]] ][0],
			   L[ Z[Z[i][j]][k] ][0],
			   L[i][0], L[j][0], L[k][0] );
			goto NOTassoc;
		    }
		}
	    }
	}
	printf("YES!\n");
	NOTassoc: ;
	printf("Lalternative? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if( Z[i][Z[i][j]] != Z[Z[i][i]][j] ){
		    printf("NO! %c.%c%c = %c != %c = %c%c.%c\n",
			   L[i][0], L[i][0], L[j][0],
                           L[ Z[i][Z[i][j]] ][0],
			   L[ Z[Z[i][i]][j] ][0],
			   L[i][0], L[i][0], L[j][0] );
		    goto NOTLalt;
		}
	    }
	}
	printf("YES!\n");
	NOTLalt: ;
	printf("Ralternative? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if( Z[i][Z[j][j]] != Z[Z[i][j]][j] ){
		    printf("NO! %c.%c%c = %c != %c = %c%c.%c\n",
			   L[i][0], L[j][0], L[j][0],
                           L[ Z[i][Z[j][j]] ][0],
			   L[ Z[Z[i][j]][j] ][0],
			   L[i][0], L[j][0], L[j][0] );
		    goto NOTRalt;
		}
	    }
	}
	printf("YES!\n");
	NOTRalt: ;
	printf("Flexible? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if( Z[i][Z[j][i]] != Z[Z[i][j]][i] ){
		    printf("NO! %c.%c%c = %c != %c = %c%c.%c\n",
			   L[i][0], L[j][0], L[i][0],
                           L[ Z[i][Z[j][i]] ][0],
			   L[ Z[Z[i][j]][i] ][0],
			   L[i][0], L[j][0], L[i][0] );
			goto NOTflex;
		}
	    }
	}
	printf("YES!\n");
	NOTflex: ;
	printf("L-Bol? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		for(k=1; k<=N; k++){
		    if( Z[Z[i][Z[j][i]]][k] != Z[i][Z[Z[j][i]][k]] ){
			printf("NO! (%c.%c%c)%c = %c != %c = %c(%c.%c%c)\n",
			   L[i][0], L[j][0], L[i][0], L[k][0],
                           L[0][ Z[Z[i][Z[j][i]]][k] ],
			   L[0][ Z[i][Z[Z[j][i]][k]] ],
                           L[i][0], L[j][0], L[i][0], L[k][0] );
			goto NOTLBol;
		    }
		}
	    }
	}
	printf("YES!\n");
	NOTLBol: ;
	printf("R-Bol? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		for(k=1; k<=N; k++){
		    if( Z[k][Z[Z[i][j]][i]] != Z[Z[Z[k][i]][j]][i] ){
			printf("NO! %c(%c%c.%c) = %c != %c = (%c%c.%c)%c\n",
			   L[k][0], L[i][0], L[j][0], L[i][0],
			       L[0][ Z[k][Z[Z[i][j]][i]] ],
			       L[0][ Z[Z[Z[k][i]][j]][i] ],
                           L[k][0], L[i][0], L[j][0], L[i][0] );
			goto NOTRBol;
		    }
		}
	    }
	}
	printf("YES!\n");
	NOTRBol: ;
	printf("MMoufang? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		for(k=1; k<=N; k++){
		    if( Z[Z[i][j]][Z[k][i]] !=  Z[Z[i][Z[j][k]]][i] ){
			printf("NO! %c%c.%c%c = %c != %c = (%c.%c%c)%c\n",
			   L[i][0], L[j][0], L[k][0], L[i][0],
			       L[0][ Z[Z[i][j]][Z[k][i]]  ],
			       L[0][ Z[Z[i][Z[j][k]]][i]  ],
                           L[i][0], L[j][0], L[k][0], L[i][0] );
			goto NOTMmouf;
		    }
		}
	    }
	}
	printf("YES!\n");
	NOTMmouf: ;
	printf("3-power-assoc? ");
        for(i=1; i<=N; i++){
	    if( Z[Z[i][i]][i] !=  Z[i][Z[i][i]] ){
		printf("NO! %c%c.%c = %c != %c = %c.%c%c\n",
		       L[i][0], L[i][0], L[i][0],
		       L[0][ Z[Z[i][i]][i] ],
		       L[0][  Z[i][Z[i][i]] ],
		       L[i][0], L[i][0], L[i][0] );
                printf("xxx=1? NO!\n");
		goto NOT3PA;
	    }
	}
	printf("YES!\n");
	printf("xxx=1? ");
        for(i=1; i<=N; i++){
	    if( Z[Z[i][i]][i] !=  Lident ){
		printf("NO! %c%c.%c = %c != %c\n",
		       L[i][0], L[i][0], L[i][0],
		       L[0][ Z[Z[i][i]][i] ], L[Lident][0] );
		goto NOTxxx1;
	    }
	}
	printf("YES!\n");
	NOTxxx1: ;
	NOT3PA: ;
	printf("xx=1? ");
        for(i=1; i<=N; i++){
	    if( Z[i][i] != Lident ){
		printf("NO! %c%c = %c != %c\n",
		       L[i][0], L[i][0],
		       L[0][ Z[i][i] ],
		       L[Lident][0] );
		goto NOTxx1;
	    }
	}
	printf("YES!\n");
	NOTxx1: ;
	printf("2-sided-inverses? ");
        for(i=1; i<=N; i++){
	    if( Linv[i] != Rinv[i] ){
		printf("NO! %c^ = %c != %c = ^%c\n",
		       L[i][0], L[ Linv[i] ][0], 
		       L[ Rinv[i] ][0], L[i][0] );
		goto NOT2si;
	    }
	}
	printf("YES!\n");
	NOT2si: ;
	printf("Antiautomorphic inverse? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if( Z[Linv[i]][Linv[j]] !=  Linv[Z[j][i]] ){
		    printf("NO! %c^.%c^ = %c != %c = (%c.%c)^\n",
			   L[i][0], L[j][0], 
			   L[0][ Z[Linv[i]][Linv[j]] ],
			   L[0][ Linv[Z[j][i]] ],
			   L[j][0], L[i][0] );
		    goto NOTaa;
		}
	    }
	}
	printf("YES!\n");
	NOTaa: ;
	printf("Lcancellative (LIP)? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if( Z[Linv[i]][Z[i][j]] != j ){
		    printf("NO! %c^.%c%c = %c != %c\n",
			   L[i][0], L[i][0], L[j][0],
                           L[  Z[Linv[i]][Z[i][j]] ][0],
			   L[j][0] );
		    goto NOTLcan;
		}
	    }
	}
	printf("YES!\n");
	NOTLcan: ;
	printf("Rcancellative (RIP)? ");
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if( Z[Z[i][j]][Rinv[j]] != i ){
		    printf("NO! %c%c.%c^ = %c != %c\n",
			   L[i][0], L[j][0], L[j][0],
                           L[   Z[Z[i][j]][Rinv[j]]  ][0],
			   L[i][0] );
		    goto NOTRcan;
		}
	    }
	}
	printf("YES!\n");
	NOTRcan: ;
	printf("Power-associative?... ");
	if(VERBOS) printf("\n");
	x=0; record = 0; recguy = 0;
        for(i=1; i<=N; i++){
	    pow[0] = Lident;
	    for(j=1; j<=N; j++){
		pow[j] = Z[i][pow[j-1]];
		if( pow[j] != pow[0] ){
		    if(VERBOS) printf("%c", L[ pow[j] ][0] );
		}else{
		    break;
		}
	    }
	    if(VERBOS) printf("%c order=%d\n", L[Lident][0], j );
	    if(j>record){ record = j; recguy = i; }
	    if(x==0){
		for(k=1; k<=j; k++){
		    for(m=1; m<k; m++){
			if(Z[pow[m]][pow[k-m]] != pow[k]){
			    x=m; y=k; z=pow[m]; r1=pow[k]; r2=pow[k-m]; si=i;
			}
		    }
		}
	    }
	}
	printf("(Biggest left-mul order %c^%d=%c) ", 
	       L[recguy][0], record, L[Lident][0] );
	if(x!=0){
	    si = L[si][0];
	    printf("...NO! %c^%d . %c^%d = %c.%c != %c = %c^%d\n",
		   si, x, si, y-x, L[z][0], L[r2][0], 
		   L[ Z[z][r2] ][0], L[r1][0], si, y );
	}else{ printf("...YES!\n"); }

        printf("Diassociative?... ");
	rec2ord=0;
	if(VERBOS) printf("\n");
	si = 1;
        for(i=2; i<=N; i++){
	    for(j=1; j<i; j++){
		for(k=1; k<=N; k++){	mark[k] = 0;    }
                mark[Lident] = 1;
		if(VERBOS) printf("<%c%c> = ", L[i][0], L[j][0] );
                for(k=1; k<=N; k++){   
		    for(x=1; x<=N; x++){   
			if(mark[x]){
			    mark[ Z[x][i] ]=1;
			    mark[ Z[x][j] ]=1;
			    mark[ Z[i][x] ]=1;
			    mark[ Z[j][x] ]=1;
			}
		    }
		}
		x=0;
                for(k=1; k<=N; k++){   
		    if(mark[k]){
			if(VERBOS) printf("%c", L[0][k] );
			x++;
		    }
		}
                if(VERBOS) printf("  order=%d  ", x);
		if(x>rec2ord){ rec2ord = x; recguy1 = i; recguy2 = j; }
                for(k=1; k<=N; k++) if(mark[k]){   
                for(x=1; x<=N; x++) if(mark[x]){   
                for(y=1; y<=N; y++) if(mark[y]){   
		    if( Z[k][Z[x][y]] != Z[Z[k][x]][y] ){
			if(VERBOS) printf("NOTdiassoc\n");
			si = 0;
			goto NOTdias;
		    }
		}}}
		if(VERBOS) printf("diassoc\n");
	    NOTdias: ;
	    }
	}
	printf("(Biggest <%c%c> order =%d) ",
	       L[recguy1][0], L[recguy2][0],  rec2ord);
	if(si==0) printf("...NO!\n");
	else printf("...YES!\n");
    }
    if(c=='L'){
	printf("\\begin{center}\n");
	printf("\\begin{tabular}{l|l}");
	for(i=0; i<=N; i++){
	    if(i != 0) printf("\\\\");
	    if(i==1){
		printf("\n\\hline");
	    }
	    printf("\n{\\tt %c} & {\\tt", L[i][0]);
	    for(j=1; j<=N; j++){ printf(" %c", L[i][j]); }
	    printf("} ");
	} 
	printf("\\\\\n\\hline\n\\end{tabular}\n");
	printf("\\end{center}\n");
	printf("\\Bcaption\n");
	printf("%d-element loop. \ncode\\# \n",N);
	printf("\\Ecaption\n");
    }
    if(c=='O'){
	OSPACE = 1-OSPACE;
	printf("toggling space-consuming output mode to %d\n", OSPACE);
    }
    if(c=='V'){
	VERBOS = 1-VERBOS;
	printf("toggling verbosity to %d\n", VERBOS);
    }
    if(c=='i'){
	printf("i not implemented\n");
    }
    if(c=='e'){
	printf("e not implemented\n");
    }
    if(c=='a'){
        printf("autobeaut  (pretty crude right now but a start)\n");
        Lident = 1;
        for(i=1; i<=N; i++){
	    for(j=1; j<=N; j++){
		if(L[i][j]!=L[0][j]) goto NOTLident2;
	    }
	    /*printf("%c ", L[i][0]);*/
	    Lident = i;
	NOTLident2: ;
	}
        if(Lident != 1) xxx( L[Lident][0], L[1][0] );

	    x=0; record = 0; recguy = 0;
	    for(i=1; i<=N; i++){
		pow[0] = Lident;
		for(j=1; j<=N; j++){
		    pow[j] = Z[i][pow[j-1]];
		    if( pow[j] == pow[0] ){
			break;
		}
		}
		if(j>record){ record = j; recguy = i; }
		if(x==0){
		    for(k=1; k<=j; k++){
			for(m=1; m<k; m++){
			    if(Z[pow[m]][pow[k-m]] != pow[k]){
				x=m; y=k; z=pow[m]; r1=pow[k];
				r2=pow[k-m]; si=i;
			    }
			}
		    }
		}
	    }
	    xxx( L[recguy][0], L[2][0] );
	    for(i=2; i<record; i++){
		xxx( L[i+1][0], L[2][i] );
	    }

        printf("\n");
    }

  }
}


