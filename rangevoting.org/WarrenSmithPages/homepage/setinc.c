/* program by WDS feb 2004 to keep track of loop-type
 * inclusions. compile with
 *      gcc setinc.c -Wall -o setinc
 * and run... 
 ******************************/

#define GRP  (1<<0)
#define MOUF (1<<1)
#define LBOL (1<<2)
#define RBOL (1<<3)
#define DIA  (1<<4)
#define PA   (1<<5)
#define IPALT  (1<<6)
#define IP   (1<<7)
#define ALT  (1<<8)
#define IPLR (1<<9)
#define LRA  (1<<10)
#define LIP  (1<<11)
#define RIP  (1<<12)
#define LALT (1<<13)
#define RALT (1<<14)
#define FLEX (1<<15)
#define AA   (1<<16)
#define PA3  (1<<17)
#define SI2  (1<<18)

/* change these modes to 1 and 0 ... */
                        /* 2^19 = 524288 */
#define INCLUDES  1     /*324*/
#define AAMIRROR  1     /*202*/
#define SOPHIST   1     /*79*/
#define FINITEINC 1     /*78*/
#define QUESTIONABLE 0  /*64*/

#define VERBOSE   1
#define NICEPRINT   1

#define MACING 0
#define REALLYMACE 0   /* these two will invoke mace4 to seek 
                        * models in all cases. Best to run in
			* empty directory since creates a lot of files.*/

/**************************************************/

#define MACEBOT 3
#define MACETOP 100

#define MACEMAXSEC 3000  /*will be 43*this seconds total runtime...
                         best to try with 1 then 5 then 25 then 125... */


#include <stdio.h>
FILE *fopen(const char *path, const char *mode);
int fclose(FILE *stream);
#include <stdlib.h>
int system(const char *string);

void spitmace4infile(int code, uint nflag, uint pflag){
    char fname[100];
    FILE *pamf;
    sprintf(fname, "si%x.in", code);
    printf("opening file %s for w\n", fname);
    pamf =    fopen(fname, "w");
    fprintf(pamf, "assign(domain_size, %d).\n", MACEBOT);
    fprintf(pamf, "assign(iterate_up_to,%d).\n", MACETOP);
    fprintf(pamf, "assign(max_seconds,%d).\n", MACEMAXSEC);
    fprintf(pamf, "assign(max_models,1).\n");
    fprintf(pamf, "op(400, infix, [+,\\,/]).\n");
    /*    fprintf(pamf, "op(300, postfix, [^]).\n"); */
    /*%set(verbose).*/
    fprintf(pamf, "clauses(loops).\n");
    fprintf(pamf, "x \\ (x + y) = y.\n"); /* defines \ */
    fprintf(pamf, "(x + y) / y = x.\n");  /* defines /  */
    fprintf(pamf, "0 + x = x.\n");  /* left ident */
    fprintf(pamf, "x + 0 = x.\n");  /* right ident */

    if(code&GRP  ){
	fprintf(pamf, "x + (y + z) = (x + y) + z. %%GRP\n");
    }else if(!(nflag&GRP)){
	fprintf(pamf, "A + (B + C) != (A + B) + C. %%!GRP\n");
    }
    if(code&LBOL ){
	fprintf(pamf, "x + (y + (x + z)) = (x + (y + x)) + z. %%LBOL\n");
    }else  if(!(nflag&LBOL)){
	fprintf(pamf, "D + (E + (D + F)) != (D + (E + D)) + F. %%!LBOL\n");
    }
    if(code&RBOL ){
	fprintf(pamf, "x + ((y + z) + y) = ((x + y) + z) + y. %%RBOL\n");
    }else  if(!(nflag&RBOL)){
	fprintf(pamf, "G + ((H + I) + H) != ((G + H) + I) + H. %%!RBOL\n");
    }
    if(code&LIP  ){
	fprintf(pamf, "(0 / x) + (x + y) = y. %%LIP\n"); 
    }else  if(!(nflag&LIP)){
	fprintf(pamf, "(0 / J) + (J + K) != K. %%!LIP\n"); 
    }
    if(code&RIP  ){
	fprintf(pamf, "(x + y) + (y \\ 0) = x. %%RIP\n");
    }else  if(!(nflag&RIP)){
	fprintf(pamf, "(L + M) + (M \\ 0) != L. %%!RIP\n");
    }
    if(code&LALT ){
	fprintf(pamf, "x + (x + y) = (x + x) + y. %%LALT\n"); 
    }else  if(!(nflag&LALT)){
	fprintf(pamf, "P + (P + Q) != (P + P) + Q. %%!LALT\n");
    }
    if(code&RALT ){
	fprintf(pamf, "(y + x) + x = y + (x + x). %%RALT\n"); 
    }else   if(!(nflag&RALT)){
	fprintf(pamf, "(R + S) + S != R + (S + S). %%!RALT\n"); 
    }
    if(code&FLEX ){
	fprintf(pamf, "x + (y + x) = (x + y) + x.  %%FLEX\n"); 
    }else   if(!(nflag&FLEX)){
	fprintf(pamf, "T + (U + T) != (T + U) + T. %%!FLEX\n");
    }
    if(code&AA   ){
	fprintf(pamf, "0 / (x + y) = (0 / y) + (0 / x). %%AA\n"); 
    }else   if(!(nflag&AA)){
	fprintf(pamf, "0 / (V + W) != (0 / W) + (0 / V). %%!AA\n"); 
    }
    if(code&PA3  ){
	fprintf(pamf, "x + (x + x) = (x + x) + x. %%PA3\n"); 
    }else    if(!(nflag&PA3)){
	fprintf(pamf, "X + (X + X) != (X + X) + X. %%!PA3\n"); 
    }
    if(code&SI2  ){
	fprintf(pamf, "(0 / x) = (x \\ 0).  %%SI2\n"); 
    }else    if(!(nflag&SI2)){
	fprintf(pamf, "(0 / Y) != (Y \\ 0). %%!SI2\n"); 
    }
    fprintf(pamf, "end_of_list.\n");
    fclose(pamf);
}

uint improveneg(int k){  /* flags unnecessary-to-say negative statements */
    uint z;
    z=k;
    if( !(k&GRP) && !(k&MOUF) )   z |=  GRP;
    if( !(k&MOUF) && !(k&LBOL) )  z |=  MOUF;
    if( !(k&MOUF) && !(k&RBOL) )  z |=  MOUF;
    if( !(k&MOUF) && !(k&DIA) )   z |=  MOUF;
    if( !(k&DIA) && !(k&PA) )     z |=  DIA;
    if( !(k&DIA) && !(k&IPALT) )  z |=  DIA;
    if( !(k&LBOL) && !(k&PA) )    z |=  LBOL;
    if( !(k&LBOL) && !(k&LALT) )  z |=  LBOL;
    if( !(k&LBOL) && !(k&LIP) )   z |=  LBOL;
    if( !(k&RBOL) && !(k&PA) )    z |=  RBOL;
    if( !(k&RBOL) && !(k&RALT) )  z |=  RBOL;
    if( !(k&RBOL) && !(k&RIP) )   z |=  RBOL;
    if( !(k&IPALT) && !(k&IP) )   z |=  IPALT;
    if( !(k&IPALT) && !(k&ALT) )  z |=  IPALT;
    if( !(k&IPALT) && !(k&IPLR) ) z |=  IPALT;
    if( !(k&IPLR) && !(k&LRA) )   z |=  IPLR;
    if( !(k&IPLR) && !(k&IP) )    z |=  IPLR;
    if( !(k&LRA) && !(k&LALT) )   z |=  LRA;
    if( !(k&ALT) && !(k&LRA) )    z |=  ALT;
    if( !(k&LRA) && !(k&RALT) )   z |=  LRA;
    if( !(k&PA) && !(k&PA3) )     z |=  PA;
    if( !(k&PA) && !(k&SI2) )     z |=  PA;
    if( !(k&IP) && !(k&LIP) )     z |=  IP;
    if( !(k&IP) && !(k&RIP) )     z |=  IP;
    if( !(k&IP) && !(k&AA) )      z |=  IP;
    if( !(k&AA) && !(k&SI2) )     z |=  AA;
    if( !(k&ALT) && !(k&RALT) )   z |=  ALT;
    if( !(k&ALT) && !(k&LALT) )   z |=  ALT;
    if( !(k&ALT) && !(k&FLEX) )   z |=  ALT;
    if( !(k&LALT) && !(k&PA3) )   z |=  LALT;
    if( !(k&RALT) && !(k&PA3) )   z |=  RALT;
    if( !(k&FLEX) && !(k&PA3) )   z |=  FLEX;
    if( !(k&FLEX) && !(k&SI2) )   z |=  FLEX;
    if( !(k&RIP) && !(k&SI2) )    z |=  RIP;
    if( !(k&LIP) && !(k&SI2) )    z |=  LIP;
    return z;
}

uint improvepos(int k){  /* flags unnecessary-to-say positive statements */
    uint z;
    z=0;
    if( (k&GRP) )  z |= MOUF;
    if( (k&MOUF) ) z |= LBOL;
    if( (k&MOUF) ) z |= RBOL;
    if( (k&MOUF) ) z |= DIA;
    if( (k&DIA) )  z |= PA;
    if( (k&DIA) )  z |= IPALT;
    if( (k&LBOL) ) z |= PA;
    if( (k&LBOL) ) z |= LALT;
    if( (k&LBOL) ) z |= LIP;
    if( (k&RBOL) ) z |= PA;
    if( (k&RBOL) ) z |= RALT;
    if( (k&RBOL) ) z |= RIP;
    if( (k&IPALT) ) z |= IP;
    if( (k&IPALT) ) z |= ALT;
    if( (k&IPALT) ) z |= IPLR;
    if( (k&IPLR) ) z |= LRA;
    if( (k&IPLR) ) z |= IP;
    if( (k&ALT) )  z |= LRA;
    if( (k&LRA) )  z |= LALT;
    if( (k&LRA) )  z |= RALT;
    if( (k&PA) )   z |= PA3;
    if( (k&PA) )   z |= SI2;
    if( (k&IP) )   z |= LIP;
    if( (k&IP) )   z |= RIP;
    if( (k&IP) )   z |= AA;
    if( (k&AA) )   z |= SI2;
    if( (k&ALT) )  z |= RALT;
    if( (k&ALT) )  z |= LALT;
    if( (k&ALT) )  z |= FLEX;
    if( (k&LALT) ) z |= PA3;
    if( (k&RALT) ) z |= PA3;
    if( (k&FLEX) ) z |= PA3;
    if( (k&FLEX) ) z |= SI2;
    if( (k&RIP) )  z |= SI2;
    if( (k&LIP) )  z |= SI2;
    return z;
}

int main(){
    int k,mct,ct;
    uint j,z,y;
    char jive[100];
    ct=mct=0;
    for(k=0; k< (1<<19); k++){
	if(INCLUDES){
	    if( (k&GRP) && !(k&MOUF) ) goto BAD;
	    if( (k&MOUF) && !(k&LBOL) ) goto BAD;
	    if( (k&MOUF) && !(k&RBOL) ) goto BAD;
	    if( (k&MOUF) && !(k&DIA) ) goto BAD;
	    if( (k&DIA) && !(k&PA) ) goto BAD;
	    if( (k&DIA) && !(k&IPALT) ) goto BAD;
	    if( (k&LBOL) && !(k&PA) ) goto BAD;
	    if( (k&LBOL) && !(k&LALT) ) goto BAD;
	    if( (k&LBOL) && !(k&LIP) ) goto BAD;
	    if( (k&RBOL) && !(k&PA) ) goto BAD;
	    if( (k&RBOL) && !(k&RALT) ) goto BAD;
	    if( (k&RBOL) && !(k&RIP) ) goto BAD;
	    if( (k&IPALT) && !(k&IP) ) goto BAD;
	    if( (k&IPALT) && !(k&ALT) ) goto BAD;
	    if( (k&IPALT) && !(k&IPLR) ) goto BAD;
	    if( (k&IPLR) && !(k&LRA) ) goto BAD;
	    if( (k&IPLR) && !(k&IP) ) goto BAD;
	    if( (k&ALT) && !(k&LRA) ) goto BAD;
	    if( (k&LRA) && !(k&LALT) ) goto BAD;
	    if( (k&LRA) && !(k&RALT) ) goto BAD;
	    if( (k&PA) && !(k&PA3) ) goto BAD;
	    if( (k&PA) && !(k&SI2) ) goto BAD;
	    if( (k&IP) && !(k&LIP) ) goto BAD;
	    if( (k&IP) && !(k&RIP) ) goto BAD;
	    if( (k&IP) && !(k&AA) ) goto BAD;
	    if( (k&AA) && !(k&SI2) ) goto BAD;
	    if( (k&ALT) && !(k&RALT) ) goto BAD;
	    if( (k&ALT) && !(k&LALT) ) goto BAD;
	    if( (k&ALT) && !(k&FLEX) ) goto BAD;
	    if( (k&LALT) && !(k&PA3) ) goto BAD;
	    if( (k&RALT) && !(k&PA3) ) goto BAD;
	    if( (k&FLEX) && !(k&PA3) ) goto BAD;
	    if( (k&FLEX) && !(k&SI2) ) goto BAD;
	    if( (k&RIP) && !(k&SI2) ) goto BAD;
	    if( (k&LIP) && !(k&SI2) ) goto BAD;
	}
	if(FINITEINC){
	    if( (k&LRA) && !(k&SI2) ) goto BAD;
	    if(QUESTIONABLE){
		if( (k&RALT) && (k&LIP) && (k&FLEX) && !(k&LALT) ) goto BAD;
		if( (k&LALT) && (k&RIP) && (k&FLEX) && !(k&RALT) ) goto BAD;
		if( (k&RALT) && (k&LIP) && (k&FLEX) && !(k&RIP) ) goto BAD;
		if( (k&LALT) && (k&RIP) && (k&FLEX) && !(k&LIP) ) goto BAD;
		if( (k&ALT) && (k&AA) && !(k&IP) ) goto BAD;
		if( (k&ALT) && (k&RIP) && !(k&IP) ) goto BAD;
		if( (k&ALT) && (k&LIP) && !(k&IP) ) goto BAD;
		if( (k&LALT) && (k&RALT) && (k&RIP) && !(k&IP) ) goto BAD;
		if( (k&LALT) && (k&RALT) && (k&LIP) && !(k&IP) ) goto BAD;
	    }
	}
	if(AAMIRROR){
	    if( (k&RBOL) && (k&AA) && !(k&LBOL) ) goto BAD;
	    if( (k&LBOL) && (k&AA) && !(k&RBOL) ) goto BAD;
	    if( (k&RIP) && (k&AA) && !(k&LIP) ) goto BAD;
	    if( (k&LIP) && (k&AA) && !(k&RIP) ) goto BAD;
	    if( (k&RALT) && (k&AA) && !(k&LALT) ) goto BAD;
	    if( (k&LALT) && (k&AA) && !(k&RALT) ) goto BAD;
	}
	if(SOPHIST){
	    if( (k&LBOL) && (k&RBOL) && !(k&MOUF) ) goto BAD;
	    if( (k&LBOL) && (k&RIP) && !(k&MOUF) ) goto BAD;
	    if( (k&RBOL) && (k&LIP) && !(k&MOUF) ) goto BAD;
	    if( (k&LBOL) && (k&FLEX) && !(k&MOUF) ) goto BAD;
	    if( (k&RBOL) && (k&FLEX) && !(k&MOUF) ) goto BAD;
	    if( (k&LBOL) && (k&RALT) && !(k&MOUF) ) goto BAD;
	    if( (k&RBOL) && (k&LALT) && !(k&MOUF) ) goto BAD;
	    if( (k&RALT) && (k&LALT) && !(k&LRA) ) goto BAD;
	    if( (k&RALT) && (k&LALT) && !(k&LRA) ) goto BAD;
	    if( (k&LRA) && (k&FLEX) && !(k&ALT) ) goto BAD;
	    if( (k&LRA) && (k&IP) && !(k&IPLR) ) goto BAD;
	    if( (k&ALT) && (k&IP) && !(k&IPALT) ) goto BAD;
	    if( (k&RIP) && (k&AA) && !(k&IP) ) goto BAD;
	    if( (k&RIP) && (k&LIP) && !(k&IP) ) goto BAD;
	    if( (k&LIP) && (k&AA) && !(k&IP) ) goto BAD;
	    if( (k&LIP) && (k&RALT) && !(k&SI2) ) goto BAD;            
	    if( (k&RIP) && (k&LALT) && !(k&SI2) ) goto BAD;            
	}
	ct++;
        if(VERBOSE){
	    printf("%d:\t",ct);
	    printf("%x\t",k);
	    if(NICEPRINT){
		z = improveneg(k);
		y = improvepos(k);

		j = k & (~y);
		if(j&GRP  ) printf("GRP, ");
		if(j&MOUF ) printf("MOUF, ");
		if(j&LBOL ) printf("LBOL, ");
		if(j&RBOL ) printf("RBOL, ");
		if(j&DIA  ) printf("DIA, ");
		if(j&PA   ) printf("PA, ");
		if(j&IPALT  ) printf("IPALT, ");
		if(j&IP   ) printf("IP, ");
		if(j&ALT  ) printf("ALT, ");
		if(j&IPLR ) printf("IPLR, ");
		if(j&LRA  ) printf("LRA, ");
		if(j&LIP  ) printf("LIP, ");
		if(j&RIP  ) printf("RIP, ");
		if(j&LALT ) printf("LALT, ");
		if(j&RALT ) printf("RALT, ");
		if(j&FLEX ) printf("FLEX, ");
		if(j&AA   ) printf("AA, ");
		if(j&PA3  ) printf("PA3, ");
		if(j&SI2  ) printf("2SI, ");

                printf("but not ");

		j = (~k) & (~z);
		if(j&GRP  ) printf("GRP, ");
		if(j&MOUF ) printf("MOUF, ");
		if(j&LBOL ) printf("LBOL, ");
		if(j&RBOL ) printf("RBOL, ");
		if(j&DIA  ) printf("DIA, ");
		if(j&PA   ) printf("PA, ");
		if(j&IPALT  ) printf("IPALT, ");
		if(j&IP   ) printf("IP, ");
		if(j&ALT  ) printf("ALT, ");
		if(j&IPLR ) printf("IPLR, ");
		if(j&LRA  ) printf("LRA, ");
		if(j&LIP  ) printf("LIP, ");
		if(j&RIP  ) printf("RIP, ");
		if(j&LALT ) printf("LALT, ");
		if(j&RALT ) printf("RALT, ");
		if(j&FLEX ) printf("FLEX, ");
		if(j&AA   ) printf("AA, ");
		if(j&PA3  ) printf("PA3, ");
		if(j&SI2  ) printf("2SI, ");
		printf("#\n");
	    }else{
		if(k&GRP  ) printf("GRP ");
		if(k&MOUF ) printf("MOUF ");
		if(k&LBOL ) printf("LBOL ");
		if(k&RBOL ) printf("RBOL ");
		if(k&DIA  ) printf("DIA ");
		if(k&PA   ) printf("PA ");
		if(k&IPALT  ) printf("IPALT ");
		if(k&IP   ) printf("IP ");
		if(k&ALT  ) printf("ALT ");
		if(k&IPLR ) printf("IPLR ");
		if(k&LRA  ) printf("LRA ");
		if(k&LIP  ) printf("LIP ");
		if(k&RIP  ) printf("RIP ");
		if(k&LALT ) printf("LALT ");
		if(k&RALT ) printf("RALT ");
		if(k&FLEX ) printf("FLEX ");
		if(k&AA   ) printf("AA ");
		if(k&PA3  ) printf("PA3 ");
		if(k&SI2  ) printf("2SI ");
		printf("\n");
	    }
	}
	if(MACING){
		mct++;
		z = improveneg(k);
		y = improvepos(k);
		spitmace4infile(k,z,y);
		sprintf(jive, "time mace4 < si%x.in > si%x.out", k,k);
		printf("Invoking %s\n", jive);
		if(REALLYMACE) system(jive);  
	}
    BAD: ;
    }
    printf("ct = %d   mct = %d\n", ct, mct);
}



