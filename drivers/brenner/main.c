/* c-routine to call brenner with input file (.inp) */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define na 2


int main(int argc,char* argv[]) { 
  char* infile;
  
  /* input to Brenner */
  int np;  /* number of atoms */
  double R[na*3];  /* coordinates */
  double cell[3];   /* unitcell boxlengths*/
  int an[na];   /* atom number */

  /* output from Brenner */
  double F[na*3];  /* force on atom */
  double ea[na];    /* energy of atom */
  double etot; /* total energy */
  

  /* Other stuff */
  int ia;
  double D;  
  int i;

  np=na;

  cell[0]=1000.;
  cell[1]=1000.;
  cell[2]=1000.;

  for(ia=0;ia<np;ia++)
    { 
      an[ia]=6;
      ea[ia]=0.;
      R[3*ia+0]=0.;
      R[3*ia+1]=0.;
      R[3*ia+2]=0.;

      F[3*ia+0]=0.;
      F[3*ia+1]=0.;
      F[3*ia+2]=0.;
    }

  /* test */
  /* BUG!!!: If we set 0. instead of 0.000001 something small we get NAN from Brenner... !!*/
  for(i=0;i<1;i++){  
    D = 0.8 + i*0.01;
    R[0]=0.;
    R[1]=0.;
    R[2]=0.;
    R[3]=D;
    R[4]=0.;
    R[5]=0.;
 


    brennerf_(&np,cell,an,R,ea,F,&etot); 

    
    printf("%f   %f ",D,etot);
    for(ia=0;ia<np;ia++)
      {printf("   %f  %f  %f",F[3*ia],F[3*ia+1],F[3*ia+2]);}
    printf("\n");


  }

}
