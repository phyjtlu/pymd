/* c-routine to call brenner with input file (.inp) */

#include "mathlink.h" 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

void Brenner(long na,double *R, long Rlen, double *uc, long ulen){ 
  char* infile;
  
  /* input to Brenner */
  int np;  /* number of atoms */
  /* double R[na*3];   coordinates */
  double cell[3];   /* unitcell boxlengths*/
  int an[na];   /* atom number */

  /* output from Brenner */
  double F[na*3];  /* force on atom */
  double ea[na];    /* energy of atom */
  double etot; /* total energy */
  
  /* output to mathematica */
  double Fea[na*4];	/* force and energy put together */

  /* Other stuff */
  int ia;
  double D;  
  int i;


  np=na;

//  cell[0]=1000.;
//  cell[1]=1000.;
//  cell[2]=1000.;

  cell[0]=uc[0];
  cell[1]=uc[1];
  cell[2]=uc[2];

for(ia=0;ia<np;ia++)
    { 
      an[ia]=6;
      ea[ia]=0.;
      /*      R[3*ia+0]=0.;
      R[3*ia+1]=0.;
      R[3*ia+2]=0.; */

      F[3*ia+0]=0.;
      F[3*ia+1]=0.;
      F[3*ia+2]=0.;
    }


    brennerf_(&np,cell,an,R,ea,F,&etot); 

    for(i=0;i<na*4;i++)
    {
	    if(i<na*3)
	    {
		    Fea[i] = F[i];
	    }else{
		    Fea[i] = ea[i-na*3];
	    }
    }

    MLPutRealList(stdlink,Fea,4*na);
}

int main(int argc, char* argv[])
{

   return MLMain(argc, argv); 
  

} 
