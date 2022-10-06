#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RND drand48()

#define NSIM 1000
#define NDIV 300
#define NPOP 50
#define MAXGC 400

int main(void)
{
  FILE *fp;
  int NGC;
  int i;
  int div, gc, sim;
  int r1, r2;
  int P[NPOP], NP[NPOP/2];
  int mcount;
  int realgc;
  int nnp;
  
  fp = fopen("sim-seg.csv", "w");
  fprintf(fp, "ngc,div,sim,realgc,h\n");
  for(NGC = 0; NGC <= MAXGC; NGC += MAXGC/2)
    {
      for(sim = 0; sim < NSIM; sim++)
	{
	  realgc = 0;
	  for(i = 0; i < NPOP; i++)
	    {
	      P[i] = (i < NPOP/2 ? 0 : 1);
	    }
	  for(div = 1; div < NDIV; div++)
	    {
	      for(gc = 0; gc < NGC; gc++)
		{
		  // choice of gene conversion protocol -- count all attempts (including ineffective) or only actual events?
		  //		  do{
		  r1 = RND*NPOP;
		  r2 = RND*NPOP;
		  if(P[r1] != P[r2]) realgc++;
		  //}while(P[r2] == P[r1]);
		  P[r2] = P[r1];
		}
	      nnp = 0;
	      for(i = 0; i < NPOP; i++)
		{
		  // choice of partitioning protocol -- realistic (hypergeometric) or binomial?
		  /* do{
		  r1 = RND*NPOP;
		  }while(P[r1] == -1);
		    
		      NP[i] = P[r1];
		      P[r1] = -1;
		  */
		  /*  r1 = RND*NPOP;
		      NP[i] = P[r1];*/
		  if(RND < 0.5)
		    {
		      NP[nnp] = P[i];
		      nnp++;
		    }
			
		}
	      mcount = 0;
	      for(i = 0; i < NPOP; i++)
		{
		  r1 = RND*nnp;
		  P[i] = NP[r1];
		  mcount += (P[i] == 1);
		}
	      fprintf(fp, "%i,%i,%i,%i,%f\n", NGC, div, sim, realgc, (double)mcount/NPOP);
	    }
	}
    }
  fclose(fp);
  
  return 0;
}
  
	
