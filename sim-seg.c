// explicitly simulate heteroplasmy distributions with divisions and gene conversion

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RND drand48()

#define NSIM 1000    // for sampling
#define NDIV 300     // maximum number of cell divisions
#define NPOP 50      // effective population size
#define MAXGC 150    // maximum gene conversion rate (AU)

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

  // open file for output
  fp = fopen("sim-seg.csv", "w");
  fprintf(fp, "ngc,div,sim,realgc,h\n");

  // loop through model parameters
  for(NGC = 0; NGC <= MAXGC; NGC += MAXGC/2)
    {
      for(sim = 0; sim < NSIM; sim++)
	{
	  realgc = 0; // tracker for actual number of gene conversion events
	  // initialise population
	  for(i = 0; i < NPOP; i++)
	    {
	      P[i] = (i < NPOP/2 ? 0 : 1);
	    }
	  // loop through divisions
	  for(div = 1; div < NDIV; div++)
	    {
	      // attempt NGC gene conversion events
	      for(gc = 0; gc < NGC; gc++)
		{
		  r1 = RND*NPOP;
		  r2 = RND*NPOP;
		  if(P[r1] != P[r2]) realgc++;
		  P[r2] = P[r1];
		}
	      nnp = 0;
	      // quick binomial division using NP as a new buffer population
	      for(i = 0; i < NPOP; i++)
		{
		  if(RND < 0.5)
		    {
		      NP[nnp] = P[i];
		      nnp++;
		    }			
		}
	      // Polya urn reamplification from NP back up to P
	      mcount = 0;
	      for(i = 0; i < NPOP; i++)
		{
		  r1 = RND*nnp;
		  P[i] = NP[r1];
		  mcount += (P[i] == 1);
		}
	      // output statistics
	      fprintf(fp, "%i,%i,%i,%i,%f\n", NGC, div, sim, realgc, (double)mcount/NPOP);
	    }
	}
    }
  fclose(fp);
  
  return 0;
}
  
	
