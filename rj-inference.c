// inference and model selection code for plant heteroplasmy v 0.1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()
#define NN 5            // number of N values in a parameter set
#define MAXS 200        // max number of different samples per family
#define MAXD 1000       // max number of individual observations
#define MAXSTR 200      // max string length for filenames

int NSTAR = 40;         // mtDNAs per cell -- set via command line
int HBIN = 20;          // number of bins for heteroplasmy -- set via command line
int MCMC = 1e6;         // length of MCMC chain
int PRIORTYPE = 0;      // type of prior to use
int RSEED = 1;          // random seed

#define NBIN 200        // number of bins for divisions
#define VERBOSE 1       // verbose output (0 or 1)
#define OUTPUT 10       // period between file outputs
#define NTERM 250       // number of terms in Kimura pdf expansion

// structure for model parameters
typedef struct {
  double h0[MAXS];
  int n[NN];
  int germline;
} Params;

// structure for data I/O
typedef struct {
  int id;
  int stage;
  int family;
  double h;
} Data;

// structure for data internally
typedef struct {
  double h[100];
  int family;
  int nbystage[4];
  int nsamp;
} RCData;

// 2F1 hypergeometric function, specialised for Kimura distribution function
double F(int a, double b, double c, double z)
{
  double f, f0, f1, aa;
  int j;
  
  f0 = 1;
  f1 = 1.-(b*z/c);
  if(a == -1) return f1; 
  if(a == 0) return f0;

  for(j = 2; j <= b-3; j++)
    {
      aa = 1-j;
      f = (aa*(1-z)*(f1-f0)+(aa+b*z-c)*f1)/(aa-c);
      f0 = f1;
      f1 = f;
      // F(a-1, b, c, z) = ( a*(1-z)*F(a, b, c, z) - a*(1-z)*F(a+1, b, c, z) + (a + b*z - c)*F(a, b, c, z) ) / (a-c);
      // F(a, b, c, z) = ( (a+1)*(1-z)*F(a+1, b, c, z) - (a+1)*(1-z)*F(a+2, b, c, z) + (a+1 + b*z - c)*F(a+1, b, c, z) ) / (a+1-c);
    }
  return f;
  //  return ( (a+1)*(1-z)*F(a+1, b, c, z, orig) - (a+1)*(1-z)*F(a+2, b, c, z, orig) + (a+1 + b*z - c)*F(a+1, b, c, z, orig) ) / (a+1-c);
  

}

// probability mass at h = 0 for Kimura(p0, b)
double f0(double p0, double b)
{
  int i;
  double sum = 0;
  double delta;

  //  printf("called f0\n");
  for(i = 1; i < NTERM; i++)
    {
      delta = (2*i+1)*p0*(1.-p0)*F(1-i, i+2, 2, 1-p0)*pow(b, i*(i+1)/2.);    
      if(i % 2 == 1) delta *= -1; 
      sum += delta;
    }

  return (1.-p0) + sum;
}

// probability mass at h = 1 for Kimura(p0, b)
double f1(double p0, double b)
{
  int i;
  double sum = 0;
  double delta;

  //printf("called F1\n");
  for(i = 1; i < NTERM; i++)
    {
      delta = (2*i+1)*p0*(1.-p0)*F(1-i, i+2, 2, p0)*pow(b, i*(i+1)/2.);    
      if(i % 2 == 1) delta *= -1; 
      sum += delta;
    }

  return p0 + sum;
}

// probability density at h = x for Kimura(p0, b)
double phi(double x, double p0, double b)
{
  int i;
  double sum = 0;
  double delta, lastdelta = 0;

  if(fabs(x) < 0.001) return f0(p0, b);
  if(fabs(x-1) < 0.001) return f1(p0, b);
  for(i = 1; i < NTERM; i++)
    {
      delta = i*(i+1)*(2*i+1)*p0*(1-p0)*F(1-i, i+2, 2, x)*F(1-i, i+2, 2, p0)*pow(b, i*(i+1)/2);
      sum += delta;
    }

  return (sum > 0 ? sum : 0);
}

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}


// a given set of (h0, n, h) values is important for looking up probabilities
// this function turns a given set into a unique integer reference, for lookup in a 1D array
int ref(double h0, int n, double h)
{
  if(n >= NBIN || n < 0 || h0 < -1e-6 || h0 > 1+1e-6 || h < -1e-6 || h > 1+1e-6)
    {
      printf("ref error %f %i %f\n", h0, n, h);
      exit(0);
    }

  return round(h0*HBIN)*NBIN*(HBIN+1) + n*(HBIN+1) + round(h*HBIN);
}

// build the database of precomputed PDF values from the Kimura distribution
// this converts physical parameters (initial h, copy number, cell divisions) into Kimura parameters p0 and b and builds a large array with lookup function above
void BuildDB(double *hist, int nstar)
{
  double h0;
  int i;
  int n;
  double h;
  double b;
  int hbin;

  // loop through initial heteroplasmies
  for(hbin = 0; hbin <= HBIN; hbin ++)
    {
      
      h0 = (double)hbin/HBIN;
      printf("%.3f\n", h0);
      // initialise cell population

      // loop through cell divisions
      for(n = 0; n < NBIN; n++)
	{
	  if(VERBOSE)
	    printf("%f %i\n", h0, n);
	  b = pow(1.-1./nstar, n);
	  for(i = 0; i <= HBIN; i++)
	    {
	      if(n == 0)
		hist[ref(h0, n, (double)i/HBIN)] = (fabs((double)i/HBIN-h0) < (0.5/HBIN) ? 1 : 0);
	      else
		{
		  hist[ref(h0, n, (double)i/HBIN)] = phi((double)i/HBIN, h0, b);
		  if(!(i == 0 || i == HBIN)) hist[ref(h0, n, (double)i/HBIN)] *= (1./HBIN);
		}
	    }
	}
    }
}

// retrieve the precomputed PDF value for a given (h0, ndiv, h)
double prob(double h, double h0, int ndiv, double *db)
{
  if(ndiv >= NBIN) return 0;
  return db[ref(h0, ndiv, h)];
}

// compute the likelihood for a given parameterisation and set of observations
// a given family's observations are stored in two 1D arrays, containing heteroplasmy values and number of values per stage
// stages: (MYL, CYL, COL, CI) = (mother young leaf, child young lead, child old leaf, child inflorescence)
// heteroplasmy sets by stage hset = { h-stage0-1, h-stage0-2, ..., h-stage1-1, h-stage1-2, ... }
// cell divisions in each stage n = {n-stage0, n-stage1, ..}
double Likelihood(double *hset, int *n, Params P, double h0, double *db, int verbose)
{
  double cs1, cs2, cs3, cs4;
  double pcs1, pcs2, pcs3, pcs4;
  double lik;
  double lik0, lik1, lik2, lik3;
  double inc;
  int ntot;
  double tlik;
  int i;
  double olik;
  int hits[100];
  double maxolik, maxpcs1, maxpcs2, maxpcs3;
  int effn[NN];
  
  // used for debugging, to find nonzero entries
  for(i = 0; i < 100; i++)
    hits[i] = 0;

  for(i = 0; i < NN; i++)
    {
      effn[i] = (P.n[i] < 0 ? 0 : P.n[i]);
    }
  
  // total number of observations, and allocated, initialised likelihood for each
  ntot = n[0]+n[1]+n[2]+n[3];

  // germline model 0
  // MS1 -> M early leaves [n0]
  // MS1 -> CS1 [n0+n1+n2+n3]
  // CS1 -> C early leaves [n0]
  // CS1 -> CS2 [n1]
  // CS2 -> C late leaves [n1]
  // CS2 -> CS3 [n1]
  // CS3 -> C inflors [n2]

  // germline model 1 
       // MS1 -> M early leaves [n0]
       // MS1 -> CS1 [n2+n3]
       // CS1 -> C early leaves [n0]
       // CS1 -> CS2 [n]
       // CS2 -> C late leaves [n1]
       // CS1 -> inflors [n2]

       // germline model 2: all separate

       // initialise
       lik = 0;
  olik = 1;
  maxolik = maxpcs1 = maxpcs2 = maxpcs3 = 0;

  // which germline model are we using?
  // loop through observations at each different stage for each, building up expected distribution by summing over all possible values of latent variables (if appropriate)
  if(P.germline == 2) // model 2: all separate lineages
    {
      lik0 = 1;
      for(i = 0; i < n[0]; i++)
	{
	  tlik = prob(hset[i], h0, effn[0], db);                       // lik: P(mother early leaf | h0 [+0])
	  if(tlik > 0) hits[i]++;                                     
	  lik0 *= tlik;
	}
      if(n[1]+n[2]+n[3] == 0) {
	lik = lik0;                                                    // no latents to consider
      }
      else
	{
	    
	  for(cs1 = 0; cs1 <= 1.0001; cs1 += 1./HBIN)                  // for a given cs1
	    {
	      pcs1 = prob(cs1, h0, effn[2]+effn[3], db);               // latent: P(CS1 | h0 [+2,3])

	      lik1 = 1;
	      for(i = n[0]; i < n[0]+n[1]; i++)
		{
		  tlik = prob(hset[i], cs1, effn[0], db);              // lik: P(child early leaf | CS1 [+0])
		  if(tlik > 0) hits[i]++;
		  lik1 *= tlik;
		}
      	      lik2 = 1;
	      for(i = n[0]+n[1]; i < n[0]+n[1]+n[2]; i++)
		{
		  tlik = prob(hset[i], cs1, effn[1], db);              // lik: P(child late leaf | CS1 [+1])
		  if(tlik > 0) hits[i]++;
		  lik2 *= tlik;
		}

	      lik3 = 1;
	      for(i = n[0]+n[1]+n[2]; i < n[0]+n[1]+n[2]+n[3]; i++)
		{
		  tlik = prob(hset[i], cs1, effn[2], db);              // lik: P(child inflor | CS1 [+2])
		  if(tlik > 0) hits[i]++;
		  lik3 *= tlik;
		}

	      olik = pcs1*lik0*lik1*lik2*lik3;                         // P(observations | latents)*P(latents)
	      if(olik > maxolik)
		{
		  maxolik = olik;
		  maxpcs1 = cs1; maxpcs2 = cs2;
		}
	      lik += olik;
	    }
	}

    }
  else if(P.germline == 1) // model 1: separate germline
    {
      lik0 = 1;
      for(i = 0; i < n[0]; i++)
	{
	  tlik = prob(hset[i], h0, effn[0], db);                      // lik: P(mother early leaf | h0 [+0])
	  if(tlik > 0) hits[i]++;                                     
	  lik0 *= tlik;
	}
      if(n[1]+n[2]+n[3] == 0) {
	lik = lik0;                                                   // no latents to consider
      }
      else
	{
	    
	  for(cs1 = 0; cs1 <= 1.0001; cs1 += 1./HBIN)                 // for a given cs1
	    {
	      pcs1 = prob(cs1, h0, effn[2]+effn[3], db);              // latent: P(CS1 | h0 [+2,3])

	      lik1 = 1;
	      for(i = n[0]; i < n[0]+n[1]; i++)
		{
		  tlik = prob(hset[i], cs1, effn[0], db);             // lik: P(child early leaf | CS1 [+0])
		  if(tlik > 0) hits[i]++;
		  lik1 *= tlik;
		}
	      lik3 = 1;
	      for(i = n[0]+n[1]+n[2]; i < n[0]+n[1]+n[2]+n[3]; i++)
		{
		  tlik = prob(hset[i], cs1, effn[2], db);             // lik: P(child inflor | CS1 [+2])
		  if(tlik > 0) hits[i]++;
		  lik3 *= tlik;
		}

	      for(cs2 = 0; cs2 <= 1.0001; cs2 += 1./HBIN)             // for a given cs2
		{
		  pcs2 = prob(cs2, cs1, effn[0], db);                 // latent: P(CS2 | CS1 [+0])

		  lik2 = 1;
		  for(i = n[0]+n[1]; i < n[0]+n[1]+n[2]; i++)
		    {
		      tlik = prob(hset[i], cs2, effn[1], db);         // lik: P(child late leaf | CS2 [+1])
		      if(tlik > 0) hits[i]++;
		      lik2 *= tlik;
		    }

		  olik = pcs1*pcs2*lik0*lik1*lik2*lik3;               // P(observations | latents)*P(latents)
		  if(olik > maxolik)
		    {
		      maxolik = olik;
		      maxpcs1 = cs1; maxpcs2 = cs2;
		    }
		  lik += olik;
		}
	    }
	}
  

    }
  else // model 0: linear development
    {
      lik0 = 1;
      for(i = 0; i < n[0]; i++)
	{
	  tlik = prob(hset[i], h0, effn[0], db);                          // lik: P(mother early leaf | h0 [+0])
	  if(tlik > 0) hits[i]++;
	  lik0 *= tlik;
	}
      if(n[1]+n[2]+n[3] == 0) {
	lik = lik0;                                                       // no latents to consider
      }
      else
	{
	  for(cs1 = 0; cs1 <= 1.0001; cs1 += 1./HBIN)                     // for a given cs1
	    {
	      pcs1 = prob(cs1, h0, effn[0]+effn[1]+effn[2]+effn[3], db);  // latent: P(CS1 | h0 [+0,1,2,3])

	      lik1 = 1;
	      for(i = n[0]; i < n[0]+n[1]; i++)
		{
		  tlik = prob(hset[i], cs1, effn[0], db);                 // lik: P(child early leaf | CS1 [+0])
		  if(tlik > 0) hits[i]++;
		  lik1 *= tlik;
		}

	      for(cs2 = 0; cs2 <= 1.0001; cs2 += 1./HBIN)                 // for a given cs2
		{
		  pcs2 = prob(cs2, cs1, effn[0], db);                     // latent: P(CS2 | CS1 [+0])

		  lik2 = 1;
		  for(i = n[0]+n[1]; i < n[0]+n[1]+n[2]; i++)
		    {
		      tlik = prob(hset[i], cs2, effn[1], db);             // lik: P(child late leaf | CS2 [+1])
		      if(tlik > 0) hits[i]++;
		      lik2 *= tlik;
		    }

		  for(cs3 = 0; cs3 <= 1.0001; cs3 += 1./HBIN)             // for a given cs3
		    {
		      pcs3 = prob(cs3, cs2, effn[1], db);                 // latent: P(CS3 | CS2 [+1])
		      lik3 = 1;
		      for(i = n[0]+n[1]+n[2]; i < n[0]+n[1]+n[2]+n[3]; i++)
			{
			  tlik = prob(hset[i], cs3, effn[2], db);         // lik: P(child inflor | CS3 [+2])
			  if(tlik > 0) hits[i]++;
			  lik3 *= tlik;
			}

		      olik = pcs1*pcs2*pcs3*lik0*lik1*lik2*lik3;          // P(observations | latents)*P(latents)
		      if(olik > maxolik)
			{
			  maxolik = olik;
			  maxpcs1 = cs1; maxpcs2 = cs2; maxpcs3 = cs3;
			}
		      lik += olik;
		    }
		}
	    }
  
	}
    }

  if(verbose == 1)
    {
      printf("max lik %.5e at cs1 = %.3f cs2 = %.3f cs3 = %.3f\n", maxolik, maxpcs1, maxpcs2, maxpcs3);
    }

  // catch issues of zero likelihood
  // these should only come up if we have a homoplasmic initial state or a zero division number
  // but if we're not doing the Kimura calculations right some big shifts may have zero probability (they should be nonzero, albeit small)
  if(lik == 0)
    {
      printf("  ! zero likelihood !\n  caused by %.3f -> ", h0);
      for(i = 0; i < n[0]; i++) printf("%.3f ", hset[i]); printf(" | ");
      for(i = n[0]; i < n[0]+n[1]; i++) printf("%.3f ", hset[i]); printf(" | ");
      for(i = n[0]+n[1]; i < n[0]+n[1]+n[2]; i++) printf("%.3f ", hset[i]); printf(" | ");
      for(i = n[0]+n[1]+n[2]; i < n[0]+n[1]+n[2]+n[3]; i++) printf("%.3f ", hset[i]); printf(" \n\n ");
    }
  
  return lik;
}



// perturb a parameterisation
void Perturb(Params *P, int nsamp, int fixmodel)
{
  int i;
  int newgermline;
  int theta[3];
  int inc;
  
  // randomly choose germline model
  newgermline = RND*3;
  // toss a coin -- if it's head and we chose a new model, apply model change step
  if(!fixmodel && RND < 0.5 && newgermline != P->germline)
    {
      P->germline = newgermline;
    }
  else
    {
      // normal perturbation to h0 vector
      for(i = 0; i < nsamp; i++)
	{
	  P->h0[i] += gsl_ran_gaussian(0.01);
	  if(P->h0[i] < 0) P->h0[i] = 0;
	  if(P->h0[i] > 1) P->h0[i] = 1;
	}
      // increment/decrements through the n vector
      for(i = 0; i < NN; i++)
	{
	  inc = (int) gsl_ran_gaussian(1);
	  P->n[i] += inc;
	  if(P->n[i] < 0) P->n[i] = 0;
	  if(P->n[i] > NBIN-1) P->n[i] = NBIN-1;
	}
    }
}

// read data from a given file
// needs [family ID] [developmental stage] [heteroplasmy] format
void ReadData(char *fname, Data *D, int *n)
{
  FILE *fp;
  int id, stage, family;
  double h;
  char ch;
  
  fp = fopen(fname, "r");
  if(fp == NULL)
    {
      printf("Couldn't find input file %s\n", fname);
      exit(-1);
    }
  // ignore header
  do{ch = fgetc(fp);}while(ch!='\n');
  // read and populate structure using CSV values
  do{
    fscanf(fp, "%i,%i,%i,%lf", &id, &family, &stage, &h);
    if(feof(fp)) break;
    D[(*n)].id = id;
    D[(*n)].family = family;
    D[(*n)].stage = stage;
    D[(*n)].h = h/100.;
    (*n)++;
  }while(!feof(fp));
  fclose(fp);
}

// main wrapper function
int main(int argc, char *argv[])
{
  double *hist;
  FILE *fp;
  int i, j, k;
  int nbystage[NN];
  Params P, newP;
  int t;
  double llik, newllik;
  int ndata;
  Data D[MAXD];
  RCData RD[MAXD];
  int count[MAXD];
  char dbfile[MAXSTR], outfile[MAXSTR], paramstr[MAXSTR], dbparamstr[MAXSTR];
  int si, stage;
  double h;
  double thislik;
  double h0;
  int ndiv;
  int nfamily;
  double fmean[MAXD], fnorm[MAXD];
  int nsamp;
  double priors;
  char ch;
  double tmp;
  int fixmodel;
  
  // use for debugging
  //double h0test[] = {0.0135,0.0716,0.6157,0.7019,0.8533};
  //double ntest[] = {5,2,1,0,14};
  
  // process command line arguments
  if(argc != 7)
    {
      printf("Need input file, rseed, nstar, hbins, prior type, and MCMC length!\n");
      printf("Suggestion: ./rj-seed.ce [filename] 1 50 100 0 1e6\n");
      printf("  Priors: (0) uniform; (1) prod 1/n; (2) prod exp(-n); (3) small germline advantage; (4) large germline advantage; (-1) only linear\n\n");
      return -1;
    }

  RSEED = atoi(argv[2]);
  NSTAR = atoi(argv[3]);
  HBIN = atoi(argv[4]);
  PRIORTYPE = atoi(argv[5]);
  MCMC = atoi(argv[6]);

  srand48(RSEED);
  fixmodel = 0;
  if(PRIORTYPE == -1) fixmodel = 1;
  
  // parameter label for output file
  sprintf(dbparamstr, "%i-%i", NSTAR, HBIN);
  sprintf(paramstr, "%s-%i-%i-%i-%i-%i", argv[1], RSEED, NSTAR, HBIN, PRIORTYPE, MCMC);

  // initialise mean heteroplasmy arrays for each family (to guess a good initial condition)
  for(i = 0; i < MAXD; i++)
    fnorm[i] = fmean[i] = 0;
  
  nfamily = 0;
  // read data from given file, and cast into useful format
  ReadData(argv[1], D, &ndata);
  for(i = 0; i < MAXD; i++)
    count[i] = 0;
  for(i = 0; i < ndata; i++)
    {
      if(D[i].id > nsamp)
	nsamp = D[i].id;
      if(D[i].family > nfamily)
	nfamily = D[i].family;
      si = D[i].id-1;
      RD[si].h[count[si]] = D[i].h;
      RD[si].family = D[i].family;
      RD[si].nbystage[D[i].stage]++;
      count[si]++;
      fmean[D[i].family] += D[i].h;
      fnorm[D[i].family] ++;
    }
  nfamily++;

  // repeat what we read for checking
  printf("Read these observations:\n");
  for(i = 0; i < nsamp; i++)
    {
      printf("ID %i, family %i:\n", i, RD[i].family);
      count[i] = 0;
      for(j = 0; j <= 3; j++)
	{
	  for(k = 0; k < RD[i].nbystage[j]; k++)
	    {
	      printf("  stage %i\t%.3f\n", j, RD[i].h[count[i]]);
	      count[i]++;
	    }
	}
    }


  // build database of PDF values
  hist = (double*)malloc(sizeof(double)*NBIN*(HBIN+1)*(HBIN+1));
  sprintf(dbfile, "dbtrj-%s.csv", dbparamstr);
  fp = fopen(dbfile, "r");
  if(fp == NULL) {
    printf("PDF database not found.\n");
    printf("Building PDF database...\n");
    BuildDB(hist, NSTAR);

    // write to file for checking
    printf("Writing database...\n");
    fp = fopen(dbfile, "w");
    fprintf(fp, "h0,n,h,Ph\n");
    for(h0 = 0; h0 <= 1.001; h0 += 0.01)
      {
	for(ndiv = 0; ndiv < NBIN; ndiv++)
	  {
	    for(h = 0; h <= 1.0001; h += 0.01)
	      {
		fprintf(fp, "%f,%i,%f,%f\n", h0, ndiv, h, hist[ref(h0, ndiv, h)]);
	      }
	  }
      }

    fclose(fp);
  } else {
    // take advantage of precomputed PDF database
    // skip header
    printf("Found existing PDF database\n");
    printf("Reading database...");
    do {ch = fgetc(fp);}while(ch != '\n');
    while(!feof(fp))
      {
	fscanf(fp, "%lf,%i,%lf,%lf", &h0, &ndiv, &h, &tmp);
	if(!feof(fp))
  	  hist[ref(h0, ndiv, h)] = tmp;
      }
    fclose(fp);
  }
      
  // initialise parameterisation
  // can use these lines to based more automatically on input data
  P.germline = 0;
  for(i = 0; i < nfamily; i++)
    P.h0[i] = (fmean[i]/fnorm[i]+0.5)/2;
  for(i = 0; i < NN; i++) P.n[i] = 10;

  // manually chosen
  //P.h0[0] = 0.05; P.h0[1] = 0.1; P.h0[2] = 0.07; P.h0[3] = 0.62; P.h0[4] = 0.65; P.h0[5] = 0.5;
  //P.n[0] = 5; P.n[1] = 5; P.n[2] = 5; P.n[3] = 5; P.n[4] = 15;

  // compute initial likelihood
  llik = 0;
  for(i = 0; i < nsamp; i++)
    {
      thislik = Likelihood(RD[i].h, RD[i].nbystage, P, P.h0[RD[i].family], hist, 0);
      printf("%e\n", thislik);
      llik += log(thislik);
    }

  // compute initial priors
  priors = 1;
  switch(PRIORTYPE) {
  case 0: break;
  case 1: for(j = 0; j < NN; j++) priors *= 1./(P.n[j]+1); break;
  case 2: for(j = 0; j < NN; j++) priors *= exp(-P.n[j]); break;
  case 3: priors *= (P.germline > 0 ? 1 : 0.1); break;
  case 4: priors *= (P.germline > 0 ? 1 : 1e-4); break;
  }
  llik += log(priors);

  // llik now stores log(likelihood*priors)
  
  printf("Start: %e\n", llik);
  
  // initialise output file
  sprintf(outfile, "outtrjswitch-%s.csv", paramstr);
  fp = fopen(outfile, "w");
  for(i=0; i<nfamily; i++)
    fprintf(fp, "h0.%i,", i);
  for(i=0; i<NN; i++)
    fprintf(fp, "n%i,", i);
  fprintf(fp, "lik,model\n");
  
  fclose(fp);

  // start MCMC
  printf("Starting MCMC...\n");
  for(t = 0; t < MCMC; t++)
    {
      printf("-- Model now %i\n", P.germline);
      // new trial parameter set
      newP = P;
      // perturb new parameter set
      Perturb(&newP, nsamp, fixmodel);
      // output params if required
      if(VERBOSE || t % 100 == 0)
	{
	  for(i = 0; i < nfamily; i++)
	    printf("%.3f,", P.h0[i]);
	  printf(" | ");
	  for(i = 0; i < NN; i++)
	    printf("%i ", P.n[i]);
	}

      printf("\n  trying... ");
      for(i = 0; i < nfamily; i++)
	printf("%.3f,", newP.h0[i]);
      printf(" | ");
      for(i = 0; i < NN; i++)
	printf("%i ", newP.n[i]);
      printf(" | %i ", newP.germline);
	    
      // compute likelihood of new param set
      newllik = 0;      
      for(i = 0; i < nsamp; i++)
	{
	  thislik = Likelihood(RD[i].h, RD[i].nbystage, newP, newP.h0[RD[i].family], hist, 0);
  	  newllik += log(thislik);
	}
      // compute priors for new param set
      priors = 1;
      switch(PRIORTYPE) {
      case 0: break;
      case 1: for(j = 0; j < NN; j++) priors *= 1./(newP.n[j]+1); break;
      case 2: for(j = 0; j < NN; j++) priors *= exp(-newP.n[j]); break;
      case 3: priors *= (newP.germline > 0 ? 1 : 0.1); break;
      case 4: priors *= (newP.germline > 0 ? 1 : 1e-4); break;
      }
      newllik += log(priors);

      // output stats if required
      if(VERBOSE || t % 100 == 0)
	{
	  printf(" | llik %e (new %e)\n", llik, newllik);
	}

      // MCMC acceptance step
      // this is actually RJMCMC. the Jacobeans of every model transition function are unity and the proposal kernel is a deterministic function of current params, so all the associated terms in the RJMCMC acceptance rule vanish and we are left with the usual acceptance prob
      if(newllik > llik || llik-newllik < RND)
	{
	  llik = newllik;
	  P = newP;
	}
    
      // output to file if required
      if(t % OUTPUT == 0)
	{
	  fp = fopen(outfile, "a");
	  for(i=0; i<nfamily; i++)
	    fprintf(fp, "%.4f,", P.h0[i]);
	  for(i=0; i<NN; i++)
	    fprintf(fp, "%i,", P.n[i]);
	  fprintf(fp, "%.4e,%i\n", llik, P.germline);
	  fclose(fp);
	}
    }
 
  return 0;
}

  

