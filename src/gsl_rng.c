/***** gsl_rng.c **********************************
 * Description: Utility functions for random number
 *   generation using the Gnu Scientific Library.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jun  3 16:39:26 2015
 **************************************************/
#include <stdio.h>
#include <time.h>
#include "gsl_rng.h"

gsl_rng *ini_gsl_rng(uint32_t useSeed){
  const gsl_rng_type *t = gsl_rng_default;
  FILE *fp = NULL;
  int seed = 0;

  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(t);
  /* seed for random number generation */
  if(useSeed){
    seed = useSeed;
  }else if((fp = fopen("randomSeed.dat","r")) != NULL){
    if(!fscanf(fp,"%d",&seed))
      printf("WARNING[gsl_rng]: Something went wrong when trying to read the the seed for the random number generator from randomSeed.dat.\n");
    fclose(fp);
  }else
    seed = -time(NULL); 

  gsl_rng_set(r,seed);
  return r;
}


void free_gsl_rng(gsl_rng *r, uint32_t useSeed){
  if(!useSeed){
    FILE *fp = fopen("randomSeed.dat","w");
    fprintf(fp,"%ld\n",gsl_rng_get(r));
    fclose(fp);
  }
  gsl_rng_free(r);
}
