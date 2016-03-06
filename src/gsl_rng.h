/***** gsl_rng.h **********************************
 * Description: Utility functions for random number
 *   generation using the Gnu Scientific Library.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jun  3 16:39:24 2015
 **************************************************/
#pragma once
#include "prelude.h"

#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

gsl_rng *ini_gsl_rng(uint32_t useSeed);
void free_gsl_rng(gsl_rng *r, uint32_t useSeed);
