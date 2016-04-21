#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <stdbool.h>  //bool, true, false
#include <inttypes.h> //(u)int(8|16|32|64)_t and others
#include <stddef.h>   //size_t, ptrdiff_t

#include <limits.h> //MAX/MIN constants for int types ...
#include <float.h>  //... and for floating point numbers

/* #include <tgmath.h> //type generic math (normal and complex) */
#include <math.h> //math

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
