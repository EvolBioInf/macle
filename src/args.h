#pragma once
#include "prelude.h"

#define PROGNAME "dnalc"
#define DESCRIPTION "Calculate DNA local neighborhood complexity"
#define VERSION "0.1"
#define COPYRIGHT "Copyright (C) 2016 Anton Pirogov, Bernhard Haubold"

typedef struct Args {
  bool h;     // help message?
  uint32_t s; // seed for random number generator

  uint32_t w; // sliding window size
  uint32_t k; // sliding interval
  bool p;     // print match length decomposition?
  bool g;     // output for plotting
  uint32_t gf; // plot output format
  bool b;     // benchmark run

  // non-parameter arguments
  size_t num_files;
  char **files;
} Args;

void parseArgs(int argc, char *argv[]);
extern Args args;
