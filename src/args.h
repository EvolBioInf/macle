#pragma once
#include <cstddef>
#include <cstdint>

#define PROGNAME "dnalc"
#define DESCRIPTION "Calculate DNA local neighborhood complexity"
#define VERSION "0.1"
#define COPYRIGHT "Copyright (C) 2016 Anton Pirogov, Bernhard Haubold"

struct Args {
  void parse(int argc, char *argv[]);

  bool h = false; // help message?
  uint32_t s = 0; // seed for random number generator

  uint32_t w = 0;  // sliding window size
  uint32_t k = 0;  // sliding interval
  bool p = false;  // print match length decomposition?
  bool g = false;  // output for plotting
  uint32_t gf = 0; // plot output format
  bool b = false;  // benchmark run

  // non-parameter arguments
  size_t num_files = 0;
  char **files = nullptr;
};

extern Args args;
