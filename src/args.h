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

  uint32_t w = 0;  // sliding window size
  uint32_t k = 0;  // sliding interval
  char m = 'b';    // complexity mode (m,r,b)
  bool j = false;  // treat one file as one single sequence

  bool i = false;  // use index (intermediate data)
  bool s = false;  // output index
  bool l = false;  // list contents of index
  uint32_t n = 0;  // number of sequence in index file to work on

  bool p = false;  // print match length decomposition?
  uint32_t g = 0;  // plot output format
  bool b = false;  // benchmark run

  // non-parameter arguments
  size_t num_files = 0;
  char **files = nullptr;
};

extern Args args;
