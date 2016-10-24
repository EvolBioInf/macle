#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>

#define PROGNAME "macle"
#define DESCRIPTION "Tool to calculate the global and local match complexity of DNA"
#define VERSION "0.1"
#define COPYRIGHT "Copyright (C) 2016 Anton Pirogov, Bernhard Haubold"

#ifdef USE_SDSL
#define BUILD_INFO "with SDSL"
#else
#define BUILD_INFO "32 bit"
#endif

struct Task {
  std::string lbl;
  int64_t idx;
  size_t start;
  size_t end;
  size_t num;
  Task(int64_t i, size_t s, size_t e);
  Task(std::string str);
};

struct Args {
  void parse(int argc, char *argv[]);

  bool h = false; // help message?

  uint32_t w = 0;  // sliding window size
  uint32_t k = 0;  // sliding interval

  bool i = false;  // use index (intermediate data)
  bool s = false;  // output index
  bool l = false;  // list contents of index
  std::vector<Task> tasks;  // number of sequence/region (+ offsets) in index file to work on
  std::vector<std::string> newnames; //new names for regions -> rename regions in index

  bool p = false;  // print match length decomposition?
  bool g = false;  // output for ./macle_plot.sh
  bool b = false;  // benchmark run

  // non-parameter arguments
  size_t num_files = 0;
  char **files = nullptr;
};

extern Args args;
