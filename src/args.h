#pragma once
#include "prelude.h"

typedef struct Args{
	bool h;     // help message?
	uint32_t s; // seed for random number generator

	uint32_t w; // sliding window size 
	uint32_t k; // sliding interval
	bool p;     // print match length decomposition?
	bool g;     // gnuplot output?
	bool b;     // benchmark run

	// non-parameter arguments
	size_t num_files;
	char **files;
} Args;

Args parseArgs(int argc, char *argv[]);
