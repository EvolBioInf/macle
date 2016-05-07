#include <iostream>
using namespace std;

#include <getopt.h>
#include "args.h"

// globally accessible arguments for convenience
Args args;

static char const opts_short[] = "hw:k:slpg:b";
static struct option const opts[] = {
    {"help", no_argument, nullptr, 'h'},
    {"window-size", required_argument, nullptr, 'w'},
    {"window-interval", required_argument, nullptr, 'k'},
    {"save-intermediate", no_argument, nullptr, 's'},
    {"load-intermediate", no_argument, nullptr, 'l'},
    {"print-factors", no_argument, nullptr, 'p'},
    {"graph", required_argument, nullptr, 'g'},
    {"benchmark", no_argument, nullptr, 'b'},
    {0, 0, 0, 0} // <- required
};

static char const usage[] =
    PROGNAME " " VERSION " - " DESCRIPTION "\n" COPYRIGHT "\n"
             "Usage: " PROGNAME " [OPTIONS] [FILES]\n"
             "OPTIONS:\n"
             "\t-h: print this help message and exit\n"
             "\t-w <NUM>: size of sliding window (default: whole sequence length)\n"
             "\t-k <NUM>: interval between sliding windows (default: w/10)\n"
             "\t-s: output intermediate data for further processing (no regular result)\n"
             "\t-l: use intermediate data from file instead of FASTA sequence file\n"
             "\t-p: print match-length and Lempel-Ziv factors and periodicities\n"
             "\t-b: print benchmarking information\n"
             "\t-g N: output pipe-ready to plot with:\n"
             "\t\tN=1 -> gnuplot -p\n"
             "\t\tN=2 -> graph -T X (part of plotutils)\n";

void Args::parse(int argc, char *argv[]) {
  int c = 0;       // getopt stores value returned (last struct component) here
  int opt_idx = 0; // getopt stores the option index here.
  while ((c = getopt_long(argc, argv, opts_short, opts, &opt_idx)) != -1) {
    switch (c) {
    case 0: // long option without a short name
      /* printf("Option: %s\n", opts[opt_idx].name); */
      break;
    case 'p':
      args.p = true;
      break;
    case 'g':
      args.gf = optarg == nullptr ? 0 : atoi(optarg);
      args.g = args.gf > 0;
      break;
    case 'b':
      args.b = true;
      break;
    case 'w':
      args.w = atoi(optarg);
      break;
    case 'k':
      args.k = atoi(optarg);
      break;
    case 's':
      args.s = true;
      break;
    case 'l':
      args.l = true;
      break;

    case 'h':
      cout << usage;
      exit(0);
      break;
    case '?': // automatic error message from getopt
      exit(1);
      break;
    default:
      cerr << "Error while parsing command line arguments. "
              "Please file a bug report."
           << endl;
      abort();
    }
  }
  // remaining arguments are not options and flags -> files
  if (optind < argc) {
    args.num_files = argc - optind;
    args.files = &argv[optind];
  }
}
