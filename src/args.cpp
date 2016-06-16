#include <iostream>
using namespace std;

#include "args.h"
#include <getopt.h>

// globally accessible arguments for convenience
Args args;

static char const opts_short[] = "hw:k:jisln:pg:b";
static struct option const opts[] = {
    {"help", no_argument, nullptr, 'h'},
    {"window-size", required_argument, nullptr, 'w'},
    {"window-interval", required_argument, nullptr, 'k'},
    {"join", required_argument, nullptr, 'j'},
    {"load-index", no_argument, nullptr, 'i'},
    {"save-index", no_argument, nullptr, 's'},
    {"list-index", no_argument, nullptr, 'l'},
    {"seq", no_argument, nullptr, 'n'},
    {"print-factors", no_argument, nullptr, 'p'},
    {"graph", required_argument, nullptr, 'g'},
    {"benchmark", no_argument, nullptr, 'b'},
    {0, 0, 0, 0} // <- required
};

static char const usage[] = PROGNAME
    " " VERSION " (" BUILD_INFO ")\n" DESCRIPTION "\n" COPYRIGHT "\n"
    "Usage: " PROGNAME " [OPTIONS] [FILES]\n"
    "OPTIONS:\n"
    "\t-h: print this help message and exit\n"
    "\t-w <NUM>: size of sliding window (default: whole sequence length)\n"
    "\t-k <NUM>: interval between sliding windows (default: w/10)\n"
    "\t-j treat all sequences in a single file as one sequence"
    " (no effect when using existing index file, default: off)\n"

    "\t-i: use index file instead of FASTA sequence file\n"
    "\t-s: output index file for further processing (no regular result)\n"
    "\t-l: list sequences stored in index file\n"
    "\t-n <NUM>: calculate for given sequence within file (default: 0=all)\n"

    "\t-p: print match-length and Lempel-Ziv factors and periodicities\n"
    "\t-b: print benchmarking information\n"
    "\t-g N: output pipe-ready to plot with:\n"
    "\t\tN=1 -> dnalc_plot.sh (-> gnuplot)\n"
    "\t\tN=2 -> graph -T X (part of plotutils)\n";

void Args::parse(int argc, char *argv[]) {
  int c = 0;       // getopt stores value returned (last struct component) here
  int opt_idx = 0; // getopt stores the option index here.
  while ((c = getopt_long(argc, argv, opts_short, opts, &opt_idx)) != -1) {
    switch (c) {
    case 0: // long option without a short name
      /* printf("Option: %s\n", opts[opt_idx].name); */
      break;
    case 'h':
      cout << usage;
      exit(0);
      break;
    case 'w':
      args.w = atoi(optarg);
      break;
    case 'k':
      args.k = atoi(optarg);
      break;
    case 'j':
      args.j = true;
      break;

    case 'i':
      args.i = true;
      break;
    case 's':
      args.s = true;
      break;
    case 'l':
      args.l = true;
      break;
    case 'n':
      args.n = atoi(optarg);
      break;

    case 'p':
      args.p = true;
      break;
    case 'g':
      args.g = optarg == nullptr ? 0 : atoi(optarg);
      break;
    case 'b':
      args.b = true;
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
