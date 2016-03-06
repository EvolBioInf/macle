#include "prelude.h"

#include <getopt.h>
#include "args.h"

// globally accessible arguments for convenience
Args args = {0, 0, 0, 0, 0, 0, 0, 0, 0};

char opts_short[] = "hs:w:k:pgb";
static struct option opts[] = {
    {"help", no_argument, NULL, 'h'},
    {"seed", required_argument, NULL, 's'},
    {"window-size", required_argument, NULL, 'w'},
    {"window-interval", required_argument, NULL, 'k'},
    {"print-factors", no_argument, NULL, 'p'},
    {"gnuplot", no_argument, NULL, 'g'},
    {"benchmark", no_argument, NULL, 'b'},
    {0, 0, 0, 0} // <- required
};

char usage[] = PROGNAME
    " " VERSION " - " DESCRIPTION "\n" COPYRIGHT "\n"
    "Usage: " PROGNAME " [OPTIONS] [FILES]\n"
    "OPTIONS:\n"
    "\t-h: print this help message and exit\n"
    "\t-s <NUM>: seed for random number generator (default: generated internally)\n"
    "\t-w <NUM>: size of sliding window\n"
    "\t-k <NUM>: interval between sliding windows\n"
    "\t-p: print match-length and Lempel-Ziv factors and periodicities\n"
    "\t-g: output with gnuplot commands (to pipe directly into 'gnuplot -p')\n"
    "\t-b: print benchmarking information\n";

void parseArgs(int argc, char *argv[]) {
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
      args.g = true;
      break;
    case 'b':
      args.b = true;
      break;
    case 's':
      args.s = atoi(optarg);
      break;
    case 'w':
      args.w = atoi(optarg);
      break;
    case 'k':
      args.k = atoi(optarg);
      break;

    case 'h':
      printf("%s", usage);
      exit(0);
      break;
    case '?': // automatic error message from getopt
      exit(1);
      break;
    default:
      fprintf(stderr, "Error while parsing command line arguments. "
                      "Please file a bug report.\n");
      abort();
    }
  }
  // remaining arguments are not options and flags -> files
  if (optind < argc) {
    args.num_files = argc - optind;
    args.files = &argv[optind];
  }
}
