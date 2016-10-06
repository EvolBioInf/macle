#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

#include "args.h"
#include "index.h"
#include "util.h"
#include <getopt.h>

// globally accessible arguments for convenience
Args args;

static char const opts_short[] = "hw:k:islr:n:f:pgb";
static struct option const opts[] = {
    {"help", no_argument, nullptr, 'h'},
    {"window-size", required_argument, nullptr, 'w'},
    {"window-interval", required_argument, nullptr, 'k'},
    {"load-index", no_argument, nullptr, 'i'},
    {"save-index", no_argument, nullptr, 's'},
    {"list-index", no_argument, nullptr, 'l'},
    {"rename-regions", required_argument, nullptr, 'r'},
    {"seq", required_argument, nullptr, 'n'},
    {"batchfile", required_argument, nullptr, 'f'},
    {"print-factors", no_argument, nullptr, 'p'},
    {"graph", required_argument, nullptr, 'g'},
    {"benchmark", no_argument, nullptr, 'b'},
    {0, 0, 0, 0} // <- required
};

static char const usage[] = PROGNAME
    " " VERSION " (" BUILD_INFO ")\n" DESCRIPTION "\n" COPYRIGHT "\n"
    "Usage: " PROGNAME " [OPTIONS] FILE\n"
    "OPTIONS:\n"
    "\t-h: print this help message and exit\n"
    "\t-w NUM: size of sliding window (default: whole sequence length)\n"
    "\t-k NUM: interval between sliding windows (default: w/10)\n"

    "\t-i: use index file instead of FASTA sequence file\n"
    "\t-s: output index file for further processing (no regular result)\n"
    "\t-l: list sequences stored in index file\n"
    "\t-r FILE: file that contains a list of regions to process\n"
    "\t-n IDX:FROM-TO: calculate for given sequence and region within file\n"
    "\t   (defaults: IDX=0, FROM=0, TO=end of whole sequence. valid syntax: IDX | IDX:FROM-TO)\n"
    "\t-f FILE: file that contains a list of regions to process\n"
    "\t   (syntax like for -n with one triple per line, not usable with -w)\n"

    "\t-p: print match factors\n"
    "\t-b: print benchmarking information\n"
    "\t-g: output to plot with dnalc_plot.sh (gnuplot wrapper)\n";

size_t stol_or_fail(string s) {
  size_t n;
  try {
    n = stol(s, nullptr, 10);
  } catch (...) {
    cerr << "ERROR: failed to parse \"" << s << "\" as number!" << endl;
    exit(1);
  }
  return n;
}

Task::Task(int64_t i, size_t s, size_t e) : lbl(""), idx(i), start(s), end(e), num(0) {}
Task::Task(string str) {
  string orig = str;
  size_t sep = str.find(":");
  lbl = "";
  idx = -1;
  start = 0;
  end = 0;
  if (sep == 0) {
    cerr << "ERROR: invalid region string \"" << orig << "\"! syntax: LBL | LBL:START-END" << endl;
    exit(1);
  }
  if (sep == string::npos) { //just label
    lbl = str;
    return;
  }
  //label:start-end
  lbl = str.substr(0, sep);
  str = str.substr(sep+1, str.size()-(sep+1));
  sep = str.find("-");
  if (sep == string::npos) {
    cerr << "ERROR: invalid region string \"" << orig << "\"! syntax: LBL | LBL:START-END" << endl;
    exit(1);
  }
  start = stol_or_fail(str.substr(0, sep)) - 1;
  str = str.substr(sep+1, str.size()-(sep+1));
  end = stol_or_fail(str) - 1;
}

void Args::parse(int argc, char *argv[]) {
  int c = 0;       // getopt stores value returned (last struct component) here
  int opt_idx = 0; // getopt stores the option index here.
  vector<string> names;
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
      if (!args.tasks.empty()) {
        cerr << "ERROR: -n incompatible with -f!" << endl;
        exit(1);
      }
      args.tasks.push_back(Task(string(optarg)));
      break;
    case 'f':
      if (!args.tasks.empty()) {
        cerr << "ERROR: -f incompatible with -n!" << endl;
        exit(1);
      }
      if (!with_file_in(optarg, [&](istream &in){
        string line;
        int num=0;
        while (in >> line) {
          Task t(line);
          t.num=num++;
          args.tasks.push_back(t);
        }
        return true;
      }))
        exit(1);
      break;
    case 'r':
      if (!with_file_in(optarg, [&](istream &in){
        string line;
        while (getline(in, line)) {
          if (line.size()>MAX_LABEL_LEN) {
            cerr << "ERROR: Name list contains too long names!" << endl;
            exit(1);
          }
          if (line.find_first_of("\t\n ")!=string::npos) {
            cerr << "ERROR: Name list contains names with whitespace!" << endl;
            exit(1);
          }
          names.push_back(line);
        }
        return true;
      })) {
        exit(1);
      }
      args.newnames = names;
      sort(names.begin(),names.end());
      for (size_t j=0; j<names.size()-1; j++) {
        if (names[i]==names[j+1]) {
          cerr << "ERROR: new names are not unique!" << endl;
          exit(1);
        }
      }

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
  if (args.tasks.empty())
    args.tasks.push_back(Task(-1,0,0));

  // remaining arguments are not options and flags -> files
  if (optind < argc) {
    args.num_files = argc - optind;
    args.files = &argv[optind];
  }

  if (args.num_files > 1)
    cerr << "WARNING: processing only first file: " << args.files[0] << endl;
  if (args.w && args.tasks.size()>1) {
    cerr << "ERROR: can not use sliding window (-w) and batch mode (-f) at the same time!" << endl;
    exit(1);
  }
  if (args.g && args.tasks.size()>1) {
    cerr << "ERROR: can not use -g and batch mode (-f) at the same time!" << endl;
    exit(1);
  }
}
