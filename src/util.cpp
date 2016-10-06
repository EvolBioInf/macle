#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstring>
#include <iostream>
#include <fstream>

//for open_or_fail
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

//for mmap stuff
#include <sys/mman.h>
#include <sys/stat.h>

#include "util.h"
#include "args.h"
using namespace std;

// generate random DNA seq of given length
string randSeq(size_t n, double gc) {
  string s;
  s.resize(n);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine gen(seed);
  bernoulli_distribution coin(0.5);
  uniform_real_distribution<double> rdistr(0, 1.0);
  for (size_t i = 0; i < n; i++) {
    if (rdistr(gen)<gc)
      s[i] = coin(gen) ? 'G' : 'C';
    else
      s[i] = coin(gen) ? 'A' : 'T';
  }
  return s;
}

// generate random sequence of given length
string randSeq(size_t n, string alphabet) {
  string s;
  s.resize(n);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine gen(seed);
  uniform_int_distribution<int> distr(0, alphabet.size() - 1);
  for (size_t i = 0; i < n; i++)
    s[i] = alphabet[distr(gen)];
  return s;
}

// calculate the GC content
double gcContent(string const &s) {
  size_t gc = 0;
  size_t at = 0;
  for (auto c : s)
    if (c == 'g' || c == 'G' || c == 'c' || c == 'C')
      gc++;
    else if (c == 'a' || c == 'A' || c == 't' || c == 'T')
      at++;
  return (double)gc/((double)gc+at);
}

// returns reverse complement DNA string
string revComp(string const &s) {
  string r = s;
  for (char &c : r)
    switch (c) {
    case 'A':
      c = 'T';
      break;
    case 'T':
      c = 'A';
      break;
    case 'G':
      c = 'C';
      break;
    case 'C':
      c = 'G';
      break;
    }
  reverse(r.begin(), r.end());
  return r;
}

std::string base_name(std::string const & path) {
  return path.substr(path.find_last_of("/\\") + 1);
}

/* open_or_fail: open file on system level and report on error */
int open_or_fail(char const *fname, int flag) {
  int fd = open(fname, flag, 0);
  if (fd < 0) {
    fprintf(stderr, "open(%s, %i) failed: %s\n", fname, flag, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return fd;
}

/* fopen_or_fail: open file on system level and report on error */
FILE *fopen_or_fail(char const *fname, char const *flags) {
  FILE *f = fopen(fname, flags);
  if (!f) {
    fprintf(stderr, "fopen(%s, %s) failed: %s\n", fname, flags, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return f;
}

// if file==nullptr, calls lambda with stdin as stream, otherwise opens file, auto-closes
template<typename T>
bool with_file(char const *file, std::function<bool(T&)> lambda, std::ios_base::openmode mode, T* def=nullptr) {
  T *streamP = def;
  std::fstream fs;
  if (file) {
    fs.open(file, mode);
    if (!fs.is_open()) {
      std::cerr << "ERROR: Could not open file: " << file << std::endl;
      return false;
    }
    streamP = &fs;
  }
  if (streamP == nullptr) {
    std::cerr << "ERROR: Invalid default file stream!" << std::endl;
    return false;
  }
  T &stream = *streamP;
  bool ret = lambda(stream);
  if (fs.is_open())
    fs.close();
  return ret;
}

bool with_file(char const *file, function<bool(fstream&)> lambda, ios_base::openmode mode) {
  return with_file(file, lambda, mode, static_cast<fstream*>(nullptr));
}
bool with_file_in(char const *file, function<bool(istream&)> lambda, ios_base::openmode mode) {
  return with_file(file, lambda, ios_base::in|mode, &cin);
}
bool with_file_out(char const *file, function<bool(ostream&)> lambda, ios_base::openmode mode) {
  return with_file(file, lambda, ios_base::out|mode, &cout);
}

size_t getFilesize(const char* filename) {
    struct stat st;
    stat(filename, &st);
    return st.st_size;
}
bool with_mmap(char const *file, function<bool(MMapReader&)> lambda) {
  int fd = open_or_fail(file, O_RDONLY);

  MMapReader r;
  if (file)
    r.sz=getFilesize(file);
  char* data = reinterpret_cast<char*>(mmap(NULL, r.sz, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0));
  if (data == MAP_FAILED) {
    cerr << "ERROR: mmap failed! Can not continue processing file! Aborting!" << endl;
    exit(EXIT_FAILURE);
  }
  r.dat = data;
  bool ret = lambda(r);
  int rc = munmap(data, r.sz);
  if (rc != 0) {
    cerr << "ERROR: Could not unmap memory! Aborting!" << endl;
    exit(EXIT_FAILURE);
  }
  close(fd);
  return ret;
}

