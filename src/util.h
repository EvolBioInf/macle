#pragma once
#include <cstdio>
#include <string>
#include <vector>
#include <list>
#include "periodicity.h"

std::string randSeq(size_t n, std::string alphabet = "ACGT");
std::string randSeq(size_t n, double gc);
std::string randRun(size_t n, size_t l, std::string alphabet = "ACGT");
double gcContent(std::string const &s);
std::string revComp(std::string const &s);
std::vector<std::list<Periodicity>> getRuns(std::string const &seq);

void fprintnf(FILE *fp, char const *str, int n);
int open_or_fail(char const *fname, int flag);
FILE *fopen_or_fail(char const *fname, char const *flags);
