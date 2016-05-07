#pragma once
#include <cstdio>
#include <string>

std::string randSeq(size_t n, std::string alphabet = "ACGT");
std::string randRun(size_t n, size_t l, std::string alphabet = "ACGT");
double gcContent(std::string const &s);
std::string revComp(std::string const &s);

void fprintnf(FILE *fp, char const *str, int n);
