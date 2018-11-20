#pragma once
#include <cstdio>
#include <string>
#include <vector>
#include <list>
#include <functional>
#include <iostream>

std::string randSeq(size_t n, std::string alphabet = "ACGT");
std::string randSeq(size_t n, double gc);
double gcContent(std::string const &s);
std::string revComp(std::string const &s);

std::string base_name(std::string const & path);
int open_or_fail(char const *fname, int flag);
FILE *fopen_or_fail(char const *fname, char const *flags);

//overloading does not work here with old g++
bool with_file(char const *file, std::function<bool(std::fstream&)> lambda, std::ios_base::openmode mode);
bool with_file_in(char const *file, std::function<bool(std::istream&)> lambda, std::ios_base::openmode mode=std::ios_base::in);
bool with_file_out(char const *file, std::function<bool(std::ostream&)> lambda, std::ios_base::openmode mode=std::ios_base::out);

/*
struct MMapReader {
  char *dat=nullptr;
  size_t sz=0;
  size_t off=0;
};
bool with_mmap(char const *file, std::function<bool(MMapReader&)> lambda);
*/
