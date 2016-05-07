/*
 * File:   minunit.h
 * Author: Zed. A. Shaw, Sam
 *
 * @see http://c.learncodethehardway.org/book/ex30.html
 *
 * Created on 27 August 2014, 22:14
 */
#pragma once

#include <iostream>
#include <string>
#include <sstream>

// 30+i foreground, 40+i background, bold->add ;1
#define ANSI_CLR_BLACK "\x1b[30m"
#define ANSI_CLR_RED "\x1b[31m"
#define ANSI_CLR_GREEN "\x1b[32m"
#define ANSI_CLR_YELLOW "\x1b[33m"
#define ANSI_CLR_BLUE "\x1b[34m"
#define ANSI_CLR_MAGENTA "\x1b[35m"
#define ANSI_CLR_CYAN "\x1b[36m"
#define ANSI_CLR_WHITE "\x1b[37m"
#define ANSI_CLR_RESET "\x1b[0m"

std::string errmsg;
int tests_run;

#define mu_assert(test, message)                                                         \
  do {                                                                                   \
    if (!(test)) {                                                                       \
      std::stringstream sstrm;                                                           \
      sstrm << __FILE__ << ":" << __LINE__ << ": " message;                              \
      std::cerr << "[ERROR] " << sstrm.str() << std::endl;                               \
      errmsg = sstrm.str();                                                              \
      return;                                                                            \
    }                                                                                    \
  } while (0)

#define mu_assert_eq(expected, observed, message)                                        \
  do {                                                                                   \
    int64_t e = (expected);                                                              \
    int64_t o = (observed);                                                              \
    if (e != o) {                                                                        \
      std::stringstream sstrm;                                                           \
      sstrm << __FILE__ << ":" << __LINE__ << ": " message;                              \
      std::cerr << "[ERROR] " << sstrm.str() << std::endl;                               \
      std::cerr << "\t(Expected: " << e << ", Observed: " << o << ")" << std::endl;      \
      errmsg = sstrm.str();                                                              \
      return;                                                                            \
    }                                                                                    \
  } while (0)

#define mu_run_test(test)                                                                \
  do {                                                                                   \
    std::cerr << "-----  " #test << std::endl;                                           \
    test();                                                                              \
    tests_run++;                                                                         \
    if (!errmsg.empty()) {                                                               \
      return;                                                                            \
    }                                                                                    \
  } while (0)

#define RUN_TESTS(name)                                                                  \
  int main(int argc, char *argv[]) {                                                     \
    tests_run = 0;                                                                       \
    argc = argc;                                                                         \
    std::cerr << "RUNNING: " << argv[0] << std::endl;                                    \
    std::cout << "----" << std::endl << "RUNNING: " << argv[0] << std::endl;             \
    name();                                                                              \
    if (!errmsg.empty()) {                                                               \
      std::cout << ANSI_CLR_RED "FAILED:" ANSI_CLR_RESET " " << errmsg << std::endl;     \
    } else {                                                                             \
      std::cout << ANSI_CLR_GREEN "ALL TESTS PASSED" << ANSI_CLR_RESET << std::endl;     \
    }                                                                                    \
    std::cout << "Tests run: " << tests_run << std::endl;                                \
    exit(errmsg.empty() ? EXIT_SUCCESS : EXIT_FAILURE);                                  \
  }
