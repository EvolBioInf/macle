/*
 * File:   minunit.h
 * Author: Zed. A. Shaw, Sam
 *
 * @see http://c.learncodethehardway.org/book/ex30.html
 *
 * Created on 27 August 2014, 22:14
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

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

#define log_err(M) fprintf(stderr, "[ERROR] %s:%d: " M "\n", __FILE__, __LINE__)

#define mu_suite_start() char *message = NULL

#define mu_assert(test, message)                                                         \
  do {                                                                                   \
    if (!(test)) {                                                                       \
      log_err(message);                                                                  \
      return message;                                                                    \
    }                                                                                    \
  } while (0)

#define mu_run_test(test)                                                                \
  do {                                                                                   \
    fprintf(stderr, "-----  " #test "\n");                                               \
    message = test();                                                                    \
    tests_run++;                                                                         \
    if (message) {                                                                       \
      return message;                                                                    \
    }                                                                                    \
  } while (0)

#define RUN_TESTS(name)                                                                  \
  int main(int argc, char *argv[]) {                                                     \
    tests_run = 0;                                                                       \
    argc = argc;                                                                         \
    fprintf(stderr, "RUNNING: %s\n", argv[0]);                                           \
    printf("----\nRUNNING: %s\n", argv[0]);                                              \
    char *result = name();                                                               \
    if (result != 0) {                                                                   \
      printf(ANSI_CLR_RED "FAILED:" ANSI_CLR_RESET " %s\n", result);                     \
    } else {                                                                             \
      printf(ANSI_CLR_GREEN "ALL TESTS PASSED\n" ANSI_CLR_RESET);                        \
    }                                                                                    \
    printf("Tests run: %d\n", tests_run);                                                \
    exit(result == 0 ? EXIT_SUCCESS : EXIT_FAILURE);                                     \
  }

int tests_run;

#ifdef __cplusplus
}
#endif
