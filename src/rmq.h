#pragma once

#include "prelude.h"

void precomputePow2RMQ(int64_t *A, size_t n, int64_t *B);
int64_t getRMQwithPow2(int64_t *A, size_t n, int64_t *B, size_t l, size_t r);
void precomputeBlockRMQ(int64_t *A, size_t n, int64_t *B);
int64_t *precomputeRMQ(int64_t *A, size_t n);
int64_t RMQ(int64_t *A, size_t n, int64_t *B, size_t l, size_t r);
