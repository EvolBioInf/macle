#include <cassert>
#include <cmath>
#include <algorithm>
using namespace std;

#include "rmq.h"

#define BLOCKSIZE(n) ((size_t)(log2(n) / 4) + 1)
#define BLOCKNUM(n) (n / BLOCKSIZE(n) + 1)

void precomputePow2RMQ(int64_t *A, size_t n, int64_t *B) {
  size_t row = log2(n) + 1;
  for (size_t i = 0; i < n; i++)
    B[i * row] = i + 1;
  for (size_t j = 1; j < row; j++) {
    for (size_t i = 0; i < n; i++) {
      size_t i2 = i + (1 << (j - 1));
      if (i2 >= n) {
        B[i * row + j] = -1;
        continue;
      }
      if (B[i2 * row + (j - 1)] - 1 < 0 ||
          A[B[i * row + (j - 1)] - 1] <= A[B[i2 * row + (j - 1)] - 1])
        B[i * row + j] = B[i * row + (j - 1)];
      else
        B[i * row + j] = B[i2 * row + (j - 1)];
    }
  }
}

int64_t getRMQwithPow2(int64_t const *A, size_t n, int64_t const *B, size_t l, size_t r) {
  assert(l <= r);
  size_t row = log2(n) + 1;
  size_t diff = r - l + 1;
  size_t range = 1;
  size_t power = 0;
  while (range <= diff) {
    power++;
    range <<= 1;
  }
  range >>= 1;
  power--;
  size_t l2 = l + (diff - range);
  int64_t cand1 = A[B[l * row + power] - 1];
  int64_t cand2 = A[B[l2 * row + power] - 1];
  return min(cand1, cand2);
}

void precomputeBlockRMQ(int64_t const *A, size_t n, int64_t *B) {
  size_t bSz = BLOCKSIZE(n);
  size_t bNum = BLOCKNUM(n);
  for (size_t i = 0; i < bNum; i++) {
    size_t next = min((i + 1) * bSz, n);
    int64_t min = INT64_MAX;
    for (size_t j = i * bSz; j < next; j++)
      if (min > A[j])
        min = A[j];
    B[i] = min;
  }
}

RMQ::RMQ(int64_t const *A, size_t len) : arr(A), n(len) {
  size_t bNum = BLOCKNUM(n);
  size_t tRow = log2(bNum) + 1;
  size_t tNum = bNum * tRow;
  this->tab = new int64_t[bNum + tNum];
  precomputeBlockRMQ(A, n, this->tab);
  precomputePow2RMQ(this->tab, bNum, &(this->tab[bNum]));
}

RMQ::~RMQ() { delete[] this->tab; }

int64_t RMQ::get(size_t l, size_t r) const {
  return RMQuery(this->arr, this->n, this->tab, l, r);
}

int64_t RMQuery(int64_t const *A, size_t n, int64_t const *B, size_t l, size_t r) {
  assert(l <= r);
  if (l == r)
    return A[l];

  size_t bSz = BLOCKSIZE(n);
  size_t bNum = BLOCKNUM(n);

  size_t lblock = l / bSz;
  size_t rblock = r / bSz;
  if (lblock == rblock) {
    int64_t min = INT64_MAX;
    for (size_t i = l; i <= r; i++)
      if (A[i] < min)
        min = A[i];
    return min;
  }

  int64_t lmin = INT64_MAX;
  int64_t rmin = INT64_MAX;
  int64_t medmin = INT64_MAX;

  size_t nextbs = min((lblock + 1) * bSz, n);
  for (size_t i = l; i < nextbs; i++)
    if (A[i] < lmin)
      lmin = A[i];
  for (size_t i = rblock * bSz; i <= r; i++)
    if (A[i] < rmin)
      rmin = A[i];
  if (rblock - lblock > 1)
    medmin = getRMQwithPow2(B, bNum, &B[bNum], lblock + 1, rblock - 1);
  return min(medmin, min(lmin, rmin));
}
