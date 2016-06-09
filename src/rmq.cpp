#include <cassert>
#include <cmath>
#include <algorithm>
using namespace std;

#include "rmq.h"

#define BLOCKSIZE(n) ((size_t)(log2(n) / 4) + 1)
#define BLOCKNUM(n) (n / BLOCKSIZE(n) + 1)

uint_vec precomputePow2RMQ(uint_vec const &A) {
  auto n = A.size();
  size_t row = log2(n) + 1;
  vector<int64_t> B(n * row);
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
  uint_vec ret(B.size());
  for (size_t i=0; i<B.size(); i++)
    ret[i] = B[i];
#ifdef USE_SDSL
  sdsl::util::bit_compress(ret);
#endif
  return ret;
}

size_t getRMQwithPow2(uint_vec const &A, uint_vec const &B, size_t l,
                       size_t r) {
  auto n = A.size();
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
  size_t cand1 = B[l * row + power] - 1;
  size_t cand2 = B[l2 * row + power] - 1;
  if (A[cand1] <= A[cand2])
    return cand1;
  return cand2;
}

uint_vec precomputeBlockRMQ(uint_vec const &A) {
  auto n = A.size();
  size_t bSz = BLOCKSIZE(n);
  size_t bNum = BLOCKNUM(n);
  uint_vec B(bNum);
  for (size_t i = 0; i < bNum; i++) {
    size_t next = min((i + 1) * bSz, n);
    size_t min = SIZE_MAX;
    for (size_t j = i * bSz; j < next; j++)
      if (min > A[j])
        min = A[j];
    B[i] = min;
  }
#ifdef USE_SDSL
  sdsl::util::bit_compress(B);
#endif
  return B;
}

RMQ_custom::RMQ_custom(uint_vec const *A) : arr(A) {
  this->btab = precomputeBlockRMQ(*A);
  this->ptab = precomputePow2RMQ(this->btab);
}

size_t RMQ_custom::operator()(size_t l, size_t r) const {
  auto &A = *arr;
  auto n = arr->size();
  assert(l <= r);
  if (l == r)
    return l;

  size_t bSz = BLOCKSIZE(n);

  size_t lblock = l / bSz;
  size_t rblock = r / bSz;
  if (lblock == rblock) {
    size_t min = SIZE_MAX;
    size_t minind=0;
    for (size_t i = l; i <= r; i++)
      if (A[i] < min) {
        min = A[i];
        minind = i;
      }
    return minind;
  }

  size_t lmin = SIZE_MAX;
  size_t rmin = SIZE_MAX;
  size_t medmin = SIZE_MAX;
  size_t lminind=0, rminind=0, medminind=0;

  size_t nextbs = min((lblock + 1) * bSz, n);
  for (size_t i = l; i < nextbs; i++)
    if (A[i] < lmin) {
      lmin = A[i];
      lminind = i;
    }
  for (size_t i = rblock * bSz; i <= r; i++)
    if (A[i] < rmin) {
      rmin = A[i];
      rminind = i;
    }
  if (rblock - lblock > 1) {
    size_t btabind = getRMQwithPow2(btab, ptab, lblock + 1, rblock - 1);
    medmin = btab[btabind];
    nextbs = min((btabind+1)*bSz, n);
    for (size_t i = btabind*bSz; i < nextbs; i++)
      if (A[i] == medmin) {
        medminind = i;
        break;
      }
  }

  if (lmin <= medmin && lmin <= rmin)
    return lminind;
  if (medmin <= rmin)
    return medminind;
  return rminind;
}
