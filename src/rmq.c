#include "prelude.h"
#include "eprintf.h"

void precomputePow2RMQ(int64_t *A, size_t n, int64_t *B) {
  size_t row = log2(n)+1;
  for (size_t i=0; i<n; i++)
    B[i*row] = i+1;
  for (size_t j=1; j<row; j++) {
    for (size_t i=0; i<n; i++) {
      size_t i2 = i+(1<<(j-1));
      if (i2 >= n) {
        B[i*row+j] = -1;
        continue;
      }
      /* printf("compare: %zu, %zu\n", A[B[i*row+(j-1)]-1], A[B[i2*row+(j-1)]-1]); */
      else if (B[i2*row+(j-1)]-1 < 0 || A[B[i*row+(j-1)]-1] <= A[B[i2*row+(j-1)]-1])
        B[i*row+j] = B[i*row+(j-1)];
      else
        B[i*row+j] = B[i2*row+(j-1)];
      /* printf("j=%zu: set %zu-%zu: %ld\n",j,i,i2, A[B[i*row+j]-1]); */
    }
  }
  /*
  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<row; j++) {
      size_t i2 = i+(1<<j)-1;
      printf("%zu-%zu: %ld, ", i, i2, B[i*row+j] == -1 ? -1 : A[B[i*row+j]-1]);
    }
    printf("\n");
  }
  */
}

int64_t getRMQwithPow2(int64_t *A, size_t n, int64_t *B, size_t l, size_t r) {
  size_t row = log2(n)+1;
  size_t diff = r-l+1;
  size_t range = 1;
  int64_t power = -1;
  while (range<=diff) {
    power++;
    range <<= 1;
  }
  range >>= 1;
  /* printf("l:%zu r:%zu range:%zu power:%zu\n",l,r, range, power); */
  size_t l2 = l+(diff-range);
  int64_t cand1 = A[B[l*row+power]-1];
  int64_t cand2 = A[B[l2*row+power]-1];
  /* printf("cand1(%zu-%zu):%ld cand2(%zu-%zu):%ld\n", l, l+range-1,
   * cand1, l2, l2+range-1, cand2);*/
  return MIN(cand1, cand2);
}

void precomputeBlockRMQ(int64_t *A, size_t n, int64_t *B) {
  size_t bSz = (size_t)log2(n)/4+1;
  size_t bNum = n/bSz+1;
  for (size_t i=0; i<bNum; i++) {
    size_t next = MIN((i+1)*bSz, n);
    int64_t min = -1;
    for (size_t j=i*bSz; j<next; j++)
      if (min==-1 || min>A[j])
        min = A[j];
    B[i] = min;
  }
}

int64_t *precomputeRMQ(int64_t *A, size_t n) {
  size_t bSz = log2(n)/4+1;
  size_t bNum = n/bSz+1;
  size_t tRow = log2(bNum)+1;
  size_t tNum = bNum*tRow;
  int64_t *B = calloc(bNum+tNum, sizeof(int64_t));
  precomputeBlockRMQ(A, n, B);
  /*
  for (size_t i=0; i<bNum; i++)
    printf("%ld ", B[i]);
  printf("\n");
  */
  precomputePow2RMQ(B, bNum, &B[bNum]);
  return B;
}

int64_t RMQ(int64_t *A, size_t n, int64_t *B, size_t l, size_t r) {
  if (l==r)
    return A[l];

  size_t bSz = log2(n)/4+1;
  size_t bNum = n/bSz+1;
  /* size_t tRow = log2(bNum)+1; */
  /* int64_t *tab = &B[bNum]; */

  size_t lblock = l/bSz;
  size_t rblock = r/bSz;
  if (lblock==rblock) {
    int64_t min=INT64_MAX;
    for (size_t i=l; i<=r; i++)
      if (A[i]<min)
        min = A[i];
    return min;
  }

  int64_t lmin = INT64_MAX;
  int64_t rmin = INT64_MAX;
  int64_t medmin = INT64_MAX;
  
  size_t nextbs = MIN((lblock+1)*bSz, n);
  for (size_t i=l; i<nextbs; i++)
    if (A[i]<lmin)
      lmin = A[i];
  for (size_t i=rblock*bSz; i<=r; i++)
    if (A[i]<rmin)
      rmin = A[i];
  /*
  for (size_t i=lblock+1; i<rblock; i++)
   if (B[i]<medmin)
     medmin = B[i];
     */
  /* printf("n:%zu bSz:%zu bNum:%zu\n", n, bSz, bNum); */
  if (rblock-lblock>1)
    medmin = getRMQwithPow2(B, bNum, &B[bNum], lblock+1, rblock-1);
  /* printf("l:%zu r:%zu lb:%zu rb:%zu lmin:%zu rmin:%zu mmin:%zu\n", l,r,lblock,rblock,lmin,rmin,medmin); */
  return MIN(medmin, MIN(lmin,rmin));
}
