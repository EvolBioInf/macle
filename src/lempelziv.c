/***** factor.c ***********************************
 * Description: Compute the longest previous factor
 *   array using a suffix array and a longest 
 *   common prefix array.
 * Reference: Crochemore, M., Ilie, L. and Smyth,
 *   W. F. (2008). A simple algorithm for com-
 *   puting the Lempel Ziv factorization. In:
 *   Data Compression Conference, p. 482-488.
 *   Computing longest previous factor in linear
 *   time and applications.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 10:29:09 2013
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "factors.h"
#include "stack.h"
#include "eprintf.h"

long minimum(long a, long b){
  if(a < b)
    return a;
  else
    return b;
}

long maximum(long a, long b){
  if(a > b)
    return a;
  else
    return b;
}

/*
 * computeLpf: Compute longest previous factor
 * Reference: M. Crochemore, L. Ilie, W.F. Smyth. 
 *   A simple algorithm for computing the Lempel-Ziv
 *   factorization, in: J.A. Storer, M.W. Marcellini
 *   (Eds.), 18th Data Compression Conference, IEEE 
 *   Computer Society, Los Alamitos, CA, 2008,
 *   pp. 482-488.
 */
size_t *computeLpf(Esa *esa){
  size_t n = esa->n;
  esa->sa = erealloc(esa->sa,(n+1)*sizeof(uint64_t));
  esa->lcp = erealloc(esa->lcp,(n+1)*sizeof(int64_t));
  size_t *lpf = (size_t *)emalloc((n+1) * sizeof(size_t));

  uint64_t *sa = esa->sa;
  int64_t *lcp = esa->lcp;
  sa[n] = -1;
  lcp[0] = 0;
  lcp[n] = 0;
  lpf[n] = 0;
  stackInit(1);
  stackPush(0);

  for(size_t i=1;i<=n;i++){
    while(!stackEmpty() &&     
	  (sa[i] < sa[stackTop()] ||  
	   (sa[i] > sa[stackTop()] && lcp[i] <= lcp[stackTop()]))){
      if(sa[i] < sa[stackTop()]){
	lpf[sa[stackTop()]] = maximum(lcp[stackTop()],lcp[i]);
	lcp[i] = minimum(lcp[stackTop()],lcp[i]);
      }else
	lpf[sa[stackTop()]] = lcp[stackTop()];
      stackPop();
    }
    if(i < n)
      stackPush(i);
  }
  freeStack();

  return lpf;
}

Fact *computeLZFact(Esa *esa){
  size_t n = esa->n;
  size_t *lpf = computeLpf(esa);

  Fact *lzf = (Fact *)emalloc(sizeof(Fact));
  lzf->fact = (size_t *)emalloc(n*sizeof(size_t));
  lzf->fact[0] = 0; 
  long i = 0;
  while(lzf->fact[i] < n){
    lzf->fact[i+1] = lzf->fact[i] + maximum(1,lpf[lzf->fact[i]]);
    i++;
  }
  lzf->lpf = lpf;
  lzf->n = i;

  lzf->str = esa->str;
  lzf->strLen = esa->n;
  return lzf;
}

