#include "prelude.h"
#include "eprintf.h"

#include "esa.h"
#include "stack.h"
#include "interval.h"

Interval *newInterval(int lcp, int lb, int rb) {
  Interval *dummy = malloc(sizeof(Interval));
  dummy->lcp = lcp;
  dummy->lb = lb;
  dummy->rb = rb;
  dummy->numChildren = 0;
  dummy->childList = NULL;
  return dummy;
}

void addChild(Interval *parent, Interval *child) {
  parent->numChildren++;
  parent->childList =
      (Interval **)erealloc(parent->childList, parent->numChildren * sizeof(Interval *));
  parent->childList[parent->numChildren - 1] = child;
}

void printInterval(Interval *in) {
  printf("%ld-[%ld..%ld]", in->lcp, in->lb, in->rb);
  if (in->numChildren)
    for (size_t i = 0; i < in->numChildren; i++) {
      Interval *child = in->childList[i];
      printf(": %ld-[%ld..%ld]", child->lcp, child->lb, child->rb);
    }
  printf("\n");
}

void printLcpTree(Interval *in) {
  printInterval(in);
  for (size_t i = 0; i < in->numChildren; i++)
    printLcpTree(in->childList[i]);
}

// returns true if an interval has been found and written to ret
bool getSubInterval(Interval *ret, Esa *esa, Interval iv, char c) {
  size_t l = iv.lb;
  size_t r = iv.rb;
  char curr;
  if (l > esa->n - 1 || r > esa->n - 1 || l > r)
    return -1;

  // get left
  size_t n=0;
  while (l <= r && n<esa->n) {
    n = (l + r) / 2;
    curr = esa->str[esa->sa[n] + iv.lcp];
    if (curr >= c)
      r = n - 1;
    else if (curr < c)
      l = n + 1;
  }
  curr = esa->str[esa->sa[l] + iv.lcp];
  int64_t nl = curr == c ? l : -1;

  // get right
  l = iv.lb;
  r = iv.rb;
  n=0;
  while (r >= l && n<esa->n) {
    n = (l + r) / 2;
    curr = esa->str[esa->sa[n] + iv.lcp];
    if (curr > c)
      r = n - 1;
    else if (curr <= c)
      l = n + 1;
  }
  curr = esa->str[esa->sa[l - 1] + iv.lcp];
  int64_t nr = curr == c ? l - 1 : -1;

  ret->lb = nl;
  ret->rb = nr;
  ret->lcp = iv.lcp + 1;

  if (nl == -1 || nr == -1)
    return false;
  return true;
}

// get interval of longest prefix matches of query using iterated binary search in ESA
// TODO: more efficient?
Interval getInterval(Esa *esa, char *query, size_t n) {
  Interval iv = {0, 0, esa->n - 1, 0, 0};
  Interval tmp = {0, 0, 0, 0, 0};
  size_t i = 0;
  for (; i < n; i++) {
    if (!getSubInterval(&tmp, esa, iv, query[i]))
      break;
    iv = tmp;
  }
  return iv;
}

// TODO: fix lcp tree stuff or implement other solution

Interval *getLcpTree(Esa *esa) {
  esa->lcp[0] = esa->lcp[esa->n] = -1;
  Interval *dummy = newInterval(-2, 0, esa->n - 1);
  Interval *root = newInterval(0, 0, esa->n - 1);
  Interval *last = NULL;

  Stack *s = newStack(esa->n);
  stackPush(s, dummy);
  stackPush(s, root);

  for (size_t i = 1; i <= esa->n; i++) {
    int lb = i - 1;
    while (esa->lcp[i] < ((Interval *)stackTop(s))->lcp) {
      ((Interval *)stackTop(s))->rb = i - 1;
      last = (Interval *)stackPop(s);
      lb = last->lb;
      if (esa->lcp[i] <= ((Interval *)stackTop(s))->lcp) {
        addChild((Interval *)stackTop(s), last);
        last = NULL;
      }
    }
    if (esa->lcp[i] > ((Interval *)stackTop(s))->lcp) {
      Interval *iv = newInterval(esa->lcp[i], lb, -1);
      if (last != NULL)
        addChild(iv, last);
      stackPush(s, iv);
    }
  }
  last = (Interval *)stackPop(s); // pop (dummy) root
  free(dummy);
  freeStack(s);
  return root;
}

// get lcp between two arbitrary suffixes
int64_t getLcp(Interval *tree, size_t i, size_t j) {
  int64_t l = MIN(i, j);
  int64_t r = MAX(i, j);
  Interval *curr = tree;
  Interval *better;
  printf("look for %zu-%zu\n", l, r);
  do {
    printf("(lb=%zu, rb=%zu, lcp=%zu)\n", curr->lb, curr->rb, curr->lcp);
    better = NULL;
    for (size_t k = 0; k < curr->numChildren; k++) {
      Interval *child = curr->childList[k];
      printf("look at (lb=%zu, rb=%zu, lcp=%zu)\n", child->lb, child->rb, child->lcp);
      if (child->lb <= l && child->rb >= r && (!better || better->lcp < child->lcp)) {
        better = child;
      }
    }
    if (better)
      curr = better;
  } while (better);
  return curr->lcp;
}
