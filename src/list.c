#include "prelude.h"
#include "eprintf.h"

#include "list.h"

List *newList() { return ecalloc(1, sizeof(List)); }

void listAppend(List **l, void *val) {
  List *item = newList();
  item->value = val;

  if (*l) {
    List *curr = *l;
    while (curr->next)
      curr = curr->next;
    curr->next = item;
  } else {
    *l = item;
  }
}

// takes pointer to pointer to list and value,
// prepends value and redirects pointer
void listPrepend(List **l, void *val) {
  List *curr = *l;
  List *item = newList();
  item->next = curr;
  item->value = val;
  *l = item;
}

// frees the list structure itself, NOT the memory at the value pointers
void freeList(List **l) {
  List *curr = *l;
  List *tmp = NULL;
  if (curr)
    do {
      tmp = curr->next;
      free(curr);
      curr = tmp;
    } while (curr);
  *l = NULL;
}
