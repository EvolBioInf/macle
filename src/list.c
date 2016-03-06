#include "list.h"
#include "eprintf.h"

List *newList() { return ecalloc(1,sizeof(List)); }

void listAppend(List **l, void *val) {
  List *item = newList();
  item->value = val;

  if (*l) {
    List *curr = *l;
    while (curr->next)
      curr = curr->next;
    curr->next = newList();
    curr->next->value = val;
  } else {
    *l = item;
  }
}

//takes pointer to pointer to list and value,
//prepends value and redirects pointer
void listPrepend(List **l, void *val) {
  List *curr = *l;
  List *item = newList();
  item->next = curr;
  item->value = val;
  *l = item;
}
