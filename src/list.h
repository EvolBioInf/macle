#pragma once

typedef struct list {
  struct list *next;
  void *value;
} List;

#define eachListItem(var, l)                                                             \
  List *(var) = l;                                                                       \
  (var);                                                                                 \
  (var) = (var)->next

List *newList();
void freeList(List **l);

void listAppend(List **l, void *val);
void listPrepend(List **l, void *val);
size_t listLength(List *l);
