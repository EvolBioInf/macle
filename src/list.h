#pragma once

typedef struct list {
  struct list *next;
  void *value;
} List;

List *newList();
void listAppend(List **l, void *val);
void listPrepend(List **l, void *val);
