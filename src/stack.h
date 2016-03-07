#pragma once

typedef void *stackel;

typedef struct stack {
  size_t currMax; // allocated
  size_t n; // used
  stackel *array;
} Stack;

Stack *newStack(size_t initialMax);
void freeStack(Stack *s);

bool stackEmpty(Stack *s);
void stackPush(Stack *s, stackel l);
stackel stackPop(Stack *s);
stackel stackTop(Stack *s);
