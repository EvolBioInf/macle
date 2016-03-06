/***** stack.h ************************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Nov 25 17:32:50 2009
 **************************************************/
#ifndef STACK_H
#define STACK_H

void stackInit(size_t maxN);
bool stackEmpty();
void stackPush(intptr_t l);
intptr_t stackPop();
intptr_t stackTop();
void freeStack();

#endif
