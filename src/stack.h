/***** stack.h ************************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Nov 25 17:32:50 2009
 **************************************************/
#ifndef STACK_H
#define STACK_H

void stackInit(int maxN);
int stackEmpty();
void stackPush(int l);
int stackPop();
int stackTop();
void freeStack();

#endif
