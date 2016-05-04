#include <chrono>
#include <stack>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace std::chrono;

#include "args.h"

static stack<high_resolution_clock::time_point> tp;

void tick() { tp.push(high_resolution_clock::now()); }

void tock(char const *str) {
  if (args.b) {
    auto const tp2 = high_resolution_clock::now();
    double span = duration_cast<duration<double>>(tp2 - tp.top()).count();
    tp.pop();
    cerr << "[BENCH:" << tp.size() << "] " << str << " " << setprecision(2) << fixed
         << span << "s" << endl;
  }
}
