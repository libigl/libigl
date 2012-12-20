#ifndef C_STR_H
#define C_STR_H
// http://stackoverflow.com/a/2433143/148668
// Suppose you have a function:
//   void func(const char * c);
// Then you can write:
//   func(C_STR("foo"<<1<<"bar"));
#include <string>
#include <sstream>
#define C_STR(X) static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << X).str().c_str()
#endif
