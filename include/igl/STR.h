#ifndef STR_H
#define STR_H
// http://stackoverflow.com/a/2433143/148668
#include <string>
#include <sstream>
// Suppose you have a function:
//   void func(std::string c);
// Then you can write:
//   func(STR("foo"<<1<<"bar"));
#define STR(X) static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << X).str()
#endif 
