// https://stackoverflow.com/a/985807/148668
#include <cassert>
#ifndef ASSERT
#ifdef NDEBUG
#define ASSERT(x) do { (void)sizeof(x);} while (0)
#else
#define ASSERT(x) assert(x)
#endif
#endif
