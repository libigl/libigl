#include <igl/example_fun.h>
using namespace igl;
using namespace std;

int main(int argc, char * argv[])
{
  double d = 4.4;
  example_fun(d);
  int i = 4;
  example_fun(i);
#ifdef IGL_HEADER_ONLY
  const char * s = "string";
  example_fun(s);
#endif
  return 0;
}
