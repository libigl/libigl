#include <igl/is_dir.h>
using namespace igl;
#include <cstdio>

int main(int argc, char * argv[])
{
  if(argc <= 1)
  {
    printf("USAGE:\n  ./example [path_1] [path_2] ... [path_n]\n");
    return 1;
  }
  // loop over arguments
  for(int i = 1; i < argc; i++)
  {
    printf("is_dir(%s) --> %s\n",argv[i],(is_dir(argv[i])?"TRUE":"FALSE"));
  }
  return 0;
}
