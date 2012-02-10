#include <igl/file_exists.h>
#include <igl/is_dir.h>
#include <igl/is_file.h>
#include <igl/is_readable.h>
#include <igl/is_writable.h>
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
    printf("file_exists(%s) --> %s\n",argv[i],(file_exists(argv[i])?"TRUE":"FALSE"));
    printf("is_dir(%s) --> %s\n",argv[i],(is_dir(argv[i])?"TRUE":"FALSE"));
    printf("is_file(%s) --> %s\n",argv[i],(is_file(argv[i])?"TRUE":"FALSE"));
    printf("is_readable(%s) --> %s\n",argv[i],(is_readable(argv[i])?"TRUE":"FALSE"));
    printf("is_writable(%s) --> %s\n",argv[i],(is_writable(argv[i])?"TRUE":"FALSE"));
  }
  return 0;
}
