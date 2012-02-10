#include <igl/file_contents_as_string.h>
using namespace igl;
#include <cstdio>
#include <string>
using namespace std;

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
    string content;
    bool success = file_contents_as_string(argv[i],content);
    if(!success)
    {
      return 1;
    }
    printf("%s",content.c_str());
  }
  return 0;
}
