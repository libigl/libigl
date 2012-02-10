#include "is_file.h"

#include <sys/stat.h>
IGL_INLINE bool igl::is_file(const char * filename)
{
  struct stat status;
  if(stat(filename,&status)!=0)
  {
    // path does not exist
    return false;
  }
  // Tests whether existing path is a regular file
  return S_ISREG(status.st_mode);
}
