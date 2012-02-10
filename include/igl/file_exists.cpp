#include "file_exists.h"

#include <sys/stat.h>

IGL_INLINE bool igl::file_exists(const char* filename)
{
  struct stat status;
  return (stat(filename,&status)==0);
}
