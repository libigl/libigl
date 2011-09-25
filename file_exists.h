#ifndef IGL_FILE_EXISTS_H
#define IGL_FILE_EXISTS_H
namespace igl
{
  // Check if a file or directory exists like PHP's file_exists function:
  // http://php.net/manual/en/function.file-exists.php
  // Input:
  //   filename  path to file
  // Returns true if file exists and is readable and false if file doesn't
  // exist or *is not readable*
  inline bool file_exists(const char * filename);
}


// Implementation
#include <sys/stat.h>

inline bool igl::file_exists(const char* filename)
{
  struct stat status;
  return (stat(filename,&status)==0);
}
#endif
