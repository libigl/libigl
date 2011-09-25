#ifndef IGL_IS_FILE_H
#define IGL_IS_FILE_H
namespace igl
{
  // Act like php's is_file function
  // http://php.net/manual/en/function.is-file.php
  // Tells whether the given filename is a regular file.
  // Input:
  //   filename  Path to the file. If filename is a relative filename, it will
  //     be checked relative to the current working directory. 
  // Returns TRUE if the filename exists and is a regular file, FALSE
  // otherwise.
  inline bool is_file(const char * filename);

}

// Implementation
#include <sys/stat.h>
inline bool igl::is_file(const char * filename)
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

#endif
