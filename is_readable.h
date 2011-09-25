#ifndef IGL_IS_READABLE_H
#define IGL_IS_READABLE_H
namespace igl
{
  // Check if a file is reabable like PHP's is_readable function:
  // http://www.php.net/manual/en/function.is-readable.php
  // Input:
  //   filename  path to file
  // Returns true if file exists and is readable and false if file doesn't
  // exist or *is not readable*
  inline bool is_readable(const char * filename);
}


// Implementation
#include <sys/stat.h>
#include <unistd.h>

inline bool igl::is_readable(const char* filename)
{
  // Check if file already exists
  struct stat status;
  if(stat(filename,&status)!=0)
  {
    return false;
  }

  // Get current users uid and gid
  uid_t this_uid = getuid();
  gid_t this_gid = getgid();

  // Dealing with owner
  if( this_uid == status.st_uid )
  {
    return S_IRUSR & status.st_mode;
  }

  // Dealing with group member
  if( this_gid == status.st_gid )
  {
    return S_IRGRP & status.st_mode;
  }

  // Dealing with other
  return S_IROTH & status.st_mode;
}
#endif
