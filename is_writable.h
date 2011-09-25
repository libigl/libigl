#ifndef IGL_IS_WRITABLE_H
#define IGL_IS_WRITABLE_H
namespace igl
{
  // Check if a file exists *and* is writable like PHP's is_writable function:
  // http://www.php.net/manual/en/function.is-writable.php
  // Input:
  //   filename  path to file
  // Returns true if file exists and is writable and false if file doesn't
  // exist or *is not writable*
  inline bool is_writable(const char * filename);
}


// Implementation
#include <sys/stat.h>
#include <unistd.h>

inline bool igl::is_writable(const char* filename)
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
    return S_IWUSR & status.st_mode;
  }

  // Dealing with group member
  if( this_gid == status.st_gid )
  {
    return S_IWGRP & status.st_mode;
  }

  // Dealing with other
  return S_IWOTH & status.st_mode;
}
#endif
