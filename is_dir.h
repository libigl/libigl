namespace igl
{
  // Act like php's is_dir function
  // http://php.net/manual/en/function.is-dir.php
  // Tells whether the given filename is a directory.
  // Input:
  //   filename  Path to the file. If filename is a relative filename, it will
  //     be checked relative to the current working directory. 
  // Returns TRUE if the filename exists and is a directory, FALSE
  // otherwise.
  bool is_dir(const char * filename);

}

// Implementation
#include <sys/stat.h>
bool igl::is_dir(const char * filename)
{
  struct stat status;
  if(stat(filename,&status)!=0)
  {
    // path does not exist
    return false;
  }
  // Tests whether existing path is a directory
  return S_ISDIR(status.st_mode);
}

