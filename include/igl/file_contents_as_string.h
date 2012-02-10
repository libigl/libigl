#ifndef IGL_FILE_CONTENTS_AS_STRING_H
#define IGL_FILE_CONTENTS_AS_STRING_H
#include "igl_inline.h"

#include <string>
namespace igl
{
  // Read a files contents as plain text into a given string
  // Inputs:
  //   file_name  path to file to be read
  // Outputs:
  //   content  output string containing contents of the given file
  // Returns true on succes, false on error
  IGL_INLINE bool file_contents_as_string(
    const std::string file_name,
    std::string & content);
}

#ifdef IGL_HEADER_ONLY
#  include "file_contents_as_string.cpp"
#endif

#endif
