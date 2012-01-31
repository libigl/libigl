#ifndef IGL_PATHINFO_H
#define IGL_PATHINFO_H
#include "igl_inline.h"

#include <string>

namespace igl
{
  //// Decided not to use these
  //const int PATHINFO_DIRNAME 01
  //const int PATHINFO_BASENAME 02
  //const int PATHINFO_EXTENSION 04
  //const int PATHINFO_FILENAME 08

  // Function like PHP's pathinfo
  //  returns information about path
  // Input:
  //  path  string containing input path
  // Outputs:
  //  dirname  string containing dirname (see dirname.h)
  //  basename  string containing basename (see basename.h)
  //  extension  string containing extension (characters after last '.')
  //  filename  string containing extension (characters of basename before last
  //    '.')
  //
  //  See also: basename, dirname
  IGL_INLINE void pathinfo(
    const std::string & path,
    std::string & dirname,
    std::string & basename,
    std::string & extension,
    std::string & filename);

}

#ifdef IGL_HEADER_ONLY
#  include "pathinfo.cpp"
#endif

#endif
