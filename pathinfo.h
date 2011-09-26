#ifndef IGL_PATHINFO_H
#define IGL_PATHINFO_H

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
  void pathinfo(
    const std::string & path,
    std::string & dirname,
    std::string & basename,
    std::string & extension,
    std::string & filename);

}

// Implementation
#include "dirname.h"
#include "basename.h"
// Verbose should be removed once everythings working correctly
#include "verbose.h"

void igl::pathinfo(
  const std::string & path,
  std::string & dirname,
  std::string & basename,
  std::string & extension,
  std::string & filename)
{
  dirname = igl::dirname(path);
  basename = igl::basename(path);
  std::string::reverse_iterator last_dot =
    std::find(
      basename.rbegin(), 
      basename.rend(), '.');
  // Was a dot found?
  if(last_dot == basename.rend())
  {
    // filename is same as basename
    filename = basename;
    // no extension
    extension = "";
  }else
  {
  // extension is substring of basename
    extension = std::string(last_dot.base(),basename.end());
    // filename is substring of basename
    filename = std::string(basename.begin(),last_dot.base()-1);
  }
}

#endif
