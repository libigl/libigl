#ifndef IGL_FILE_CONTENTS_AS_STRING_H
#define IGL_FILE_CONTENTS_AS_STRING_H

#include <string>
namespace igl
{
  // Read a files contents as plain text into a given string
  // Inputs:
  //   file_name  path to file to be read
  // Outputs:
  //   content  output string containing contents of the given file
  // Returns true on succes, false on error
  inline bool file_contents_as_string(
    const std::string file_name,
    std::string & content);
}

// Implementation
#include <fstream>
#include <cstdio>

inline bool igl::file_contents_as_string(
  const std::string file_name,
  std::string & content)
{
  std::ifstream ifs(file_name.c_str());
  // Check that opening the stream worked successfully
  if(!ifs.good())
  {
    fprintf(
      stderr,
      "IOError: file_contents_as_string() cannot open %s\n",
      file_name.c_str());
    return false;
  }
  // Stream file contents into string
  content = std::string(
    (std::istreambuf_iterator<char>(ifs)),
    (std::istreambuf_iterator<char>()));
  return true;
}
#endif
