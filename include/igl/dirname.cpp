#include "dirname.h"

#include <algorithm>
#include "verbose.h"

IGL_INLINE std::string igl::dirname(const std::string & path)
{
  if(path == "")
  {
    return std::string("");
  }
  // http://stackoverflow.com/questions/5077693/dirnamephp-similar-function-in-c
  std::string::const_reverse_iterator last_slash =
    std::find(
      path.rbegin(), 
      path.rend(), '/');
  if( last_slash == path.rend() )
  {
    // No slashes found
    return std::string(".");
  }else if(1 == (last_slash.base() - path.begin()))
  {
    // Slash is first char
    return std::string("/");
  }else if(path.end() == last_slash.base() )
  {
    // Slash is last char
    std::string redo = std::string(path.begin(),path.end()-1);
    return igl::dirname(redo);
  }
  return std::string(path.begin(),last_slash.base()-1);
}
