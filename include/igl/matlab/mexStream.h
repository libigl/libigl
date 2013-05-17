#ifndef IGL_MEX_STREAM_H
#define IGL_MEX_STREAM_H
#include <iostream>
namespace igl
{
  // http://stackoverflow.com/a/249008/148668
  
  // Class to implement "cout" for mex files to print to the matlab terminal
  // window.
  //
  // Insert at the beginning of mexFunction():
  //  mexStream mout;
  //  std::streambuf *outbuf = std::cout.rdbuf(&mout); 
  //  ...
  //  ALWAYS restore original buffer to avoid memory leak problems in matlab
  //  std::cout.rdbuf(outbuf);
  //
  class mexStream : public std::streambuf
  {
    public:
    protected:
      virtual std::streamsize xsputn(const char *s, std::streamsize n); 
      virtual int overflow(int c = EOF);
  }; 
}
#ifdef IGL_HEADER_ONLY
#  include "mexStream.cpp"
#endif
#endif
