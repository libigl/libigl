#include "mexStream.h"
#include "mex.h"

std::streamsize igl::mexStream::xsputn(const char *s, std::streamsize n) 
{
  mexPrintf("%.*s",n,s);
  mexEvalString("drawnow;"); // to dump string.
  return n;
}

int igl::mexStream::overflow(int c) 
{
    if (c != EOF) {
      mexPrintf("%.1s",&c);
      mexEvalString("drawnow;"); // to dump string.
    }
    return 1;
}
