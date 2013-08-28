#include "mexErrMsgTxt.h"

// Overload mexErrMsgTxt to check an assertion then print text only if
// assertion fails
#include "mex.h"
IGL_INLINE void igl::mexErrMsgTxt(bool assertion, const char * text)
{
  if(!assertion)
  {
    ::mexErrMsgTxt(text);
  }
}

