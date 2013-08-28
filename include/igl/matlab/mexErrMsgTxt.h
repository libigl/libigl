#include "../igl_inline.h"
namespace igl
{
  // Wrapper for mexErrMsgTxt that only calls error if test fails
  IGL_INLINE void mexErrMsgTxt(bool test, const char * message);
}
