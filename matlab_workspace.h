#ifndef IGL_WRITE_MATLAB_WORKSPACE
#define IGL_WRITE_MATLAB_WORKSPACE
#include "igl_inline.h"

namespace igl
{
  // Class which contains data of a matlab workspace which can be written to a
  // .mat file and loaded from matlab
  // 
  // This depends on matlab at compile time (though it shouldn't necessarily
  // have to) but it does not depend on running the matlab engine at run-time.
  class matlab_workspace
  {
  };
}

#ifdef IGL_HEADER_ONLY
#  include "matlab_workspace.cpp"
#endif

#endif

