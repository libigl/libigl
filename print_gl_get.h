#ifndef PRINT_gl_get_H
#define PRINT_gl_get_H

#if __APPLE__
#  include <OpenGL/gl.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#endif

namespace igl
{
  // Prints the value of pname found by issuing glGet*(pname,*)
  // Inputs:
  //   pname  enum key to gl parameter
  inline void print_gl_get(GLenum pname);
}


// implementation
#include <cstdio>
inline void igl::print_gl_get(GLenum pname)
{
  double dM[16];

  int rows = 4;
  int cols = 4;
  switch(pname)
  {
    case GL_MODELVIEW_MATRIX:
    case GL_PROJECTION_MATRIX:
    {
      rows = 4;
      cols = 4;
      glGetDoublev(pname,dM);
      for(int i = 0;i<rows;i++)
      {
        for(int j = 0;j<cols;j++)
        {
          printf("%lg ",dM[j*rows+i]);
        }
        printf("\n");
      }
      break;
    }
    default:
      fprintf(stderr,"ERROR in print_gl_get(), gl enum not recognized.\n");
  }
}

#endif
