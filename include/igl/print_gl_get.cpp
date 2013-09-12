#include "print_gl_get.h"
#ifndef IGL_NO_OPENGL

#include <cstdio>
IGL_INLINE void igl::print_gl_get(GLenum pname)
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
