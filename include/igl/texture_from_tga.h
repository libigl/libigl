#ifndef IGL_TEXTURE_FROM_TGA_H
#define IGL_TEXTURE_FROM_TGA_H
#include "igl_inline.h"

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

#include <string>

namespace igl
{
  // Read an image from a .tga file and use it as a texture
  //
  // Input:
  //  tga_file  path to .tga file
  // Output:
  //  id  of generated openGL texture
  // Returns true on success, false on failure
  IGL_INLINE bool texture_from_tga(const std::string tga_file, GLuint & id);
}

#ifdef IGL_HEADER_ONLY
#  include "texture_from_tga.cpp"
#endif

#endif
