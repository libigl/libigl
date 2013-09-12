#ifndef IGL_TEXTURE_FROM_TGA_H
#define IGL_TEXTURE_FROM_TGA_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include "OpenGL_convenience.h"

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
#endif
