#ifndef IGL_TEXTURE_FROM_PNG_H
#define IGL_TEXTURE_FROM_PNG_H
#ifndef IGL_NO_OPENGL
#include "../igl_inline.h"

#include "../OpenGL_convenience.h"

#include <string>

namespace igl
{
  // Read an image from a .png file and use it as a texture
  //
  // Input:
  //  png_file  path to .png file
  // Output:
  //  id  of generated openGL texture
  // Returns true on success, false on failure
  IGL_INLINE bool texture_from_png(const std::string png_file, GLuint & id);
}

#ifdef IGL_HEADER_ONLY
#  include "texture_from_png.cpp"
#endif

#endif
#endif

