#ifndef IGL_TEXTURE_FROM_FILE_H
#define IGL_TEXTURE_FROM_FILE_H
#ifndef IGL_NO_OPENGL
#include "../igl_inline.h"

#include "../OpenGL_convenience.h"

#include <string>

namespace igl
{
  // Read an image from an image file and use it as a texture. Officially, only
  // .tga and .png are supported. Any filetype read by ImageMagick's `convert`
  // will work via an unsafe system call.
  //
  // Input:
  //  filename  path to image file
  // Output:
  //  id  of generated openGL texture
  // Returns true on success, false on failure
  IGL_INLINE bool texture_from_file(const std::string filename, GLuint & id);
}

#ifdef IGL_HEADER_ONLY
#  include "texture_from_file.cpp"
#endif

#endif
#endif



