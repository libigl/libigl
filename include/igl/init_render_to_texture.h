#ifndef IGL_INIT_RENDER_TO_TEXTURE_H
#define IGL_INIT_RENDER_TO_TEXTURE_H
#include "igl_inline.h"
#include "OpenGL_convenience.h"
#include <cstdlib>
namespace igl
{
  // Create a texture+framebuffer+depthbuffer triplet bound for rendering into
  // the texture;
  //
  // Inputs:
  //   width  image width
  //   height image height
  // Outputs:
  //   tex_id  id of the texture
  //   fbo_id  id of the frame buffer object
  //   dfbo_id  id of the depth frame buffer object
  IGL_INLINE void init_render_to_texture(
    const size_t width,
    const size_t height,
    GLuint & tex_id,
    GLuint & fbo_id,
    GLuint & dfbo_id);
}
#ifndef IGL_STATIC_LIBRARY
#  include "init_render_to_texture.cpp"
#endif
#endif
