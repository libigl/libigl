#include "render_to_png.h"
#include <YImage.hpp>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#endif

IGL_INLINE bool igl::render_to_png(
  const std::string png_file,
  const int width,
  const int height,
  const bool alpha,
  const bool fast)
{
  YImage img;
  img.resize(width,height);
  glReadPixels(
    0,
    0,
    width,
    height,
    GL_RGBA,
    GL_UNSIGNED_BYTE,
    img.data());
  img.flip();
  if(!alpha)
  {
    for(int i = 0;i<width;i++)
    for(int j = 0;j<height;j++)
    {
      img.at(i,j).a = 255;
    }
  }
  return img.save(png_file.c_str(),fast);
}
