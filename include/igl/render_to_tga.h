#ifndef IGL_RENDER_TO_TGA_H
#define IGL_RENDER_TO_TGA_H
#include "igl_inline.h"

#include <string>
namespace igl
{
  // Render current open GL image to .tga file
  // Inputs:
  //   tga_file  path to output .tga file
  //   width  width of scene and resulting image
  //   height height of scene and resulting image
  ///  alpha  whether to include alpha channel
  // Returns true only if no errors occured
  //
  // See also: png/render_to_png which is slower but writes .png files
  IGL_INLINE bool render_to_tga(
    const std::string tga_file,
    const int width,
    const int height,
    const bool alpha);
}

#ifdef IGL_HEADER_ONLY
#  include "render_to_tga.cpp"
#endif

#endif
