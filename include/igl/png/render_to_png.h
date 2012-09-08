#ifndef IGL_RENDER_TO_PNG_H
#define IGL_RENDER_TO_PNG_H
#include <igl/igl_inline.h>

#include <string>
namespace igl
{
  //
  // Render current open GL image to .png file
  // Inputs:
  //   png_file  path to output .png file
  //   width  width of scene and resulting image
  //   height height of scene and resulting image
  //   alpha  whether to include alpha channel
  //   fast  sacrifice compression ratio for speed
  // Returns true only if no errors occured
  //
  // See also: igl/render_to_tga which is faster but writes .tga files
  IGL_INLINE bool render_to_png(
    const std::string png_file,
    const int width,
    const int height,
    const bool alpha = true,
    const bool fast = false);
}

#ifdef IGL_HEADER_ONLY
#  include "render_to_png.cpp"
#endif

#endif
