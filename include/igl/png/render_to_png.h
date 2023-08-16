// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PNG_RENDER_TO_PNG_H
#define IGL_PNG_RENDER_TO_PNG_H
#include "../igl_inline.h"

#include <string>
namespace igl
{
  namespace png
  {
    /// Render current open GL image to .png file
    ///
    /// @param[in] png_file  path to output .png file
    /// @param[in] width  width of scene and resulting image
    /// @param[in] height height of scene and resulting image
    /// @param[in] alpha  whether to include alpha channel
    /// @param[in] fast  sacrifice compression ratio for speed
    /// @return true only if no errors occurred
    ///
    /// \see render_to_tga which is faster but writes .tga files
    IGL_INLINE bool render_to_png(
      const std::string png_file,
      const int width,
      const int height,
      const bool alpha = true,
      const bool fast = false);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "render_to_png.cpp"
#endif

#endif
