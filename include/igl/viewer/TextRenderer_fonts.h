// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Wenzel Jacob <wenzel@inf.ethz.ch>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_TEXT_RENDERER_FONTS_H
#define IGL_TEXT_RENDERER_FONTS_H

#include <stdint.h>

#ifndef IGL_STATIC_LIBRARY
namespace
{
#endif
  extern uint8_t igl_entypo_ttf[];
  extern uint32_t igl_entypo_ttf_size;

  extern uint8_t igl_roboto_bold_ttf[];
  extern uint32_t igl_roboto_bold_ttf_size;

  extern uint8_t igl_roboto_regular_ttf[];
  extern uint32_t igl_roboto_regular_ttf_size;

#ifndef IGL_STATIC_LIBRARY
}
#endif

#ifndef IGL_STATIC_LIBRARY
namespace
{
  #include "TextRenderer_fonts.cpp"
}
#endif

#endif
