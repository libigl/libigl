// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "hsv_to_rgb.h"
#include <cmath>


template <typename T>
IGL_INLINE void igl::hsv_to_rgb(const T * hsv, T * rgb)
{
  igl::hsv_to_rgb( 
      hsv[0],hsv[1],hsv[2],
      rgb[0],rgb[1],rgb[2]);
}

template <typename T>
IGL_INLINE void igl::hsv_to_rgb( 
  const T & h, const T & s, const T & v, 
  T & r, T & g, T & b)
{
  // From medit
  double f,p,q,t,hh;
  int    i;
  hh = ((int)h % 360) / 60.;
  i = (int)floor(hh);    /* largest int <= h     */
  f = hh - i;                    /* fractional part of h */
  p = v * (1.0 - s);
  q = v * (1.0 - (s * f));
  t = v * (1.0 - (s * (1.0 - f)));

  switch(i) {
  case 0: r = v; g = t; b = p; break;
  case 1: r = q; g = v; b = p; break;
  case 2: r = p; g = v; b = t; break;
  case 3: r = p; g = q; b = v; break;
  case 4: r = t; g = p; b = v; break;
  case 5: r = v; g = p; b = q; break;
  }
}
