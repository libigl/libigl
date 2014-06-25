// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "rgb_to_hsv.h"

template <typename R,typename H>
IGL_INLINE void igl::rgb_to_hsv(const R * rgb, H * hsv)
{
  // http://en.literateprograms.org/RGB_to_HSV_color_space_conversion_%28C%29
  R rgb_max = 0.0;
  R rgb_min = 1.0;
  rgb_max = (rgb[0]>rgb_max?rgb[0]:rgb_max);
  rgb_max = (rgb[1]>rgb_max?rgb[1]:rgb_max);
  rgb_max = (rgb[2]>rgb_max?rgb[2]:rgb_max);
  rgb_min = (rgb[0]<rgb_min?rgb[0]:rgb_min);
  rgb_min = (rgb[1]<rgb_min?rgb[1]:rgb_min);
  rgb_min = (rgb[2]<rgb_min?rgb[2]:rgb_min);
  //hsv[2] = rgb_max;
  hsv[2] = rgb_max;
  if(hsv[2] == 0)
  {
    hsv[0]=hsv[1]=0;
    return;
  }
  // normalize
  R rgb_n[3];
  rgb_n[0] = rgb[0]/hsv[2];
  rgb_n[1] = rgb[1]/hsv[2];
  rgb_n[2] = rgb[2]/hsv[2];
  // Recomput max min?
  rgb_max = 0;
  rgb_max = (rgb_n[0]>rgb_max?rgb_n[0]:rgb_max);
  rgb_max = (rgb_n[1]>rgb_max?rgb_n[1]:rgb_max);
  rgb_max = (rgb_n[2]>rgb_max?rgb_n[2]:rgb_max);
  rgb_min = 1;
  rgb_min = (rgb_n[0]<rgb_min?rgb_n[0]:rgb_min);
  rgb_min = (rgb_n[1]<rgb_min?rgb_n[1]:rgb_min);
  rgb_min = (rgb_n[2]<rgb_min?rgb_n[2]:rgb_min);
  hsv[1] = rgb_max - rgb_min;
  if(hsv[1] == 0)
  {
    hsv[0] = 0;
    return;
  }
  rgb_n[0] = (rgb_n[0] - rgb_min)/(rgb_max - rgb_min);
  rgb_n[1] = (rgb_n[1] - rgb_min)/(rgb_max - rgb_min);
  rgb_n[2] = (rgb_n[2] - rgb_min)/(rgb_max - rgb_min);
  // Recomput max min?
  rgb_max = 0;
  rgb_max = (rgb_n[0]>rgb_max?rgb_n[0]:rgb_max);
  rgb_max = (rgb_n[1]>rgb_max?rgb_n[1]:rgb_max);
  rgb_max = (rgb_n[2]>rgb_max?rgb_n[2]:rgb_max);
  rgb_min = 1;
  rgb_min = (rgb_n[0]<rgb_min?rgb_n[0]:rgb_min);
  rgb_min = (rgb_n[1]<rgb_min?rgb_n[1]:rgb_min);
  rgb_min = (rgb_n[2]<rgb_min?rgb_n[2]:rgb_min);
  if (rgb_max == rgb_n[0]) {
    hsv[0] = 0.0 + 60.0*(rgb_n[1] - rgb_n[2]);
    if (hsv[0] < 0.0) {
      hsv[0] += 360.0;
    }
  } else if (rgb_max == rgb_n[1]) {
    hsv[0] = 120.0 + 60.0*(rgb_n[2] - rgb_n[0]);
  } else /* rgb_max == rgb_n[2] */ {
    hsv[0] = 240.0 + 60.0*(rgb_n[0] - rgb_n[1]);
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit instanciation
template void igl::rgb_to_hsv<float, double>(float const*, double*);
template void igl::rgb_to_hsv<double, double>(double const*, double*);
#endif
