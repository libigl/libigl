// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#ifndef __EMBREE_BSPLINE_FILTER_H__
#define __EMBREE_BSPLINE_FILTER_H__

#include "filters/filter.h"

namespace embree
{
  /*! Implements a B-Spline filter. */
  class BSplineFilter : public Filter
  {
  public:
    BSplineFilter () : Filter(2.0f,4.0f,4.0f,2) { init(); }

    float eval(const Vec2f distanceToCenter) const
    {
      float d = length(distanceToCenter);
      if (d > 2.0f) return 0.0f;
      else if (d < 1.0f) {
        float t = 1.0f - d;
        return ((((-3.0f*t)+3.0f)*t+3.0f)*t+1.0f)/6.0f;
      }
      else {
        float t = 2.0f - d;
        return t*t*t/6.0f;
      }
    }
  };
}

#endif
