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

#ifndef __EMBREE_BOX_FILTER_H__
#define __EMBREE_BOX_FILTER_H__

#include "filters/filter.h"

namespace embree
{
  /*! Implements a simple box filter. */
  class BoxFilter : public Filter
  {
  public:

    /*! Constructs a box filter of specified half width. */
    BoxFilter (const float halfWidth = 0.5f)
      : Filter(1.414213562f * halfWidth, 2.0f * halfWidth, 2.0f * halfWidth, uint32(ceil(halfWidth - 0.5f))), halfWidth(halfWidth) {}

    float eval(const Vec2f distanceToCenter) const {
      if (fabsf(distanceToCenter.x) <= halfWidth && fabsf(distanceToCenter.y) <= halfWidth)  return 1.0f;
      else return 0.0f;
    }
  private:
    float halfWidth; //!< Half the width of the box filter
  };
}

#endif
