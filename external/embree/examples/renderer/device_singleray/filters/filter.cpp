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

#include "filters/filter.h"
#include "sys/stl/array2d.h"

namespace embree
{
  void Filter::init()
  {
    float invTableSize = 1.0f / tableSize;
    Array2D<float> absoluteValues(tableSize, tableSize);

    for (uint32 x = 0; x < tableSize; ++x) {
      for (uint32 y = 0; y < tableSize; ++y) {
        Vec2f pos((x+0.5f)*invTableSize*width - width*0.5f,
                  (y+0.5f)*invTableSize*height - height*0.5f);
        absoluteValues.set(x, y, fabsf(eval(pos)));
      }
    }
    distribution.init(absoluteValues, tableSize, tableSize);
  }

  Vec2f Filter::sample(const Vec2f uv) const
  {
    Sample2f result = distribution.sample(uv);
    result.value.x = result.value.x/tableSize*width  - width*0.5f;
    result.value.y = result.value.y/tableSize*height - height*0.5f;
    return result;
  }
}
