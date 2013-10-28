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

#ifndef __EMBREE_FILTER_H__
#define __EMBREE_FILTER_H__

#include "samplers/distribution2d.h"

namespace embree
{
  /*! Interface to pixel filters. */
  class Filter : public RefCount
  {
  public:

    /*! Constructs a filter interface. */
    Filter(float radius, float width, float height, int32 border, uint32 tableSize = 256)
      : radius(radius), width(width), height(height), border(border), tableSize(tableSize) {}

    /*! Initializes the automatic importance sampling feature. */
    void init();

    /*! Filters need a virtual destructor. */
    virtual ~Filter() {}

    /*! Evaluates the filter at some location relative to the pixel center. */
    virtual float eval(const Vec2f distanceToCenter) const = 0;

    /*! Draws a filter sample. */
    virtual Vec2f sample(const Vec2f uv) const;

  public:
    const float radius;   //!< Radius of the filter
    const float width;    //!< Width of the filter
    const float height;   //!< Height of the filter
    const int32 border;   //!< Number of border pixels of the filter

  private:
    const uint32 tableSize;      //!< Size of the table used to approximate sampling.
    Distribution2D distribution; //!< Distribution for automatic sampling.
  };
}

#endif
