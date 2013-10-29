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

#ifndef __EMBREE_TEXTURE_H__
#define __EMBREE_TEXTURE_H__

#include "../default.h"

namespace embree
{
  /*! Interface for different textures. Textures implement a mapping
   *  from a surface point to a color. */
  class Texture : public RefCount {
    ALIGNED_CLASS
  public:

    /*! Texture virtual destructor. */
    virtual ~Texture() {};

    /*! Returns the color for a surface point p. \param p is the
     *  location to query the color for. The range is 0 to 1. */
    virtual Color4 get(const Vec2f& p) const = 0;
  };
}

#endif
