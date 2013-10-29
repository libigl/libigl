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

#ifndef __EMBREE_NEAREST_NEIGHBOR_H__
#define __EMBREE_NEAREST_NEIGHBOR_H__

#include "image/image.h"
#include "../textures/texture.h"

namespace embree
{
  /*! Implements an image mapped texture with nearest neighbor
   *  lookup. */
  class NearestNeighbor : public Texture
  {
  public:

    /*! Construction from image. */
    NearestNeighbor (const Ref<Image>& image)
      : image(image) {
    }

    /*! Construction from parameters. */
    NearestNeighbor (const Parms& parms) {
      image = parms.getImage("image");
    }

    Color4 get(const Vec2f& p) const {
      float   s1 = p.x-floor(p.x), t1 = p.y-floor(p.y);
      int si = (int)(s1*float(image->width)), ti = (int)(t1*float(image->height));
      int ix = clamp(si, int(0), int(image->width-1));
      int iy = clamp(ti, int(0), int(image->height-1));
      return image->get(ix,iy);
    }

  protected:
    Ref<Image> image; //!< Image mapped to surface.
  };
}

#endif
