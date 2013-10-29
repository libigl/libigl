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

#ifndef __EMBREE_DEFAULT_TONEMAPPER_H__
#define __EMBREE_DEFAULT_TONEMAPPER_H__

#include "../tonemappers/tonemapper.h"

namespace embree
{
  /*! Default tonemapper. Performs gamma correction and vignette effect. */
  class DefaultToneMapper : public ToneMapper {
    ALIGNED_CLASS
  public:

    /*! Default tonemapper constructor. */
    DefaultToneMapper (const Parms& parms) 
    {
      gamma = parms.getFloat("gamma",1.0f);
      rcpGamma   = rcp(gamma);
      vignetting = parms.getBool("vignetting",false);
    }

    /*! Evaluates the tonemapper, */
    virtual Color eval (const Color& color, const int x, const int y, const Ref<SwapChain>& swapchain) const
    {
      Color color0 = color;
      if (unlikely(gamma != 1.0f)) 
        color0 = pow(color0,rcpGamma);
      
      if (unlikely(vignetting)) {
        const int width = swapchain->getWidth();
        const int height = swapchain->getHeight();
        float d = length((Vec2f(float(x),float(y)) - 0.5f*Vec2f(float(width),float(height))) * rcp(0.5f*float(width)));
        color0 *= powf(cosf(d*0.5f),3.0f);
      }
      return color0;
    }

  protected:
    float gamma;     //!< Gamma value
    float rcpGamma;  //!< Reciprocal gamma value.
    bool vignetting; //!< Add a vignetting effect.
  };
}

#endif
