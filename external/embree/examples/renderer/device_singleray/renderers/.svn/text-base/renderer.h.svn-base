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

#ifndef __EMBREE_RENDERER_H_
#define __EMBREE_RENDERER_H_

#include "../api/parms.h"
#include "../api/scene.h"
#include "../api/swapchain.h"
#include "../tonemappers/tonemapper.h"
#include "../cameras/camera.h"

namespace embree
{
  /*! Renderer interface definition. */
  class Renderer : public RefCount {
    ALIGNED_CLASS
  public:

    /*! tile size to use for rendering */
    enum { TILE_SIZE = 16 };

    /*! Renderers need a virtual destructor. */
    virtual ~Renderer() {}

    /*! Renders a single frame. */
    virtual void renderFrame(const Ref<Camera>&       camera,   /*!< Camera to render from.      */
                             const Ref<BackendScene>& scene,    /*!< Scene to render.            */
                             const Ref<ToneMapper>&   toneMapper, /*!< Tonemapper to use.          */
                             Ref<SwapChain>           film,  /*!< Framebuffer to render into. */
                             int accumulate) = 0;               /*!< Accumulation mode.          */
    
  };
}

#endif
