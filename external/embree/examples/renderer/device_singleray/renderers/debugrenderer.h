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

#ifndef __EMBREE_DEBUG_RENDERER_H__
#define __EMBREE_DEBUG_RENDERER_H__

#include "../renderers/renderer.h"

namespace embree
{
  /*! Simple renderer for testing the ray shooting performance. The
   *  renderer performs a series of diffuse bounces when given a
   *  recursion depth greater than 1. */
  class DebugRenderer : public Renderer
  {
  public:

    /*! Construction from parameters. */
    DebugRenderer (const Parms& parms);

    /*! Renders a single frame. */
    void renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain > film, int accumulate);

  private:

    class RenderJob
    {
    public:
      RenderJob (Ref<DebugRenderer> renderer, const Ref<Camera>& camera, const Ref<BackendScene>& scene, 
                 const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate);
       
    private:

      /*! start functon */
      TASK_RUN_FUNCTION(RenderJob,renderTile);
      
      /*! finish function */
      TASK_COMPLETE_FUNCTION(RenderJob,finish);

      /*! Arguments of renderFrame function */
    private:
      Ref<DebugRenderer> renderer;
      Ref<Camera> camera;            //!< Camera to render from
      Ref<BackendScene> scene;       //!< Scene to render
      Ref<ToneMapper> toneMapper;    //!< Tonemapper to use.
      Ref<FrameBuffer> framebuffer;  //!< Framebuffer to render into
      Ref<SwapChain > swapchain;   //!< Swapchain to render into
      int accumulate;                //!< Accumulation mode

      /*! Precomputations. */
    private:
      float rcpWidth;                //!< Reciprocal width of framebuffer.
      float rcpHeight;               //!< Reciprocal height of framebuffer.
      size_t numTilesX;              //!< Number of tiles in x direction.
      size_t numTilesY;              //!< Number of tiles in y direction.
      
    private:
      double t0;                     //!< start time of rendering
      Atomic tileID;                 //!< ID of current tile
      Atomic atomicNumRays;          //!< for counting number of shoot rays
      TaskScheduler::Task task;
    };

    /*! Configuration */
  private:
    size_t maxDepth;                  //!< Maximal recursion depth
    size_t spp;
  };
}

#endif
