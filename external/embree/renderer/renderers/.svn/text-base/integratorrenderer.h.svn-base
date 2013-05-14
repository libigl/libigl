// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_INTEGRATOR_RENDERER_H__
#define __EMBREE_INTEGRATOR_RENDERER_H__

#include "renderers/renderer.h"
#include "integrators/integrator.h"
#include "samplers/sampler.h"
#include "filters/filter.h"
#include "renderers/progress.h"

namespace embree
{
  /*! Renderer that uses a given integrator, sampler, and pixel
   *  filter. */
  class IntegratorRenderer : public Renderer
  {
  public:

    /*! Construction from parameters. */
    IntegratorRenderer (const Parms& parms);

    /*! Renders a single frame. */
    void renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain<> > film, int accumulate);

  private:

    /*! Render function called once for each thread and frame. */
    void renderThread();
    static void run_renderThread(const TaskScheduler::ThreadInfo& thread, IntegratorRenderer* This, size_t) { This->renderThread(); }

    /*! Configuration */
  private:
    int maxDepth;                  //!< Maximal recursion depth.
    float gamma;                   //!< Gamma to use for framebuffer writeback.
    
  private:
    Ref<Integrator> integrator;    //!< Integrator to use.
    Ref<SamplerFactory> samplers;  //!< Sampler to use.
    Ref<Filter> filter;            //!< Pixel filter to use.

    /*! Arguments of renderFrame function */
  private:
    Ref<Camera> camera;            //!< Camera to render from.
    Ref<BackendScene> scene;       //!< Scene to render.
    Ref<ToneMapper> toneMapper;    //!< Tonemapper to use.
    Ref<FrameBuffer> framebuffer;  //!< Framebuffer to render into.
    Ref<SwapChain<> > swapchain;   //!< Swapchain to render into

    /*! Precomputations. */
  private:
    float rcpWidth;                //!< Reciprocal width of framebuffer.
    float rcpHeight;               //!< Reciprocal height of framebuffer.
    size_t numTilesX;              //!< Number of tiles in x direction.
    size_t numTilesY;              //!< Number of tiles in y direction.
    int iteration;                 //!< Accumulation iteration of framebuffer.

  private:
    int accumulate;                //!< Accumulation mode
    Atomic tileID;                 //!< ID of current tile
    Atomic finishedTiles;          //!< number of finished tiles
    Atomic atomicNumRays;          //!< for counting number of shoot rays
  private:
    bool showProgress;             //!< Set to true if user wants rendering progress shown
    Progress progress;             //!< Progress printer
  };
}

#endif
