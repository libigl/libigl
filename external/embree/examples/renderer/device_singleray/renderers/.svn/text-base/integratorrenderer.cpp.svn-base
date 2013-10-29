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

#include "renderers/integratorrenderer.h"

/* include all integrators */
#include "integrators/pathtraceintegrator.h"

/* include all samplers */
#include "samplers/sampler.h"

/* include all image filters */
#include "filters/boxfilter.h"
#include "filters/bsplinefilter.h"

namespace embree
{
  IntegratorRenderer::IntegratorRenderer(const Parms& parms)
    : iteration(0)
  {
    /*! create integrator to use */
    std::string _integrator = parms.getString("integrator","pathtracer");
    if (_integrator == "pathtracer") integrator = new PathTraceIntegrator(parms);
    else throw std::runtime_error("unknown integrator type: "+_integrator);

    /*! create sampler to use */
    std::string _samplers = parms.getString("sampler","multijittered");
    if (_samplers == "multijittered"   ) samplers = new SamplerFactory(parms);
    else throw std::runtime_error("unknown sampler type: "+_samplers);

    /*! create pixel filter to use */
    std::string _filter = parms.getString("filter","bspline");
    if      (_filter == "none"   ) filter = NULL;
    else if (_filter == "box"    ) filter = new BoxFilter;
    else if (_filter == "bspline") filter = new BSplineFilter;
    else throw std::runtime_error("unknown filter type: "+_filter);

    /*! get framebuffer configuration */
    gamma = parms.getFloat("gamma",1.0f);

    /*! show progress to the user */
    showProgress = parms.getInt("showprogress",0);
  }

  void IntegratorRenderer::renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate) 
  {
    if (accumulate == 0) iteration = 0;
    new RenderJob(this,camera,scene,toneMapper,swapchain,accumulate,iteration);
    iteration++;
  }

  IntegratorRenderer::RenderJob::RenderJob (Ref<IntegratorRenderer> renderer, const Ref<Camera>& camera, const Ref<BackendScene>& scene, 
                                            const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate, int iteration)
    : renderer(renderer), camera(camera), scene(scene), toneMapper(toneMapper), swapchain(swapchain), 
      accumulate(accumulate), iteration(iteration), tileID(0), atomicNumRays(0)
  {
    numTilesX = ((int)swapchain->getWidth() +TILE_SIZE-1)/TILE_SIZE;
    numTilesY = ((int)swapchain->getHeight()+TILE_SIZE-1)/TILE_SIZE;
    rcpWidth  = rcp(float(swapchain->getWidth()));
    rcpHeight = rcp(float(swapchain->getHeight()));
    this->framebuffer = swapchain->buffer();
    this->framebuffer->startRendering(numTilesX*numTilesY);
    if (renderer->showProgress) new (&progress) Progress(numTilesX*numTilesY);

    if (renderer->showProgress) progress.start();
    renderer->samplers->reset();
    renderer->integrator->requestSamples(renderer->samplers, scene);
    renderer->samplers->init(iteration,renderer->filter);

#if 1
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_renderTile,this,TaskScheduler::getNumThreads(),_finish,this,"render::tile");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
    rtcDebug();
#else
    new (&task) TaskScheduler::Task (NULL,_renderTile,this,TaskScheduler::getNumThreads(),_finish,this,"render::tile");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
#endif
  }

  void IntegratorRenderer::RenderJob::finish(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    if (renderer->showProgress) progress.end();
    double dt = getSeconds()-t0;

     /*! print fps, render time, and rays per second */
    std::ostringstream stream;
    stream << "render  ";
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(2);
    stream << 1.0f/dt << " fps, ";
    stream.precision(0);
    stream << dt*1000.0f << " ms, ";
    stream.precision(3);
    stream << atomicNumRays/dt*1E-6 << " mrps";
    std::cout << stream.str() << std::endl;

    delete this;
  }

  void IntegratorRenderer::RenderJob::renderTile(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    /*! create a new sampler */
    IntegratorState state;
    if (taskIndex == taskCount-1) t0 = getSeconds();
    
    /*! tile pick loop */
    while (true)
    {
      /*! pick a new tile */
      size_t tile = tileID++;
      if (tile >= numTilesX*numTilesY) break;

      /*! process all tile samples */
      const int tile_x = (tile%numTilesX)*TILE_SIZE;
      const int tile_y = (tile/numTilesX)*TILE_SIZE;
      Random randomNumberGenerator(tile_x * 91711 + tile_y * 81551 + 3433*swapchain->firstActiveLine());
      
      for (size_t dy=0; dy<TILE_SIZE; dy++)
      {
        size_t y = tile_y+dy;
        if (y >= swapchain->getHeight()) continue;

        if (!swapchain->activeLine(y)) continue;
        size_t _y = swapchain->raster2buffer(y);
        
        for (size_t dx=0; dx<TILE_SIZE; dx++)
        {
          size_t x = tile_x+dx;
          if (x >= swapchain->getWidth()) continue;

          const int set = randomNumberGenerator.getInt(renderer->samplers->sampleSets);

          Color L = zero;
          size_t spp = renderer->samplers->samplesPerPixel;
          for (size_t s=0; s<spp; s++)
          {
            PrecomputedSample& sample = renderer->samplers->samples[set][s];
            const float fx = (float(x) + sample.pixel.x)*rcpWidth;
            const float fy = (float(y) + sample.pixel.y)*rcpHeight;

            Ray primary; camera->ray(Vec2f(fx,fy), sample.getLens(), primary);
            primary.time = sample.getTime();
            
            state.sample = &sample;
            state.pixel = Vec2f(fx,fy);
            L += renderer->integrator->Li(primary, scene, state);
          }
          const Color L0 = swapchain->update(x, _y, L, spp, accumulate);
          const Color L1 = toneMapper->eval(L0,x,y,swapchain);
          framebuffer->set(x, _y, L1);
        }
      }
      
      /*! print progress bar */
      if (renderer->showProgress) progress.next();

      /*! mark one more tile as finished */
      framebuffer->finishTile();
    }

    /*! we access the atomic ray counter only once per tile */
    atomicNumRays += state.numRays;
  }
}
