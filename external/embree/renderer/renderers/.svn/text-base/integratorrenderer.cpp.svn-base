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
    : iteration(0), accumulate(false), tileID(0), finishedTiles(0)
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
  
  void IntegratorRenderer::renderThread()
  {
    /*! create a new sampler */
    size_t numRays = 0;
    Sampler* sampler = samplers->create();

    /*! tile pick loop */
    while (true)
    {
      /*! pick a new tile */
      size_t tile = tileID++;
      if (tile >= numTilesX*numTilesY) break;

      /*! compute tile pixel range */
      Vec2i start((int)((tile%numTilesX)*TILE_SIZE),(int)((tile/numTilesX)*TILE_SIZE));
      Vec2i end (min(int(framebuffer->getWidth()),int(start.x+TILE_SIZE))-1,min(int(framebuffer->getHeight()),int(start.y+TILE_SIZE))-1);

      /*! configure the sampler with the tile pixels */
      sampler->init(Vec2i((int)framebuffer->getWidth(), (int)framebuffer->getHeight()), start, end, iteration);
      framebuffer->clear(start,end);

      /*! process all tile samples */
      while (!sampler->finished()) 
      {
        PrecomputedSample sample; sampler->proceed(sample);
        Ray primary; camera->ray(sample.raster*Vec2f(rcpWidth,rcpHeight), sample.getLens(), primary);
        primary.time = sample.getTime();
        Col3f L = integrator->Li(primary, scene, sample, numRays);
        if (!finite(L.r+L.g+L.b) || L.r < 0 || L.g < 0 || L.b < 0) L = zero;
        framebuffer->accumulate(sample.getIntegerRaster(), start, end, L, 1.0f);
      }
    
      /*! set or accumulate next buffer */
      swapchain->accumulate(framebuffer,start,end,accumulate);
      //swapchain->accumulateAndNormalize(framebuffer,start,end,accumulate);
 
      /*! run tonemapper */
      toneMapper->eval(framebuffer,framebuffer,start,end);

      /*! print progress bar */
      if (showProgress) progress.next();

      /*! mark one more tile as finished */
      framebuffer->finishTile();
    }

    /*! we access the atomic ray counter only once per tile */
    atomicNumRays += numRays;
    delete sampler;
  }

  void IntegratorRenderer::renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain<> > swapchain, int accumulate)
  {
    /*! precompute some values */
    numTilesX = ((int)swapchain->getWidth() +TILE_SIZE-1)/TILE_SIZE;
    numTilesY = ((int)swapchain->getHeight()+TILE_SIZE-1)/TILE_SIZE;
    rcpWidth  = 1.0f/float(swapchain->getWidth());
    rcpHeight = 1.0f/float(swapchain->getHeight());
    if (accumulate == 0) iteration = 0;
    if (showProgress) new (&progress) Progress(numTilesX*numTilesY);

    /*! render frame */
    if (showProgress) progress.start();
    double t = getSeconds();
    this->accumulate = accumulate;
    this->tileID = 0;
    this->finishedTiles = 0;
    this->atomicNumRays = 0;
    this->samplers->reset();
    this->integrator->requestSamples(this->samplers, scene);
    this->samplers->init(iteration,filter);
    this->camera = camera;
    this->scene = scene;
    this->toneMapper = toneMapper;
    this->framebuffer = swapchain->buffer();
    this->framebuffer->startRendering(numTilesX*numTilesY);
    this->swapchain = swapchain;
    scheduler->start();
    scheduler->addTask(TaskScheduler::ThreadInfo(),TaskScheduler::GLOBAL_FRONT,
                       (TaskScheduler::runFunction)&run_renderThread,this,scheduler->getNumThreads(),NULL,NULL,"render::tile");
    scheduler->stop();
    iteration++;
    this->camera = null;
    this->scene = null;
    this->toneMapper = null;
    this->framebuffer = null;
    this->swapchain = null;
    double dt = getSeconds()-t;
    if (showProgress) progress.end();

     /*! print fps, render time, and rays per second */
    std::ostringstream stream;
    stream << "render  ";
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(2);
    stream << 1.0f/dt << " fps, ";
    stream.precision(0);
    stream << dt*1000.0f << " ms, ";
    stream.precision(3);
    stream << atomicNumRays/dt*1E-6 << " Mrps";
    std::cout << stream.str() << std::endl;
  }
}
