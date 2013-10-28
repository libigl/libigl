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

#include "debugrenderer.h"
#include "math/random.h"

namespace embree
{
  DebugRenderer::DebugRenderer(const Parms& parms)  
  {
    maxDepth = parms.getInt("maxDepth",1);
    spp      = parms.getInt("sampler.spp",1);
  }

  void DebugRenderer::renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate) 
  {
    new RenderJob(this,camera,scene,toneMapper,swapchain,accumulate);
  }
  
  DebugRenderer::RenderJob::RenderJob (Ref<DebugRenderer> renderer, const Ref<Camera>& camera, const Ref<BackendScene>& scene, 
                                       const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate)
    : renderer(renderer), camera(camera), scene(scene), toneMapper(toneMapper), swapchain(swapchain), accumulate(accumulate)
  {
    this->tileID = 0;
    this->atomicNumRays = 0;
    numTilesX = ((int)swapchain->getWidth() +TILE_SIZE-1)/TILE_SIZE;
    numTilesY = ((int)swapchain->getHeight()+TILE_SIZE-1)/TILE_SIZE;
    rcpWidth  = rcp(float(swapchain->getWidth()));
    rcpHeight = rcp(float(swapchain->getHeight()));
    this->framebuffer = swapchain->buffer();
    this->framebuffer->startRendering(numTilesX*numTilesY);

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

  void DebugRenderer::RenderJob::finish(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    double dt = getSeconds()-t0;
    std::cout.precision(3);
    std::cout << "render  " << rcp(dt) << " fps, " << dt*1000.0f << " ms, " << atomicNumRays/dt*1E-6 << " Mrps" << std::endl;
    delete this;
  }

  void DebugRenderer::RenderJob::renderTile(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    size_t numRays = 0;
    if (taskIndex == taskCount-1) t0 = getSeconds();

    /*! tile pick loop */
    while (true)
    {
      /*! pick a new tile */
      size_t tile = tileID++;
      if (tile >= numTilesX*numTilesY) break;

      /*! compute tile location */
      Random rand(int(tile)*1024);
      size_t x0 = (tile%numTilesX)*TILE_SIZE;
      size_t y0 = (tile/numTilesX)*TILE_SIZE;

      Vec2i start((int)x0,(int)y0);
      Vec2i end (min(int(framebuffer->getWidth()),int(start.x+TILE_SIZE))-1,min(int(framebuffer->getHeight()),int(start.y+TILE_SIZE))-1);

      /*! loop over all pixels of the tile */
      for (size_t dy=0; dy<TILE_SIZE; dy++)
      {
        size_t iy = y0+dy; float fy = iy*rcpHeight;
        if (iy >= framebuffer->getHeight()) continue;
        if (!swapchain->activeLine(iy)) continue;
        iy = swapchain->raster2buffer(iy);

        for (size_t dx=0; dx<TILE_SIZE; dx++)
        {
          /*! ignore tile pixels outside framebuffer */
          size_t ix = x0+dx; float fx = ix*rcpWidth;
          if (ix >= framebuffer->getWidth()) continue;

	  for (size_t i=0; i<renderer->spp; i++)
	  {
	    /*! create primary ray */
	    Ray ray; camera->ray(Vec2f(fx,fy), Vec2f(rand.getFloat(),rand.getFloat()), ray);
	    
	    for (size_t depth=0; depth<renderer->maxDepth; depth++)
	    {
	      /*! shoot current ray */
              scene->intersector->intersect(ray);
              numRays++;
	      if (!ray) break;
	      
	      /*! compute new ray through diffuse bounce */
	      Vector3f Nf = normalize(ray.Ng);
	      if (dot(-ray.dir,Nf) < 0) Nf = -Nf;
	      if (depth+1<renderer->maxDepth) new (&ray) Ray(ray.org+0.999f*ray.tfar*ray.dir,cosineSampleHemisphere(rand.getFloat(),rand.getFloat(),Nf),4.0f*float(ulp)/**hit.error*/);
	    }  

	    /*! update framebuffer */
	    if (!ray) framebuffer->set(ix,iy,one);
	    //else      framebuffer->get(ix,iy) = Vec4f(abs(hit.Ns.x),abs(hit.Ns.y),abs(hit.Ns.z),1.0f);
            /*else      {
              hit.Ng = normalize(hit.Ng);
              framebuffer->set(ix,iy,Color(abs(hit.Ng.x),abs(hit.Ng.y),abs(hit.Ng.z)));
              }*/
	    //else      framebuffer->set(ix,iy,Color(clamp(dot(ray.dir,normalize(hit.Ng))),0,0));
	    //else      framebuffer->set(ix,iy,Color(hit.u,hit.v,1.0f-hit.u-hit.v));
	    //else      framebuffer->get(ix,iy) = Vec4f(hit.st.x,hit.st.y,1.0f-hit.st.x-hit.st.y,1.0f);
	    //else      framebuffer->get(ix,iy) = Vec4f(hit.st.x,0,0,1.0f);
	    //else      framebuffer->get(ix,iy) = Vec4f(abs(dot(ray.dir,normalize(hit.dPds))),0,0,1.0f);
	    else framebuffer->set(ix,iy,Color(((3434553*((unsigned)(ray.id0+ray.id1+3243)))%255)/255.0f,
                                              ((7342453*((unsigned)(ray.id0+ray.id1+8237)))%255)/255.0f,
                                              ((9234454*((unsigned)(ray.id0+ray.id1+2343)))%255)/255.0f));
	  }
        }
      }
      framebuffer->finishTile();
    }

    /*! we access the atomic ray counter only once per tile */
    atomicNumRays += numRays;
  }
}
