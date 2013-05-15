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

#include "debugrenderer.h"
#include "math/random.h"

namespace embree
{
  DebugRenderer::DebugRenderer(const Parms& parms) : accumulate(false) {
    maxDepth = parms.getInt("maxDepth",1);
    spp = parms.getInt("sampler.spp",1);
  }

  void DebugRenderer::renderThread()
  {
    size_t numRays = 0;

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
        for (size_t dx=0; dx<TILE_SIZE; dx++)
        {
          /*! ignore tile pixels outside framebuffer */
          size_t ix = x0+dx, iy = y0+dy;
          if (ix >= framebuffer->getWidth() || iy >= framebuffer->getHeight()) continue;

	  for (int s=0; s<spp; s++)
	  {	
	    /*! create primary ray */
	    DifferentialGeometry hit;
	    Ray ray; camera->ray(Vec2f(ix*rcpWidth,iy*rcpHeight), Vec2f(rand.getFloat(),rand.getFloat()), ray);
	    
	    for (ssize_t depth=0; depth<maxDepth; depth++)
	    {
	      /*! shoot current ray */
	      new (&hit) DifferentialGeometry;
	      scene->intersector->intersect(ray,hit);
	      scene->postIntersect(ray,hit);
	      numRays++;
	      if (!hit) break;
	      
	      /*! compute new ray through diffuse bounce */
	      Vec3f Nf = hit.Ng;
	      if (dot(-ray.dir,Nf) < 0) Nf = -Nf;
	      if (depth+1<maxDepth) {
		ray = Ray(ray.org+0.999f*hit.t*ray.dir,cosineSampleHemisphere(rand.getFloat(),rand.getFloat(),Nf),4.0f*float(ulp)*hit.error);
	      }
	    }
	      
	    /*! update framebuffer */
	    if (!hit) framebuffer->get(ix,iy) = zero;
	    //else      framebuffer->get(ix,iy) = Vec4f(-dot(ray.dir,hit.Ns),0,0,1.0f);
	    else      framebuffer->get(ix,iy) = Vec4f(hit.u,hit.v,1.0f-hit.u-hit.v,1.0f);
	    //else      framebuffer->get(ix,iy) = Vec4f(hit.st.x,hit.st.y,1.0f-hit.st.x-hit.st.y,1.0f);
	    //else      framebuffer->get(ix,iy) = Vec4f(hit.st.x,0,0,1.0f);
	    //else      framebuffer->get(ix,iy) = Vec4f(abs(dot(ray.dir,normalize(hit.dPds))),0,0,1.0f);
	    //framebuffer->get(ix,iy) = Vec4f(ix*rcpWidth,iy*rcpHeight,0.0f,0.0f);
	  }
        }
      }
      framebuffer->finishTile();
    }

    /*! we access the atomic ray counter only once per tile */
    atomicNumRays += numRays;
  }

  void DebugRenderer::renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain<> > swapchain, int accumulate)
  {
    /*! precompute some values */
    numTilesX = ((int)swapchain->getWidth() +TILE_SIZE-1)/TILE_SIZE;
    numTilesY = ((int)swapchain->getHeight()+TILE_SIZE-1)/TILE_SIZE;
    rcpWidth  = rcp(float(swapchain->getWidth()));
    rcpHeight = rcp(float(swapchain->getHeight()));

    /*! render frame */
    double t = getSeconds();
    this->accumulate = accumulate;
    this->tileID = 0;
    this->atomicNumRays = 0;
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
    this->camera = null;
    this->scene = null;
    this->toneMapper = null;
    this->framebuffer = null;
    this->swapchain = null;
    double dt = getSeconds()-t;

    /*! print framerate */
    std::cout << "render  " << rcp(dt) << " fps, " << dt*1000.0f << " ms, " << atomicNumRays/dt*1E-6 << " Mrps" << std::endl;
  }
}
