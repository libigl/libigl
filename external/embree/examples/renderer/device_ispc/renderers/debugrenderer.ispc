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

#include "renderer.isph"

#define TILE_SIZE 16

#if __IVL_SSE4 
#  define WARP_WIDTH 2
#else
#  define WARP_WIDTH 4
#endif

uniform int DebugRenderer__frameID = 0;

struct DebugRenderer {
  uniform Renderer base;
  uniform int maxDepth;
  uniform int spp;
};

vec3f randomID(const int ID) {
  int r = ((ID+13)*17*23) & 255;
  int g = ((ID+15)*11*13) & 255;
  int b = ((ID+17)* 7*19) & 255;
  const float oneOver255f = 1.f/255.f;
  return make_vec3f(r*oneOver255f,g*oneOver255f,b*oneOver255f);
}

task void DebugRenderer__renderTile(uniform DebugRenderer* uniform this,
                                    const uniform int numTiles_x,
                                    const uniform Camera      *uniform camera,
                                    const uniform Scene       *uniform scene,
                                    uniform FrameBuffer *uniform fb)
{
  const uniform int tileID = taskIndex;
  const uniform int tile_y = tileID / numTiles_x;
  const uniform int tile_x = tileID - tile_y * numTiles_x;
  const uniform int x0 = tile_x * TILE_SIZE;
  const uniform int x1 = min(x0+TILE_SIZE,(int)fb->size.x);
  const uniform int y0 = tile_y * TILE_SIZE;
  const uniform int y1 = min(y0+TILE_SIZE,(int)fb->size.y);
  for (uniform int y=y0;y<y1;y++) {
    if (!activeLine(y)) continue;
    size_t _y = raster2buffer(y);
    for (varying int x=x0+programIndex;x<x1;x+=programCount) {
     
      for (uniform int s=0; s<this->spp; s++) 
      {
      vec3f col;
#if 0
      // -------------------------------------------------------
      // ray cast, display eyelight(gnormal)
      // -------------------------------------------------------
      const float fx = (x+.5f)*fb->invSize.x;
      const float fy = (y+.5f)*fb->invSize.y;
      Ray ray;
      camera->initRay(camera,ray,
                            make_vec2f(fx,fy),
                            make_vec2f(0.f,0.f));
      scene->intersector->intersect(scene->intersector,ray,hit);
      vec3f Ng = make_vec3f(0,0,0);
      if (hadHit(ray)) Ng = normalize(ray.Ng);
      col = make_vec3f(minf(1.f,.3f+.8f*absf(dot(Ng,normalize(ray.dir)))));
      // col = randomID(ray.meshID+ray.triID);
#elif 1
      // -------------------------------------------------------
      // ray cast, but only display triangle ID
      // -------------------------------------------------------
      const float fx = (x+.5f)*fb->invSize.x;
      const float fy = (y+.5f)*fb->invSize.y;
      Ray ray; 
      camera->initRay(camera,ray,
                            make_vec2f(fx,fy),
                            make_vec2f(0.f,0.f));
      ray.u = ray.v = 0.0f;
      scene->intersector->intersect(scene->intersector,ray);

      col = make_vec3f(ray.u,ray.v,1.0f-ray.u-ray.v);
      //col = randomID(ray.geomID+ray.primID);
#elif 0
      // -------------------------------------------------------
      // ray cast, but only display shading normal
      // -------------------------------------------------------
      const float fx = (x+.5f)*fb->invSize.x;
      const float fy = (y+.5f)*fb->invSize.y;
      Ray ray;
      camera->initRay(camera,ray,
                            make_vec2f(fx,fy),
                            make_vec2f(0.f,0.f));
      DifferentialGeometry dg;
      scene->intersector->intersect(scene->intersector,ray); 
      //postIntersect(scene,ray,dg);
      vec3f Ng = make_vec3f(0,0,0);
      if (hadHit(ray)) Ng = normalize(ray.Ng);
      col = make_vec3f(abs(Ng.x),abs(Ng.y),abs(Ng.z));
#elif 0
      // -------------------------------------------------------
      // visualize ray direction
      // -------------------------------------------------------
      const float fx = (x+.5f)*fb->invSize.x;
      const float fy = (y+.5f)*fb->invSize.y;
      Ray ray;
      camera->initRay(camera,
			    ray,
                            make_vec2f(fx,fy),
                            make_vec2f(0.f,0.f));
      col = abs(normalize(ray.dir));
#else
       // -------------------------------------------------------
      // simplest: just a simple pattern
      // -------------------------------------------------------
      col.x = ((x+DebugRenderer__frameID)&255)/255.f;
      col.y = ((y+DebugRenderer__frameID)&255)/255.f;
      col.z = ((x+y+DebugRenderer__frameID)&255)/255.f;
#endif
      fb->set(fb,x,_y,col);
      }
    }
  }
}

uniform int DebugRenderer_renderFrame(uniform Renderer* uniform _this,
                                      const uniform Camera* uniform camera,
                                      const uniform Scene* uniform scene,
                                      const uniform ToneMapper* uniform toneMapper,
                                      uniform SwapChain *uniform swapchain,
                                      const uniform int accuMode)
{
  uniform DebugRenderer* uniform this = (uniform DebugRenderer* uniform) _this;
  uniform int numTiles_x = (swapchain->width +(TILE_SIZE-1))/TILE_SIZE;
  uniform int numTiles_y = (swapchain->height+(TILE_SIZE-1))/TILE_SIZE;
  uniform int numTiles = numTiles_x * numTiles_y;
  
  launch[numTiles] DebugRenderer__renderTile(this,numTiles_x,camera,scene,SwapChain__get_buffer(swapchain));
  sync;
  rtcDebug();
  return swapchain->width*swapchain->height*this->spp;
}

void DebugRenderer__Constructor(uniform DebugRenderer* uniform this,
                                const uniform int& maxDepth, 
                                const uniform int& spp)
{
  Renderer__Constructor(&this->base,Renderer__Destructor,NULL,DebugRenderer_renderFrame);
  this->maxDepth = maxDepth;
  this->spp = spp;
}

export void* uniform DebugRenderer__new(const uniform int& maxDepth, 
                                        const uniform int& spp)
{
  uniform DebugRenderer *uniform this = uniform new uniform DebugRenderer;
  DebugRenderer__Constructor(this,maxDepth,spp);
  return this;
}
