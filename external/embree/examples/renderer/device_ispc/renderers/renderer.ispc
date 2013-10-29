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

export void Renderer__renderFrameInit(void* uniform _this,
                                      void* uniform scene)
{
  uniform Renderer* uniform this = (uniform Renderer* uniform) _this;
  if (!this->renderFrameInit) return;
  this->renderFrameInit(this, (const uniform Scene* uniform) scene);
}

export uniform int Renderer__renderFrame(void* uniform _this,
                                         void* uniform camera,
                                         void* uniform scene,
                                         void* uniform toneMapper,
                                         void* uniform swapchain,
                                         const uniform int accuMode)
{
  uniform Renderer* uniform this = (uniform Renderer* uniform) _this;
  return this->renderFrame(this,
                           (const uniform Camera* uniform) camera,
                           (const uniform Scene* uniform) scene,
                           (const uniform ToneMapper* uniform) toneMapper,
                           (uniform SwapChain *uniform) swapchain,
                           accuMode);
}

export uniform bool Renderer__pick(void* uniform _camera,
                                   uniform float x,
                                   uniform float y,
                                   void* uniform _scene,
                                   uniform vec3f& p)
{
  const uniform Camera* uniform camera = (const uniform Camera* uniform) _camera;
  const uniform Scene*  uniform scene  = (const uniform Scene*  uniform) _scene;
  
  Ray ray;
  camera->initRay(camera, ray,make_vec2f(x,y), make_vec2f(0.0f,0.0f));
  scene->intersector->intersect(scene->intersector,ray);

  const vec3f hit_p = add(ray.org,mul(ray.dir,ray.tfar));
  p = make_vec3f(extract(hit_p.x,0),extract(hit_p.y,0),extract(hit_p.z,0));
  return extract((int)hadHit(ray),0);
}
