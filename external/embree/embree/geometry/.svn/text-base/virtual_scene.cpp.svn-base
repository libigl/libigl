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

#include "virtual_scene.h"

namespace embree
{
  AccelRegistry VirtualScene::accels; 

  BuilderRegistry VirtualScene::builders; 

  IntersectorRegistry<Intersector1> VirtualScene::intersectors1;
    
#if defined(__SSE__)
  IntersectorRegistry<Intersector4> VirtualScene::intersectors4;
#endif
    
#if defined(__AVX__)
  IntersectorRegistry<Intersector8> VirtualScene::intersectors8;
#endif
    
#if defined(__MIC__)
  IntersectorRegistry<Intersector16> VirtualScene::intersectors16;
#endif
    
  VirtualScene::VirtualScene (size_t numObjects, const char* accelTy) 
    : RTCGeometry(empty), numObjects(numObjects) 
  {
    objects = (Object*) alignedMalloc(numObjects*sizeof(Object));
    for (size_t i=0; i<numObjects; i++) new (&objects[i]) Object;
    accel = accels.create(accelTy,this);
  }
  
  VirtualScene::~VirtualScene () {
    alignedFree(objects);
  }

  void VirtualScene::clearRegistry ()
  {
    accels.clear();
    builders.clear();
    intersectors1.clear();
#if defined(__SSE__)
    intersectors4.clear();
#endif
#if defined(__AVX__)
    intersectors8.clear();
#endif
#if defined(__MIC__)
    intersectors16.clear();
#endif
  }

  size_t VirtualScene::size() const {
    return numObjects;
  }
    
  BBox3f VirtualScene::bounds(size_t i) const {
    return objects[i].worldBounds;
  }

  void VirtualScene::build(TaskScheduler::Event* event, std::string builderName) {
    builders.build(builderName,event,accel);
  }

  void VirtualScene::freeze() {
  }

  Intersector1* VirtualScene::intersector1(std::string travName) const {
    assert(accel);
    return intersectors1.get(accel->name(),travName,accel);
  }

#if defined(__SSE__)
  Intersector4* VirtualScene::intersector4(std::string travName) const {
    assert(accel);
    return intersectors4.get(accel->name(),travName,accel);
  }
#endif

#if defined(__AVX__)
  Intersector8* VirtualScene::intersector8(std::string travName) const {
    assert(accel);
    return intersectors8.get(accel->name(),travName,accel);
  }
#endif

#if defined(__MIC__)
  Intersector16* VirtualScene::intersector16(std::string travName) const {
    assert(accel);
    return intersectors16.get(accel->name(),travName,accel);
  }
#endif
}
