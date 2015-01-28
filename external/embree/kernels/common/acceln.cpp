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

#include "acceln.h"
#include "embree2/rtcore_ray.h"

namespace embree
{
  AccelN::AccelN () : N(0) { 
    memset(accels,0,sizeof(accels));
  }

  AccelN::~AccelN() 
  {
    for (size_t i=0; i<N; i++)
      delete accels[i];
  }

  void AccelN::add(Accel* accel) 
  {
    assert(N<16);
    assert(accel);
    if (N<16) accels[N++] = accel;
  }
  
  void AccelN::intersect (void* ptr, RTCRay& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++)
      This->accels[i]->intersect(ray);
  }

  void AccelN::intersect4 (const void* valid, void* ptr, RTCRay4& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++)
      This->accels[i]->intersect4(valid,ray);
  }

  void AccelN::intersect8 (const void* valid, void* ptr, RTCRay8& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++)
      This->accels[i]->intersect8(valid,ray);
  }

  void AccelN::intersect16 (const void* valid, void* ptr, RTCRay16& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++)
      This->accels[i]->intersect16(valid,ray);
  }

  void AccelN::occluded (void* ptr, RTCRay& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++) {
      This->accels[i]->occluded(ray); 
      if (ray.geomID == 0) break;
    }
  }

  void AccelN::occluded4 (const void* valid, void* ptr, RTCRay4& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++) {
      This->accels[i]->occluded4(valid,ray);
      // FIXME: exit if already hit something
    }
  }

  void AccelN::occluded8 (const void* valid, void* ptr, RTCRay8& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++) {
      This->accels[i]->occluded8(valid,ray);
      // FIXME: exit if already hit something
    }
  }

  void AccelN::occluded16 (const void* valid, void* ptr, RTCRay16& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->N; i++) {
      This->accels[i]->occluded16(valid,ray);
      // FIXME: exit if already hit something
    }
  }

  void AccelN::print(size_t ident)
  {
    for (size_t i=0; i<N; i++)
    {
      for (size_t j=0; j<ident; j++) std::cout << " "; 
      std::cout << "accels[" << i << "]" << std::endl;
      accels[i]->intersectors.print(ident+2);
    }
  }

  void AccelN::immutable()
  {
    for (size_t i=0; i<N; i++)
      accels[i]->immutable();
  }

  void AccelN::build (size_t threadIndex, size_t threadCount) 
  {
    size_t validAccelIndex = 0;
    size_t validAccelCount = 0;
    
    /* build all acceleration structures */
    for (size_t i=0; i<N; i++) 
    {
      accels[i]->build(threadIndex,threadCount);
      if (accels[i]->bounds.empty()) continue;
      validAccelIndex = i;
      validAccelCount++;
    }

    if (validAccelCount == 1) {
      intersectors = accels[validAccelIndex]->intersectors;
    }
    else 
    {
      intersectors.ptr = this;
      intersectors.intersector1 = Intersector1(&intersect,&occluded,"AccelN::intersector1");
      intersectors.intersector4 = Intersector4(&intersect4,&occluded4,"AccelN::intersector4");
      intersectors.intersector8 = Intersector8(&intersect8,&occluded8,"AccelN::intersector8");
      intersectors.intersector16= Intersector16(&intersect16,&occluded16,"AccelN::intersector16");
    }
    
    /*! calculate bounds */
    bounds = empty;
    for (size_t i=0; i<N; i++) 
      bounds.extend(accels[i]->bounds);
  }
}
